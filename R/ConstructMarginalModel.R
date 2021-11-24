#' ConstructMarginalModel -- the main function that will estimate the significance of the marginal model fitted to a CIFTI dscalar dataset
#'
#' This function wraps all of the other functions in the MarginalModelCifti package. The output is a significance map, which is either uncorrected or cluster detection was performed.
#' @param external_df A data frame comprising non-brain measures to model. Can be specified as a string to a csv file with appropriate headers.
#' @param concfile A character string denoting a single column text file that lists the dscalars in the same order as the external_df and wave frames.
#' @param structtype A character string denoting whether the map is volumetric ('volume'), surface-based ('surface'), a pconn ('pconn'), or a NIFTI connectivity matrix ('niiconn').
#' @param structfile A character string denoting the structural map file, used for cluster detection on surfaces only.
#' @param matlab_path A character string denoting the path to the matlab compiler. Please use v91.
#' @param surf_command A character string denoting the path to the SurfConnectivity command.
#' @param wave A data frame comprising the order of the waves (longitudinal measures) for each independent case.
#' @param notation A character string representing the model formula. Uses Wilkinson notation.
#' @param zcor A matrix defining the correlation structure set by the user, can be used to specify custom dependencies (e.g. site or kinship matrices).
#' @param corstr A character string that defines the correlation structure in a predetermined way. "Independence" is recommend for most cases.
#' @param family_dist A character string denoting the assumed underlying distribution of all input data.
#' @param dist_type A character string denoting the distribution to use for the wild bootstrap.
#' @param z_thresh A numeric that represents the threshold for Z statistics when performing cluster detection.
#' @param nboot A numeric that represents the number of wild bootstraps to perform.
#' @param p_thresh A numeric that represents the p value threshold for assessing significance.
#' @param sigtype A character string denoting cluster ('cluster'),  point ('point'), or enrichment ('enrichment') comparisons.
#' @param id_subjects A character string denoting the column for the subject ID. Only needed when FastSwE is FALSE
#' @param output_directory A character string denoting the path to the output MRI statistical maps
#' @param ncores An integer denoting the number of CPU cores to use when conducting permutation tests
#' @param fastSWE A boolean that determine the sandwhich estimator approach. If set to FALSE, will use standard R package geeglm. If set to TRUE will use custom-built estimator using rfast.
#' @param adjustment A character string denoting the small sample size adjustment to use when fastSwE is set to TRUE. Is NULL by default.
#' @param norm_external_data A boolean. If set to true, external data will be normed prior to analysis.
#' @param norm_internal_data A boolean. If set to true, MRI data will be normed per datapoint prior to analysis.
#' @param marginal_outputs A boolean. If set to true, marginal values will be output as statistical maps
#' @param marginal_matrix A numeric matrix depicting how to draw the map. Only needed if marginal_outputs is set to TRUE
#' @param enrichment_path A string depicting the path to the enrichment code. Used for enrichment analysis only (i.e. when sigtype is set to 'enrichment').
#' @param modules A csv file or array that depicts the modules for the enrichment analysis.
#' @param wb_command A character string denoting the path to the wb_command file
#' @param subsetfile A character string denoting the path to a subset selection file, to select only a subset of the participants
#' @param permutation_directory A character string denoting the path to the directory containing ONLY permutation tests. If set to NULL, permutation tests will be generated on the fly
#' @param analysismode A character string denoting how much of the analysis to perform. Can be set to "full" to include cluster detection and permutation testing, or set to "statmaponly" to only output statistical maps
#' @#' @keywords wild bootstrap sandwich estimator marginal model CIFTI scalar
#' @export
#' @examples
#' max_cc <- ComputeMM_WB(cifti_map,zscore_map,resid_map,fit_map,type,external_df,
#' notation,family_dist,structtype,thresh,structfile,matlab_path,surf_command,correctiontype)
ConstructMarginalModel <- function(external_df,
                                   concfile,
                                   structtype,
                                   structfile,
                                   matlab_path,
                                   surf_command,
                                   wave = NULL,
                                   notation=formula(y~x),
                                   zcor = NULL,
                                   corstr = NULL,
                                   family_dist = NULL,
                                   dist_type="radenbacher",
                                   z_thresh = 2.3,
                                   nboot=1000,
                                   p_thresh = 0.5,
                                   sigtype = 'cluster',
                                   id_subjects='subid',
                                   output_directory='~/',
                                   ncores=1,
                                   fastSwE = TRUE,
                                   adjustment = NULL,
                                   norm_external_data = TRUE,
                                   norm_internal_data = TRUE,
                                   marginal_outputs = FALSE,
                                   marginal_matrix = NULL,
                                   enrichment_path = NULL,
                                   modules = NULL,
                                   wb_command = 'wb_command',
                                   subsetfile = NULL,
                                   permutation_directory = NULL,
                                   analysismode = 'full',
                                   checkorder = TRUE,
                                   roisubset=NULL){
  initial_time = proc.time()
  require(purrr)
  require(cifti)
  require(gifti)
  require(geepack)
  require(Matrix)
  require(oro.nifti)
  require(mmand)
  require(parallel)
  require(Rfast)
  require(stringr)
  require(R.matlab)
  require(SparseM)
  require(raveio)
  cifti_firstsub= NULL
  zeros_array = NULL
  curr_directory = getwd()
  setwd(output_directory)
  print("reading imaging list")
  if (structtype == 'dataframe') {
    ciftilist = list()
    ciftilist$file <- concfile
  } else
  {
    ciftilist <- read.csv(concfile,header=FALSE,col.names="file")
    if (is.null(subsetfile) == FALSE){
      subset_list <- read.csv(subsetfile,header=FALSE,col.names="subset")
      print("parsing subset list")
      list_bin_expanded <- sapply(1:length(subset_list$subset),
                      function(x) {
                        SubSelection(as.character(ciftilist$file),
                                     as.character(subset_list$subset[x])
                                     )
                        }
                      )
      list_bin = as.logical(rowsums(list_bin_expanded))
      new_ciftilist = list()
      new_ciftilist$file = ciftilist$file[list_bin]
      ciftilist = as.data.frame(new_ciftilist)
    }
  }
  print("loading non-imaging data")
  external_df_file = FALSE
  if (is.character(external_df)) {
    external_df_file = TRUE
    external_df <- read.csv(external_df,header=TRUE)
    if (is.null(subsetfile) == FALSE){
      external_df <- external_df[list_bin,]
    }
  }
  wave_file = FALSE
  if (is.character(wave)) {
    print("parsing longitudinal data")
    wave_file = TRUE
    wave <- read.csv(wave,header=TRUE)
    if (is.null(subsetfile) == FALSE){
      wave = wave[list_bin,]
    }
  }
  if (checkorder == TRUE) {
    print('checking conc file order against external data file')
    list_order_expanded <- sapply(1:length(external_df[[id_subjects]]),
                                function(x) {
                                  SubSelection(as.character(ciftilist$file),
                                               as.character(external_df[[id_subjects]][x])
                                  )
                                }
    )
    list_order = as.logical(rowsums(list_order_expanded))
    if (sum(list_order) == 0) {
      print('no subjects found -- perhaps your conc file path does not contain a matching subject id pattern?')
    }
    if (sum(list_order) > 0) {
      new_ciftilist = list()
      new_ciftilist$file = ciftilist$file[list_order]
      ciftilist = as.data.frame(new_ciftilist)
    }
    print('conc order check completed')
    if (wave_file) {
      print('comparing external data file and nested wave file')
      if (is.null(wave[[id_subjects]]) == TRUE){
        print('could not find matching subject id column in longitudinal file -- will assume order is correct but please check!!!!')
        new_wave <- wave
      }
      if (is.null(wave[[id_subjects]]) == FALSE){
        print('matching id column found in longitudinal file -- proceeding with matching to external data for analysis')
        new_wave_temp <- wave[wave[[id_subjects]] %in% external_df[[id_subjects]],]
        new_wave <- new_wave_temp[,!names(new_wave_temp) %in% id_subjects]
      }
    }
  }
  if (checkorder == FALSE){
    print('checking orders has been disabled -- please check your files!')
    new_wave <- wave    
  }
  print("loading imaging data")
  start_load_time = proc.time()
  CiftiInputs <- LoadBrainMetrics(ciftilist,structtype,roisubset)
  Nelm = CiftiInputs$Nelm
  cifti_alldata =  CiftiInputs$cifti_alldata
  cifti_scalarmap = CiftiInputs$cifti_scalarmap
  cifti_firstsub = CiftiInputs$cifti_firstsub 
  zeros_array = CiftiInputs$zeros_array
  cifti_nonans = CiftiInputs$nonans
  cifti_dim = CiftiInputs$cifti_dim
  cifti_splitfile <- strsplit(as.character(ciftilist$file[1]))
  surf_template_filename <- cifti_splitfile[[1]][length(cifti_splitfile[[1]])]
  surf_template_file=paste(output_directory,'/',surf_template_filename)
  surf_template_copycommand = paste('cp ',as.character(ciftilist$file[1]), ' ', surf_template_filename)
  system(surf_template_copycommand)
  print("parse non-imaging data")
  external_df <- external_df[cifti_nonans,]
  df_nan <- FilterDFNA(external_df = external_df,notation = notation)
  external_df = external_df[df_nan,]
  if (structtype != "surface"){
    cifti_scalarmap = cifti_scalarmap[df_nan,]
  }
  if (structtype == "surface"){
    if (fastSwE == FALSE) {
      cifti_scalarmap = cifti_scalarmap[df_nan,]
    }
  }
  cifti_alldata = cifti_alldata[df_nan,]
  predictors <- attr(terms(notation),"term.labels")
  measnames <- c("intercept",predictors)    
  if (fastSwE == TRUE){
    external_df <- ParseDf(external_df = external_df,notation = notation,norm_data=norm_external_data)
    nmeas <- dim(external_df)[2]
  }
  if (wave_file) {
    wave <- DetermineNestedGroups(new_wave)
    wave = wave[cifti_nonans]
    wave = wave[df_nan]
  }
  finish_load_time = proc.time()-start_load_time
  cat("loading data complete. Time elapsed: ",finish_load_time[3],"s")  
  print("running marginal model on observed data")
  start_model_time = proc.time()
  if (fastSwE == FALSE){
    cifti_map <- lapply(cifti_scalarmap,ComputeMM,external_df=external_df,notation=notation,family_dist=family_dist,corstr=corstr,zcor=zcor,wave=wave,id_subjects=id_subjects)
    finish_model_time = proc.time() - start_model_time
    cat("modeling complete. Time elapsed: ",finish_model_time[3],"s")
    varlist <- all.vars(notation)
    nmeas <- length(varlist)
    start_normthresh_time = proc.time()
    print("Normalizing observed marginal model estimates")
    zscore_map <- map(cifti_map,ComputeZscores,nmeas)
    save(zscore_map,file = "zscore_observed.Rdata")
    resid_map <- map(cifti_map,ComputeResiduals,nmeas)
    fit_map <- map(cifti_map,ComputeFits,nmeas)
    print("thresholding observed z scores")
    thresh_map <- map(zscore_map,ThreshMap,zthresh=z_thresh)
    save(thresh_map,file = "zscore_thresh_observed.Rdata")
    finish_normthresh_time = proc.time() - start_normthresh_time
    cat("thresholding complete. Time elapsed: ", finish_normthresh_time[3],"s")
  } else
  {
    cifti_map <- lm.fit(external_df,cifti_alldata)
    beta_map <- cifti_map$coefficients
    SaveBWASfile(BWAS_statmap=beta_map,structtype=structtype,
                 surf_template_file=surf_template_file,
                 output_prefix='beta_map_',measnames=measnames,
                 surf_command=surf_command,wb_command=wb_command,
                 zeros_array=zeros_array,matlab_path=matlab_path,
                 output_directory=output_directory)
    if (marginal_outputs == TRUE){
      marginal_map = array(data = 0, dim = c(dim(marginal_matrix)[1],dim(beta_map)[2]))
      for (curr_marginal in 1:dim(marginal_matrix)[1]){
        marginal_map[curr_marginal,] <- sapply(1:dim(beta_map)[2],function(x){CalculateMarginalValue(beta_map[,x],marginal_matrix[curr_marginal,])})
      }
      SaveBWASfile(BWAS_statmap=marginal_map,structtype=structtype,
                   surf_template_file=surf_template_file,
                   output_prefix='marginal_map_',measnames=measnames,
                   surf_command=surf_command,wb_command=wb_command,
                   zeros_array=zeros_array,matlab_path=matlab_path,
                   output_directory=output_directory)
    }
    print("Finshed beta map outputs")
    resid_map <- cifti_map$residuals
    fit_map <- cifti_map$fitted.values
    Nelm <- dim(beta_map)[2]
    print("running fast sandwich estimator")
    t_map <- ComputeFastSwE(X=external_df,nested=wave,Nelm=Nelm,resid_map=resid_map,npredictors=nmeas,beta_map=beta_map,adjustment=adjustment)
    SaveBWASfile(BWAS_statmap=t_map,structtype=structtype,
                 surf_template_file=surf_template_file,
                 output_prefix='t_map_',measnames=measnames,
                 surf_command=surf_command,wb_command=wb_command,
                 zeros_array=zeros_array,matlab_path=matlab_path,
                 output_directory=output_directory)
    finish_model_time = proc.time() - start_model_time
    cat("modeling complete. Time elapsed: ",finish_model_time[3],"s")
    start_normthresh_time = proc.time()
    print("Normalizing observed marginal model estimates")
    zscore_map <- t(sapply(1:nmeas,function(x) {(t_map[x,] - mean(t_map[x,is.finite(t_map[x,])]))/sd(t_map[x,is.finite(t_map[x,])])}))
    SaveBWASfile(BWAS_statmap=zscore_map,structtype=structtype,
                 surf_template_file=surf_template_file,
                 output_prefix='zscore_map_',measnames=measnames,
                 surf_command=surf_command,wb_command=wb_command,
                 zeros_array=zeros_array,matlab_path=matlab_path,
                 output_directory=output_directory)
    print("thresholding observed z scores")
    thresh_map <- t(sapply(1:nmeas,function(x) abs(zscore_map[x,]) > z_thresh))
    thresh_map[is.na(thresh_map)] <- NaN
    SaveBWASfile(BWAS_statmap=thresh_map,structtype=structtype,
                 surf_template_file=surf_template_file,
                 output_prefix='thresh_map_',measnames=measnames,
                 surf_command=surf_command,wb_command=wb_command,
                 zeros_array=zeros_array,matlab_path=matlab_path,
                 output_directory=output_directory)      
    finish_normthresh_time = proc.time() - start_normthresh_time    
    cat("thresholding complete. Time elapsed: ", finish_normthresh_time[3],"s")
  }
  if (analysismode == 'statmaps'){
    all_maps = t_map
    return(all_maps)
  }
  if (analysismode == 'full'){    
    print("performing cluster detection")
    start_clust_time = proc.time()
    if (sigtype == 'cluster'){
      all_cc = matrix(data=NA,nrow=Nelm,ncol=nmeas)     
      for (curr_meas in 1:nmeas){
        thresh_array = unlist(thresh_map)
        mask_vector = 1:nmeas == curr_meas
        thresh_array <- thresh_array[mask_vector]
        thresh_array <- as.numeric(thresh_array)
        thresh_array[is.na(thresh_array)] <-  0
        if (structtype == 'volume'){
          thresh_vol <- RevertVolume(thresh_array,cifti_dim)
          observed_cc = GetVolAreas(thresh_vol)
          all_cc[,curr_meas] = ConvertVolume(observed_cc,cifti_dim)
        } else
        {
          observed_cc <- GetSurfAreas(thresh_array,structfile,matlab_path,surf_command)
          all_cc[,curr_meas] = observed_cc
        }
      }
      SaveBWASfile(BWAS_statmap=all_cc,structtype=structtype,
                   surf_template_file=surf_template_file,
                   output_prefix='observed_clusters_',measnames=measnames,
                   surf_command=surf_command,wb_command=wb_command,
                   zeros_array=zeros_array,matlab_path=matlab_path,
                   output_directory=output_directory)
    }
      if (sigtype == 'enrichment'){
        for (curr_meas in 1:nmeas){
          thresh_array = unlist(thresh_map)
          mask_vector = 1:nmeas == curr_meas
          thresh_array <- thresh_array[mask_vector]
          thresh_array <- as.numeric(thresh_array)
          thresh_array[is.na(thresh_array)] <-  0
          thresh_mat <- Matrix2Vector(pconn_data = zeros_array,
                                      pconn_vector = thresh_array,
                                      direction = "to_matrix")
          observed_cc <- EnrichmentAnalysis(metric_data = thresh_mat,
                                            ncols = dim(thresh_mat)[1],
                                            modules = modules, 
                                            enrichment_path = enrichment_path,
                                            matlab_path = matlab_path,
                                            output_file = paste(output_directory,'/','observed_chisqrd_',measnames[curr_meas],sep=""),
                                            tempname=paste(output_directory,'/','observed_',measnames[curr_meas],sep=""))
        if (curr_meas == 1) {
          all_cc = array(data=NA,dim=c(dim(observed_cc)[1],dim(observed_cc)[2],nmeas)) 
        }
        all_cc[,,curr_meas] <- unlist(observed_cc)
        }
      }
      finish_clust_time = proc.time() - start_clust_time
      cat("cluster detection complete. Time elapsed", finish_clust_time[3],"s")
      if(is.null(permutation_directory)){
        print("performing permutation test")
        start_perm_time = proc.time()
        seeds <- sample(.Machine$integer.max,size=nboot)
        cl <- makeForkCluster(nnodes = ncores)
        WB_cc <- parSapply(cl,1:nboot,ComputeMM_WB,resid_map=resid_map,
                           fit_map=fit_map,
                           type=dist_type,
                           external_df=external_df,
                           notation=notation,
                           family_dist=family_dist,
                           structtype=structtype,
                           thresh = z_thresh,
                           structfile=structfile,
                           matlab_path=matlab_path,
                           surf_command=surf_command,
                           corstr=corstr,
                           wave=wave,
                           zcor=zcor,
                           correctiontype = sigtype,
                           id_subjects=id_subjects,
                           cifti_dim=cifti_dim,
                           nmeas=nmeas,
                           seeds=seeds,
                           fastSwE=fastSwE,
                           adjustment=adjustment,
                           enrichment_path = enrichment_path,
                           modules = modules,
                           cifti_firstsub = zeros_array,
                           output_directory = output_directory)
        stopCluster(cl)
        finish_perm_time = proc.time() - start_perm_time
        cat("permutation testing complete. Time elapsed",finish_perm_time[3],"s")
      } else
        {
        print("aggregating p values from already performed permutation test")
        start_perm_time = proc.time()
        permutation_file_list <- dir(permutation_directory)
        count = 1
        for (permutation_file in permutation_file_list){
          if (count == 1){
            WB_cc = read.csv(paste(permutation_directory,permutation_file,sep='/'),header=FALSE)
          }
          WB_cc_temp = read.csv(paste(permutation_directory,permutation_file,sep='/'),header=FALSE)
          WB_cc = rbind(WB_cc,WB_cc_temp)      
        }
        finish_perm_time = proc.time() - start_perm_time
        cat("aggregation complete. Time elapsed",finish_perm_time[3],"s")
      }
      print("calculating p values for observed data using null distribution(s)")
      if (is.null(sigtype)) {
        for (curr_meas in 1:nmeas){
          pval_map <- map(zscore_map[curr_meas,],CalculatePvalue,WB_cc=NaN,nboot=NaN,sigtype=sigtype)
          if (curr_meas == 1){
            all_maps = list(pval_map)
          } else
          {
            all_maps <- append(all_maps,list(pval_map))
          }
        }        
      }      
      if (sigtype == 'cluster'){
        for (curr_meas in 1:nmeas){
          pval_map <- map(all_cc[,curr_meas],CalculatePvalue,WB_cc=WB_cc[curr_meas,],nboot=nboot,sigtype=sigtype)
          if (curr_meas == 1){
            all_maps = list(pval_map)
          } else
          {
            all_maps <- append(all_maps,list(pval_map))
          }
        }
      }
      if (sigtype == 'point'){
        for (curr_meas in 1:nmeas){
          pval_map <- map(zscore_map[curr_meas,],CalculatePvalue,WB_cc=WB_cc[curr_meas,],nboot=nboot,sigtype=sigtype)
          if (curr_meas == 1){
            all_maps = list(pval_map)
          } else
          {
            all_maps <- append(all_maps,list(pval_map))
          }
        }
      }
      if (sigtype == 'enrichment'){
        for (curr_meas in 1:nmeas){
          pval_map <- map(all_cc[,,curr_meas],CalculatePvalue,WB_cc=WB_cc[curr_meas,],nboot=nboot,sigtype=sigtype)
          if (curr_meas == 1){
            all_maps = list(pval_map)
          } else
          {
            all_maps <- append(all_maps,list(pval_map))
          }
        }
      }
    SaveBWASfile(BWAS_statmap=all_maps,structtype=structtype,
                 surf_template_file=surf_template_file,
                 output_prefix='pval_map_',measnames=measnames,
                 surf_command=surf_command,wb_command=wb_command,
                 zeros_array=zeros_array,matlab_path=matlab_path,sigtype=sigtype,
                 output_directory=output_directory)
    setwd(curr_directory)
    all_time = proc.time() - initial_time
    cat("ConstructMarginalModel complete. Time elapsed", all_time[3],"s")
    return(all_maps)
    surf_template_rmcommand = paste('rm ', surf_template_filename)
    system(surf_template_rmcommand)
  }
}

