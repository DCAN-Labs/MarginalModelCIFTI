#' ConstructMarginalModel -- the main function that will estimate the significance of the marginal model fitted to a CIFTI dscalar dataset
#'
#' This function wraps all of the other functions in the MarginalModelCifti package. The output is a significance map, which is either uncorrected or cluster detection was performed.
#' @param external_df A data frame comprising non-brain measures to model. Can be specified as a string to a csv file with appropriate headers.
#' @param concfile A character string denoting a single column text file that lists the dscalars in the same order as the external_df and wave frames.
#' @param structtype A character string denoting whether the map is volumetric ('volume') or surface-based.
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
#' @param sigtype A character string denoting cluster ('cluster') or point ('point') comparisons.
#' @param id_subjects A character string denoting the column for the subject ID. Only needed when FastSwE is FALSE
#' @param output_directory A character string denoting the path to the output MRI statistical maps
#' @param ncores An integer denoting the number of CPU cores to use when conducting permutation tests
#' @param fastSWE A boolean that determine the sandwhich estimator approach. If set to FALSE, will use standard R package geeglm. If set to TRUE will use custom-built estimator using rfast.
#' @param adjustment A character string denoting the small sample size adjustment to use when fastSwE is set to TRUE. Is NULL by default.
#' @param norm_external_data A boolean. If set to true, external data will be normed prior to analysis.
#' @param norm_internal_data A boolean. If set to true, MRI data will be normed per datapoint prior to analysis.
#' @keywords wild bootstrap sandwich estimator marginal model CIFTI scalar
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
                                   notation,
                                   zcor = NULL,
                                   corstr,
                                   family_dist,
                                   dist_type,
                                   z_thresh,
                                   nboot,
                                   p_thresh,
                                   sigtype,
                                   id_subjects='subid',
                                   output_directory='~/',
                                   ncores=1,
                                   fastSwE = TRUE,
                                   adjustment = NULL,
                                   norm_external_data = FALSE,
                                   norm_internal_data = FALSE){
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
  curr_directory = getwd()
  setwd(output_directory)
  ciftilist <- read.csv(concfile,header=FALSE,col.names="file")
  print("loading imaging data")
  start_load_time = proc.time()
  if (structtype == 'volume'){
    cifti_file <- PrepVolMetric(as.character(ciftilist$file[1]))
    cifti_dim <- dim(cifti_file)
    cifti_alldata <- array(dim = c(length(ciftilist$file),cifti_dim[1]*cifti_dim[2]*cifti_dim[3]))
    count = 1
    for (filename in ciftilist$file) {
      cifti_alldata[count,] = ConvertVolume(PrepVolMetric(as.character(ciftilist$file[count])),cifti_dim)
      count = count + 1
    }      
    remove(count)
    Nelm <- dim(cifti_alldata)[2]
    if (norm_internal_data == TRUE){
      for (count in 1:Nelm){
        if (var(cifti_alldata[,count]) != 0){ 
          cifti_alldata[,count] <- ((cifti_alldata[,count] - mean(cifti_alldata[,count]))/sd(cifti_alldata[,count]))
        }
      }
    }
    cifti_scalarmap <- as.data.frame(cifti_alldata)
  } else
  {
    cifti_alldata <- t(data.frame((lapply(as.character(ciftilist$file),PrepSurfMetric))))
    Nelm <- dim(cifti_alldata)[2]
    cifti_dim <- Nelm
    if (norm_internal_data == TRUE){
      for (count in 1:Nelm){
        if (var(cifti_alldata[,count]) != 0){ 
          cifti_alldata[,count] <- ((cifti_alldata[,count] - mean(cifti_alldata[,count]))/sd(cifti_alldata[,count]))
        }     
      }
    }
    if (fastSwE == FALSE){
      cifti_index <- 1:length(cifti_alldata)
      cifti_scalarmap <- map(cifti_index,ReframeCIFTIdata,cifti_rawmeas=cifti_alldata)
      cifti_dim=NULL
    }
  }
  print("loading non-imaging data")
  if (is.character(external_df)) {
    external_df <- read.csv(external_df,header=TRUE)
    if (fastSwE == TRUE){
      external_df <- ParseDf(external_df = external_df,notation = notation,norm_data=norm_external_data)
      nmeas <- dim(external_df)[2]
    }
  }
  if (is.character(wave)) {
    print("loading longitudinal data")
    wave <- read.csv(wave,header=TRUE)
    if (dim(wave)[2] > 1) wave <- DetermineNestedGroups(wave)
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
    resid_map <- cifti_map$residuals
    fit_map <- cifti_map$fitted.values
    Nelm <- dim(beta_map)[2]
    t_map <- ComputeFastSwE(X=external_df,nested=wave,Nelm=Nelm,resid_map=resid_map,npredictors=nmeas,beta_map=beta_map,adjustment=adjustment)
    if (structtype == 'surface'){
      for (curr_map in 1:dim(t_map)[1]){
        WriteVectorToGifti(metric_data = t_map[curr_map,],
                           surf_template_file = as.character(ciftilist[1,1]),
                           surf_command = surf_command,
                           matlab_path = matlab_path,
                           output_file = paste(output_directory,'/','t_map',curr_map,'.func.gii',sep=""))      
      }       
    } else
    {
      for (curr_map in 1:dim(t_map)[1]){
        temp_map <- RevertVolume(t_map[curr_map,],cifti_dim)
        cifti_file[] <- temp_map
        writeNIfTI(nim = cifti_file,filename = paste(output_directory,'/','t_map',curr_map,sep=""))
      }
    }
    finish_model_time = proc.time() - start_model_time
    cat("modeling complete. Time elapsed: ",finish_model_time[3],"s")
    start_normthresh_time = proc.time()
    print("Normalizing observed marginal model estimates")
    zscore_map <- t(sapply(1:nmeas,function(x) {(t_map[x,] - mean(t_map[x,is.finite(t_map[x,])]))/sd(t_map[x,is.finite(t_map[x,])])}))
    if (structtype == 'surface'){
      for (curr_map in 1:dim(zscore_map)[1]){
        WriteVectorToGifti(metric_data = zscore_map[curr_map,],
                           surf_template_file = as.character(ciftilist[1,1]),
                           surf_command = surf_command,
                           matlab_path = matlab_path,
                           output_file = paste(output_directory,'/','zscore_map',curr_map,'.func.gii',sep=""))      
      }       
    } else
    {
      for (curr_map in 1:dim(zscore_map)[1]){
        temp_map <- RevertVolume(zscore_map[curr_map,],cifti_dim)
        cifti_file[] <- temp_map
        writeNIfTI(nim = cifti_file,filename = paste(output_directory,'/','zscore_map',curr_map,sep=""))
      }      
    }
    print("thresholding observed z scores")
    thresh_map <- t(sapply(1:nmeas,function(x) abs(zscore_map[x,]) > z_thresh))
    thresh_map[is.na(thresh_map)] <- NaN
    if (structtype == 'surface'){
      for (curr_map in 1:dim(thresh_map)[1]){
        WriteVectorToGifti(metric_data = thresh_map[curr_map,],
                           surf_template_file = as.character(ciftilist[1,1]),
                           surf_command = surf_command,
                           matlab_path = matlab_path,
                           output_file = paste(output_directory,'/','thresh_map',curr_map,'.func.gii',sep=""))      
      }       
    } else
    {
      for (curr_map in 1:dim(thresh_map)[1]){
        temp_map <- RevertVolume(thresh_map[curr_map,],cifti_dim)
        cifti_file[] <- temp_map
        writeNIfTI(nim = cifti_file,filename = paste(output_directory,'/','thresh_map',curr_map,sep=""))
      }
    }    
    finish_normthresh_time = proc.time() - start_normthresh_time    
    cat("thresholding complete. Time elapsed: ", finish_normthresh_time[3],"s")
  }
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
    if (structtype == 'surface'){
      for (curr_map in 1:dim(all_cc)[2]){
        WriteVectorToGifti(metric_data = all_cc[,curr_map],
                           surf_template_file = as.character(ciftilist[1,1]),
                           surf_command = surf_command,
                           matlab_path = matlab_path,
                           output_file = paste(output_directory,'/','observed_clusters',curr_map,'.func.gii',sep=""))      
      }       
    } else
    {
      for (curr_map in 1:dim(all_cc)[2]){
        temp_map <- RevertVolume(all_cc[curr_map,],cifti_dim)
        cifti_file[] <- temp_map
        writeNIfTI(nim = cifti_file,filename = paste(output_directory,'/','observed_clusters',curr_map,sep=""))
      }
    }      
    finish_clust_time = proc.time() - start_clust_time
    cat("cluster detection complete. Time elapsed", finish_clust_time[3],"s")
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
                       adjustment=adjustment)
    stopCluster(cl)
    finish_perm_time = proc.time() - start_perm_time
    cat("permutation testing complete. Time elapsed",finish_perm_time[3],"s")
    print("calculating p values for observed data using null distribution(s)")
    for (curr_meas in 1:nmeas){
      pval_map <- map(all_cc[,curr_meas],CalculatePvalue,WB_cc=WB_cc[curr_meas,],nboot=nboot,sigtype=sigtype)
      if (curr_meas == 1){
        all_maps = pval_map
      } else
      {
        all_maps <- list(all_maps,pval_map)
      }
    }
    if (structtype == 'surface'){
      for (curr_map in 1:length(all_maps)){
        WriteVectorToGifti(metric_data = unlist(all_maps[curr_map]),
                           surf_template_file = as.character(ciftilist[1,1]),
                           surf_command = surf_command,
                           matlab_path = matlab_path,
                           output_file = paste(output_directory,'/','observed_cluster_pval_',curr_map,'.func.gii',sep=""))      
      }       
    } else
    {
      for (curr_map in 1:length(all_maps)){
        temp_map <- RevertVolume(unlist(all_maps[curr_map]),cifti_dim)
        cifti_file[] <- temp_map
        writeNIfTI(nim = cifti_file,filename = paste(output_directory,'/','observed_cluster_pval_',curr_map,sep=""))
      }
    }  
  } else
  {
    all_maps <- map(zscore_map,CalculatePvalue,WB_cc=NaN,nboot=NaN,sigtype=sigtype)
    if (structtype == 'surface'){
      for (curr_map in 1:length(all_maps)){
        WriteVectorToGifti(metric_data = unlist(all_maps[curr_map]),
                           surf_template_file = as.character(ciftilist[1,1]),
                           surf_command = surf_command,
                           matlab_path = matlab_path,
                           output_file = paste(output_directory,'/','observed_cluster_pval_',curr_map,'.func.gii',sep=""))      
      }       
    } else
    {
      for (curr_map in 1:length(all_maps)){
        temp_map <- RevertVolume(unlist(all_maps[curr_map]),cifti_dim)
        cifti_file[] <- temp_map
        writeNIfTI(nim = cifti_file,filename = paste(output_directory,'/','observed_cluster_pval_',curr_map,sep=""))
      }
    }  
  }
  setwd(curr_directory)
  all_time = proc.time() - initial_time
  cat("ConstructMarginalModel complete. Time elapsed", all_time[3],"s")
  return(all_maps)
}
