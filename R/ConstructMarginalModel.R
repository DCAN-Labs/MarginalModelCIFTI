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
                                   ncores=1){
  initial_time = proc.time()
  library("purrr")
  library("cifti")
  library("gifti")
  library("geepack")
  library("Matrix")
  library("oro.nifti")
  library("mmand")
  library("parallel")
  curr_directory = getwd()
  setwd(output_directory)
  ciftilist <- read.csv(concfile,header=FALSE,col.names="file")
  print("loading imaging data")
  start_load_time = proc.time()
  if (structtype == 'volume'){
    cifti_file <- PrepVolMetric(as.character(ciftilist$file[1]))
    cifti_dim <- dim(cifti_file)
    cifti_alldata <- array(dim = c(length(ciftilist$file),cifti_dim))
    remove(cifti_file)
    count = 1
    for (filename in ciftilist$file) {
      cifti_alldata[count,,,] = PrepVolMetric(as.character(ciftilist$file[count]))
      count = count + 1
    }
    remove(count)
    cifti_scalarmap <- as.data.frame(cifti_alldata)
  } else
  {
    cifti_alldata <- transpose(data.frame((lapply(as.character(ciftilist$file),PrepSurfMetric))))
    cifti_index <- 1:length(cifti_alldata)
    cifti_scalarmap <- map(cifti_index,ReframeCIFTIdata,cifti_rawmeas=cifti_alldata)
    cifti_dim=NULL
  }
  print("loading non-imaging data")
  if (is.character(external_df)) {
    external_df <- read.csv(external_df,header=TRUE)
  }
  if (is.character(wave)) {
    print("loading longitudinal data")
    wave <- read.csv(wave,header=TRUE)
  }
  finish_load_time = proc.time()-start_load_time
  cat("loading data complete. Time elapsed: ",finish_load_time[3],"s")  
  print("running marginal model on observed data")
  start_model_time = proc.time()
  cifti_map <- lapply(cifti_scalarmap,ComputeMM,external_df=external_df,notation=notation,family_dist=family_dist,corstr=corstr,zcor=zcor,wave=wave,id_subjects=id_subjects)
  ncomps <- length(coef(cifti_map[min(which((lengths(cifti_map) > 1) == TRUE))]$V))
  finish_model_time = proc.time() - start_model_time
  cat("modeling complete. Time elapsed: ",finish_model_time[3],"s")
  start_normthresh_time = proc.time()
  print("Normalizing observed marginal model estimates")
  zscore_map <- map(cifti_map,ComputeZscores,ncomps)
  save(zscore_map,file = "zscore_observed.Rdata")
  resid_map <- map(cifti_map,ComputeResiduals,ncomps)
  fit_map <- map(cifti_map,ComputeFits,ncomps)
  print("thresholding observed z scores")
  thresh_map <- map(zscore_map,ThreshMap,zthresh=z_thresh)
  save(thresh_map,file = "zscore_thresh_observed.Rdata")
  finish_normthresh_time = proc.time() - start_normthresh_time
  cat("thresholding complete. Time elapsed: ", finish_normthresh_time[3],"s")
  print("performing cluster detection")
  start_clust_time = proc.time()
  if (sigtype == 'cluster'){
    varlist <- all.vars(notation)
    nmeas <- length(varlist)
    if (sigtype == 'surface'){
      all_cc = matrix(data=NA,nrow=length(thresh_map),ncol=nmeas)     
    } else 
    {
      all_cc = array(data=NA,dim=c(nmeas,cifti_dim))    
      }
    for (curr_meas in 1:nmeas){
      thresh_array = unlist(thresh_map)
      mask_vector = 1:nmeas == curr_meas
      thresh_array <- thresh_array[mask_vector]
      thresh_array <- as.numeric(thresh_array)
      thresh_array[is.na(thresh_array)] <-  0
      if (structtype == 'volume'){
        thresh_vol <- RevertVolume(thresh_array,cifti_dim)
        observed_cc = GetVolAreas(thresh_vol)
        all_cc[curr_meas,,,] = observed_cc
      } else
      {
        observed_cc <- GetSurfAreas(thresh_array,structfile,matlab_path,surf_command)
        all_cc[,curr_meas] = observed_cc
      }
      save(all_cc,file = "connected_components_observed.Rdata")
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
                                        nmeas=ncomps,
                                        seeds=seeds)
    stopCluster(cl)
    finish_perm_time = proc.time() - start_perm_time
    cat("permutation testing complete. Time elapsed",finish_perm_time[3],"s")
    print("calculating p values for observed data using null distribution(s)")
    if (structtype=='surface') {
      for (curr_meas in 1:nmeas){
        pval_map <- map(all_cc[,curr_meas],CalculatePvalue,WB_cc=WB_cc[curr_meas,],nboot=nboot,sigtype=sigtype)
        if (curr_meas == 1){
          all_maps = pval_map
        } else
        {
          all_maps <- list(all_maps,pval_map)
        }
      }
    } else 
    {
      all_maps = array(data=0,dim = c(nmeas,cifti_dim))
      for (curr_meas in 1:nmeas){
        pval_map <- apply(all_cc[curr_meas,,,],c(1,2,3),CalculatePvalue,WB_cc=WB_cc[curr_meas,],nboot=nboot,sigtype=sigtype)
        all_maps[curr_meas,,,] <- pval_map
      }
    }
    save(all_maps,file = "pval_maps.Rdata")    
  } else
  {
    all_maps <- map(zscore_map,CalculatePvalue,WB_cc=NaN,nboot=NaN,sigtype=sigtype)
    save(all_maps,file = "pval_maps.Rdata")
  }
  setwd(curr_directory)
  all_time = proc.time() - initial_time
  cat("ConstructMarginalModel complete. Time elapsed", all_time[3],"s")
  return(all_maps)
}
