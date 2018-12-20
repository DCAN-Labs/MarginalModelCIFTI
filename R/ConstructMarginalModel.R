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
                                   id_subjects='subid'){
  library("purrr")
  library("cifti")
  library("gifti")
  library("geepack")
  library("Matrix")
  library("oro.nifti")
  library("mmand")
  ciftilist <- read.csv(concfile,header=FALSE,col.names="file")
  if (structtype == 'volume'){
    cifti_alldata <- data.frame((lapply(as.character(ciftilist$file),PrepVolMetric)))
    cifti_scalarmap <- cifti_alldata
  } else
  {
    cifti_alldata <- transpose(data.frame((lapply(as.character(ciftilist$file),PrepSurfMetric))))
    cifti_index <- 1:length(cifti_alldata)
    cifti_scalarmap <- map(cifti_index,ReframeCIFTIdata,cifti_rawmeas=cifti_alldata)
  }
  if (is.character(external_df)) {
    external_df <- read.csv(external_df,header=TRUE)
  }
  if (is.character(wave)) {
    wave <- read.csv(wave,header=TRUE)
  }
#  cifti_map = list()
#  for (curr_cifti_meas in 1:length(cifti_scalarmap)){
#    cifti_meas <- cifti_scalarmap[curr_cifti_meas]
#    if (sum(is.na(unlist(cifti_meas))) > 0){
#      cifti_measb = data.frame(y=numeric(length(cifti_meas)))
#      data_to_fit <-  cbind(cifti_measb,external_df)
#    } else
#    {
#      data_to_fit <-  cbind(cifti_meas,external_df)
#    }
#    message(data_to_fit)
#    cifti_map[curr_cifti_meas] <- geeglm(notation, data=data_to_fit, id=data_to_fit[[id_subjects]], family=family_dist,
#                         corstr=corstr, waves=wave,zcor=zcor)
#  }
#  data_to_fit <- external_df 
  cifti_map <- lapply(cifti_scalarmap,ComputeMM,external_df=external_df,notation=notation,family_dist=family_dist,corstr=corstr,zcor=zcor,wave=wave,id_subjects=id_subjects)
  zscore_map <- map(cifti_map,ComputeZscores)
  resid_map <- map(cifti_map,ComputeResiduals)
  fit_map <- map(cifti_map,ComputeFits)
  thresh_map <- map(zscore_map,ThreshMap,zthresh=z_thresh)
  if (sigtype == 'cluster'){
    varlist <- all.vars(notation)
    nmeas <- length(varlist)
    all_cc = matrix(data=NA,nrow=length(thresh_map),ncol=nmeas)
    for (curr_meas in 1:nmeas){
      thresh_array = unlist(thresh_map)
      mask_vector = 1:nmeas == curr_meas
      thresh_array <- thresh_array[mask_vector]
      thresh_array <- as.numeric(thresh_array)
      thresh_array[is.na(thresh_array)] <-  0
      if (structtype == 'volume'){
        observed_cc <- GetVolAreas(thresh_array)
      } else
      {
        observed_cc <- GetSurfAreas(thresh_array,structfile,matlab_path,surf_command)
      }
      all_cc[,curr_meas] = observed_cc
    }
    WB_cc <- replicate(nboot,ComputeMM_WB(resid_map,
                                        fit_map,
                                        dist_type,
                                        external_df,
                                        notation,
                                        family_dist,
                                        structtype,
                                        thresh = z_thresh,
                                        structfile,
                                        matlab_path,
                                        surf_command,
                                        corstr,
                                        wave,
                                        zcor,
                                        correctiontype = sigtype,
                                        id_subjects))
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
    all_maps <- map(zscore_map,CalculatePvalue,WB_cc=NaN,nboot=NaN,sigtype=sigtype)
  }
  return(all_maps)
}
