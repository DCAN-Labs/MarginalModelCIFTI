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
  ciftilist <- read.csv(concfile,header=FALSE,col.names="file")
  if (structtype == 'volume'){
    cifti_alldata <- data.frame((lapply(as.character(ciftilist$file),PrepVolMetric)))
    cifti_scalarmap <- cifti_alldata
  }
  else{
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
  cifti_map <- map(cifti_scalarmap,ComputeMM,external_df=external_df,notation=notation,family_dist=family_dist,corstr=corstr,zcor=zcor,wave=wave,id_subjects=id_subjects)
  zscore_map <- map(cifti_map,ComputeZscores)
  resid_map <- map(cifti_map,ComputeResiduals)
  fit_map <- map(cifti_map,ComputeFits)
  thresh_map <- zscore_map > z_thresh
  if (sigtype == 'cluster'){
    if (structtype == 'volume'){
      observed_cc <- GetVolAreas(thresh_map)
    }
    else{
      observed_cc <- GetSurfAreas(thresh_map,structfile,matlab_path,surf_command)
    }
    WB_cc <- replicate(nboot,ComputeMM_WB(cifti_map = cifti_map,
                                          zscore_map,
                                          resid_map,
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
                                          correction_type = sigtype,
                                          id_subjects))
    pval_map <- map(observed_cc,CalculatePvalue,WB_cc=WB_cc,nboot=nboot,sigtype=sigtype)
  }
  else{
    pval_map <- map(zscore_map,CalculatePvalue,WB_cc=NaN,nboot=NaN,sigtype=sigtype)
  }
  return(pval_map)
}
