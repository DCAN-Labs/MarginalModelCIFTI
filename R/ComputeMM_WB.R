#' ComputeMM_WB -- calculate the maximum cluster extent from a wild bootstrap
#'
#' This function wraps many of the other functions to perform cluster detection analysis using
#' the wild bootstrap.
#' @param cifti_map surface array or nifti containing the scalar values -- just used to determine the number of participants
#' @param zscore_map The z-scored map from the geeglm analysis
#' @param resid_map The residual map from the geeglm analysis
#' @param fit_map The fitted values map from the geeglm analysis
#' @param type A character string denoting the distribution to use for the wild bootstrap
#' @param external_df A data frame comprising non-brain measures to model.
#' @param notation A character string representing the model formula. Uses Wilkinson notation.
#' @param family_dist A character string denoting the assumed underlying distribution of all input data
#' @param structtype A character string denoting whether the map is volumetric ('volume') or surface-based
#' @param matlab_path A character string denoting the path to the matlab compiler. Please use v91.
#' @param surf_command A character string denoting the path to the SurfConnectivity command.
#' @param correctiontype A character string denoting cluster ('cluster') or point ('point') comparisons
#' @keywords wild bootstrap
#' @export
#' @examples
#' max_cc <- ComputeMM_WB(cifti_map,zscore_map,resid_map,fit_map,type,external_df,
#' notation,family_dist,structtype,thresh,structfile,matlab_path,surf_command,correctiontype)
ComputeMM_WB <- function(cifti_map,
                         zscore_map,
                         resid_map,
                         fit_map,
                         type,
                         external_df,
                         notation,
                         family_dist,
                         structtype,
                         thresh,
                         structfile,
                         matlab_path,
                         surf_command,
                         corstr,
                         wave,
                         zcor,
                         correctiontype,
                         subid) {
  library("purrr")
  library("cifti")
  library("gifti")
  library("geepack")
  library("Matrix")
  library("oro.nifti")
  library("mmand")
  ncomps <- nrow(cifti_map)
  if (type=='radenbacher'){
    bootstrap_mult <- sample(c(-1,1),ncomps,replace=TRUE)
  } else if (type=='mammen'){
    bootstrap_mult <- sample(c((1+sqrt(5))/2,(1-sqrt(5))/2),ncomps,prob=c((sqrt(5)-1)/(2*sqrt(5)),(sqrt(5)+1)/(2*sqrt(5))),replace=TRUE)
  } else if (type=='webb4') {
    bootstrap_mult <- sample(c(sqrt(3/2)*-1,sqrt(1/2)*-1,sqrt(1/2),sqrt(3/2)),ncomps,replace=TRUE)
  } else if (type=='webb6') {
    bootstrap_mult <- sample(c(sqrt(3/2)*-1,-1,sqrt(1/2)*-1,sqrt(1/2),1,sqrt(3/2)),ncomps,replace=TRUE)
  } else {
    return('error: invalid argument for distribution type, please specify one of: radenbacher,mammen,webb4,webb6')
  }
  df <- data.frame(
    fit_val = fit_map,
    resid_val = resid_map,
    sample_val = bootstrap_mult
  )
  cifti_bootmap <- pmap(df,ApplyWB_toData)
  MM_bootmap <- map(cifti_bootmap,ComputeMM,external_df=external_df,notation=notation,family_dist=family_dist,corstr=corstr,zcor=zcor,waves=wave,subid=subid)
  zscore_bootmap <- map(MM_bootmap,ComputeZscores)
  thresh_bootmap <- zscore_bootmap > thresh
  if (correctiontype=='point') {
    return(max(zscore_bootmap))
  }
  else if (correctiontype=='cluster'){
    if (structtype=='volume'){
      cc <- GetVolAreas(thresh_bootmap)
      }
    else{
      cc <- GetSurfAreas(thresh_bootmap,structfile,matlab_path,surf_command)
      }
    return(max(cc))
    }
}