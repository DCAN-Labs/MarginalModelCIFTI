#' ComputeMM -- compute the marginal model at a single vertex/voxel
#' 
#' This function will compute the marginal model on a single vertex or voxel. 
#' The function can be called in parallel to operate on a map. Such usage can be found in the main function
#  or in ComputeMM_WB.R. See below for examples.
#' @param cifti_meas An array of N numeric measures, where N is the number of subjects
#' @param external_df A data frame with NxM measures. Should include all features that the user wants to model.
#' @param notation A character string representing the modelling formula. Uses the wilkinson notation.
#' @param family_dist A character string representing the type of distribution to use for modelling. 
#' @keywords marginal model sandwich estimator
#' @export
#' @examples 
#' geeglm_obj <- ComputeMM(cifti_meas,external_df,notation,family_dist)
ComputeMM <- function(cifti_meas,external_df,notation,family_dist,corstr,zcor,waves) {
  library("geepack")
  data_to_fit <-  cbind(cifti_meas,external_df)
  geeglm_obj <- geeglm(notation, id=subid, family=family_dist,
                       corstr=corstr, data=data_to_fit, waves=wave,zcor=zcor)
  return(geeglm_obj)
}