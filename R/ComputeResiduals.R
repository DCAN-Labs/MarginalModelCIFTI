#' ComputeResiduals -- get residuals from a geeglm_object
#' 
#' This function will compute residuals from a marginal model, which can be used to generate wild bootstraps.
#' The function can be called in parallel to operate on a map. Such usage can be found in the main function
#  or in ComputeMM_WB.R. See below for examples.
#' @param geeglm_obj An GLM object that is produced from computing the marginal model
#' @keywords residuals bootstrapping
#' @export
#' @examples 
#' resids <- ComputeResiduals(geeglm_obj)
ComputeResiduals <- function(geeglm_obj) {
  residuals_model <- resid(geeglm_obj)
}