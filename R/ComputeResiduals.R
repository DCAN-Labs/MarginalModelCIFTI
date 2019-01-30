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
ComputeResiduals <- function(geeglm_obj,nmaps) {
  if (is.numeric(unlist(geeglm_obj))) {
    residual_model = numeric(nmaps)
  } else
  {
    if (is.object(geeglm_obj)) {
      gee_obj = geeglm_obj
    } else
    {
      gee_obj <- geeglm_obj$V
    }
    residual_model <- resid(gee_obj)    
  }
  return(residual_model)
}