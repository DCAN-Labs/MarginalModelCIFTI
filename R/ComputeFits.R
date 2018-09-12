#' ComputeFits -- get estimated fits from a geeglm_object
#' 
#' This function will compute fits from a marginal model, which can be used to generate wild bootstraps.
#' The function can be called in parallel to operate on a map. Such usage can be found in the main function
#  or in ComputeMM_WB.R. See below for examples.
#' @param geeglm_obj An GLM object that is produced from computing the marginal model
#' @keywords estimates bootstrapping
#' @export
#' @examples 
#' fits <- ComputeFits(geeglm_obj)
ComputeFits <- function(geeglm_obj) {
  model_fits <- fitted(geeglm_obj)
}