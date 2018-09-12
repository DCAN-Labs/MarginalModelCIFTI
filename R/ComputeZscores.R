#' ComputeZscores -- convert estimated Wald statistics to z scores
#' 
#' This function will compute z scores from Wald statistics, which can be used to infer statistical significance.
#' The function can be called in parallel to operate on a map. Such usage can be found in the main function
#  or in ComputeMM_WB.R. See below for examples.
#' @param geeglm_obj An GLM object that is produced from computing the marginal model
#' @keywords normalization
#' @export
#' @examples 
#' zscores <- ComputeZscores(geeglm_obj)
ComputeZscores <- function(geeglm_obj) {
  zscores_model <- qnorm(0.99999999999999-pchisq((coef(geeglm_obj)/sqrt(diag(geeglm_obj$geese$vbeta)))^2,1, lower.tail=F))
}