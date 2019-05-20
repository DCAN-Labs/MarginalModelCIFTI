#' CalculateMarginalValue -- a  function that will calculate a fitted value for a set of parameters
#'
#' This function is used in CalculateMarginalModel to calculate an intermediate map for visualizations
#' @param beta_measures The map of parameters for the fitted formula, including an intercept as the first parameter.
#' @param marginal_measures A numeric vector of measures to input per parameter. Ordered in the same order as parameters. Nothing should be specified for the first parameter
#' @keywords sandwich estimator marginal model fast outputs
#' @export
#' @examples
#' marginal_value <- CalculateMarginalValue(beta_measures=beta_measures,marginal_measures=as.array(c(8,9,10)))
CalculateMarginalValue <- function(beta_measures,marginal_measures){
  #first determine number of parameters
  nparams <- length(marginal_measures)
  #initialize marginal_value as the intercept value
  marginal_value = beta_measures[1]
  for (curr_param in 1:nparams) {
    marginal_value = marginal_value + beta_measures[curr_param+1]*marginal_measures[curr_param]
  }
  return(marginal_value)
}