#' CalculatePvalue -- calculate the p value for cluster detection or via a lookup table
#'
#' This function wraps many of the other functions to perform cluster detection analysis using
#' the wild bootstrap.
#' @param observed_cc the observed statistical value, can be the area of a connected component or the Z transformed Wald statistic
#' @param WB_cc The distribution of areal extents derived from the wild bootstrap -- used when performing cluster detection
#' @param nboot The number of wild bootstraps performed, used to calculate the p value for cluster detection
#' @param sigtype A character string representing the type of significance test, will perform cluster detection if "cluster" is specified.
#' @keywords pvalues cluster
#' @export
#' @examples
#' p <- CalculatePvalue(observed_cc,WB_cc,nboot,sigtype)
CalculatePvalue <- function(observed_cc,WB_cc,nboot,sigtype){
  if(sigtype == 'cluster') {
    p <- sum(WB_cc > observed_cc)/nboot
    if (p == 0) {
      p <- 1/nboot
    }
  } 
  if(sigtype == 'enrichment'){
    p <- sum(WB_cc > observed_cc)/nboot
    if (p == 0) {
      p <- 1/nboot
    }
  }     
  if(sigtype == 'point'){
    p <- sum(WB_cc > abs(observed_cc))/nboot
    if (p == 0) {
      p <- 1/nboot
    }   
  }
  if(is.null(sigtype)){
      p <- 2*pnorm(abs(observed_cc)*-1)
    }
  return(p)
}