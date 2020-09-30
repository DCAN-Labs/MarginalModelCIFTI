#' DetermineNestedGroups -- takes mulitiple within/between subject measures to create unique blocks for calculating SwE
#' 
#' This function will identify unique groups based on the input wave dataset. 
#' These groups will have their covariance structure calculated separately for the SwE
#' See below for examples.
#' @param wave A data frame typically read from a csv file.
#' @keywords blocks wave
#' @export
#' @examples 
#' final_wave <- DetermineNestedGroups(wave)
DetermineNestedGroups <- function(wave) {
  nvars <- dim(wave)[2]
  if (is.null(nvars)){
    final_wave <- as.numeric(factor(wave))
  } else if (nvars == 2) {
    final_wave <- as.numeric(factor(paste0(wave[,1],wave[,2])))
  } else if (nvars > 2){
    pasted_wave <- paste0(wave[,1],wave[,2])
    for (curr_col in 3:nvars){
      pasted_wave <- paste0(pasted_wave,wave[,curr_col])
    }
    final_wave <- as.numeric(factor(pasted_wave))
  } else
  {
    final_wave <- as.numeric(factor(wave))
  }
  return(final_wave)
}