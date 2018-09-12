#' PrepSurfMetric -- A function for loading gifti metric data
#' 
#' This function is a simple wrapper for reading gifti metric data
#' @param funcgii_file A character string representing the path to a func.gii file.
#' @keywords gifti metric
#' @export
#' @examples 
#' gifti_data <- PrepSurfMetric(giifunc_file)
#' gifti_dataset <- lapply(giifunc_filelist,PrepSurfMetric)
PrepSurfMetric <- function(funcgii_file) {
  library("gifti")
  gifti_datafile <- read_gifti(funcgii_file)
  return(gifti_datafile$data)
}