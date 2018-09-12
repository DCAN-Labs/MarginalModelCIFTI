#' PrepCIFTI -- A function for laoding cifti data directly. Not useful for cluster detection
#' 
#' This function is a simple wrapper for reading cifti data
#' @param cifti_file A character string representing the path to a cifti file.
#' @keywords cifti scalar
#' @export
#' @examples 
#' cifti_data <- PrepCIFTI(cifti_file)
#' cifti_dataset <- lapply(cifti_file,PrepCIFTI)
PrepCIFTI <- function(cifti_file) {
  library("cifti")
  cifti_datafile <-  read_cifti(cifti_file)
  return(cifti_datafile$data)
}