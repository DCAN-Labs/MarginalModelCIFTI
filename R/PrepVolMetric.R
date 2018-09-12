#' PrepVolMetric -- A function for loading nifti volume metric data
#' 
#' This function is a simple wrapper for reading nifti metric volumes
#' @param vol_file A character string representing the path to a gifti file.
#' @keywords nifti volume
#' @export
#' @examples 
#' niftifile <- PrepVolMetric(nifti_file)
#' niftifiles <- lapply(nifti_filelist,PrepVolMetric)
PrepVolMetric <- function(vol_file){
  library("oro.nifti")
  niftfile <- readNIfTI(vol_file)
}