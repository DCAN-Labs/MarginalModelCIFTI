#' PrepNiiConnMetric -- A function for loading nifti connectivity metric data
#' 
#' This function is a simple wrapper for reading connectivity matrices from niftis. 
#' The function will check to see whether the qform is needed before loading in the data.
#' Generally, the qform shouldn't be needed, if it is used, there may be a warning when using the dataset.
#' @param niiconn_file A character string representing the path to a nifti connectivity matrix file.
#' @param flip A boolean representing TRUE or FALSE. If set to TRUE, will rotate the matrix 180 degrees so that it is ordered properly with a parcel file
#' @keywords nifti volume conn
#' @export
#' @examples 
#' niftifile <- PrepNiiConnMetric(nifti_file)
#' niftifiles <- lapply(nifti_filelist,PrepNiiConnMetric)
PrepNiiConnMetric <- function(niiconn_file,flip = TRUE){
  library("oro.nifti")
  niiconnmat <- readNIfTI(niiconn_file,reorient = FALSE)
  if (niiconnmat[1,1] != 1){
    print("Warning: default load did not orient matrix properly -- using qform to flip matrix")
    niiconnmat <- readNIfTI(niiconn_file,reorient = TRUE)
  }
  if (flip == TRUE){
    niiconnmat_t <- t(apply(niiconnmat,2,rev))
    niiconnmat_t2 <- t(apply(niiconnmat_t,2,rev))
    niiconnmat[,] = niiconnmat_t2
  }
  return(niiconnmat)
}