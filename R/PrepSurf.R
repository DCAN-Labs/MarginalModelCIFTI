#' PrepSurf -- A function for loading gifti surfaces
#' 
#' This function is a simple wrapper for read_gifti.
#' @param surf_file A character string representing the path to a gifti file.
#' @keywords gifti surface
#' @export
#' @examples 
#' surface_data <- PrepSurf(gifti_file)
#' surfaces_data <- lapply(gifti_filelist,PrepSurf)

PrepSurf <- function(surf_file){
  library("gifti")
  surface_data <- read_gifti(surf_file)
}