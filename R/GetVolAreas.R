#' GetVolAreas -- Get computed size of connected volumetric components
#' 
#' This function will compute the size of connected components in statistical volumes
#' @param bin_image A nifti array containing only binary values (e.g. 0;1)
#' @keywords volume cluster
#' @export
#' @examples 
#' component_sizes <- GetVolAreas(bin_image)
GetVolAreas <- function(bin_image) {
  library("mmand")
  components <- components(bin_image,c(1,1,1))
  W <- as.vector(na.exclude(unique(components)))
  return(sapply(1:length(W),function(x) sum(components==x,na.rm=TRUE)))
}