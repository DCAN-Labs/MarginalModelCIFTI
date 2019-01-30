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
  component_mat <- components(bin_image,c(1,1,1))
  components_list <- unique(as.list(component_mat))
  components_vol <- array(data = 0, dim = dim(bin_image))
  for (curr_comp in 1:length(components_list)) {
    components_vol[component_mat == curr_comp] <- as.numeric(components_list[curr_comp])
    }
  return(components_vol)
}