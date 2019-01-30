#' RevertVolume -- convert a list of values into a 3D MRI volume. 
#' 
#' This function will return a 3D array for items within a list. The dimensions of the image must be
#' specified and match the same number of items as the list. Useful when converting an array to a list
#' See below for examples.
#' @param vol_list A list item that contains values per voxel
#' @param vol_dim A numeric of three values to use to identify the size of each dimension (X,Y,Z)
#' @keywords volume conversion
#' @export
#' @examples 
#' thresh_vol <- RevertVolume(thresh_map,vol_dim)
RevertVolume <- function(vol_list,vol_dim){
  vol_array = array(data = 0,dim = vol_dim)
  i = 1
  j = 1
  k = 1
  vol_dim = dim(vol_array)
  for (currframe in 1:length(vol_list)) { 
    vol_array[i,j,k] = as.numeric(vol_list[currframe])
    i = i + 1
    if (i > vol_dim[1]) {
      i = 1
      j = j + 1
      if (j > vol_dim[2]) {
        j = 1
        k = k + 1
      }
    }
  }
  return(vol_array)
}