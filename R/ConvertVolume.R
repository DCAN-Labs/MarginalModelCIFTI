#' ConvertVolume -- convert a 3D MRI volume into a list of values
#' 
#' This function will return a 3D array for items within a list. The dimensions of the image must be
#' specified and match the same number of items as the list. Useful when converting an array to a list
#' See below for examples.
#' @param vol_array A 3D array that contains values per voxel
#' @param vol_dim A numeric of three values to use to identify the size of each dimension (X,Y,Z)
#' @keywords volume conversion
#' @export
#' @examples 
#' thresh_vol <- ConvertVolume(thresh_map,vol_dim)
ConvertVolume <- function(vol_array,vol_dim){
  vol_list = array(dim = c(1,vol_dim[1]*vol_dim[2]*vol_dim[3]))
  i = 1
  j = 1
  k = 1
  vol_dim = dim(vol_array)
  for (currframe in 1:dim(vol_list)[2]) { 
    vol_list[1,currframe] <- vol_array[i,j,k] 
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
  return(vol_list)
}