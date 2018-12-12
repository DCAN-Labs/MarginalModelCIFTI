#' ThreshMap -- threshold a list of values and return a binary list, can handle multiple values within each item of a list
#' 
#' This function will return a logical value for items within a list, used with purrr's map function, can be used to
#' parallelize thresholding on a list of many values (e.g. from CIFTIs). See below for examples.
#' @param zscore_val A list item that contains values
#' @param z_thresh A value to use to threshold the image
#' @keywords thresholding
#' @export
#' @examples 
#' thresh <- ThreshMap(zscore_val,z_thresh)
ThreshMap = function(zscore_val,zthresh){
  thresh_val <- zscore_val > zthresh
  return(thresh_val)
}
