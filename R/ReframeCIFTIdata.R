#' ReframeCIFTIdata -- converts CIFTI scalar data from multiple subjects into a greyordinate x subject data frame that is easier to manipulate
#'
#' This function wraps many of the other functions to perform cluster detection analysis using
#' the wild bootstrap.
#' @param index the greyordinate to extract and convert into a data frame
#' @param cifti_rawmeas The object created by PrepCIFTI, when multiple subjects are loaded (e.g.  "transpose(data.frame((lapply(as.character(ciftilist$file),PrepCIFTI))))")
#' @keywords CIFTI scalar import
#' @export
#' @examples
#' cifti_meas <- ReframeCIFTIData(index,cifti_rawmeas)
ReframeCIFTIdata <- function(index,cifti_rawmeas) {
  cifti_meas <- data.frame(y = as.numeric(unlist(cifti_rawmeas[index])))
  return(cifti_meas)
}