#' SubSelection -- A function that returns indices for subselecting datasets
#'
#' This function wraps all of the other functions in the MarginalModelCifti package. The output is a significance map, which is either uncorrected or cluster detection was performed.
#' @param sublist A data frame comprising non-brain measures to model. Can be specified as a string to a csv file with appropriate headers.
#' @param subsetlist A character string denoting a single column text file that lists the dscalars in the same order as the external_df and wave frames.
#' @#' @keywords wild bootstrap sandwich estimator marginal model CIFTI scalar
#' @export
#' @examples
#' list_bin <- SubSelection(sublist,subsetlist)

SubSelection <- function(sublist,subsetlist){
  masklist <- str_detect(sublist,subsetlist)
  return(masklist)
}
