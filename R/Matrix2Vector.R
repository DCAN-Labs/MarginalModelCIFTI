#' Matrix2Vector -- a function that will map a matrix to a vector and vice versa
#'
#' This function is used when loading pconns in ConstructMarginalModel. 
#' @param pconn_data The pconn matrix represented as a 2D matrix
#' @param pconn_vector The pconn vector to map to the matrix. Only used when "to_matrix" is specified
#' @param direction A string denoting the direction of the transform. Can be either "to_matrix" or "to_vector". Default is "to_vector"
#' @keywords pconn cifti marginal model
#' @export
#' @examples
#' data <- Matrix2Vector(pconn_data=cifti_data$data,pconn_vector=NULL,direction="to_vector")
#' data <- Matrix2Vector(pconn_data=cifti_data$data,pconn_vector=pvalue_vector,direction="to_matrix")
Matrix2Vector <- function(pconn_data=NULL,pconn_vector=NULL,direction="to_vector") {
  if(direction=="to_matrix"){
    data_out <- pconn_data
    data_out[upper.tri(data_out)] <- pconn_vector
    data_out[lower.tri(data_out)] <- t(data_out[upper.tri(data_out)])
  } else
  {
    data_out <- pconn_data[upper.tri(pconn_data)]
  }
  return(data_out)
}