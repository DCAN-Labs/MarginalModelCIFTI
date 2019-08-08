#' FilterDFNA -- a fast function to check and filter cases based on bad behavioral data. Used when fastSwE is set to true
#'
#' This function is used to extract data from a data frame and convert to a design matrix. The matrix is then input into SwE for fast computations.
#' @param external_df The data frame used to extract data
#' @param notation The formula specified by the user
#' @keywords sandwich estimator marginal model fast matrix
#' @export
#' @examples
#' X <- FilterDFNA(external_df,notation)
FilterDFNA <- function(external_df,notation){
  ncases <- dim(external_df)[1]
  predictors <- attr(terms(notation),"term.labels")
  X <- cbind(rep(1,ncases),external_df[[predictors[1]]])
  if (length(predictors) > 1){
    for (curr_pred in 2:length(predictors)){
      X <- cbind(X,external_df[[predictors[curr_pred]]])
    }
  }
  X_bin = rowsums(is.na(X)) == 0
  return(X_bin)
}