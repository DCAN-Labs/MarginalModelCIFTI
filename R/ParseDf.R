#' ParseDf -- a fast function to extract data from a data frame and convert to a design matrix. Used when fastSwE is set to true
#'
#' This function is used to extract data from a data frame and convert to a design matrix. The matrix is then input into SwE for fast computations.
#' @param external_df The data frame used to extract data
#' @param notation The formula specified by the user
#' @param norm_data A TRUE or FALSE flag. If set to TRUE will norm the data prior to processing.
#' @keywords sandwich estimator marginal model fast matrix
#' @export
#' @examples
#' X <- ParseDf(external_df,notation)
ParseDf <- function(external_df,notation,norm_data){
  ncases <- dim(external_df)[1]
  predictors <- attr(terms(notation),"term.labels")
  if (norm_data==TRUE){
    X <- cbind(rep(1,ncases),(external_df[[predictors[1]]] - mean(external_df[[predictors[1]]]))/var(external_df[[predictors[1]]]))
    if (length(predictors) > 1){
      for (curr_pred in 2:length(predictors)){
        X <- cbind(X,(external_df[[predictors[curr_pred]]] - mean(external_df[[predictors[curr_pred]]]))/var(external_df[[predictors[curr_pred]]]))
      }
    }    
  } else 
  {
    X <- cbind(rep(1,ncases),external_df[[predictors[1]]])
    if (length(predictors) > 1){
      for (curr_pred in 2:length(predictors)){
        X <- cbind(X,external_df[[predictors[curr_pred]]])
      }
    }
  }
  return(X)
}