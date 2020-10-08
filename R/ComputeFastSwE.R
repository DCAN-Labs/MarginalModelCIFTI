#' ComputeFastSwE -- a fast function to calculate errors via the sandwich estimator. Modified from Tom Nichols, Feb 1, 2019
#'
#' This function is used to quickly calculate error from a marginal model via a sandwich estimator. The function can handle a single block (i.e. repeated measure)
#' @param X The predictor data recast as a matrix
#' @param nested A vector denoting the blocks (i.e. groups) which are used to calculate variance separately if set to NULL will compute as if a single group
#' @param Nelm An integer denoting the number of voxels/vertices
#' @param resid_map The residuals from the estimated data fits, usually by lm.fit or lmfit
#' @param npredictors The number of predictors in X
#' @param beta_map The beta map from the estimated data fits, usually by lm.fit or lmfit
#' @param adjustment The residuals will be adjusted according to the small sample size adjustment. Acceptable values are "HC2", "HC3", and NULL.
#' @keywords sandwich estimator marginal model fast matrix
#' @export
#' @examples
#' T_map <- ComputeFastSwE(X=external_df,nested=nested,Nelm=Nelm,resid_map=resid_map,npredictors=npredictors,beta_map=beta_map,adjustment=adjustment)
ComputeFastSwE <- function(X,nested,Nelm,resid_map,npredictors,beta_map,adjustment){
  if (is.null(adjustment)) {
    hat_adjust <- NULL
  } else
  {
    hat_adjust <- diag(X%*%((t(X)%*%X)^-1)%*%t(X))
  }
  # Computation of SwE standard errors
  start_time <- Sys.time()
  SE.swe = matrix(0,nrow=npredictors,ncol=Nelm)
  Bread  = solve(t(X)%*%X)
  BreadX = Bread%*%t(X)
  S      = array(0,dim=c(npredictors,npredictors,Nelm))
  S0     = array(0,dim=c(1,npredictors,Nelm))
  if (is.null(nested)){
    if (is.null(adjustment)){
    } else if (adjustment == "HC2"){
      resid_map <- resid_map/sqrt(1 - hat_adjust)
    } else if (adjustment == "HC3"){
      resid_map <- resid_map/(1 - hat_adjust)
    }
    Ns = dim(X)[1]
    e = array(resid_map,c(1,Ns,Nelm))
    if (npredictors==1){
      S0[] = apply(e,3,function(x)x%*%(BreadX))
    } else
    {
      S0[] = apply(e,3,function(x)x%*%t(BreadX))
    }
    # Full `Bread*Meat*Bread' contribution
    S = S + array(apply(S0,3,function(x)t(x)%*%x),dim=c(npredictors,npredictors,Nelm))
  } else
  {
    Nnest <- max(nested)
    for (s in 1:Nnest) {
      I=(s==nested)
      Ns=sum(I)
      if (is.null(adjustment)){
        e = array(resid_map[I,],c(1,Ns,Nelm))
      } else if (adjustment == "HC2"){
        e = array(resid_map[I,]/sqrt(1 - hat_adjust),c(1,Ns,Nelm))
      } else if (adjustment == "HC3"){
        e = array(resid_map[I,]/(1 - hat_adjust),c(1,Ns,Nelm))
      }
      # half of Meat times t(BreadX)
      if (npredictors==1){
        S0[] = apply(e,3,function(x)x%*%(BreadX[,I]))
      }else 
      {
        S0[] = apply(e,3,function(x)x%*%t(BreadX[,I]))
      }
      # Full `Bread*Meat*Bread' contribution for nest s
      S = S + array(apply(S0,3,function(x)t(x)%*%x),dim=c(npredictors,npredictors,Nelm))
    }
  }
  SE.swe[]=sqrt(apply(S,3,diag))
  T.swe = beta_map/SE.swe
  return(T.swe)
}