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
  # Compute the pseudoinverse of the design matrix.
  # In the future we may compute this with pracma's pinv(),
  # or it may be passed into the function directly.
  Xpinv = solve(X%*%t(X))%*%t(X);

  # Apply small sample size adjustment to residuals.
  if (is.null(adjustment)){
  } else if (adjustment == "HC2"){
    resid_map <- resid_map/sqrt(1 - diag(Xpinv))
  } else if (adjustment == "HC3"){
    resid_map <- resid_map/(1 - diag(Xpinv))
  }

  # Allocate memory into which we will store the computed standard errors.
  SE.swe = matrix(0,nrow=npredictors,ncol=Nelm)

  # Iterate over voxels/vertices using the apply() function.
  if (is.null(nested)){
    # Assume error terms are independently but *not* identically distributed.
    # That is, the errors are heteroskedastic without any blocking/nesting.

    # This is optimized code.  Recall Xpinv is the Moore-Penrose
    # pseudoinverse of the design matrix X.  What we actually want to compute
    # separately for each voxel/vertice's residuals (as a column vector e) is:
    # SE = sqrt( diag( Xpinv %*% diag(e^2) %*% t(Xpinv) ) )
    # The first optimization is to replace raising to a power with:
    # SE = sqrt( diag( Xpinv %*% diag(e*e) %*% t(Xpinv) ) )
    # Since the middle is a diagonal matrix, rather than actually allocate
    # the entire matrix in memory, we can broadcast the diagonal using the sweep
    # function and then do simple element-wise multiplication:
    # SE = sqrt( diag( sweep(Xpinv,MARGIN=2,e*e,`*`) %*% t(Xpinv) ) )
    # Next, we realize that we only care about the diagonal of the product of
    # the Xpinv %*% t(Xpinv) (the element-wise multiplication by sweep() is
    # factored out).  We can compute this by taking the sum of each row of
    # Xpinv * t(Xpinv).  Putting it all together:
    SE.swe = apply(resid_map,2,function(e)sqrt(rowSums(sweep(Xpinv,MARGIN=2,e*e,`*`)*Xpinv))
  } else
  {
    # Iterate over each nesting and compute partitions of diag(e^2).
    block_resid_map = array(0, dim=c(dim(X)[1], Nelm))
    Nnest <- max(nested)
    for (s in 1:Nnest) {
      I=(s==nested)
      Ns=sum(I)
      block_resid_map[I,] = apply(resid_map[I,],2,function(e)array(mean(e.^2), dim=c(Ns)))
    }

    # Same computation as non-nested case using blocked residual map.
    SE.swe = apply(block_resid_map,2,function(e)sqrt(rowSums(sweep(Xpinv,MARGIN=2,e*e,`*`)*Xpinv))
  }
  # T-values = estimated parameters divided by their standard errors
  T.swe = beta_map/SE.swe
  return(T.swe)
}
