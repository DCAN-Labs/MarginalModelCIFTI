#' RecursiveNested -- a recursive function that will calculate the sandwich estimator errors controlling for a blocked variable (e.g. within family/school/age)
#'
#' This function is used in the fastSWEcalculation.R function to quickly calculate SwE for a nested variable. 
#' @param iter The current iteration, when called it should be set to 1. 
#' @param nests The range of blocks within the nested design should be ordered in same order as the set of whole numbers (e.g. 1,2,3,4...)
#' @param nested An integer vector denoting the blocks ordered by the input data
#' @param datavec An integer denoting the number of vertices/voxels to be tested
#' @param Sinit The initial covariance structure used by the sandwich estimator, default is a 1xPxdatavec numeric array set to zeros
#' @param residarray The residuals from the initial fit of the model to the voxel/vertex-wide dataset, should be a numeric array of num_casesXnum_voxel/vertices
#' @param Breadvec The "bread" of the sandwich estimator, calculated from the predictor variables P, should be a numeric array of Pxnum_cases
#' @param npred The number of predictors, P. Should be an integer.
#' @param S The covariance structure estimated from the sandwich estimator, updated upon each recursive estimation.
#' @keywords sandwich estimator marginal model fast recursive
#' @export
#' @examples
#' SwE <- RecursiveNested(iter=1,nests=1:Nschool,nested=Ischool,datavec=Nelm,Sinit=S0,residarray=res,Breadvec=BreadX,npred=P,adjustment=adjustment,hat_adjust=hat_adjust,S)
RecursiveNested <- function(iter,nests,nested,datavec,Sinit,residarray,Breadvec,npred,adjustment,X,Sb) {
  if(iter > length(nests)) return(Sb)
  else {
    I=(nests[iter]==nested)
    Ns=sum(I)
    # half of Meat times t(Breadvec)
    if (is.null(adjustment)){
    } else if (adjustment == "HC2"){
      hat_adjust <- diag(X[I,]%*%((t(X[I,])%*%X[I,])^-1)%*%t(X[I,]))
      residarray[I,] <- residarray[I,]/sqrt(1-hat_adjust)
    } else if (adjustment == "HC3"){
      hat_adjust <- diag(X[I,]%*%((t(X[I,])%*%X[I,])^-1)%*%t(X[I,]))
      residarray[I,] <- residarray[I,]/(1-hat_adjust)      
    }
    e = array(residarray[I,],c(1,Ns,datavec))
    Sinit[] = apply(e,3,function(x)x%*%t(Breadvec[,I]))
    Sb = Sb + array(apply(Sinit,3,function(x)t(x)%*%x),dim=c(npred,npred,datavec))
    return(RecursiveNested(iter+1,nests,nested,datavec,Sinit,residarray,Breadvec,npred,Sb))
  }
}