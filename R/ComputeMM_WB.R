#' ComputeMM_WB -- calculate the maximum cluster extent from a wild bootstrap
#'
#' This function wraps many of the other functions to perform cluster detection analysis using
#' the wild bootstrap.
#' @param cifti_map surface array or nifti containing the scalar values -- just used to determine the number of participants
#' @param zscore_map The z-scored map from the geeglm analysis
#' @param resid_map The residual map from the geeglm analysis
#' @param fit_map The fitted values map from the geeglm analysis
#' @param type A character string denoting the distribution to use for the wild bootstrap
#' @param external_df A data frame comprising non-brain measures to model.
#' @param notation A character string representing the model formula. Uses Wilkinson notation.
#' @param family_dist A character string denoting the assumed underlying distribution of all input data
#' @param structtype A character string denoting whether the map is volumetric ('volume') or surface-based
#' @param matlab_path A character string denoting the path to the matlab compiler. Please use v91.
#' @param surf_command A character string denoting the path to the SurfConnectivity command.
#' @param correctiontype A character string denoting cluster ('cluster') or point ('point') comparisons
#' @keywords wild bootstrap
#' @export
#' @examples
#' max_cc <- ComputeMM_WB(cifti_map,zscore_map,resid_map,fit_map,type,external_df,
#' notation,family_dist,structtype,thresh,structfile,matlab_path,surf_command,correctiontype)
ComputeMM_WB <- function(resid_map,
                         fit_map,
                         type,
                         external_df,
                         notation,
                         family_dist,
                         structtype,
                         thresh,
                         structfile,
                         matlab_path,
                         surf_command,
                         corstr,
                         wave,
                         zcor,
                         correctiontype,
                         id_subjects) {
  library("purrr")
  library("cifti")
  library("gifti")
  library("geepack")
  library("Matrix")
  library("oro.nifti")
  library("mmand")
  ncomps <- length(unlist(resid_map[1]))
  if (type=='radenbacher'){
    bootstrap_mult <- sample(c(-1,1),ncomps,replace=TRUE)
  } else if (type=='mammen'){
    bootstrap_mult <- sample(c((1+sqrt(5))/2,(1-sqrt(5))/2),ncomps,prob=c((sqrt(5)-1)/(2*sqrt(5)),(sqrt(5)+1)/(2*sqrt(5))),replace=TRUE)
  } else if (type=='webb4') {
    bootstrap_mult <- sample(c(sqrt(3/2)*-1,sqrt(1/2)*-1,sqrt(1/2),sqrt(3/2)),ncomps,replace=TRUE)
  } else if (type=='webb6') {
    bootstrap_mult <- sample(c(sqrt(3/2)*-1,-1,sqrt(1/2)*-1,sqrt(1/2),1,sqrt(3/2)),ncomps,replace=TRUE)
  } else 
  {
    return('error: invalid argument for distribution type, please specify one of: radenbacher,mammen,webb4,webb6')
  }
  fit_mat= matrix(unlist(fit_map),nrow=ncomps,ncol=length(fit_map))
  resid_mat = matrix(unlist(resid_map),nrow=ncomps,ncol=length(resid_map))
  cifti_mat = fit_mat + resid_mat*bootstrap_mult
  cifti_bootmap = relist(cifti_mat,skeleton=fit_map)
  cifti_bootindex <- 1:length(cifti_bootmap)
  cifti_bootscalarmap <- map(cifti_bootindex,ReframeCIFTIdata,cifti_rawmeas=cifti_bootmap)  
#  MM_bootmap = list()
#  for (curr_cifti_bootmeas in 1:length(cifti_bootscalarmap)){
#    cifti_bootmeas <- cifti_bootscalarmap[curr_cifti_bootmeas]
#    if (sum(is.na(unlist(cifti_bootmeas))) > 0){
#      cifti_bootmeasb = data.frame(y=numeric(length(cifti_bootmeas)))
#      data_to_fit <-  cbind(cifti_bootmeasb,external_df)
#    } else
#    {
#      data_to_fit <-  cbind(cifti_bootmeas,external_df)
#    }
#    MM_bootmap[curr_cifti_meas] <- geeglm(notation, data=data_to_fit, id=data_to_fit[[id_subjects]], family=family_dist,
#                                         corstr=corstr, waves=wave,zcor=zcor)
#  }
#  data_to_fit <- external_df 
  MM_bootmap <- map(cifti_bootscalarmap,ComputeMM,external_df=external_df,notation=notation,family_dist=family_dist,corstr=corstr,zcor=zcor,wave=wave,id_subjects=id_subjects)
  zscore_bootmap <- map(MM_bootmap,ComputeZscores)
  thresh_bootmap <- map(zscore_bootmap,ThreshMap,zthresh=thresh)
  if (correctiontype=='point') {
    return(max(zscore_bootmap))
  }
  else if (correctiontype=='cluster'){
    varlist <- all.vars(notation)
    nmeas <- length(varlist)
    all_cc = list(1:nmeas)
    for (curr_meas in 1:nmeas){
      thresh_bootarray = unlist(thresh_bootmap)
      mask_bootvector = 1:nmeas == curr_meas
      thresh_bootarray <- thresh_bootarray[mask_bootvector]
      thresh_bootarray <- as.numeric(thresh_bootarray)
      thresh_bootarray[is.na(thresh_bootarray)] <-  0
      if (structtype=='volume'){
        cc <- GetVolAreas(thresh_bootarray)
        } else
        {
        cc <- GetSurfAreas(thresh_bootarray,structfile,matlab_path,surf_command)
        }
      all_cc[curr_meas] = max(cc)
    }
      return(all_cc)
  }
}