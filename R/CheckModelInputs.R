#' CheckModelInputs -- A Modularlized function for MarginalModelCifti to check loaded data
#'
#' This function checks the imaging data for ConstructMarginalModel. Created to further modularlize the R package.
#' @param external_df A data frame comprising non-brain measures to model. Can be specified as a string to a csv file with appropriate headers.
#' @param concfile A character string denoting a single column text file that lists the dscalars in the same order as the external_df and wave frames.
#' @param structtype A character string denoting whether the map is volumetric ('volume'), surface-based ('surface'), a pconn ('pconn'), or a NIFTI connectivity matrix ('niiconn').
#' @param structfile A character string denoting the structural map file, used for cluster detection on surfaces only.
#' @param matlab_path A character string denoting the path to the matlab compiler. Please use v91.
#' @param surf_command A character string denoting the path to the SurfConnectivity command.
#' @keywords wild bootstrap sandwich estimator marginal model CIFTI scalar
#' @export
#' @examples
#' cifti_data <- CheckModelInputs(concfile,)
