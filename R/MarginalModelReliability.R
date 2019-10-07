#' MarginalModelReliability -- this function will perform a marginal model analysis and calculate reliability between two subsets of data
#'
#' This function wraps all of the other functions in the MarginalModelCifti package. The output is a significance map, which is either uncorrected or cluster detection was performed.
#' @param subset_file A csv file with two columns listing the subsets per group
#' @param group1_external_df A data frame comprising non-brain measures to model. Can be specified as a string to a csv file with appropriate headers.
#' @param group2_external_df A data frame comprising non-brain measures to model. Can be specified as a string to a csv file with appropriate headers.
#' @param group1_concfile A character string denoting a single column text file that lists the dscalars in the same order as the external_df and wave frames.
#' @param group2_concfile A character string denoting a single column text file that lists the dscalars in the same order as the external_df and wave frames.
#' @param structtype A character string denoting whether the map is volumetric ('volume'), surface-based ('surface'), or a pconn ('pconn').
#' @param structfile A character string denoting the structural map file, used for cluster detection on surfaces only.
#' @param matlab_path A character string denoting the path to the matlab compiler. Please use v91.
#' @param surf_command A character string denoting the path to the SurfConnectivity command.
#' @param group1_wave A data frame comprising the order of the waves (longitudinal measures) for each independent case.
#' @param group2_wave A data frame comprising the order of the waves (longitudinal measures) for each independent case.
#' @param notation A character string representing the model formula. Uses Wilkinson notation.
#' @param zcor A matrix defining the correlation structure set by the user, can be used to specify custom dependencies (e.g. site or kinship matrices).
#' @param corstr A character string that defines the correlation structure in a predetermined way. "Independence" is recommend for most cases.
#' @param family_dist A character string denoting the assumed underlying distribution of all input data.
#' @param dist_type A character string denoting the distribution to use for the wild bootstrap.
#' @param z_thresh A numeric that represents the threshold for Z statistics when performing cluster detection.
#' @param nboot A numeric that represents the number of wild bootstraps to perform.
#' @param p_thresh A numeric that represents the p value threshold for assessing significance.
#' @param sigtype A character string denoting cluster ('cluster'),  point ('point'), or enrichment ('enrichment') comparisons.
#' @param id_subjects A character string denoting the column for the subject ID. Only needed when FastSwE is FALSE
#' @param output_directory A character string denoting the path to the output MRI statistical maps
#' @param ncores An integer denoting the number of CPU cores to use when conducting permutation tests
#' @param fastSWE A boolean that determine the sandwhich estimator approach. If set to FALSE, will use standard R package geeglm. If set to TRUE will use custom-built estimator using rfast.
#' @param adjustment A character string denoting the small sample size adjustment to use when fastSwE is set to TRUE. Is NULL by default.
#' @param norm_external_data A boolean. If set to true, external data will be normed prior to analysis.
#' @param norm_internal_data A boolean. If set to true, MRI data will be normed per datapoint prior to analysis.
#' @param marginal_outputs A boolean. If set to true, marginal values will be output as statistical maps
#' @param marginal_matrix A numeric matrix depicting how to draw the map. Only needed if marginal_outputs is set to TRUE
#' @param enrichment_path A string depicting the path to the enrichment code. Used for enrichment analysis only (i.e. when sigtype is set to 'enrichment').
#' @param modules A csv file or array that depicts the modules for the enrichment analysis.
#' @param wb_command A character string denoting the path to the wb_command file
#' @#' @keywords wild bootstrap sandwich estimator marginal model CIFTI scalar
#' @export
#' @examples
#' max_cc <- ComputeMM_WB(cifti_map,zscore_map,resid_map,fit_map,type,external_df,
#' notation,family_dist,structtype,thresh,structfile,matlab_path,surf_command,correctiontype)


