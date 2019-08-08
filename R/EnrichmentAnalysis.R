#' WriteMatrixToCifti -- Output an R array matrix into an appropriate CIFTI format
#' 
#' This function will write an R data matrix as a CIFTI pconn file, properly indexed.
#' Unfortunately, this function requires a compiled matlab package, SurfConnectivity.
#' Contact Eric Feczko to get the package.
#' @param metric_data A binarized matrix
#' @param ncols The number of columns in the matrix -- used to write the array to a text file
#' @param modules A numeric vector representing the modules for the binarized modules
#' @param matlab_path A character string representing the path to the matlab compiler
#' @param enrichment_path A character string representing the path to the compiled enrichment folder
#' @param output_file A character string representing the path to the output chisquared file
#' @keywords surface cluster
#' @export
#' @examples 
#' WriteMatrixToCifti(metric_data,ncols,surf_template_file,matlab_path,enrichment_path,output_file)
EnrichmentAnalysis <- function(metric_data,ncols,modules,matlab_path,enrichment_path,output_file) {
  run_gifti_command = paste(enrichment_path,'/',"run_CountSignificantEffectsByModules.sh",sep="")
  temp_text_file_m = paste(getwd(),'/temp_m.txt',sep='')
  modules_data = read.csv(file = modules,header=TRUE,sep = '')
  temp_text_file_modules = paste(getwd(),'/temp_modules.txt',sep='')
  write(as.numeric(metric_data),temp_text_file_m,ncolumns=ncols)
  write(as.numeric(unlist(modules_data)),temp_text_file_modules,ncolumns = 1)
  poss_args= c(matlab_path,
               temp_text_file_m,
               temp_text_file_modules,
               'ExportText',
               output_file)
  system2(run_gifti_command,args=poss_args)
  chi2_mat <- read.csv(output_file,sep = " ",header=FALSE)
  return(chi2_mat)
}