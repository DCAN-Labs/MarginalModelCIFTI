#' WriteMatrixToCifti -- Output an R array matrix into an appropriate CIFTI format
#' 
#' This function will write an R data matrix as a CIFTI pconn file, properly indexed.
#' Unfortunately, this function requires a compiled matlab package, SurfConnectivity.
#' Contact Eric Feczko to get the package.
#' @param metric_data_matrix A numeric vector containing any values. Indexing must match the template file
#' @param ncols The number of columns in the matrix -- used to write the array to a text file
#' @param surf_template_file A character string representing the path to the surface template file, only used for indexing.
#' @param matlab_path A character string representing the path to the matlab compiler
#' @param surf_command A character string representing the path to the compiled SurfConnectivity folder
#' @param output_file A character string representing the path to the output gifti file
#' @keywords surface cluster
#' @export
#' @examples 
#' WriteMatrixToCifti(metric_data,ncols,surf_template_file,matlab_path,surf_command,output_file)
WriteMatrixToCifti <- function(metric_data,ncols,surf_template_file,matlab_path,surf_command,output_file,wb_command) {
  run_gifti_command = paste(surf_command,'/',"run_MakeCiftiFromText.sh",sep="")
  temp_text_file = paste(getwd(),'/temp.txt',sep='')
  write(as.numeric(metric_data),temp_text_file,ncolumns=ncols)
  poss_args= c(matlab_path,
               surf_template_file,
               temp_text_file,
               output_file,
               wb_command)
  system2(run_gifti_command,args=poss_args)
}