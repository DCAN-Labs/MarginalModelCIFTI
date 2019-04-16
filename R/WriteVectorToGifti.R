#' WriteVectorToGifti -- Output an R data vector into an appropriate GIFTI format
#' 
#' This function will write an R data vector as a GIFTI file, properly indexed.
#' Unfortunately, this function requires a compiled matlab package, SurfConnectivity.
#' Contact Eric Feczko to get the package.
#' @param metric_data A numeric vector containing any values. Indexing must match the template file
#' @param surf_template_file A character string representing the path to the surface template file, only used for indexing.
#' @param matlab_path A character string representing the path to the matlab compiler
#' @param surf_command A character string representing the path to the compiled SurfConnectivity folder
#' @param output_file A character string representing the path to the output gifti file
#' @keywords surface cluster
#' @export
#' @examples 
#' WriteVectorToGifti(metric_data,surf_template_file,matlab_path,surf_command,output_file)
WriteVectorToGifti <- function(metric_data,surf_template_file,matlab_path,surf_command,output_file) {
  run_gifti_command = paste(surf_command,'/',"run_MakeGiftiFromText.sh",sep="")
  temp_text_file = paste(getwd(),'/temp.txt',sep='')
  write(as.numeric(metric_data),temp_text_file,ncolumns=1)
  poss_args= c(matlab_path,
               surf_template_file,
               temp_text_file,
               output_file)
  system2(run_gifti_command,args=poss_args)
}