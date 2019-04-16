#' GetSurfAreas -- Get computed size of connected surface components
#' 
#' This function will compute the size of connected components in statistical surfaces.
#' Unfortunately, this function requires a compiled matlab package, SurfConnectivity.
#' Contact Eric Feczko to get the package.
#' @param metric_data A numeric array containing only binary values (e.g. 0;1)
#' @param surf_file A character string representing the path to the surface file
#' @param matlab_path A character string representing the path to the matlab compiler
#' @param surf_command A character string representing the path to the compiled surface command
#' @keywords surface cluster
#' @export
#' @examples 
#' component_sizes <- GetSurfAreas(metric_data,surf_file,matlab_path,surf_command)
GetSurfAreas <- function(metric_data,surf_file,matlab_path,surf_command) {
  run_surf_command = paste(surf_command,'/',"run_ComputeComponents.sh",sep="")
  metricchar <- shQuote(paste(c('[',metric_data,']'),collapse=" "))
  poss_args= c(matlab_path,
               surf_file,
               metricchar)
  cluster_surf <- system2(run_surf_command,args=poss_args,stdout=TRUE)
  cluster_num <- as.numeric(cluster_surf[(which(regexpr("dp_val",cluster_surf) > 0)+2):
                                           (which(regexpr("dp_val",cluster_surf) > 0)+32493)])
  return(cluster_num)
}