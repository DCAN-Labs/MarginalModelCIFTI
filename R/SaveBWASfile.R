#' SaveBWASfile -- A Modularlized function for MarginalModelCifti to save statistical outputs
#'
#' This function saves stat map data for ConstructMarginalModel. Created to further modularlize the R package.
#' @param BWAS_statmap A matrix containing a brain-wide assocation statistical map, can be anything as long as it matches the dimensions of the inputs
#' @param structtype A character string denoting whether the map is volumetric ('volume'), surface-based ('surface'), a pconn ('pconn'), or a NIFTI connectivity matrix ('niiconn').
#' @param output_prefix A character string depicting the prefix of the output filenames, usually depicting the map metric
#' @param measnames A list of measures (covariates) being tested by the model, used to label the maps appropriately.
#' @param surf_command A path to the SurfConnectivity toolbox used to perform CIFTI/GIFTI writing
#' @param wb_command A full path/filename to thte wb_command file, used for ciftis
#' @param zeros_array A 2D array filled with zeros representing the initial structure for pconns and niicons, used to speed up file-saving
#' @param matlab_path A full path/filename to the matlab runtime compiler engine, must be V91
#' @param sigtype A character string denoting whether the map is volumetric ('volume'), surface-based ('surface'), a pconn ('pconn'), or a NIFTI connectivity matrix ('niiconn').
#' @keywords wild bootstrap sandwich estimator marginal model CIFTI scalar
#' @export
#' @examples
#' SaveBWASfile(BWAS_statmap,structtype,output_prefix,measnames)
#'
SaveBWASfile <- function(BWAS_statmap,
                      structtype,
                      output_prefix,
                      measnames,
                      surf_command=NULL,
                      wb_command=NULL,
                      zeros_array=NULL,
                      matlab_path=NULL,
                      sigtype=NULL) {
  require(cifti)
  require(gifti)
  require(geepack)
  require(Matrix)
  require(oro.nifti)    
  if (structtype == 'surface'){
    for (curr_map in 1:dim(BWAS_statmap)[1]){
      WriteVectorToGifti(metric_data = BWAS_statmap[curr_map,],
                         surf_template_file = as.character(ciftilist[1,1]),
                         surf_command = surf_command,
                         matlab_path = matlab_path,
                         output_file = paste(output_directory,'/',output_prefix,measnames[curr_map],'.func.gii',sep=""))      
    }       
  } else
    if (structtype == 'volume'){
      for (curr_map in 1:dim(BWAS_statmap)[1]){
        temp_map <- RevertVolume(BWAS_statmap[curr_map,],cifti_dim)
        cifti_file[] <- temp_map
        writeNIfTI(nim = cifti_file,filename = paste(output_directory,'/',output_prefix,measnames[curr_map],sep=""))
      }
    }
  if (structtype == 'pconn'){
    if (sigtype == 'enrichment') {
      for (curr_map in 1:length(BWAS_statmap)) {
        write.table(array(unlist(BWAS_statmap[curr_map]),dim=c(sqrt(length(unlist(BWAS_statmap[curr_map]))),sqrt(length(unlist(BWAS_statmap[curr_map]))))),file = paste(output_directory,'/',output_prefix,measnames[curr_map],sep=""))
      }
    } else
    {
      for (curr_map in 1:dim(BWAS_statmap)[1]){
        temp_mat <- Matrix2Vector(pconn_data = zeros_array,
                                  pconn_vector = BWAS_statmap[curr_map,],
                                  direction = "to_matrix")
        WriteMatrixToCifti(metric_data = temp_mat,
                           ncols = dim(temp_mat)[1],
                           surf_template_file = as.character(ciftilist$file[1]),
                           matlab_path = matlab_path,
                           surf_command = surf_command,
                           output_file = paste(output_directory,'/',output_prefix,measnames[curr_map],sep=""),
                           wb_command = wb_command)
      }
    }
  }
  if (structtype == 'niiconn'){
    if (sigtype == 'enrichment') {
      for (curr_map in 1:length(BWAS_statmap)) {
        write.table(array(unlist(BWAS_statmap[curr_map]),dim=c(sqrt(length(unlist(BWAS_statmap[curr_map]))),sqrt(length(unlist(BWAS_statmap[curr_map]))))),file = paste(output_directory,'/',output_prefix,measnames[curr_map],sep=""))
      }
    } else
    {    
      for (curr_map in 1:dim(BWAS_statmap)[1]){
        temp_map <- Matrix2Vector(pconn_data = zeros_array,
                                  pconn_vector = BWAS_statmap[curr_map,],
                                  direction = "to_matrix")
        cifti_firstsub[] = temp_map
        writeNIfTI(nim=cifti_firstsub,filename=paste(output_directory,'/',output_prefix,measnames[curr_map],sep=""))
      }
    }
  }
  if (structtype == 'subjectcsv'){
    for (curr_map in 1:dim(BWAS_statmap)[1]){
      write.csv(x = curr_map,file = paste(output_directory,'/',output_prefix,measnames[curr_map],sep="",header=FALSE))
    }
  }
  if (structtype == 'dataframe'){
    for (curr_map in 1:dim(BWAS_statmap)[1]){
      write.csv(x = curr_map,file = paste(output_directory,'/',output_prefix,measnames[curr_map],sep="",header=FALSE))
    }
  }
}