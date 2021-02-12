#' LoadBrainMetrics -- A Modularlized function for MarginalModelCifti to load MR data
#'
#' This function loads the imaging data for ConstructMarginalModel. Created to further modularlize the R package.
#' @param ciftilist A data frame comprising a list of path/filenames to load
#' @param structtype A character string denoting whether the map is volumetric ('volume'), surface-based ('surface'), a pconn ('pconn'), or a NIFTI connectivity matrix ('niiconn').
#' @keywords wild bootstrap sandwich estimator marginal model CIFTI scalar
#' @export
#' @examples
#' cifti_data <- LoadBrainMetrics(cifitilist,structtype)
#'
   
LoadBrainMetrics <- function(ciftilist,structtype) {
  CiftiInputs = list()
  require(cifti)
  require(gifti)
  require(geepack)
  require(Matrix)
  require(oro.nifti)
  require(R.matlab)
  require(SparseM)
  if (structtype == 'volume'){
    zeros_array = NULL
    cifti_file <- PrepVolMetric(as.character(ciftilist$file[1]))
    cifti_dim <- dim(cifti_file)
    cifti_alldata <- array(dim = c(length(ciftilist$file),cifti_dim[1]*cifti_dim[2]*cifti_dim[3]))
    count = 1
    for (filename in ciftilist$file) {
      cifti_alldata[count,] = ConvertVolume(PrepVolMetric(as.character(ciftilist$file[count])),cifti_dim)
      count = count + 1
    }      
    remove(count)
    Nelm <- dim(cifti_alldata)[2]
    cifti_nonans <- sapply(1:dim(cifti_alldata)[1],function(x){ sum(is.na(cifti_alldata[x,]))==0})
    cifti_alldata <- cifti_alldata[cifti_nonans,]    
    if (norm_internal_data == TRUE){
      for (count in 1:Nelm){
        if (var(cifti_alldata[,count]) != 0){ 
          cifti_alldata[,count] <- ((cifti_alldata[,count] - mean(cifti_alldata[,count]))/sd(cifti_alldata[,count]))
        }
      }
    }
    cifti_scalarmap <- as.data.frame(cifti_alldata)
  } 
  if (structtype == 'surface') {
    zeros_array = NULL
    cifti_alldata <- t(data.frame((lapply(as.character(ciftilist$file),PrepSurfMetric))))
    Nelm <- dim(cifti_alldata)[2]
    cifti_dim <- Nelm
    cifti_nonans <- sapply(1:dim(cifti_alldata)[1],function(x){ sum(is.na(cifti_alldata[x,]))==0})
    cifti_alldata <- cifti_alldata[cifti_nonans,]    
    if (norm_internal_data == TRUE){
      for (count in 1:Nelm){
        if (var(cifti_alldata[,count]) != 0){ 
          cifti_alldata[,count] <- ((cifti_alldata[,count] - mean(cifti_alldata[,count]))/sd(cifti_alldata[,count]))
        }     
      }
    }
    if (fastSwE == FALSE){
      cifti_index <- 1:length(cifti_alldata)
      cifti_scalarmap <- map(cifti_index,ReframeCIFTIdata,cifti_rawmeas=cifti_alldata)
      cifti_dim=NULL
    } else
    {
     cifti_scalarmap <- cifti_alldata 
    }
  }
  if (structtype == 'pconn') {
    cifti_firstsub <- PrepCIFTI(as.character(ciftilist$file[1]))
    cifti_alldata <- array(data = NA, dim = c(length(ciftilist$file),dim(cifti_firstsub)[1]*(dim(cifti_firstsub)[1]-1)/2))
    zeros_array <- array(data=0,dim = c(dim(cifti_firstsub)[1],dim(cifti_firstsub)[2]))
    count = 1
    for (filename in ciftilist$file) { 
      cifti_temp <- PrepCIFTI(as.character(filename))
      cifti_alldata[count,] <- Matrix2Vector(cifti_temp)
      count = count + 1
    }
    Nelm <- dim(cifti_alldata)[2]
    cifti_dim <- Nelm
    cifti_nonans <- sapply(1:dim(cifti_alldata)[1],function(x){ sum(is.na(cifti_alldata[x,]))==0})
    cifti_alldata <- cifti_alldata[cifti_nonans,]
    if (norm_internal_data == TRUE){
      for (count in 1:Nelm){
        if (var(cifti_alldata[,count]) != 0) { 
          cifti_alldata[,count] <- ((cifti_alldata[,count] - mean(cifti_alldata[,count]))/sd(cifti_alldata[,count]))
        }     
      }
    }
    cifti_scalarmap <- as.data.frame(cifti_alldata)  
  }
  if (structtype == 'niiconn'){
    cifti_firstsub <- PrepNiiConnMetric(as.character(ciftilist$file[1]))
    cifti_dim <- dim(cifti_firstsub)
    cifti_alldata <- array(data = NA, dim = c(length(ciftilist$file),dim(cifti_firstsub)[1]*(dim(cifti_firstsub)[1]-1)/2))
    zeros_array <- array(data=0,dim = c(dim(cifti_firstsub)[1],dim(cifti_firstsub)[2]))
    count = 1
    for (filename in ciftilist$file) { 
      cifti_temp <- PrepNiiConnMetric(filename)
      cifti_alldata[count,] <- Matrix2Vector(cifti_temp)
      count = count + 1
    }
    Nelm <- dim(cifti_alldata)[2]
    cifti_dim <- Nelm
    cifti_nonans <- sapply(1:dim(cifti_alldata)[1],function(x){ sum(is.na(cifti_alldata[x,]))==0})
    cifti_alldata <- cifti_alldata[cifti_nonans,]
    if (norm_internal_data == TRUE){
      for (count in 1:Nelm){
        if (var(cifti_alldata[,count]) != 0) { 
          cifti_alldata[,count] <- ((cifti_alldata[,count] - mean(cifti_alldata[,count]))/sd(cifti_alldata[,count]))
        }     
      }
    }
    cifti_scalarmap <- as.data.frame(cifti_alldata)
  }
  if (structtype == 'subjectcsv'){
    cifti_firstsub = read.csv(as.character(ciftilist$file[1]),header = FALSE,sep = ',')
    cifti_dim = dim(cifti_firstsub)
    cifti_alldata <- array(data=0,dim = c(length(ciftilist$file),dim(cifti_firstsub)[1]))
    zeros_array <- array(data=0,dim = c(dim(cifti_firstsub)[1],dim(cifti_firstsub[2])))
    count = 1
    for (filename in ciftilist$file) { 
      cifti_temp <- read.csv(filename,header = FALSE,sep = ',')
      cifti_alldata[count,] <- as.array(cifti_temp)
      count = count + 1
    }
    Nelm <- dim(cifti_alldata)[2]
    cifti_nonans <- sapply(1:dim(cifti_alldata)[1],function(x){ sum(is.na(cifti_alldata[x,]))==0})
    cifti_alldata <- cifti_alldata[cifti_nonans,]
    if (norm_internal_data == TRUE){
      for (count in 1:Nelm){
        if (var(cifti_alldata[,count]) != 0) { 
          cifti_alldata[,count] <- ((cifti_alldata[,count] - mean(cifti_alldata[,count]))/sd(cifti_alldata[,count]))
        }     
      }
    }
    cifti_scalarmap <- as.data.frame(cifti_alldata)    
  }
  if (structtype == 'dataframe'){
    cifti_alldata = read.csv(as.character(ciftilist$file),header = TRUE, sep = ',')
    cifti_firstsub = cifti_alldata[1,]
    cifti_dim = dim(cifti_firstsub)
    zeros_array <- array(data=0,dim = c(dim(cifti_firstsub)[1],dim(cifti_firstsub[2])))
    Nelm <- dim(cifti_alldata)[2]
    cifti_nonans <- sapply(1:dim(cifti_alldata)[1],function(x){ sum(is.na(cifti_alldata[x,]))==0})
    cifti_alldata <- cifti_alldata[cifti_nonans,]
    if (norm_internal_data == TRUE){
      for (count in 1:Nelm){
        if (var(cifti_alldata[,count]) != 0) { 
          cifti_alldata[,count] <- ((cifti_alldata[,count] - mean(cifti_alldata[,count]))/sd(cifti_alldata[,count]))
        }     
      }
    }
    cifti_scalarmap <- as.data.frame(cifti_alldata)
  }
  CiftiInputs$Nelm = Nelm
  CiftiInputs$cifti_alldata = cifti_alldata
  CiftiInputs$cifti_scalarmap = cifti_scalarmap
  CiftiInputs$cifti_firstsub = cifti_firstsub
  CiftiInputs$zeros_array = zeros_array
  CiftiInputs$nonans = cifti_nonans
  CiftiInputs$cifti_dim = cifti_dim
  return(CiftiInputs)
}