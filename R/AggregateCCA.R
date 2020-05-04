#' AggregateCCA -- a function that pulls the weights for each canonical model and variable and returns a data frame 
#'
#' This function is used as an input for wordcloudCCA.R to visualize outputs
#' @param CCAfold CCA variable output from CrossFoldCCA.R
#' @param x_data data frame containing the X data, if provided, will create weights for each CCA mode correlated against the X domain
#' @param y_data data frame containing the Y data, if provided, will create weights for each CCA mode correlated against the Y domain
#' @param threshold a positive number indicating the threshold to determine whether a variable is contributing
#' @param filename character vector representing the filename for the csv output. If provided will create a csv with the output data frame.
#' @keywords CCA, wordcloud
#' @export
#' @usage
#' AggregateCCA(CCAfold = CCAfold10, x_data = NULL, y_data = NULL, filename = NULL)

AggregateCCA <- function(CCAfold, 
                         x_data = NULL,
                         y_data = NULL,
                         sig_threshold = NULL,
                         filename = NULL){
# pull relevant informations from CCAfold object. As of current version each mode represents either 5 or 7 lists:
  # 1) scores for the X domain organized as a 2D matrix (case X mode) -- training set
  # 1a)  scores for the X domain organized as a 2D matrix (case X mode) -- test set
  # 2) weights for the X domain organized as a 2d matrix (var X mode)
  # 3) scores for the Y domain organized as a 2D matrix (case X mode) -- training set
  # 3a) scores for the Y domain organized as a 2D matrix (case X mode) -- test set
  # 4) weights for the Y domain organized as a 2d matrix (var X mode)
  # 5) Canonical correlation per mode
  n_subs <- dim(CCAfold[[1]])[2]
  n_xvars <- dim(CCAfold[[2]])[1]
  n_yvars <- dim(CCAfold[[4]])[1]
  n_modes <- dim(CCAfold[[2]])[2]
  x_names <- row.names(CCAfolds10[[2]])
  y_names <- row.names(CCAfolds10[[4]])
  # default aggregates
  x_weights <- array(data = 0, dim = c(n_xvars,n_modes))
  y_weights <-  array(data = 0, dim = c(n_yvars,n_modes))
  nloops = length(CCAfold)/5
  if (is.null(sig_threshold)){
    for (curr_run in seq(1,length(CCAfold),5)) {
      x_weights = x_weights + CCAfold[[curr_run+1]]
      y_weights = y_weights + CCAfold[[curr_run+3]]
    } 
    x_weights = x_weights/nloops
    y_weights = y_weights/nloops
  }
  else {
      for (curr_run in seq(1,length(CCAfold),5)) {
        x_weights = x_weights + as.numeric(CCAfold[[curr_run+1]] >= sig_threshold)
        y_weights = y_weights + as.numeric(CCAfold[[curr_run+3]] >= sig_threshold)
      }
    x_weights = x_weights/nloops
    y_weights = y_weights/nloops
    }
  if (is.null(x_data)){
    if (is.null(y_data)){
      CCA_weights = list(x_weights = as.data.frame(x_weights,row.names = x_names), y_weights = as.data.frame(y_weights,row.names = y_names))      
    }
  }
  
  # pull aggregates for weights correlated with x domain
  if (is.null(x_data) == FALSE){
  }
  # pull aggregates for weights correlated with y domain
  if (is.null(y_data) == FALSE){
    
  }
  if (is.null(x_data) == FALSE){
    if (is.null(y_data) == FALSE){
      CCA_weights = list(x_weights = as.data.frame(x_weights,row.names = x_names), y_weights = as.data.frame(y_weights,row.names = y_names))  
    }
    else {
      CCA_weights = list(x_weights = as.data.frame(x_weights,row.names = x_names), y_weights = as.data.frame(y_weights,row.names = y_names))  
    }
  }
  else if(is.null(y_data) == FALSE)  {
    CCA_weights = list(x_weights = as.data.frame(x_weights,row.names = x_names), y_weights = as.data.frame(y_weights,row.names = y_names))  
  }
  #write outputs to files
  if (is.null(filename) == FALSE){
    write.csv(CCA_weights$x_weights,file = paste(filename,'_x_weights.txt',sep=''))
    write.csv(CCA_weights$y_weights,file = paste(filename,'_y_weights.txt',sep=''))
  }
return(CCA_weights)
}