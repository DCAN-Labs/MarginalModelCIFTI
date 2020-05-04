#' AggregateCCA -- a function that pulls the weights for each canonical model and variable and returns a data frame 
#'
#' This function is used as an input for wordcloudCCA.R to visualize outputs
#' @param CCAfold CCA variable output from CrossFoldCCA.R
#' @param x_data data frame containing the X data, if provided, will create weights for each CCA mode correlated against the X domain
#' @param y_data data frame containing the Y data, if provided, will create weights for each CCA mode correlated against the Y domain
#' @param threshold a positive number indicating the threshold to determine whether a variable is contributing
#' @param filename character vector representing the filename for the csv output. If provided will create a csv with the output data frame.
#' @param crossfold If set to TRUE, will inform that the CCAfold variable contains train and test folds.
#' @keywords CCA, wordcloud
#' @export
#' @usage
#' AggregateCCA(CCAfold = CCAfold10, x_data = NULL, y_data = NULL, filename = NULL)

AggregateCCA <- function(CCAfold, 
                         x_data = NULL,
                         y_data = NULL,
                         sig_threshold = NULL,
                         filename = NULL,
                         crossfold=FALSE){
# pull relevant informations from CCAfold object. As of current version each mode represents either 5 or 7 lists:
  # 1) scores for the X domain organized as a 2D matrix (case X mode) -- training set
  # 1a)  scores for the X domain organized as a 2D matrix (case X mode) -- test set
  # 2) weights for the X domain organized as a 2d matrix (var X mode)
  # 3) scores for the Y domain organized as a 2D matrix (case X mode) -- training set
  # 3a) scores for the Y domain organized as a 2D matrix (case X mode) -- test set
  # 4) weights for the Y domain organized as a 2d matrix (var X mode)
  # 5) Canonical correlation per mode
  listsize=5
  skip=0
  if (crossfold){
    listsize=7
    skip=1
  }
  n_subs <- dim(CCAfold[[1]])[2]
  n_xvars <- dim(CCAfold[[2]])[1]
  n_yvars <- dim(CCAfold[[4]])[1]
  n_modes <- dim(CCAfold[[2]])[2]
  x_names <- row.names(CCAfolds10[[2]])
  y_names <- row.names(CCAfolds10[[4]])
  # default aggregates
  x_weights <- array(data = 0, dim = c(n_xvars,n_modes))
  y_weights <-  array(data = 0, dim = c(n_yvars,n_modes))
  nloops = length(CCAfold)/listsize
  if (is.null(sig_threshold)){
    for (curr_run in seq(1,length(CCAfold),listsize)) {
      x_weights = x_weights + CCAfold[[curr_run+1]]
      y_weights = y_weights + CCAfold[[curr_run+3+skip*2]]
    } 
    x_weights = x_weights/nloops
    y_weights = y_weights/nloops
  }
  else {
      for (curr_run in seq(1,length(CCAfold),listsize)) {
        x_weights = x_weights + as.numeric(abs(CCAfold[[curr_run+1]]) >= sig_threshold)
        y_weights = y_weights + as.numeric(abs(CCAfold[[curr_run+3+skip*2]]) >= sig_threshold)
      }
    x_weights = x_weights/nloops
    y_weights = y_weights/nloops
    }
  if (is.null(x_data)){
    if (is.null(y_data)){
      CCA_weights = list(x_weights = as.data.frame(x_weights,row.names = x_names),
                         y_weights = as.data.frame(y_weights,row.names = y_names))      
    }
  }
  
  # pull aggregates for weights correlated with x domain
  if (is.null(x_data) == FALSE){
    xCCA_corr <- array(data = 0,dim = c(n_xvars,n_modes))
    if(is.null(sig_threshold)){
      for (curr_run in seq(1,length(CCAfold),listsize)) {
        xCCA_corr <- xCCA_corr + t(cor(CCAfold[[curr_run]]*CCAfold[[curr_run+2+skip]],x_data,use="complete.obs"))
      }
    }
    else {
      for (curr_run in seq(1,length(CCAfold),listsize)) {
        xCCA_corr <- xCCA_corr + as.numeric(abs(t(cor(CCAfold[[curr_run]]*CCAfold[[curr_run+2+skip]],x_data,use="complete.obs"))) >= sig_threshold)
      }      
    }
    xCCA_corr = xCCA_corr/nloops
  }
  # pull aggregates for weights correlated with y domain
  if (is.null(y_data) == FALSE){
    yCCA_corr <- array(data = 0,dim = c(n_yvars,n_modes))
    if(is.null(sig_threshold)){
      for (curr_run in seq(1,length(CCAfold),listsize)) {
        yCCA_corr <- yCCA_corr + t(cor(CCAfold[[curr_run]]*CCAfold[[curr_run+2+skip]],y_data,use="complete.obs"))
      }
    }
    else {
      for (curr_run in seq(1,length(CCAfold),listsize)) {
        yCCA_corr <- yCCA_corr + as.numeric(abs(t(cor(CCAfold[[curr_run]]*CCAfold[[curr_run+2+skip]],y_data,use="complete.obs"))) >= sig_threshold)
      }      
    }
    yCCA_corr = yCCA_corr/nloops
  }
  if (is.null(x_data) == FALSE){
    if (is.null(y_data) == FALSE){
      CCA_weights = list(x_weights = as.data.frame(x_weights,row.names = x_names),
                         y_weights = as.data.frame(y_weights,row.names = y_names), 
                         xCCA_corr = as.data.frame(xCCA_corr,row.names = x_names),
                         yCCA_corr = as.data.frame(yCCA_corr,row.names = y_names))  
    }
    else {
      CCA_weights = list(x_weights = as.data.frame(x_weights,row.names = x_names),
                         y_weights = as.data.frame(y_weights,row.names = y_names), 
                         xCCA_corr = as.data.frame(xCCA_corr,row.names = x_names))
      }
  }
  else if(is.null(y_data) == FALSE)  {
    CCA_weights = list(x_weights = as.data.frame(x_weights,row.names = x_names),
                       y_weights = as.data.frame(y_weights,row.names = y_names), 
                       yCCA_corr = as.data.frame(yCCA_corr,row.names = y_names))
    }
  #write outputs to files
  if (is.null(filename) == FALSE){
    write.csv(CCA_weights$x_weights,file = paste(filename,'_x_weights.txt',sep=''))
    write.csv(CCA_weights$y_weights,file = paste(filename,'_y_weights.txt',sep=''))
    if (is.null(x_data) == FALSE){
      write.csv(CCA_weights$xCCA_corr, file = paste(filename,'_x_corr.txt',sep=''))
    }
    if (is.null(y_data) == FALSE){
      write.csv(CCA_weights$yCCA_corr, file = paste(filename,'_y_corr.txt',sep=''))
    }
  }
return(CCA_weights)
}