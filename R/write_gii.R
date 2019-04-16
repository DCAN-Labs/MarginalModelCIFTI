write_gii <- function(gifti_object,outputfile,encoding = 'GZipBase64Binary'){
  requires(stringr)
  requires(base64enc)
  requires(mmap)
  #set defaults for data info
  def_endian = .Platform$endian
  def_encoding = encoding
  def_intent = 'NIFTI_INTENT_NONE'
  def_datatype = "NIFTI_TYPE_FLOAT32"
  def_externalfilename = ""
  def_externalfileoffset = ""
  def_offset = 0
  #edit object data attributes
  for (curr_array in 1:length(gifti_object$data_info[,1])) {
    #data storage is normal enforcing conventions here
    gifti_object$data_info[curr_array,]$ArrayIndexingOrder <- "ColumnMajorOrder"
    if (gifti_object$data_info[curr_array,]$DataType == ""){
      gifti_object$data_info[curr_array,]$DataType <-  def_datatype
      cat("Warning: data type is empty. It will be set to default: ", def_datatype) 
    } 
    if (gifti_object$data_info[curr_array,]$Intent == ""){
      gifti_object$data_info[curr_array,]$Intent = def_intent
      cat("Warning: Intent is empty. It will be set to default: ", def_intent)
    }
    gifti_object$data_info[curr_array,]$Encoding = def_encoding
    gifti_object$data_info[curr_array,]$Endian = def_endian
    gifti_object$data_info[curr_array,]$ExternalFileName = def_externalfilename
    gifti_object$data_info[curr_array,]$ExternalFileOffset = def_externalfileoffset
    if (gifti_object$data_info[curr_array,]$Encoding == 'ExternalFileBinary'){
      extfilename = gifti_object$data_info[curr_array,]$ExternalFileName
      if (extfilename == ""){
        outputfilename <- basename(outputfile)
        extfilename = paste(basename(outputfile),'.dat',sep="")
      }
       gifti_object$data_info[curr_array,]$ExternalFileName = paste(dirname(outputfile),extfilename,sep="/")
       gifti_object$data_info[curr_array,]$ExternalFileOffset = as.character(def_offset)
    } else if (gifti_object$data_info[curr_array,]$Encoding == 'ExternalFileBinary') {
      
    } else if (gifti_object$data_info[curr_array,]$Encoding == 'ExternalFileBinary') {
      
    } else if (gifti_object$data_info[curr_array,]$Encoding == 'ExternalFileBinary') {
      
    } else 
    {
      stop(paste('[GIFTI] Uknown data encoding: %s.',gifti_object$data_info[curr_array,]$Encoding))
    }
  }
#open outputfile and start adding XML
  sink(outputfile)
  cat('<?xml version="1.0" encoding="UTF-8"?>\n')
  cat('<!DOCTYPE GIFTI SYSTEM "http://www.nitrc.org/frs/download.php/115/gifti.dtd">\n')
  cat('<GIFTI Version="1.0"  NumberOfDataArrays="',length(gifti_datafile$data_info[,1]),'">\n',sep="")
#write MetaData
  cat(str_pad('<MetaData',width = nchar('<MetaData')+3,side="left"))
  if (length(gifti_object$meta)==0){
    cat('/>\n')
  } else
  {
    cat('>\n')
    for (curr_meta in 1:length(gifti_object$meta))
    {
      cat(str_pad('<MD>\n',width = nchar('<MD>\n')+6,side="left"))
      cat(str_pad(paste('<Name><! [CDATA[',names(gifti_object$meta[curr_meta]),']]></Name>\n',sep=""),
              width = nchar(paste('<Name><! [CDATA[',names(gifti_object$meta[curr_meta]),']]></Name>\n',sep=""))+9,
              side="left"))
      cat(str_pad(paste('<Value><! [CDATA[',gifti_object$meta[curr_meta],']]></Value>\n',sep=""),
              width = nchar(paste('<Value><! [CDATA[',gifti_object$meta[curr_meta],']]></Value>\n',sep=""))+9,
              side="left"))
      cat(str_pad('</MD>\n',width = nchar('</MD>\n')+6,side="left"))
    }
    cat(str_pad('</MetaData>\n',width = nchar('</MetaData>\n')+3,side="left"))
  }
#write labelfile
  cat(str_pad('<LabelTable',width = nchar('<LabelTable')+3,side="left"))
  if (length(gifti_object$label[,1])==0){
    cat('/>\n')
  } else
  {
    cat('>\n')
    for (curr_label in 1:length(gifti_object$label[,1]))
    {
      if (sum(is.na(gifti_object$label[curr_label,])) < length(gifti_object$label[curr_label,])){
        rgb_label = paste(' Red="',gifti_object$label[curr_label,"Red"],
                          '" Green="',gifti_object$label[curr_label,"Green"],
                          '" Blue="',gifti_object$label[curr_label,"Blue"],
                          '" Alpha="',gifti_object$label[curr_label,"Alpha"],
                          '"',
                          sep="")
      } else 
        {
        rgb_label = ''
        }
      cat(str_pad(paste('<Label Key="',gifti_object$label[curr_label,"Key"],rgb_label,'><! [CDATA[LabelKey',
                        gifti_object$label[curr_label,"Key"],']]></Label>\n',
                        sep=""),
                  width=nchar(paste('<Label Key="',gifti_object$label[curr_label,"Key"],rgb_label,'><! [CDATA[LabelKey',
                                                  gifti_object$label[curr_label,"Key"],']]></Label>\n',
                                                  sep=""))+6,
                  side="left"))
      cat(str_pad('</LabelTable>\n',width=nchar('</LabelTable>\n')+3,side="left"))      
    }
  }
#prepare data array
    #write attributes
    for (curr_array in 1:length(gifti_object$data_info[,1])){ 
      cat(str_pad(str_pad('<DataArray',width=nchar('<DataArray')+3,side="left"),
                  width=nchar(str_pad('<DataArray',width=nchar('<DataArray')+3,side="left"))+3,
                  side="right"))
      if (def.offset != 0) {
        gifti_object$data_info[curr_array,]$ExternalFileOffset = as.character(def.offset)
      }
      array_names <- sort(names(gifti_object$data_info[curr_array,]))
      for (curr_name in 1:length(array_names)){
        if (array_names[curr_name] == 'ExternalFileName'){
          curr_value <- basename(gifti_object$data_info[curr_array,array_names[curr_name]])
        } else
        {
          curr_value <- gifti_object$data_info[curr_array,array_names[curr_name]]
        }
        cat(array_names[curr_name],'="',curr_value,'"')
        if (curr_name == length(array_names)) {
          cat('>\n')
        } else
        {
          cat('\n')
          cat(str_pad("",width=15,side="left"))
        }
      }
    #write MetaData
      cat(str_pad('<MetaData>\n',width=nchar('<MetaData>\n')+6,side="left"))
      if (length(gifti_object$data[curr_array]$meta) > 0) {
        for (curr_meta in 1:length(gifti_object$data[curr_array]$meta)) {
          cat(str_pad('<MD>\n',width=nchar('<MD>\n')+9,side='left'))
          cat(str_pad(cat('<Name><![CDATA[',names(gifti_object$data[curr_array]$meta[curr_meta]),']]></Name>\n'),
                      width=nchar(paste('<Name><![CDATA[',names(gifti_object$data[curr_array]$meta[curr_meta]),']]></Name>\n'))+12,
                      side="left"))
          cat(str_pad(cat('<Value><![CDATA[',gifti_object$data[curr_array]$meta[curr_meta],']]></Value>\n'),
                      width=nchar(paste('<Value><![CDATA[',gifti_object$data[curr_array]$meta[curr_meta],']]></Value>\n'))+12,
                      side="left"))
          cat(str_pad('</MD>\n',width=nchar('</MD>\n')+9,side="left"))
        }
      }
      cat(str_pad,'</MetaData>\n',width=nchar('</MetaData>\n'),side="left")
    #Write any coordinate system transformation matrices
      if (length(gifti_object$transformations[curr_array]) > 0){
        for (curr_trans in 1:length(gifti_object$transformations[curr_array])){
          cat(str_pad('<CoordinateSystemTransformMatrix>\n',width=nchar('<CoordinateSystemTransformMatrix>\n')+6,side="left"))  
          cat(str_pad(cat('<DataSpace><! [CDATA[',gifti_object$transformations[curr_array]$DataSpace,']]></DataSpace>\n'),
              width=nchar(paste('<DataSpace><! [CDATA[',gifti_object$transformations[curr_array]$DataSpace,']]></DataSpace>\n'))+9,
              side="left"))
          cat(str_pad(cat('<TransformedSpace><! [CDATA[',gifti_object$transformations[curr_array]$TransformedSpace,']]></TransformedSpace>\n'),
              width=nchar(paste('<TransformedSpace><! [CDATA[',gifti_object$transformations[curr_array]$TransformedSpace,']]></TransformedSpace>\n'))+9,
              side="left"))
          cat(str_pad(cat('<MatrixData>',gifti_object$transformations[curr_array]$MatrixData,'</MatrixData>\n'),
              width=nchar(paste('<MatrixData>',gifti_object$transformations[curr_array]$MatrixData,'</MatrixData>\n'))+9,
              side="left"))
          cat(str_pad('</CoordinateSystemTransformMatrix>\n',width=nchar('</CoordinateSystemTransformMatrix>\n')+6,side="left"))    
        }
      }
    #saving binary data
      cat(str_pad('<Data>',width=nchar('<Data>')+6,side="left"))
      sink()
      
      fid <- gzfile(outputfile,open="w")
      close(fid)
    }
}