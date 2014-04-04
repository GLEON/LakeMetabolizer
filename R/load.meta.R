#load.meta

load.meta = function(fPath){
  
  fid = file(fPath, 'rt', warn=FALSE)
  tmp = readLines(fid, 100)
  close(fid)
  
  tmp = strsplit(tmp, '[\t,]') #splits on tabs *and* commas
  
  header = tolower(tmp[[1]])
  
  id.col = which(header=='id')
  val.col = which(header=='value')
  
  if(length(id.col) == 0 || length(val.col) == 0){
    stop('Metadata file must contain "ID" and "Value" header elements.')
  }
  
  output = list()
  
  for(i in 2:length(tmp)){
    id = tolower(tmp[[i]][id.col])
    val = as.numeric(tmp[[i]][val.col])
    
    if(is.na(val)){
      val = tmp[[i]][val.col]
    }
    output[id] = val
  }
  
  return(output)  
  
}