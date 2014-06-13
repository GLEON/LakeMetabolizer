
load.all.data <- function(lake.name, data.path, checkMerge=TRUE){
 
	files <- dir(data.path)
	files <- files[grepl(paste(lake.name, '\\.', sep=''), files, ignore.case=TRUE)]
	files <- file.path(data.path, files)
 
	# error out if no matching files found
	if(length(files) < 1){
		stop('No files for site:', lake.name, ' found in path:', data.path)
	}
 
	# Treat the metadata file differently
	meta.indx <- grep('.meta', files)
	if(length(meta.indx) == 1){
		metadata <- load.meta(files[meta.indx])
		files <- files[meta.indx*-1]
	}else{ 
		if(length(meta.indx) > 1){
			stop('Metadata file pattern matched more than one file')
		}else{
					metadata <- NA
		}
	}
 
	 #Drop bathy if we find it (could maybe handle this if we wanted)
	 files <- files[!grepl('.bth', files)]
 
	 #Load and merge all Timeseries files
	 data <- load.ts(files[1])
  
	for(i in 2:length(files)){
		tmp <- load.ts(files[i])
		# Hmm, if datetime in a file is screwed up, this blows up
		# I wonder if I could devise a good check. Hmm.
		# See pred.Merge in helper.functions.R ~RDB
		merge.out <- pred.merge(data, tmp, all=TRUE)
		if(merge.out$merge > 1.5*merge.out$desired){
			if(checkMerge){
				stop("a data frame to be merged contains many duplicated values in merge column")
			}else{
				warning("a data frame to be merged contains many duplicated values in merge column")
			}
		}

		data <- merge(data, tmp, all=TRUE, by='datetime')
	}
  
	return(list(data=data, metadata=metadata))
}
