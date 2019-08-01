# Runs univariate cox regression gene-by-gene
#
# Parameters:
# project - a project name (file name without the file endings). Example: Breast cancer.GSE3494.HGU133A_EntrezCDF
#
# Example calls:
# Run_Cox("Breast cancer.GSE3494.HGU133A_EntrezCDF")
#
# To use a validation dataset, uncomment the lines relating to valid
normalize <- function(project) {

# Require the necessary libraries
require(fdrtool)
require(qvalue) 
require(preprocessCore)

	# Read in the gene expression (data) and clinical annotation (info) datasets
    data=readRDS(paste("PRECOG_DMFS/RemovedNA/",project,".data.RData",sep=""))
    info=read.csv(paste("PRECOG_DMFS/RemovedNA/",project,".info.tsv",sep=""),sep="\t",header=TRUE)

	# Set the first and second column names of data to be "Name" and "Description", respectively
	colnames(data)[1] <- "Name"
	colnames(data)[2] <- "Description"
	
	# Set the first column name of info to be "Array"
	colnames(info)[1] <- "Array"
	
	# TCGA dataset have a different format for storing sample names in the data and info files
	# This regex expression discards everything after and including the third period for the column names
	# of the data file.
	if (grepl("TCGA", project)) {
		lapply(colnames(data), function(y) gsub("^([^.]*.[^.]*.[^.]*)..*$", "\\1", y))
	}
	
	#first_row_info <- info[,1]
	#first_row_info[1] <- "Name"
	#first_row_info <- prepend(first_row_info, "Name", before = 1)
	#colnames(info) <- first_row_info
	
	# Note: TCGA Array Error in eval(substitute(select), nl, parent.frame()) : 
    # object 'Array' not found


	# Save the names of data
	dhead=names(data)
	
	# Annot is the first two columns of data, "Name" and "Description"
	annot=data[,1:2]

	# Datan is the rest of data, the actual gene expression values
	datan=as.matrix(data[,3:ncol(data)]) 

	# Convert to log2 space if it is not already in it
	maxex=max(datan,na.rm=T)
	if (maxex > 100) {
	datan=log2(datan+1)
	datan=ifelse(is.finite(datan),datan,NA)
	}
	
	# Normalize quantiles if it is not a TCGA dataset
	if(grepl("TCGA",project)){datan <- datan}else{datan=normalize.quantiles(datan,copy=T)}


	# Calculate the standard deviation, mean, and median
	sds=apply(datan,1,sd,na.rm=T)
	mns= apply(datan,1,mean,na.rm=T)
	mds= apply(datan,1,median,na.rm=T)

	# Filter out genes with low variation (standard deviation less than 0.00001)
	lo=which(sds<0.00001)
	if (length(lo)>0) {
		annot=annot[-lo,]
		datan=datan[-lo,]
		sds=sds[-lo]
		mns=mns[-lo]
		mds=mds[-lo]
	}

	# Filter out rows/cols with more than 80% NAs
	xy=dim(datan)
	rownas=rowSums(is.na(datan))
	colnas=colSums(is.na(datan))
	rrow=which(rownas>=0.8*xy[2])
	rcol=which(colnas>=0.8*xy[1])
	if (length(rrow)>0) {
		annot=annot[-rrow,]
		datan=datan[-rrow,]
		sds = sds[-rrow]
		mns = mns[-rrow]
		mds = mds[-rrow]
	}
	if (length(rcol)>0) {
		datan=datan[,-rcol]
		info=info[-rcol,]
	}

	# Standardize genes so they have mean of 0 and variance of 1
	datan = (datan-mns)/sds

	# Impute missing values using k-nearest-neighbors imputation
	if (length(which(is.na(datan)))>0) {
	library(impute) #modify myself to get xx in correct format
	xx=impute.knn(datan)
	datan=xx$data
	cat(paste("# IMPUTEDMISSING: TRUE\n"),file=project)
	} else {
	cat(paste("# IMPUTEDMISSING: FALSE\n"),file=project)
	}

	# Call split to split datan
	my_split <- split(t(datan))
	
	# Print what's happening
	print("Shuffled and split datan")
	
	# Extract the train part
	train <- t(my_split[[1]])
	
	#valid_orig <- t(my_split[[3]])
	
	# Extract the test part
	test <- t(my_split[[2]])
	print("Extracted gene expression split")
	
	# Convert all to numeric
	class(train) <- "numeric"
	#class(valid) <- "numeric"
	class(test) <- "numeric"
	
	# Print what's happening
	print("Converted to numeric")

	# Split the info file
	my_split_info <- split(info)
	
	# Print what's happening
	print("Shuffled and split info")
	
	# Extract the train prat
	train_info <- my_split_info[[1]]
	
	#valid_info <- my_split_info[[3]]
	
	# Extract the test part
	test_info <- my_split_info[[2]]
	
	# Cbind annot back to train
	train <- cbind(annot, train)

	# Write the train and test gene expression (data) files
	write.table(train, file=paste("PRECOG_DMFS/Split/",project,".","DMFS",".train.data.tsv",sep=""),sep="\t",row.names=F,quote=F,col.names=T)
	#write.table(valid, file=paste("PRECOG_DMFS/Split/",project,".","DMFS",".valid.data.tsv",sep=""),sep="\t",row.names=F,quote=F,col.names=T)
	write.table(test, file=paste("PRECOG_DMFS/Split/",project,".","DMFS",".test.data.tsv",sep=""),sep="\t",row.names=F,quote=F,col.names=T)
	
	# Write the train and test clinical annotation (info) files
	write.table(train_info, file=paste("PRECOG_DMFS/Split/",project,".","DMFS",".train.info.tsv",sep=""),sep="\t",row.names=F,quote=F,col.names=T)
	#write.table(valid_info, file=paste("PRECOG_DMFS/Split/",project,".","DMFS",".valid.info.tsv",sep=""),sep="\t",row.names=F,quote=F,col.names=T)
	write.table(test_info, file=paste("PRECOG_DMFS/Split/",project,".","DMFS",".test.info.tsv",sep=""),sep="\t",row.names=F,quote=F,col.names=T)
	
	# Also save the train gen expression file as a .data.RData object for convenience
	saveRDS(train, file=paste("PRECOG_DMFS/Split/",project,".","DMFS",".train.data.RData",sep=""))
}	
