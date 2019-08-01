# Removes samples with have NA in the output columns from both gene and annotation dataset
#
# Parameters:
# project - a project name (file name without the file endings). Example: Breast cancer.GSE3494.HGU133A_EntrezCDF
#
# Example calls:
# removena("Breast cancer.GSE3494.HGU133A_EntrezCDF")
removena <- function(project) {
    # Read in the clinical annotation (info) file
    info=read.csv(paste("PRECOG_DMFS/Original/",project,".info.tsv",sep=""),sep="\t",header=TRUE)
    
    # Set the first column name of the info file to be "Array"
    colnames(info)[1] <- "Array"
    
    # Read in the gene expresion file
    data=readRDS(paste("PRECOG_DMFS/Original/",project,".data.RData",sep=""))  
	
	data_cols <- colnames(data)
	
	# TCGA dataset have a different format for storing sample names in the data and info files
	# This regex expression discards everything after and including the third period for the column names
	# of the data file.
	if (grepl("TCGA", project)) {
		data_cols <- lapply(data_cols, function(y) gsub("^([^.]*.[^.]*.[^.]*)..*$", "\\1", y))
	}
	
	# Set the first two columns to be named Name and Description, respectively
	data_cols[1] <- "Name"
	data_cols[2] <- "Description"
	
	# Convert data_cols to a vector
	data_cols <- unlist(data_cols)
	
	# Set the column names of data to be data_cols
	names(data) <- data_cols

	
	info[info=="NA"] <- NA

    # Select samples which don't have NA in the DMFS_Time and DMFS_Status columns
    notnasamples <- unlist(subset(info, (!(is.na(DMFS_Status)) & !(is.na(DMFS_Time))), select = (Array)))
	
	
	# Subset the info dataset to just the acceptable samples
    new_info <- info[info$Array %in% notnasamples,]
	
	
	# Convert notnasamples to correct format (remove factors)
    notnasamples <- as.character(unname(notnasamples))
	
	# Add "Description" and "Name" to beginning of notnasamples
    notnasamples <- c("Description", notnasamples)
    notnasamples <- c("Name", notnasamples)
	
	# Subset gene expression dataset to just the acceptable samples (and the "Name" and "Description" columns)
    new_data <- data[notnasamples]

    # Write both new datasets
    write.table(new_info, paste("PRECOG_DMFS/RemovedNA/",project,".info.tsv",sep=""), sep="\t", row.names=FALSE)
    saveRDS(new_data, paste("PRECOG_DMFS/RemovedNA/",project,".data.RData",sep=""))

    # Print what's happening
    print("Saving new datasets")
}