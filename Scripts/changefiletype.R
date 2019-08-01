# Converts objects in a directory between RData and .tsv files so that they can easily be visually inspected
#
# Example calls:
# rdatatotsv()

changefiletype <- function() {
require("tm")
	# Get a list of all the files in the path
    file_list <- list.files(path = "PRECOG_DMFS/Original")
	
	# Iterate over the file_list
	for (i in 1:length(file_list)){
		# Get the name of a file in file_list
        temp_data <- file_list[i]
		
		# If the file contains a certain ending (choose either .data.RData or .data.tsv)
        if (grepl(".data.tsv", temp_data)) {
			# Read the file appropriately using one of the two lines below
			temp_file <- read.csv(paste("PRECOG_DMFS/",temp_data,sep=""),sep="\t")
            #temp_file <- readRDS(paste("PRECOG_DMFS/",temp_data,sep=""))
			
			# Remove characters to extract just the project name
            project <- substr(temp_data, 1, nchar(temp_data)-4)
			
			# Save/write the file appropriately (choose either .RData or .tsv)
			saveRDS(temp_file, file = paste("PRECOG_DMFS/",project,".RData",sep=""))
		    #write.table(temp_file, paste("PRECOG_DMFS/",project,".tsv",sep=""), sep = "\t", row.names = FALSE)

        }
    }
}