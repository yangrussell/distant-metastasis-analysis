# Find significant genes
#
# Parameters:
# projects - the name of the file without file extensions. Example: Breast cancer.GSE3494.HGU133A_EntrezCDF
# method - either threshold or nlowest. Using threshold will select genes which have a significance less than the given threshold.
#  		    Using nlowest will select the n genes with the lowest significance value.
# n - if method is set to threshold, then n represents the significance threshold
# 	  if method is set to nlowest, then n represents the number of genes that will be selected
# 
# Example calls:
# findsignificant("Breast cancer.GSE3494.HGU133A_EntrezCDF", method = "threshold", n = 0.05)
# findsignificant("Breast cancer.GSE3494.HGU133A_EntrezCDF", method = "nlowest", n = 10)
#
# To use a validation dataset, uncomment the lines relating to valid

findsignificant <- function(projects, method, sigtype, n) {
	
    # Read in the gene expresion files
	data=readRDS(paste("PRECOG_DMFS/Split/",projects,".DMFS.train.data.RData",sep=""))
	data_test=read.csv(paste("PRECOG_DMFS/Split/",projects,".DMFS.test.data.tsv",sep=""), sep = "\t")
	#data_valid=read.csv(paste("PRECOG_DMFS/Split/",projects,".DMFS.valid.data.tsv",sep=""),sep="\t")

    # Read in the survres file
    survres=read.csv(paste("C:/Users/russe/Downloads/PRECOG_DMFS/Results/",projects,".DMFS.survres.tsv",sep=""),sep="\t",header=TRUE)

	# Set the first column name in survres to be "ID_REF"
    colnames(survres)[1] <- "ID_REF"
    
	# If the method being used is threshold, then we subset genes for which the sigtype (either pval or FDRtool_qval) is less than n
	if (method == "s") {
		print("Correct")
		if (sigtype == "pval") {
			sortedsurvres <- subset(survres, pval<n, select = (ID_REF))
		}
		if (sigtype == "FDRtool_qval") {
			sortedsurvres <- subset(survres, FDRtool_qval<as.numeric(n), select = (ID_REF))
		}
	}
	
	# If the method being used is nlowest, then we subset the n genes with lowest sigtype
	if (method == "n") {
		if (sigtype == "pval") {
			sortedsurvres <- survres[order(survres$pval),]
			sortedsurvres <- sortedsurvres[1:n,]
		}
		if (sigtype == "FDRtool_qval") {
			sortedsurvres <- survres[order(survres$FDRtool_qval),]
			sortedsurvres <- sortedsurvres[1:n,]
		}
	}
    
	# Significant is a list of the selected genes
    significant <- subset(sortedsurvres, select = (ID_REF))
	    
    # Print what's happening
    print(paste(nrow(significant), "genes found ..."))
	
	num <- nrow(significant)

    # Convert to correct format (remove factors)
    significant <- as.character(unname(unlist(significant)))
	
	#significant <- c("4609_at", significant)
    
	# Print the list of genes that have been subsetted
    print(significant)

	dir.create(paste("C:/Users/russe/Downloads/PRECOG_DMFS/Significant/",n,"qvaluetrain", sep = ""))
	dir.create(paste("C:/Users/russe/Downloads/PRECOG_DMFS/Significant/",n,"qvaluetest", sep = ""))

    # Subset the train gene expression dataset and write it
	# To change the location where the file is saved, change the file path (the first argument in the paste function)
    new_data <- data[data$Name %in% significant,]
    write.table(new_data, paste("PRECOG_DMFS/Significant/",n,"qvaluetrain/",projects,".covar.qvalue.data.tsv",sep=""), sep="\t", row.names=FALSE)

	# Subset the test gene expression dataset and write it
	# To change the location where the file is saved, change the file path (the first argument in the paste function)
	new_test = data_test[data$Name %in% significant,]
	write.table(new_test, paste("PRECOG_DMFS/Significant/",n,"qvaluetest/",projects,".model.qvalue.data.tsv",sep=""), sep="\t", row.names=FALSE)
	
	#new_valid = data_valid[data$Name %in% significant,]
	#write.table(new_valid, paste("PRECOG_DMFS/Significant/allvalidlowest10qvalue/",projects,".qvalue.data.tsv",sep=""), sep="\t", row.names=FALSE)
	

    # Print what's happening
    print("Saving new .tsv files ...")
	
    return(num)
}