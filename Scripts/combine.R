# Combines gene expression and clinical annotation (info) datasets into one file
# Note: uncomment the lines relating to valid if you want a validation dataset as well
combine <- function(projects, n) {
	# Read the gene expression datasets
	covar_gene <- read.csv(paste("PRECOG_DMFS/Significant/",n,"qvaluetrain/",projects,".covar.qvalue.data.tsv",sep=""),sep="\t",header=TRUE)
	#valid_gene <- read.csv(paste("PRECOG_DMFS/Significant/bootlowest/",projects,".qvalue.data.tsv",sep=""),sep="\t",header=TRUE)
	model_gene <- read.csv(paste("PRECOG_DMFS/Significant/",n,"qvaluetest/",projects,".model.qvalue.data.tsv",sep=""),sep="\t",header=TRUE)
	# Read the info datasets
	covar_info <- read.csv(paste("PRECOG_DMFS/Split/",projects,".DMFS.train.info.tsv",sep=""),sep="\t", header = TRUE)
	#valid_info <- read.csv(paste("PRECOG_DMFS/Split/",projects,".DMFS.valid.info.tsv",sep=""),sep="\t", header = TRUE)
	model_info <- read.csv(paste("PRECOG_DMFS/Split/",projects,".DMFS.test.info.tsv",sep=""),sep="\t", header = TRUE)
   
	# Transpose the gene expression datasets
	covar_gene <- t(covar_gene)
	model_gene <- t(model_gene)
	
	# Save the column names of the train gene expression dataset
	cnames <- covar_gene[1,]
	
	# Subset the train gene expression dataset to not include the first two columns
	covar_gene <- covar_gene[3:nrow(covar_gene),]
	
	# Set the column names of train gene expression dataset to cnames
	colnames(covar_gene) <- cnames
	
	#colnames(valid_gene) <- cnames
	
	# The test gene expression dataset has the same names
	colnames(model_gene) <- cnames
		
	# Cbind the covar_info and covar_gene datasets
	covar_comb <- cbind(covar_info, covar_gene)

	#valid_comb <- cbind(valid_info, valid_gene)
	
	# Cbind the model_info and model_gene datasets
	model_comb <- cbind(model_info, model_gene)

   # Write the new datasets
   write.table(covar_comb, paste("PRECOG_DMFS/Combined/",projects,".covar.combined.tsv",sep=""), sep="\t", row.names=FALSE)
   #write.table(valid_comb, paste("PRECOG_DMFS/Combined/",projects,".valid.combined.tsv",sep=""), sep="\t", row.names=FALSE)
   write.table(model_comb, paste("PRECOG_DMFS/Combined/",projects,".model.combined.tsv",sep=""), sep="\t", row.names=FALSE)
}