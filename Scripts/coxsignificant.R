# A function that find significant genes and combines two datasets
#
# Example calls:
# coxsignifiant()
coxsignificant <- function(project, method, sigtype, n) {
		
	# Print that the dataset is being processed
	print(paste("Now processing",project))	
	
	# Find significant with given values for the parameters
	num <- findsignificant(project = project, method = method, sigtype = sigtype, n = n)
	
	# Combine the gene expression and clinical annotation (info) files so that they can be
	# used for multivariate cox regression, random survival forests, and coxnet
	combine(project, n)
	
	return(num)
}