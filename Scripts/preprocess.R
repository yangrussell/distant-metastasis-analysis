# A master function to perform all preprocessing steps on all the datasets
#
# Example calls:
# preprocess()
preprocess <- function(project, method, sigtype, n) {
		
	# Print that the dataset is being processed
	print(paste("Now processing",project))

	# Call removena
	removena(project)
	
	# Run the Cox function
	#Run_Cox("Breast cancer.Vijver_BreastCancer.Vijver_Agilent")	
	
	# Find significant with given values for the parameters
	num <- findsignificant(project = project, method = method, sigtype = sigtype, n = n)
	
	# Combine the gene expression and clinical annotation (info) files so that they can be
	# used for multivariate cox regression, random survival forests, and coxnet
	combine(project, n)
	
	return(num)
}