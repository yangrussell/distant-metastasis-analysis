# Splits a given dataset into training/testing parts
#
# Parameters:
# data - a dataset to be split
#
# Example calls:
# split("Breast cancer.GSE3494.HGU133A_EntrezCDF")
#
# To use a validation dataset, uncomment the lines relating to valid

split <- function(data) {
	# Set a seed for reproducibility and so that the gene expression and clinical annotation files shuffle the samples the same way
	set.seed(0)
	
	# Randomly permute the rows in the data
	rows <- sample(nrow(data))
	
	# Reorder data by the permuted rows
	data <- data[rows,]
	
	# Create split points in the data
	split1 <- round(nrow(data)*0.30)
	split2 <- round(nrow(data)*1)
	
	# Covar is everything from the beginning to split1
	covar <- data[1:split1,]
		
	# Model is everything after split1
	model <-  data[(split1+1):split2,]
	
	#valid <- data[split2:nrow(data),]
	
	# Create a list of covar and model datasets, which is automatically returned by R
	my_split <- list(covar, model)
}