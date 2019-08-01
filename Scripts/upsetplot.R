# Generates a UpSet plot to visualize the intersections of statistically significant gene expressions
#
# Example calls:
# upsetplot()

upsetplot <- function() {
	# Get list of files
    file_list <- list.files(path = "PRECOG_DMFS/Significant/")
	
	# Initialize an empty list
    lt = list()
	
	# Iterate through the list of files
    for (i in 1:length(file_list)){
		# Get the name of a file
        temp_data <- file_list[i]
		
		# Split the filename by period
        chopped = strsplit(temp_data, "\\.")
		
		# Convert to vector
        chopped <- unlist(chopped)
		
		# The name is the second part
        name = chopped[[2]]
		
		# Read the .tsv file
        temp_data_read <- read.csv(paste("PRECOG_DMFS/Significant/pvalue0.001/",temp_data,sep=""), sep = "\t")
		
		# Extract the genes (the first column)
        temp_genes = temp_data_read[,1]
		
		# Convert to vector
        temp_genes_vector <- unlist(temp_genes)
		
		# Add a new entry in the list, which holds the vector of genes for that file
        lt[[name]] = temp_genes_vector
    }
	
	# Make a combination matrix for intersection using the list lt
    m1 = make_comb_mat(lt, mode = "intersect", top_n_sets = 10)

	# Create the UpSet plot for
    UpSet(m1, pt_size = unit(2, "mm"), lwd = 1)
}