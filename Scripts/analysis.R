# Asks user questions, then performs preprocessing steps and bootstrap analysis.
#
# Example calls:
# analysis()

analysis <- function() {
	print("You will be asked a series of questions. NOTE: there is no robust input validation")
	options(warn=-1) # Suppress warnings related to dir.create()
	
	# Ask the user for a project name
	project <- readline(prompt = "Enter a project name: ")
	
	# Regardless of the project, we will always need to do certain preprocessing setps
	print("Removing samples with NA")
	# Call removena
	removena(project)
	print("Normalizing, transforming, and imputing missing values")
	# Call normalize to do basic preprocessing
	normalize(project)
	
	# Ask the user for a type of feature selection/dimensionality reduction
	covariate_selection <- readline(prompt = "Do you want to do univariate feature selection (u), PCA (p), or diffusion maps (d)? Enter u, p, or d: ")
	
	# Ask the user if they want to do Multivariate Cox Regression, do basic processing
	do_cox <- readline(prompt = "Do you want to do Multivariate Cox Regression? Enter y or n: ")
	if (do_cox == "y") {
		do_cox <- TRUE
	} else if (do_cox == "n") {
		do_cox <- FALSE
	}
	else {
		stop("Invalid input")
	}
	
	# Ask the user if they want to do Random Survival Forests, do basic processing
	do_rsf <- readline(prompt = "Do you want to do Random Survival Forests? Enter y or n: ")
	if (do_rsf == "y") {
		do_rsf <- TRUE
	} else if (do_rsf == "n") {
		do_rsf <- FALSE
	}
	else {
		stop("Invalid input")
	}
	
	# Ask the user if they want to do Ridge Regularized Cox Regression, do basic processing
	do_coxnet <- readline(prompt = "Do you want to do Ridge Regularized Cox Regression? Enter y or n: ")
	if (do_coxnet == "y") {
		do_coxnet <- TRUE
	} else if (do_coxnet == "n") {
		do_coxnet <- FALSE
	}
	else {
		stop("Invalid input")
	}
	
	# If the user wants to do univariate feature selection, then they have a choice between using a significance threshold or trying the 5, 10, 15, 20, 25, and 30
	# most significant genes
	if (covariate_selection=="u") {
		# Ask which method they want to use
		method <- readline(prompt = "Do you want to select genes by a significance threshold (s) or try the 5, 10, 15, 20, 25, & 30 most significant genes (n)? Enter s or n: ")
		if (method=="s") {
			# If a significance threshold is used, then we need to know if the threshold is in terms of a p-value or q-value
			sigtype <- readline("Do you want to select genes lower than a specific p-value (p) or q-value (q)? Enter p or q: ")
			# Asks for the numeric value of the threshold
			n <- as.numeric(readline("Enter a significance threshold: "))
		} else if (method=="n") {
			# We still need to know p-value or q-value if we are trying the 5, 10, 15, 20, 25, and 30 most significant genes
			sigtype <- readline("Do you want to select a certain number of most significant genes by p-value (p) or q-value (q)? Enter p or q: ")
			print("I will try 5, 10, 15, 20, 25, and 30 genes")
		}
	# If instead, the user wants to do PCA, then PCA will be conducted on the covariate selection set, and the same transformations will be applied to
	# the bootstrapping (model) set
	} else if (covariate_selection=="p") {
		print("Performing PCA")
		
		# Read the gene expression covariate selection set
		covar <- readRDS(paste("PRECOG_DMFS/Split/",project,".DMFS.train.data.RData",sep=""))
		
		# Remove the first two columns, as they are not relevant to the PCA
		covar$Name <- NULL
		covar$Description <- NULL
		
		# Transpose the matrix so that samples are the rows and genes are the columns
		covar <- t(covar)
		
		# Run PCA
		covar.pca <- prcomp(covar, center = TRUE, scale. = TRUE)
		
		# Get the standard deviation
		std_dev <- covar.pca$sdev
		
		# Square it to get the variance
		pr_var <- std_dev^2
		
		# Get the proportion of variance explained
		prop_varex <- pr_var/sum(pr_var)
		
		# Create a scree plot and save it
		png("PRECOG_DMFS/Pictures/scree.png", width = 1500, height = 1000 )
		plot(prop_varex, xlab = "Principal Component",
             ylab = "Proportion of Variance Explained",
             type = "b")
		dev.off()
		
		# Create a cumulative scree plot and save it
		png("PRECOG_DMFS/Pictures/cumscree.png", width = 1500, height = 1000 )
		plot(cumsum(prop_varex), xlab = "Principal Component",
              ylab = "Cumulative Proportion of Variance Explained",
              type = "b")
		dev.off()
		print("Saved scree plot and cumulative scree plot for PCA to Pictures directory")
		
		# Read the bootstrapping gene expression dataset
		model <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.test.data.tsv",sep=""), sep = "\t")
		
		# Transpose it so that samples are rows and genes are columns
		model <- t(model)
		
		# Apply the PCA computed earlier to the new dataset
		model_pca <- predict(covar.pca, newdata = model)
		
		# Convert to data frame
		model_pca <- as.data.frame(model_pca)
		
		# Read the bootstrapping clinical annotation dataset
		model_info <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.test.info.tsv",sep=""),sep="\t")
		
		# Create a new dataset with the clinical annotation on the left, and the new PCA components on the right
		new_model <- cbind(model_info, model_pca)
		
		# Convert the covariate selection set PCA componenets to a dataframe
		covar_pca <- as.data.frame(covar.pca$x)
		
		# Read the covariate selection set clinical annotatoin
		covar_info <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.train.info.tsv",sep=""),sep="\t")
		
		# Create a new dataset with the clinical annotation on the left, and the new PCA components on the right
		new_covar <- cbind(covar, covar_pca)
		
		# Write all the new combined datasets
		write.table(new_covar, paste("PRECOG_DMFS/Combined/",project,".covar.combined.tsv",sep=""), sep="\t", row.names=FALSE)
		write.table(new_model, paste("PRECOG_DMFS/Combined/",project,".model.combined.tsv",sep=""), sep="\t", row.names=FALSE)
		
		# We don't have a choice for method. We can only do the 5, 10, 15, 20, 25, 30 first PCs because it doesn't make sense to choose PCs by significance
		method <- "n"
		
		# This doesn't matter, it is just required as a parameter
		sigtype <- "p"
		
		print("I will try 5, 10, 15, 20, 25, and 30 principal components")
		
	# The user can also choose diffusion maps as the dimensionality reduction step
	} else if (covariate_selection=="d") {
		print("Performing Diffusion Mapping")
		
		# Read the gene expression covariate selection set
		covar <- readRDS(paste("PRECOG_DMFS/Split/",project,".DMFS.train.data.RData",sep=""))
		
		# Remove the first two columns, as they are not relevant to the PCA
		covar$Name <- NULL
		covar$Description <- NULL
		
		# Transpose the matrix so that samples are the rows and genes are the columns
		covar <- t(covar)

		# Read the bootstrapping gene expression dataset
		model <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.test.data.tsv",sep=""), sep = "\t")
		
		# Transpose it so that samples are rows and genes are columns
		model <- t(model)
		
		# Bind the covar and model together, then compute the pairwise Euclidean Distance Matrix
		D <- as.matrix(dist(rbind(covar, model)))
		
		# Extract the number of samples in the covariate selection set
		num_covar <- nrow(covar)
		
		# Extract the number of samples in the bootstrapping set
		num_model <- nrow(model)
		
		# Subset the distance matrix int otwo parts
		Dorig <- D[1:num_covar, 1:num_covar]
		DExt <- D[(num_covar+1):(num_covar+num_model), 1:num_covar]
		
		# Run diffusion process on the original datset
		dmap <- diffuse(Dorig)
		
		# Use the Nystrom extension to get approximate diffusion coordinates on the bootstrapping set
		dmapExt <- nystrom(dmap, DExt)
		
		# Read the bootstrapping clinical annotation dataset
		model_info <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.test.info.tsv",sep=""),sep="\t")
		# Create a new dataset with the clinical annotation on the left, and the new diffusion coordinates on the right
		new_model <- cbind(model_info, dmapExt)
		
		# Read the covariate selection set clinical annotation dataset
		covar_info <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.train.info.tsv",sep=""),sep="\t")
		# Create a new dataset with the clinical annotation on the left, and the new diffusion coordinates on the right
		new_covar <- cbind(covar_info, dmap$X)
		
		# Save diffusion map plots
		png("PRECOG_DMFS/Pictures/diffmapnewold.png", width = 1500, height = 1000 )
		plot(dmapExt[,1:2],pch=8,col=2,
			main="Diffusion map, black = original, red = new data",
			xlab="1st diffusion coefficient",ylab="2nd diffusion coefficient")
		points(dmap$X[,1:2],pch=19,cex=.5)
		dev.off()
		png("PRECOG_DMFS/Pictures/diffmap.png", width = 1500, height = 1000 )
		plot(dmap)
		dev.off()
		print("Saved diffusion map images to Pictures directory")
		
		# Write the new combined datasets
		write.table(new_covar, paste("PRECOG_DMFS/Combined/",project,".covar.combined.tsv",sep=""), sep = "\t", row.names = FALSE)
		write.table(new_model, paste("PRECOG_DMFS/Combined/",project,".model.combined.tsv",sep=""), sep = "\t", row.names = FALSE)
		
		# We don't have a choice for method. We can only do the 5, 10, 15, 20, 25, 30 first PCs because it doesn't make sense to choose PCs by significance
		method <- "n"
		
		# This doesn't matter, it is just required as a parameter
		sigtype <- "p"
		
		print("I will try 5, 10, 15, 20, 25, and 30 diffusion coordinates")
	}
	
	# Convert p to pval and q to FDRtool_qval
	if (sigtype == "p") {
		sigtype <- "pval"
	} else if (sigtype == "q") {
		sigtype <- "FDRtool_qval"
	}
	else {
		stop("Invalid input")
	}

	# If we are trying 5, 10, 15, 20, 25, 30 covariates
	if (method == "n") {
		# Will go through 5, 10, 15, 20, 25, 30
		for (i in seq(from = 5, to = 30, by = 5)) {
			# Print the number
			print(paste(i, "covariates"))
			if (covariate_selection == "u") {
				# Call preprocess with that number if doing univariate cox regression for feature selection
				coxsignificant(project = project, method = method, sigtype = sigtype, n = i)
			}
			# Read in the project file
			my_data <- read.csv(paste("PRECOG_DMFS/Combined/",project,".model.combined.tsv", sep = ""), sep = "\t")
						
			# Perform bootstrapped analysis
			bootanalysis(my_data = my_data, method = method, n = i, covariate_selection = covariate_selection, do_cox = do_cox, do_rsf = do_rsf, do_coxnet = do_coxnet)
		}
	# If we are using a significance threshold
	} else if (method == "s") {
		num <- coxsignificant(project = project, method = method, sigtype = sigtype, n = n)
		my_data <- read.csv(paste("PRECOG_DMFS/Combined/",project,".model.combined.tsv", sep = ""), sep = "\t")
		bootanalysis(my_data = my_data, method = method, n = num, covariate_selection = covariate_selection, do_cox = do_cox, do_rsf = do_rsf, do_coxnet = do_coxnet)
	}
}