# Performs preprocessing steps and bootstrap analysis for 5, 10, 15, 20, 25, 30 covariates, using the project of choice
#
# Parameters:
# pca - TRUE means PCA should be applied to the dataset first
# do_cox - TRUE means that the cox bootstrapping analysis will be done
# do_rsf - TRUE means that the random survival forest bootstrapping analysis will be done
# do_coxnet - TRUE means that the ridge penalized cox regression bootstrapping analysis will be done
#
# Example calls:
# analysis("Breast cancer.GSE3494.HGU133A_EntrezCDF", pca = TRUE, do_cox = TRUE, do_rsf = TRUE, do_coxnet = FALSE)

analysis <- function() {
	print("You will be asked a series of questions. NOTE: there is no robust input validation")
	options(warn=-1)
	project <- readline(prompt = "Enter a project name: ")
	print("Removing samples with NA")
	removena(project)
	print("Normalizing, transforming, and imputing missing values")
	normalize(project)
	covariate_selection <- readline(prompt = "Do you want to do univariate feature selection (u), PCA (p), or diffusion maps (d)? Enter u, p, or d: ")
	do_cox <- readline(prompt = "Do you want to do Multivariate Cox Regression? Enter y or n: ")
	if (do_cox == "y") {
		do_cox <- TRUE
	} else if (do_cox == "n") {
		do_cox <- FALSE
	}
	else {
		stop("Invalid input")
	}
	
	do_rsf <- readline(prompt = "Do you want to do Random Survival Forests? Enter y or n: ")
	
	if (do_rsf == "y") {
		do_rsf <- TRUE
	} else if (do_rsf == "n") {
		do_rsf <- FALSE
	}
	else {
		stop("Invalid input")
	}
	
	do_coxnet <- readline(prompt = "Do you want to do Ridge Regularized Cox Regression? Enter y or n: ")
	
	if (do_coxnet == "y") {
		do_coxnet <- TRUE
	} else if (do_coxnet == "n") {
		do_coxnet <- FALSE
	}
	else {
		stop("Invalid input")
	}
	
	if (covariate_selection=="u") {
		method <- readline(prompt = "Do you want to select genes by a significance threshold (s) or try the 5, 10, 15, 20, 25, & 30 most significant genes (n)? Enter s or n: ")
		if (method=="s") {
			sigtype <- readline("Do you want to select genes lower than a specific p-value (p) or q-value (q)? Enter p or q: ")
			n <- as.numeric(readline("Enter a significance threshold: "))
		} else if (method=="n") {
			sigtype <- readline("Do you want to select a certain number of most significant genes by p-value (p) or q-value (q)? Enter p or q: ")
			print("I will try 5, 10, 15, 20, 25, and 30 genes")
		}
	} else if (covariate_selection=="p") {
		print("Performing PCA")
		covar <- readRDS(paste("PRECOG_DMFS/Split/",project,".DMFS.train.data.RData",sep=""))
		covar$Name <- NULL
		covar$Description <- NULL
		covar <- t(covar)
		covar.pca <- prcomp(covar, center = TRUE, scale. = TRUE)
		screeplot(covar.pca)
		std_dev <- covar.pca$sdev
		pr_var <- std_dev^2
		prop_varex <- pr_var/sum(pr_var)
		png("PRECOG_DMFS/Pictures/scree.png", width = 1500, height = 1000 )
		plot(prop_varex, xlab = "Principal Component",
             ylab = "Proportion of Variance Explained",
             type = "b")
		dev.off()
		png("PRECOG_DMFS/Pictures/cumscree.png", width = 1500, height = 1000 )
		plot(cumsum(prop_varex), xlab = "Principal Component",
              ylab = "Cumulative Proportion of Variance Explained",
              type = "b")
		dev.off()
		print("Saved scree plot and cumulative scree plot for PCA to Pictures directory")
		model <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.test.data.tsv",sep=""), sep = "\t")
		model <- t(model)
		model_pca <- predict(covar.pca, newdata = model)
		model_pca <- as.data.frame(model_pca)
		model_info <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.test.info.tsv",sep=""),sep="\t")
		new_model <- cbind(model_info, model_pca)
		covar_pca <- as.data.frame(covar.pca$x)
		covar_info <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.train.info.tsv",sep=""),sep="\t")
		new_covar <- cbind(covar, covar_pca)
		write.table(new_covar, paste("PRECOG_DMFS/Combined/",project,".covar.combined.tsv",sep=""), sep="\t", row.names=FALSE)
		#write.table(valid_comb, paste("PRECOG_DMFS/Combined/",projects,".valid.combined.tsv",sep=""), sep="\t", row.names=FALSE)
		write.table(new_model, paste("PRECOG_DMFS/Combined/",project,".model.combined.tsv",sep=""), sep="\t", row.names=FALSE)
		method <- "n"
		sigtype <- "p" # Doesn't matter
		print("I will try 5, 10, 15, 20, 25, and 30 principal components")
	} else if (covariate_selection=="d") {
		print("Performing Diffusion Mapping")
		covar <- readRDS(paste("PRECOG_DMFS/Split/",project,".DMFS.train.data.RData",sep=""))
		covar$Name <- NULL
		covar$Description <- NULL
		covar <- t(covar)
		model <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.test.data.tsv",sep=""), sep = "\t")
		model <- t(model)
		D <- as.matrix(dist(rbind(covar, model)))
		num_covar <- nrow(covar)
		num_model <- nrow(model)
		Dorig <- D[1:num_covar, 1:num_covar]
		DExt <- D[(num_covar+1):(num_covar+num_model), 1:num_covar]
		dmap <- diffuse(Dorig)
		dmapExt <- nystrom(dmap, DExt)
		
		model_info <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.test.info.tsv",sep=""),sep="\t")
		new_model <- cbind(model_info, dmapExt)
		
		covar_info <- read.csv(paste("PRECOG_DMFS/Split/",project,".DMFS.train.info.tsv",sep=""),sep="\t")
		new_covar <- cbind(covar_info, dmap$X)
		
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
		write.table(new_covar, paste("PRECOG_DMFS/Combined/",project,".covar.combined.tsv",sep=""), sep = "\t", row.names = FALSE)
		write.table(new_model, paste("PRECOG_DMFS/Combined/",project,".model.combined.tsv",sep=""), sep = "\t", row.names = FALSE)
		method <- "n"
		sigtype <- "p" # Doesn't matter
		print("I will try 5, 10, 15, 20, 25, and 30 diffusion coordinates")
	}
	
	
	if (sigtype == "p") {
		sigtype <- "pval"
	} else if (sigtype == "q") {
		sigtype <- "FDRtool_qval"
	}
	else {
		stop("Invalid input")
	}

	if (method == "n") {
		# Will go through 5, 10, 15, 20, 25, 30
		for (i in seq(from = 5, to = 30, by = 5)) {
			# Print the number
			print(paste(i, "covariates"))
			if (covariate_selection == "u") {
				# Call preprocess with that number
				preprocess(project = project, method = method, sigtype = sigtype, n = i)
			}
			# Read in the project file
			my_data <- read.csv(paste("PRECOG_DMFS/Combined/",project,".model.combined.tsv", sep = ""), sep = "\t")
						
			# Perform bootstrapped analysis
			bootanalysis(my_data = my_data, method = method, n = i, covariate_selection = covariate_selection, do_cox = do_cox, do_rsf = do_rsf, do_coxnet = do_coxnet)
		}
	} else if (method == "s") {
		num <- preprocess(project = project, method = method, sigtype = sigtype, n = n)
		my_data <- read.csv(paste("PRECOG_DMFS/Combined/",project,".model.combined.tsv", sep = ""), sep = "\t")
		bootanalysis(my_data = my_data, method = method, n = num, covariate_selection = covariate_selection, do_cox = do_cox, do_rsf = do_rsf, do_coxnet = do_coxnet)
	}
}