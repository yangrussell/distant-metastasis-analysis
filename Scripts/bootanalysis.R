# Performs bootstrapped multivariate cox regression, random survival forest analysis, and coxnet
# The given data is subjected to sampling with replacement n times (where n is the number of samples
# in the dataset). Probabilistically, there will be about 1/e = 37% of the dataset that is not drawn
# (out-of-bag). This is used to compute an oob concordance index. The process outlined above is repeated
# 100 times, and all the train & oob concordance indices are reported. Histograms and qq (quantiles-quantiles)
# plots of the train and oob concordance indices are saved, and the mean scores are printed.
# 
# Parameters:
# my_data - a combined data file (contains the clinical annotation with a certain number of covariates IN THE SAME FILE)
# n - the number of covariates that are being used
# pca - TRUE if pca is being used
# do_cox - TRUE if cox analysis is to be conducted
# do_rsf - TRUE if random survival forest analysis is to be conducted
# do_coxnet - TRUE if ridge penalized cox regression is to be conducted
#
# Example calls:
# data <- read.csv(...)
# bootanalysis(data, n = 10, pca = FALSE, do_cox = TRUE, do_rsf = FALSE, do_coxnet = TRUE)

bootanalysis <- function(my_data, method, n, covariate_selection, do_cox, do_rsf, do_coxnet) {
	# Get the column names of the data
	cols <- colnames(my_data)
	if (method == "n") {
		# If PCA is not being used, then covariates that will be used appear at the end (right)
		# of the combined data file. Get the column numbers of the first and last columns
		if (covariate_selection == "u") {
			lastcol <- length(cols)
			firstcol <- lastcol - n + 1
		}
		
		# If PCA is being used, then all 50 PCA componenets are in the combined data file. They
		# are numbered PC1, PC2, PC3, ... , PC50. Get the column numbers of the first and last columns
		# that are necessary to get the specified number of components (n)
		if (covariate_selection == "d") {
			firstcol <- length(cols) - 49
			lastcol <- firstcol + n - 1
		}
		
		if (covariate_selection == "p") {
			firstcol <- length(cols) - 74
			lastcol <- firstcol + n - 1
		}
	} else if (method == "s") {
		lastcol <- length(cols)
		firstcol <- lastcol - n + 1
	}
		
		
	# Slice the column names to extract only the column names of the covariates we care about
	gene_cols <- cols[firstcol:lastcol]
	# Print them out
	print(gene_cols)
	
	# Coxnet cannot handle 0s in the time column. Replace
	# all 0s in the DMFS_Time column with 1e-10
	for (i in 1:length(my_data$DMFS_Time)) {
		if (my_data$DMFS_Time[i] == 0) {
			my_data$DMFS_Time[i] = 1e-10
		}
	}
	
	# Cox analysis
	if (do_cox) {
		
		# Use a boostrapping object, passing in (1) the data, (2) the cox function as a statistic (defined below),
		# (3) the number of repetitions, which is 100 in this case, and (4) the columns we care about, which are passed
		# as a parameter of the cox function
		cox_results <- boot(data = my_data, statistic = cox, R = 100, genes = gene_cols)

		# Save a histogram of the cox training concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxtrain",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		hist(cox_results$t[,1], xlab = "Train concordance index", col = "#7FBF7F", cex.lab = 3, cex.axis = 3)
		dev.off()
		
		# Save a qq plot of the cox training concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxtrainq",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		qqnorm(cox_results$t[,1], pch = 21, col = "blue", bg = "blue", cex = 2, xlab = "Theoretical quantiles", ylab = "Ordered values", cex.lab = 3, cex.axis = 3)
		qqline(cox_results$t[,1], col = "red")
		dev.off()
		
		# Save a histogram of the cox OOB concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxoob",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		hist(cox_results$t[,2], xlab = "OOB concordance index", col = "#7FBF7F", cex.lab = 3, cex.axis = 3)
		dev.off()
		
		# Save a qq plot of the cox oob concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxoobq",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		qqnorm(cox_results$t[,2], pch = 21, col = "blue", bg = "blue", cex = 2, xlab = "Theoretical quantiles", ylab = "Ordered values", cex.lab = 3, cex.axis = 3)
		qqline(cox_results$t[,2], col = "red")
		dev.off()

		# Print the final mean train and test concordance indices, rounded to 3 decimal places
		print("Cox train & test")
		print(round(mean(cox_results$t[,1]), digits = 3))
		print(round(mean(cox_results$t[,2]), digits = 3))
	}
	
	if (do_rsf) {
	
		# Use a boostrapping object, passing in (1) the data, (2) the rsf function as a statistic (defined below),
		# (3) the number of repetitions, which is 100 in this case, and (4) the columns we care about, which are passed
		# as a parameter of the rsf function
		rsf_results <- boot(data = my_data, statistic = rsf, R = 100, genes = gene_cols)

		# Save a histogram of the rsf training concordance indices
		png(paste("PRECOG_DMFS/Pictures/rsftrain",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		hist(rsf_results$t[,1], xlab = "Train concordance index", col = "#7FBF7F", cex.lab = 3, cex.axis = 3)
		dev.off()
		
		# Save a qq plot of the rsf train concordance indices
		png(paste("PRECOG_DMFS/Pictures/rsftrainq",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		qqnorm(rsf_results$t[,1], pch = 21, col = "blue", bg = "blue", cex = 2, xlab = "Theoretical quantiles", ylab = "Ordered values", cex.lab = 3, cex.axis = 3)
		qqline(rsf_results$t[,1], col = "red")
		dev.off()
		
		# Save a histogram of the rsf oob concordance indices
		png(paste("PRECOG_DMFS/Pictures/rsfoob",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		hist(rsf_results$t[,2], xlab = "OOB concordance index", col = "#7FBF7F", cex.lab = 3, cex.axis = 3)
		dev.off()
		
		# Save a qq plot of the rsf oob concordance indices
		png(paste("PRECOG_DMFS/Pictures/rsfoobq",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		qqnorm(rsf_results$t[,2], pch = 21, col = "blue", bg = "blue", cex = 2, xlab = "Theoretical quantiles", ylab = "Ordered values", cex.lab = 3, cex.axis = 3)
		qqline(rsf_results$t[,2], col = "red")
		dev.off()

		# Print the final mean train and oob rsf concordance indices, rounded to 3 digits
		print("RSF train & test")
		print(round(mean(rsf_results$t[,1]), digits = 3))
		print(round(mean(rsf_results$t[,2]), digits = 3))
	}
	
	if (do_coxnet) {
		# Use a boostrapping object, passing in (1) the data, (2) the coxnet function as a statistic (defined below),
		# (3) the number of repetitions, which is 100 in this case, and (4) the columns we care about, which are passed
		# as a parameter of the coxnet function
		coxnet_results <- boot(data = my_data, statistic = coxnet, R = 100, genes = gene_cols)
		
		# Save a histogram of the coxnet 1se training concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxnettrain1se",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		hist(coxnet_results$t[,1], xlab = "Train concordance index", col = "#7FBF7F", cex.lab = 3, cex.axis = 3)
		dev.off()
		
		# Save a qq plot of the coxnet 1se training concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxnettrain1seq",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		qqnorm(coxnet_results$t[,1], pch = 21, col = "blue", bg = "blue", cex = 2, xlab = "Theoretical quantiles", ylab = "Ordered values", cex.lab = 3, cex.axis = 3)
		qqline(coxnet_results$t[,1], col = "red")
		dev.off()
		
		# Save a histogram of the coxnet 1se oob concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxnetoob1se",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		hist(coxnet_results$t[,2], xlab = "OOB concordance index", col = "#7FBF7F", cex.lab = 3, cex.axis = 3)
		dev.off()
		
		# Save a qq plot of the coxnet 1se oob concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxnetoob1seq",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		qqnorm(coxnet_results$t[,2], pch = 21, col = "blue", bg = "blue", cex = 2, xlab = "Theoretical quantiles", ylab = "Ordered values", cex.lab = 3, cex.axis = 3)
		qqline(coxnet_results$t[,2], col = "red")
		dev.off()
		
		# Save a histogram of the coxnet min train concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxnettrainmin",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		hist(coxnet_results$t[,3], xlab = "Train concordance index", col = "#7FBF7F", cex.lab = 3, cex.axis = 3)
		dev.off()
		
		# Save a qq plot of the coxnet min train concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxnettrainminq",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		qqnorm(coxnet_results$t[,3], pch = 21, col = "blue", bg = "blue", cex = 2, xlab = "Theoretical quantiles", ylab = "Ordered values", cex.lab = 3, cex.axis = 3)
		qqline(coxnet_results$t[,3], col = "red")
		dev.off()
		
		# Save a histogram of the coxnet min oob concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxnetoobmin",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		hist(coxnet_results$t[,4], xlab = "OOB concordance index", col = "#7FBF7F", cex.lab = 3, cex.axis = 3)
		dev.off()
		
		# Save a qq plot of the coxnet min oob concordance indices
		png(paste("PRECOG_DMFS/Pictures/coxnetoobminq",n,".png",sep=""), width = 1500, height = 1000 )
		par(mar = c(10, 10, 10, 10), mgp = c(7, 2, 0))
		qqnorm(coxnet_results$t[,4], pch = 21, col = "blue", bg = "blue", cex = 2, xlab = "Theoretical quantiles", ylab = "Ordered values", cex.lab = 3, cex.axis = 3)
		qqline(coxnet_results$t[,4], col = "red")
		dev.off()

		# Print the final mean train/test concordance indices for 1se and min, rounded to 3 decimal places
		print("Coxnet train1se, test1se, trainmin, testmin")
		print(round(mean(coxnet_results$t[,1]), digits = 3))
		print(round(mean(coxnet_results$t[,2]), digits = 3))
		print(round(mean(coxnet_results$t[,3]), digits = 3))
		print(round(mean(coxnet_results$t[,4]), digits = 3))
	}
	
}

# The cox function computes a train/test concordance index for a certain bootstrapping sample
# Parameters:
# genes - which COVARIATES (not necessarily genes) should be used in the analysis
# data - the data object
# indices - selected indices (there may be repeated numbers in this) for the samples that the model should be trained on
cox <- function(genes, data, indices) {
	# The implementation of the boot function is that the first time, it will give the original indices (in sorted order)
	# This messes up the analysis, since there will be an empty out-of-bag set, causing the models to fail. Therefore, we only
	# Run the bootstrapping if the list of indices is UNSORTED.
	if (is.unsorted(indices)) {
		train <- data[indices,] # Get the training samples
		test <- data[-indices,] # Everything else is out-of-bag samples, and will be used for testing
		
		# Create a cox model using the specified covariates
		cox_model <- coxph(as.formula(paste("Surv(DMFS_Time, DMFS_Status) ~", paste(genes, collapse = "+"), sep = "")), data = train)
		
		# Get the train concordance index
		my_coxtrain <- summary(cox_model)$concordance[[1]]
		
		# Run predictions using the trained model on the oob data
		cox_predictions <- predict(cox_model, newdata = test)
		
		# Calculate the concordance index using a function from the survcomp package
		cox_cindex_test <- concordance.index(cox_predictions, surv.time = test$DMFS_Time, surv.event = test$DMFS_Status)
		
		# Extract the test concordance index
		my_coxtest <- cox_cindex_test$c.index
		
		# Return a vector of the train and test concordance indices
		return(c(my_coxtrain, my_coxtest))
	}
	# If instead, the list of indices is sorted, then we return an arbitrary value. In this case, 0.5 is returned for the train
	# and oob concordance indices, but it does not really matter for our purposes, because I do not include the scores for the 
	# sorted (original) list of indices when calculating the mean train and oob concordance indices. This can be any 2-length numeric vector
	# and the results of the analysis will not be affected.
	else {
		return(c(0.5, 0.5))
	}
}

# The rsf function computes a train/test concordance index for a certain bootstrapping sample
# Parameters:
# genes - which COVARIATES (not necessarily genes) should be used in the analysis
# data - the data object
# indices - selected indices (there may be repeated numbers in this) for the samples that the model should be trained on
rsf <- function(genes, data, indices) {
	# The implementation of the boot function is that the first time, it will give the original indices (in sorted order)
	# This messes up the analysis, since there will be an empty out-of-bag set, causing the models to fail. Therefore, we only
	# Run the bootstrapping if the list of indices is UNSORTED.
	if (is.unsorted(indices)) {
		train <- data[indices,] # Get the training samples
		test <- data[-indices,] # Everything else is out-of-bag samples, and will be used for testing
		
		# Create a rsf model using the specified covariates
		rsf_model <- rfsrc(as.formula(paste("Surv(DMFS_Time, DMFS_Status) ~", paste(genes, collapse = "+"), sep = "")), data = train)
		
		# Get the train concordance index (calculated as 1 minus the error rate)
		my_rsftrain <- 1-rsf_model$err.rate[rsf_model$ntree]
		
		# Run predictions using the trained model on the oob data
		rsf_predictions <- predict(rsf_model, newdata = test)
		
		# Get the test concordance index (calculated as 1 minus the error rate)
		my_rsftest <- 1-rsf_predictions$err.rate[rsf_predictions$ntree]
		
		# Return a vector of the train and test concordance indices
		return(c(my_rsftrain, my_rsftest))
	}
	else {
	# If instead, the list of indices is sorted, then we return an arbitrary value. In this case, 0.5 is returned for the train
	# and oob concordance indices, but it does not really matter for our purposes, because I do not include the scores for the 
	# sorted (original) list of indices when calculating the mean train and oob concordance indices. This can be any 2-length numeric vector
	# and the results of the analysis will not be affected.
		return(c(0.5, 0.5))
	}
}

# The coxnet function computes a train/test concordance index for a certain bootstrapping sample
# Parameters:
# genes - which COVARIATES (not necessarily genes) should be used in the analysis
# data - the data object
# indices - selected indices (there may be repeated numbers in this) for the samples that the model should be trained on
coxnet <- function(genes, data, indices) {
	# The implementation of the boot function is that the first time, it will give the original indices (in sorted order)
	# This messes up the analysis, since there will be an empty out-of-bag set, causing the models to fail. Therefore, we only
	# Run the bootstrapping if the list of indices is UNSORTED.
	if (is.unsorted(indices)) {
		train <- data[indices,] # Get the training samples
		test <- data[-indices,] # Everything else is out-of-bag samples, and will be used for testing
		
		# Create a string that contains the covariates with + symbols in between, and a - 1 at the end
		# Example: covar1+covar2+covar3 - 1
		# The coxnet model is created differently than cox and rsf. We first create an object called x,
		# which is a matrix of the covariates. The - 1 is used to remove the intercept.
		sub <- paste("~", paste(genes, collapse = "+"), "-1")
						
						
		# Create the x matrix
		x <- model.matrix(as.formula(sub), train)
		
		# Create the Surv object
		y <- Surv(train$DMFS_Time, train$DMFS_Status)
		# alpha = 0, so we only do Ridge regression (no Lasso). This is because Lasso sometimes fails and returns NA
		coxnet_model <- cv.glmnet(x, y, family = "cox", alpha = 0)
		
		# Create trainng and testing gene datasets
		newx_train <- as.matrix(train[, unlist(genes)])
		newx_test <- as.matrix(test[, unlist(genes)])
	
		# Run predictions. We test two values of lambda, the regularization parameter: lambda.min and lambda.1se.
		# lambda.min is the one that results in the lowest cross validation error. lambda.1se is the one which results
		# in a more parsimonious model and has a error one standard error away from the error in lamda.min
		# Run predictions on train and test for both lambda.1se and lambda.min so we can get the train and test concordance indices
		# in both cases.
		predictions_train_1se <- predict(coxnet_model, s = coxnet_model$lambda.1se, newx = newx_train)
		predictions_test_1se <- predict(coxnet_model, s = coxnet_model$lambda.1se, newx = newx_test)
		predictions_train_min <- predict(coxnet_model, s = coxnet_model$lambda.min, newx = newx_train)
		predictions_test_min <- predict(coxnet_model, s = coxnet_model$lambda.min, newx = newx_test)
		
		# Calculate concordance indices
		cindex_train_1se <- concordance.index(predictions_train_1se, surv.time = train$DMFS_Time, surv.event = train$DMFS_Status)
		cindex_test_1se <- concordance.index(predictions_test_1se, surv.time = test$DMFS_Time, surv.event = test$DMFS_Status)
		cindex_train_min <- concordance.index(predictions_train_min, surv.time = train$DMFS_Time, surv.event = train$DMFS_Status)
		cindex_test_min <- concordance.index(predictions_test_min, surv.time = test$DMFS_Time, surv.event = test$DMFS_Status)
		
		# Extract the actual concordance indices from the objects above
		my_coxnettrain1se <- cindex_train_1se$c.index
		my_coxnetoob1se <- cindex_test_1se$c.index
		my_coxnettrainmin <- cindex_train_min$c.index
		my_coxnetoobmin <- cindex_test_min$c.index

		return(c(my_coxnettrain1se, my_coxnetoob1se, my_coxnettrainmin, my_coxnetoobmin))
	}
	else {
	# If instead, the list of indices is sorted, then we return an arbitrary value. In this case, 0.5 is returned for the train
	# and oob concordance indices, but it does not really matter for our purposes, because I do not include the scores for the 
	# sorted (original) list of indices when calculating the mean train and oob concordance indices. This can be any 2-length numeric vector
	# and the results of the analysis will not be affected.
		return(c(0.5, 0.5, 0.5, 0.5))
	}
}


