## Datasets

# Implementation for Greedy FS

#' @title Greedy forward selection for features in Cox PH models using Condordance index
#'
#' @description x Input data matrix
#'
#' @param x Input data matrix
#' @param y Input Surv-response (survival::Surv)
#' @param maxk Maximum number of variables to add (Default: 30)
#'
#' @examples
#'
#'
#' @noRd
#' @keywords internal
greedyfw <- function(
	x,
	y,
	maxk = 50,
	verb = 1, # Level of verbosity
	fold = 10, # Cross-validation folds
	seed, # If there should be a set seed
	runCV = TRUE, # Should CV be run
	...
){
	# Extract names of the variables/columns while replacing the problematic '-', '/' etc with '.' as per R conventions
	vars <- gsub("-|/", ".", colnames(x))
	colnames(x) <- gsub("-|/", ".", colnames(x))

	# Variables selected thus far
	selected <- c()
	cindex	<- c()
	if(runCV) cvs <- list()
	
	# Iterate over candidate variables until desired var count 
	for(p in 1:min(maxk, ncol(x))){
		if(verb>=1) cat(paste0("\np: ", p, " / ", min(maxk, ncol(x)), "\n"))
		# Loop over remaining variables and selected one which adds most to C-index
		cs <- unlist(lapply(vars, FUN=\(v){
			tryCatch(
				expr = {
					# Construct appropriate data matrix
					feats <- as.data.frame(x[,c(selected, v)])
					colnames(feats) <- c(selected, v)
					f <- as.formula(paste0("y ~ ", paste(c(selected, v), collapse="+")))
					if(verb>=2){
						print(f)
						print(head(feats))
					}
					fit <- survival::coxph(f, data = feats)
					ci <- fit[["concordance"]]["concordance"]
					ci
				},
				# Error when trying to extract fit/c-index
				error = function(e){
					if(verb>=2) print(paste("Error occurred with variable", v, "; returning ci=0.5"))
					ci <- 0.5
					ci
				}
			)
		}))
		w <- which.max(cs)
		# Add C-index maximizing variable into the selected variables, remove it from candidate variables
		selected <- c(selected, vars[w])
		vars <- vars[-w]
		# Store concordance of original fit
		cindex <- c(cindex, cs[w])

		if(verb>=2){
			print("Selected:")
			print(selected)
			print("C-indices:")
			print(cindex)
		}
		
		# Run CV for thus far selected variables
		if(runCV) cvs[[length(cvs)+1]] <- greedyfw.cv(x = x[,selected,drop=FALSE], y = y, fold = fold, seed = seed)

		if(verb>=1 & runCV){
			print("CV-folds:")
			print(cvs)
		}
	}
	
	# Return either just fits or fits together with CVs	
	if(runCV){
		cvs <- do.call("cbind", cvs)
		rownames(cvs) <- paste0("fold", 1:fold)
		colnames(cvs) <- paste0("k", 1:ncol(cvs))
		list(selected = selected, cindex = cindex, cvs = cvs)
	}else{
		list(selected = selected, cindex = cindex)
	}
}

#' Cross-validation for greedy forward selection with Cox PH and c-index
greedyfw.cv <- function(
	x,
	y,
	fold = 10,
	verb = TRUE,
	seed,
	...
){
	# Internal cv-sample allocation function, modified from the ePCR-package
	cv <- function(
		# Original x data frame
		x,
		# Number of CV-folds
		fold = fold,
		# Should some strata be balanced over the bins? By default no, should be a vector equal to the number of rows in x
		strata = rep(1, times=nrow(x)),
		# Should data be randomized when creating cv-folds in addition to allocating to bins
		shuffle = TRUE,
		# Seed number for reproducibility
		# If NULL then seed is not set
		seed = seed
	){
		# Seed number
		if(!is.null(seed)) set.seed(seed)
	
		# Allocate folds
		uniqs <- unique(strata)
		folds <- rep(1:fold, times=ceiling(nrow(x)/fold))[1:nrow(x)]
		ifelse(shuffle,
			ord <- unlist(lapply(uniqs, FUN=function(z) sample(which(strata==z)))),
			ord <- unlist(lapply(uniqs, FUN=function(z) which(strata==z))))
		
		whichfold <- vector(length=nrow(x))
		whichfold[ord] <- folds
	
		# Construct test and train sets based on the folds
		test <- lapply(1:fold, FUN=function(z) which(z==whichfold))
		train <- lapply(1:fold, FUN=function(z) do.call("c", test[-z]))
		
		# Return cv sample indices per each fold for train and test sets respectively
		list(train = train, test = test)
	}

	# Generate cv-folds
	cvsets <- cv(x = x, fold = fold, seed = seed)
	
	if(verb) print("CV sets generated")
	
	# Loop over the CV folds, make a list of predictions and the real values
	cvs <- lapply(1:fold, FUN=function(z){
		if(verb>=1) cat(paste0("\nCV fold ", z, " of ", fold, "\n"))

		# Train-test split
		x_train <- as.data.frame(x[cvsets$train[[z]],,drop=FALSE])
		y_train <- y[cvsets$train[[z]]]
		x_test <- as.data.frame(x[cvsets$test[[z]],,drop=FALSE])
		y_test <- y[cvsets$test[[z]]]
		
		colnames(x_train) <- colnames(x)
		colnames(x_test) <- colnames(x)

		# Fit using train data and predict for test		
		fit_train <- survival::coxph(y_train ~ ., data = x_train)		
		preds_test <- predict(fit_train, newdata = x_test)
		
		# Calculate test-prediction c-index
		tryCatch(
			expr = {
				ci_pred <- survival::coxph(y_test ~ preds_test)[["concordance"]]["concordance"]
				ci_pred
			}, 
			e = {
				NA_real_
			}
		)
	})	
	do.call("c", cvs)
}
