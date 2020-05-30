####
#
# Relevant helper functions
#
####

#' Cross-validation for casso-fitted model objects over k-range
#'
#' TODO
cv.casso <- function(
	# casso-object
	fit,
	# k-fold
	fold = 10,
	# RNG seed (integer) that can be set for exact reproducibility
	seed = NULL,
	...
){
	# Internal cv-sample allocation function, modified from the ePCR-package
	cv <- function(
		# Original x data frame
		x,
		# Number of CV-folds
		fold = 10,
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
	cvsets <- cv(fit@x, fold = fold, seed = seed)
	cvsets
}



#' Bootstrapping for casso-fitted model objects
#'
#' TODO
bs.casso <- function(
	# casso-object
	fit,
	# How many bootstrapped datasets to generate
	bootstrap = 100,
	# RNG seed (integer) that can be set for exact reproducibility
	seed = NULL,
	# Level of verbosity (<1 supresses everything in R; parameter also passed to Fortran subroutine)
	verb = 1,
	...
){
	# Seed number for RNG reproducibility
	if(!is.null(seed)) set.seed(seed)
	# Generate and fit bootstrapped datasets
	bperks <- lapply(1:bootstrap, FUN=function(i){
		if(verb>=1) print(paste("Bootstrap iteration", i))
		# Sampling with replacement from rows
		samps <- sample(1:nrow(fit@x), replace=TRUE)
		xtemp <- fit@x[samps,]
		ytemp <- fit@y[samps,]
		ftemp <- casso(x = xtemp, y = ytemp, k = fit@k, w = fit@w, family = fit@family, print = verb, start = fit@start)
		
		ftemp@bperk # Return bootstrapped beta per ks
	})
	# Return bootstrapped beta per ks
	bperks
}
	
	