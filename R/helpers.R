####
#
# Relevant helper functions
#
####

#' @title Cross-validation for oscar-fitted model objects over k-range
#' @description Create a cross-validation matrix with the chosen goodness metric with n-folds. Based on the goodness metric, one ought to pick optimal cardinality (parameter 'k').
#' @param fit oscar-model object
#' @param fold Number of cross-validation folds, Default: 10
#' @param seed Random seed for reproducibility with NULL indicating that it is not set, Default: NULL
#' @param verb Level of verbosity with higher integer giving more information, Default: 0
#' @return A matrix with goodness of fit over folds and k-values
#' @details TODO
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{predict.glm}}
#'  \code{\link[survival]{predict.coxph}},\code{\link[survival]{coxph}},\code{\link[survival]{Surv}}
#' @rdname cv
#' @export 
#' @importFrom stats predict.glm
#' @importFrom survival coxph Surv
cv.oscar <- function(
	# oscar-object
	fit,
	# k-fold
	fold = 10,
	# RNG seed (integer) that can be set for exact reproducibility
	seed = NULL,
	# Level of verbosity mainly for debugging
	verb = 0
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
	cvsets <- cv(fit@x, fold = fold, seed = seed)
	cvsets
	# Verbose
	if(verb>=1) print(str(cvsets))
	
	# Fit family to use
	#family <- fit@family
	# K-steps
	ks <- 1:nrow(fit@k)
	# Loop over the CV folds, make a list of predictions and the real values
	cvs <- lapply(1:fold, FUN=function(z){
		if(verb>=1) print(paste("CV fold", z))
		# Constructing appropriate model object
		if(fit@family == "cox"){
			# Cox model is 2-column in y-response
			fittmp <- oscar::oscar(x=fit@x[cvsets$train[[z]],], y=fit@y[cvsets$train[[z]],], family=fit@family, k=fit@k, w=fit@w, verb=verb, start=fit@start)
		}else if(fit@family %in% c("mse", "gaussian", "logistic")){
			# All other models have a y-vector
			fittmp <- oscar::oscar(x=fit@x[cvsets$train[[z]],], y=c(fit@y)[cvsets$train[[z]]], family=fit@family, k=fit@k, w=fit@w, verb=verb, start=fit@start)
		}else{
			stop(paste("Incorrect family-parameter fit@family:", fit@family))
		}
		if(verb>=1) print("CV fit, predicting...")
		# Perform predictions over all k-value fits
		# Model specificity in predictions (?)
		pred <- lapply(1:fit@kmax, FUN=function(ki){
			#if(verb>=2) print(f)
			x <- fit@x[cvsets$test[[z]],drop=FALSE]
			#colnames(x) <- colnames(fit@x)
			#x <- as.matrix(fit@x[cvsets$test[[z]],])
			# MSE/Gaussian
			if(fit@family %in% c("mse", "gaussian")){
				predict.oscar(fit, type = "response", k = ki, newdata = x)
				#as.vector(unlist(stats::predict.glm(f, type="response", newdata=x)))
				
			# Logistic	
			}else if(fit@family %in% c("logistic")){
				predict.oscar(fit, type = "response", k = ki, newdata = x)
				#as.vector(unlist(stats::predict.glm(f, type="response", newdata=x)))
			# Cox
			}else if(fit@family %in% c("cox")){
				predict.oscar(fit, type = "response", k = ki, newdata = x)
				#as.vector(unlist(survival:::predict.coxph(f, type="risk", newdata=x)))
			}
		})
		# True values; vectorization if not Cox ph model
		if(fit@family=="cox"){
			true <- fit@y[cvsets$test[[z]],]
		}else if(fit@family %in% c("mse", "gaussian", "logistic")){
			true <- c(fit@y)[cvsets$test[[z]]]
		}else{
			stop(paste("Incorrect family-parameter fit@family:", fit@family))
		}
		# Return predictions vs. real y-values
		list(
			# Predicted values
			pred = pred, 
			# True values; vectorization if not Cox ph model
			true = true
		)
	})
	
	if(verb>=2){
		print("cvs prior to goodness measure")
		print(class(cvs))
		print(cvs)
	}
	
	# Construct goodness as a function of k
	# Loop over cv folds
	# $pred slot holds k-iterations
	# $true slot holds the real answer
	cvs <- lapply(cvs, FUN=function(z){
		if(verb>=2){
			print("true:")
			print(z$true)
			print("z$pred:")
			print(z$pred)
		}
		lapply(z$pred, FUN=function(q){
			# MSE/Gaussian, mean squared error
			if(fit@family %in% c("mse", "gaussian")){
				mean((q-z$true)^2)
			# Logistic; ROC-AUC by default
			}else if(fit@family %in% c("logistic") & fit@metric == "auc"){
				pROC::auc(z$true ~ q)
				#sum(as.integer(q>0.5)==z$true)/length(q)
			# Logistic; correct classification rate if metric desired is accuracy
			}else if(fit@family %in% c("logistic") & fit@metric == "accuracy"){
				sum(as.integer(q>0.5)==z$true)/length(q)
			# Cox proportional hazards model; concordance-index
			}else if(fit@family %in% c("cox")){
				# NOTE: Assuming first colmn is the survival time, second is the observed event status
				survival::coxph(survival::Surv(time = z$true[,1], event = z$true[,2]) ~ q)$concordance["concordance"]
			}else{
				stop(paste("Invalid 'family'", fit@family, "or 'metric':", fit@metric))
			}
		})
	})
	
	# Rows: cv-folds, cols: k-values
	cvs <- do.call("rbind", lapply(cvs, FUN=function(z) do.call("c", z)))
	rownames(cvs) <- paste("cv_", 1:nrow(cvs), sep="")
	# Return cv results
	cvs
}

#' @title Bootstrapping for oscar-fitted model objects
#' @description This model bootstraps the fitting of a given oscar object (re-fits the model for data that is equal in size but sampled with replacement). The output objects give insight into robustness of the oscar-coefficient path, as well as relative importance of model objects.
#' @param fit oscar-model object
#' @param bootstrap Number of bootstrapped datasets, Default: 100
#' @param seed Random seed for reproducibility with NULL indicating that it is not set, Default: NULL
#' @param verb Level of verbosity with higher integer giving more information, Default: 0
#' @return 3-dimensional array with dimensions corresponding to k-steps, beta coefficients, and bootstrap runs
#' @details TODO
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname bs
#' @export 
bs.oscar <- function(
	# oscar-object
	fit,
	# How many bootstrapped datasets to generate
	bootstrap = 100,
	# RNG seed (integer) that can be set for exact reproducibility
	seed = NULL,
	# Level of verbosity (<1 supresses everything in R; parameter also passed to Fortran subroutine)
	verb = 0
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
		# Wrap expression inside try for catching errors
		try({
			ftemp <- oscar::oscar(x = xtemp, y = ytemp, k = fit@k, w = fit@w, family = fit@family, kmax = fit@kmax, print = verb, start = fit@start, verb = verb)		
		})
		# Return successfully fitted model
		if(!class(ftemp)=="try-error"){
			ftemp@bperk # Return bootstrapped beta per ks
		}else{
			NA # Model fitting issues, return NA-bperk
		}
	})
	# Extract bootstrapped beta coefficient values and create a 3-dim array of regular size
	dat <- array(as.numeric(unlist(lapply(bperks, FUN=c))), 
		dim=c(fit@kmax, ncol(bperks[[1]]), bootstrap), 
		dimnames=list(paste0("k_",1:fit@kmax), colnames(bperks[[1]]), paste0("bs_",1:bootstrap))
	)	
	# Return 3-dim array with first dim as k, second as beta coef, third as bootstrap runs
	dat
}
	
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param bs PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname bs.k
#' @export
bs.k <- function(
	bs	# Bootstrapped list from bs.oscar
){
	# Omit entries with try-errors
	if(any(unlist(lapply(bs, FUN=class))=="try-error")){
		bs <- bs[-which(unlist(lapply(bs, FUN=class))=="try-error")]
	}
	
	# Choices of variables as a function of k
	bs <- lapply(bs, FUN=function(z){
		apply(z, MARGIN=1, FUN=function(q){
			colnames(z)[which(!q==0)]
		})
	})
	
	# Concatenate all together into a data.frame
	bs <- as.data.frame(
		do.call("rbind", 
			lapply(1:length(bs), FUN=function(z){
				do.call("rbind", lapply(1:length(bs[[z]]), FUN=function(q){
					cbind(
						kth = rep(paste("k_", q, sep=""), times=length(bs[[z]][[q]])), # k:th var
						bsn = z, # z:th bootstrap run
						var = bs[[z]][[q]] # The variables in order
					)
				}))
			})
		)
	)
	
	bs
}
