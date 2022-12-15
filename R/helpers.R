####
#
# Relevant helper functions
#
####

#' @title Cross-validation for oscar-fitted model objects over k-range
#'
#' @description Create a cross-validation matrix with the chosen goodness metric with n-folds. Based on the goodness metric, one ought to pick optimal cardinality (parameter 'k').
#'
#' @param fit oscar-model object
#' @param fold Number of cross-validation folds, Default: 10
#' @param seed Random seed for reproducibility with NULL indicating that it is not set, Default: NULL
#' @param strata Should stratified cross-validation be used; separate values indicate balanced strata. Default: Unit vector, which will treat all observations equally.
#' @param verb Level of verbosity with higher integer giving more information, Default: 0
#' @param ... Additional parameters passed to oscar-function
#'
#' @return A matrix with goodness of fit over folds and k-values
#'
#' @details A k-fold cross-validation is run by mimicking the parameters contained in the original oscar S4-object. This requires the original data at slots @x and @y.
#*
#' @examples 
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_cv <- oscar.cv(fit, fold=10, seed=123)
#'   fit_cv
#' }
#' @rdname oscar.cv
#' @export 
#' @importFrom survival coxph Surv
oscar.cv <- function(
	# oscar-object
	fit,
	# k-fold
	fold = 10,
	# RNG seed (integer) that can be set for exact reproducibility
	seed = NULL,
	# Should some strata be balanced over the bins? By default no, should be a vector equal to the number of rows in x
        strata = rep(1, times=nrow(fit@x)),
	# Level of verbosity mainly for debugging
	verb = 0,
	# Additional parameters passed to oscar-function
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
	cvsets <- cv(fit@x, fold = fold, seed = seed, strata = strata)
	cvsets
	# Verbose
	if(verb>=2) print(str(cvsets))
	
	# K-steps
	ks <- 1:nrow(fit@k)
	# Loop over the CV folds, make a list of predictions and the real values
	cvs <- lapply(1:fold, FUN=function(z){
		if(verb>=1) cat(paste0("\nCV fold ", z, " of ", fold, "\n"))
		
		# Remove redundant columns but give out a warning
		x <- fit@x[cvsets$train[[z]],]
		k <- fit@k
		w <- fit@w
		if(any(apply(x, MARGIN=2, FUN=\(q){ all(q==unique(q)[1]) }))){
			omits <- which(apply(x, MARGIN=2, FUN=\(q){ all(q==unique(q)[1]) }))
			x <- x[,-omits]
			k <- k[-omits,-omits]
			w <- w[-omits]
			warning(paste("CV matrix 'x' contains redundant data columns in fold", z, " and these are removed in the cross-validation"))
		}
		# Constructing appropriate model object
		if(fit@family == "cox"){
			# Cox model is 2-column in y-response
			y <- survival::Surv(fit@y[cvsets$train[[z]],])
		}else if(fit@family %in% c("mse", "gaussian", "logistic")){
			# All other models have a y-vector
			y <- c(fit@y)[cvsets$train[[z]]]
		}else{
			stop(paste("Incorrect family-parameter fit@family:", fit@family))
		}
		# Exhaustively use same set of fitting parameters in CV fits as in original fit
		fittmp <- oscar::oscar(
			x = x, 
			y = y, 
			family = fit@family, 
			k = k, 
			w = w, 
			kmax = fit@kmax,
			verb = verb, 
			start = fit@start, 
			#rho = fit@rho,
			solver = fit@solver,
			in_selection = fit@in_selection, 
			percentage = fit@percentage, 
			...
		)
		
		if(verb>=1) print("CV fit, predicting...")
		# Perform predictions over all k-value fits
		# Model specificity in predictions (?)
		pred <- lapply(1:fit@kmax, FUN=function(ki){
			x <- fit@x[cvsets$test[[z]],,drop=FALSE]
			# MSE/Gaussian
			if(fit@family %in% c("mse", "gaussian")){
				oscar::predict(fit, type = "response", k = ki, newdata = x)
			# Logistic	
			}else if(fit@family %in% c("logistic")){
				oscar::predict(fit, type = "response", k = ki, newdata = x)
			# Cox
			}else if(fit@family %in% c("cox")){
				oscar::predict(fit, type = "response", k = ki, newdata = x)
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
		
		if(verb>=1) print(paste("Finished fold", z))
		
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
				#pROC::auc(response = z$true, predictor = c(q))
				# Less 'cat' output
				invisible(as.numeric(pROC::auc(pROC::roc(response=z$true, predictor=c(q), levels=c(0,1), direction="<"))))
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
#'
#' @description This model bootstraps the fitting of a given oscar object (re-fits the model for data that is equal in size but sampled with replacement). The output objects give insight into robustness of the oscar-coefficient path, as well as relative importance of model objects.
#'
#' @param fit oscar-model object
#' @param bootstrap Number of bootstrapped datasets, Default: 100
#' @param seed Random seed for reproducibility with NULL indicating that it is not set, Default: NULL
#' @param verb Level of verbosity with higher integer giving more information, Default: 0
#' @param ... Additional parameters passed to oscar-function
#'
#' @return 3-dimensional array with dimensions corresponding to k-steps, beta coefficients, and bootstrap runs
#'
#' @details The function provides a fail-safe try-catch in an event of non-convergence of the model fitting procedure. This may occur for example if a bootstrapped data matrix has a column consist of a single value only over all observations.
#'
#' @examples 
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_bs <- oscar.cv(fit, bootstrap = 20, seed = 123)
#'   fit_bs
#' }
#' @rdname oscar.bs
#' @export 
oscar.bs <- function(
	# oscar-object
	fit,
	# How many bootstrapped datasets to generate
	bootstrap = 100,
	# RNG seed (integer) that can be set for exact reproducibility
	seed = NULL,
	# Level of verbosity (<1 supresses everything in R; parameter also passed to Fortran subroutine)
	verb = 0,
	# Additional parameters passed to oscar-function
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
		# Wrap expression inside try for catching errors
		try({
			ftemp <- oscar::oscar(x = xtemp, y = ytemp, k = fit@k, w = fit@w, family = fit@family, kmax = fit@kmax, print = verb, start = fit@start, verb = verb, in_selection=fit@in_selection, percentage = fit@percentage,...)		
		})
		# Return successfully fitted model
		if(!inherits(ftemp,"try-error")){
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
	
#' @title Reformatting bootstrap output for cardinality k rows
#'
#' @description The function reformats bootstrapped runs to a single long data.frame, where all bootstrapped runs are covered along with the choices for the variables at each cardinality 'k'.
#'
#' @param bs Bootstrapped list from oscar.bs
#'
#' @return Reformatted data.frame
#'
#' @examples 
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_bs <- oscar.bs(fit, bootstrap = 20, seed = 123)
#'   ll <- oscar.bs.k(fit_bs)
#'   head(ll)
#'   tail(ll)
#' }
#' @rdname oscar.bs.k
#' @export
oscar.bs.k <- function(
	bs	# Bootstrapped list from oscar.bs
){
	# Omit entries with try-errors
	if(any(unlist(lapply(bs, FUN=function(x) { inherits(x, "try-error") })))){
		bs <- bs[-which(unlist(lapply(bs, FUN=function(x) { inherits(x, "try-error") })))]
		warning(paste("try-errors detected in some bootstrap runs; failed bootstrap count:", sum(unlist(lapply(bs, FUN=function(x) { inherits(x, "try-error") })))))
	}
	
	# Choices of variables as a function of k
	bs <- apply(bs, MARGIN=3, FUN=function(z){
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

#' @title Create a sparse matrix representation of betas as a function of k
#'
#' @description Variable estimates (rows) as a function of cardinality (k, columns). Since a model can drop out variables in favor of two better ones as k increases, this sparse representation helps visualize which variables are included at what cardinality.
#'
#' @param fit oscar-model object
#' @param kmax Create matrix until kmax-value; by default same as for fit object, but for high dimensional tasks one may wish to reduce this
#'
#' @return A sparse matrix of variables (rows) as a function of cardinality k (columns), where elements are the beta estimates.
#'
#' @details Uses sparseMatrix-class from Matrix-package
#'
#' @examples 
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   oscar.sparsify(fit, kmax=5)
#' }
#'
#' @rdname oscar.sparsify
#' @importFrom Matrix sparseMatrix
#' @export 
oscar.sparsify <- function(
	fit,
	kmax = fit@kmax
){
	if(!inherits(fit, "oscar")){
		stop("'fit' should be a fit oscar object")
	}
	
	# Order in which variables are first observed as non-zero
	#varorder <- unique(unlist(fit@kperk))
	# For models fit kits, a more sophisticated approach is required
	varorder <- unique(unlist(apply(fit@bperk, MARGIN=1, FUN=function(z) which(!z==0))))

	# Ordered beta matrix, order variables and transpose
	bkorder <- t(fit@bperk[,varorder])
	# Extract pairwise indices for non-zero elements
	nonzeroes <- which(!bkorder==0, arr.ind=TRUE)
	

	smat <- Matrix::sparseMatrix(
		i = nonzeroes[,1], # row indices of non-zero elements
		j = nonzeroes[,2], # col indices of non-zero elements
		x = c(bkorder[!bkorder==0]) # beta estimates at {i,j} running indices
	)
	dimnames(smat) <- dimnames(bkorder)	
	
	smat[,1:kmax]
}

#' @title Binary logical indicator matrix representation of an oscar object's coefficients (zero vs. non-zero, i.e. feature inclusion)
#'
#' @description Create a sparse matrix with binary indicator 1 indicating that a coefficient was non-zero, and value 0 (or . in sparse matrix) indicating that a coefficient was zero (i.e. feature not included)
#'
#' @param fit Fit oscar-model object
#' @param kmax Create matrix until kmax-value; by default same as for fit object, but for high dimensional tasks one may wish to reduce this
#'
#' @return A binary logical indicator matrix of variables (rows) as a function of cardinality k (columns), where elements are binary indicators for 1 as non-zero and 0 as zero.
#'
#' @details The matrix consists of TRUE/FALSE values, and is very similar to the oscar.sparsify, where the function provides estimate values in a sparse matrix format.
#'
#' @examples 
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   oscar.binarize(fit, kmax=5)
#' }
#'
#' @rdname oscar.binarize
#' @export 
oscar.binarize <- function(
	fit, # oscar model object
	kmax = fit@kmax # limit to kmax
){
	if(!inherits(fit, "oscar")){
		stop("'fit' should be a fit oscar object")
	}
	
	# Full sparse matrix representation
	binmat <- apply(oscar::oscar.sparsify(fit), MARGIN=2, FUN=function(z) { ifelse(z==0, FALSE, TRUE) })
	
	binmat[,1:kmax]
}

#' @title Retrieve a set of pareto-optimal points for an oscar-model based on model goodness-of-fit or cross-validation
#'
#' @description This function retrieves the set of pareto optimal points for an oscar model fit in n-proportional time as cardinality axis is readily sorted. It is advisable to optimize model generalization (via cross-validation) rather than mere goodness-of-fit.
#'
#' @param fit Fit oscar S4-object
#' @param cv A cross-validation matrix as produced by oscar.cv; if CV is not provided, then goodness-of-fit from fit object itself is used rather than cross-validation generalization metric
#' @param xval The x-axis to construct pareto front based on; by default 'cost' vector for features/kits, can also be 'cardinality'/'k'
#' @param weak If weak pareto-optimality is allowed; by default FALSE.
#' @param summarize Function that summarizes over cross-validation folds; by default, this is the mean over the k-folds.
#' 
#' @return A data.frame containing points and indices at which pareto optimal points exist
#'
#' @examples 
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_cv <- oscar.cv(fit, fold=10)
#'   oscar.pareto(fit, cv=fit_cv)
#' }
#'
#' @rdname oscar.pareto
#' @export
oscar.pareto <- function(
	fit,
	cv,
	xval = "cost",
	weak = FALSE,
	summarize = mean
){
	# 
	if(xval == "cost"){
		xs <- unlist(lapply(1:fit@kmax, FUN=function(k) { oscar::cost(fit, k) }))
		ord <- order(xs)
	}else if(xval %in% c("cardinality", "k")){
		xs <- 1:fit@kmax
		ord <- 1:fit@kmax
	}else{
		stop(paste("Invalid xval value:",xval))
	}
	# x-axis is always the cardinality k
	# y-axis is either the goodness of fit or cross-validation performance. Latter is preferred.
	# No CV provided
	if(missing(cv)){
		ys <- fit@goodness
	# CV provided
	}else{
		# Apply summarization function over columns (the cardinality values), as rows are the folds in CV
		ys <- apply(cv, MARGIN=2, FUN=summarize)
	}
	# For certain metrics, test in the other direction
	if(fit@metric %in% c("mse")){
		ys <- -ys
	}
	# Loop over the x-values in ascending order; first one is always a pareto point, so we'll store it
	# Assuming x-axis is cardinality, we want lower cardinality models i.e. conservative in terms of their complexity
	# xs are actually cardinalities starting from 1, so we can just use integer indices
	df <- data.frame(y = ys, x = xs, ord = ord)
	df <- df[ord,]	
	ymax <- df$y[1]
	paretos <- 1
	for(i in 1:nrow(df)){
		# Always testing for higher values, as the metrics in other direction were already flipped
		if(df$y[i] > ymax){
			paretos <- c(paretos, i)
			ymax <- df$y[i]
		}
	}
	# Result for {x,y,pareto yes/no}
	res <- data.frame(k = df$x, y = df$y, paretofront = FALSE, ord = ord)
	# Replace y with the actual model metric
	names(res)[1:2] <- c(xval, fit@metric)
	res[paretos,"paretofront"] <- TRUE
	# If weak pareto-optimal points should be removed as per default
	if(!weak){
		res <- do.call("rbind", by(res, INDICES=res[,1], FUN=function(x) { x[which(!x[,2]==max(x[,2])),3]<-FALSE; x }))
		rownames(res) <- NULL
	}
	# Return to original point ordering
	res <- res[order(res$ord),]
	res
}

#' @title Return total cost of model fits if the cost is not included in the oscar object
#'
#' @description If at least one measurement from a kit is included in the model, the kit cost is added.
#'
#' @param object Fit oscar S4-object 
#'
#' @return A vector for numeric values of total kit costs at different cardinalities.
#'
#' @examples 
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   oscar.cost.after(fit)
#' }
#'
#' @rdname oscar.cost.after
#' @export
oscar.cost.after <- function(object){

	# The values are readily stored in oscar S4-objects
	kit.matrix <- object@k
	cost.vector <- object@w

  # Assume that the features are in the same order in x-matrix of oscar object and kit matrix
  # Assume that the kits are in the same order in kit.matrix and cost.vector
  
  #Sanity checks
  if(nrow(kit.matrix)!=length(cost.vector)){
    stop("Number of kits in the kit matrix is different from the length of cost vector.")
  }
  if(ncol(kit.matrix)!=ncol(object@x)){
    stop("Number of predictors in the kit.matrix differs from the number of predictors in x-matrix of oscar object.")
  }
  
  costs <- c()
  for(i in 1:object@kmax){ # Go through each cardinality
    cost.tmp <-0
    nzero <- which(object@bperk[i,]!=0)

    for(j in 1:ncol(kit.matrix)){ #Go through kits
      if(any(kit.matrix[j,nzero]!=0)){
        cost.tmp <- cost.tmp+cost.vector[j] # Add kit price if any feature is incl.
      }
    }
    costs<-c(costs,as.numeric(cost.tmp))
  }
  return(costs)
}


		
