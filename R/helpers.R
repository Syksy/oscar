####
#
# Relevant helper functions
#
####

#' @title Cross-validation for casso-fitted model objects over k-range
#' @description FUNCTION_DESCRIPTION
#' @param fit PARAM_DESCRIPTION
#' @param fold PARAM_DESCRIPTION, Default: 10
#' @param seed PARAM_DESCRIPTION, Default: NULL
#' @param verb PARAM_DESCRIPTION, Default: 0
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[casso]{character(0)}}
#'  \code{\link[stats]{predict.glm}}
#'  \code{\link[survival]{predict.coxph}},\code{\link[survival]{coxph}},\code{\link[survival]{Surv}}
#' @rdname cv.casso
#' @export 
#' @importFrom casso casso
#' @importFrom stats predict.glm
#' @importFrom survival predict.coxph coxph Surv
cv.casso <- function(
	# casso-object
	fit,
	# k-fold
	fold = 10,
	# RNG seed (integer) that can be set for exact reproducibility
	seed = NULL,
	# Level of verbosity mainly for debugging
	verb = 0,
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
			fittmp <- casso::casso(x=fit@x[cvsets$train[[z]],], y=fit@y[cvsets$train[[z]],], family=fit@family, k=fit@k, w=fit@w, verb=verb)
		}else if(fit@family %in% c("mse", "gaussian", "logistic")){
			# All other models have a y-vector
			fittmp <- casso::casso(x=fit@x[cvsets$train[[z]],], y=c(fit@y)[cvsets$train[[z]]], family=fit@family, k=fit@k, w=fit@w, verb=verb)
		}else{
			stop(paste("Incorrect family-parameter fit@family:", fit@family))
		}
		if(verb>=1) print("CV fit, predicting...")
		# Perform predictions over all k-value fits
		# Model specificity in predictions (?)
		pred <- lapply(fittmp@fits, FUN=function(f){
			#if(verb>=2) print(f)
			x <- as.data.frame(fit@x[cvsets$test[[z]],])
			colnames(x) <- colnames(fit@x)
			# MSE/Gaussian
			if(fit@family %in% c("mse", "gaussian")){
				as.vector(unlist(stats::predict.glm(f, type="response", newdata=x)))
			# Logistic	
			}else if(fit@family %in% c("logistic")){
				as.vector(unlist(stats::predict.glm(f, type="response", newdata=x)))
			# Cox
			}else if(fit@family %in% c("cox")){
				as.vector(unlist(survival:::predict.coxph(f, type="risk", newdata=x)))
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
			# Logistic placeholder; correct classification rate
			}else if(fit@family %in% c("logistic")){
				sum(as.integer(q>0.5)==z$true)/length(q)
			# Cox proportional hazards model; concordance-index
			}else if(fit@family %in% c("cox")){
				# NOTE: Assuming first colmn is the survival time, second is the observed event status
				survival::coxph(survival::Surv(time = z$true[,1], event = z$true[,2]) ~ q)$concordance["concordance"]
			}
		})
	})

	#if(verb>=2){
	#	print("cvs prior to wrapping up")	
	#	print(class(cvs))
	#	print(cvs)
	#}
		
	# Rows: cv-folds, cols: k-values
	cvs <- do.call("rbind", lapply(cvs, FUN=function(z) do.call("c", z)))
	rownames(cvs) <- paste("cv_", 1:nrow(cvs), sep="")
	# Return cv results
	cvs
}

#' @title Bootstrapping for casso-fitted model objects
#' @description FUNCTION_DESCRIPTION
#' @param fit PARAM_DESCRIPTION
#' @param bootstrap PARAM_DESCRIPTION, Default: 100
#' @param seed PARAM_DESCRIPTION, Default: NULL
#' @param verb PARAM_DESCRIPTION, Default: 0
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[casso]{character(0)}}
#' @rdname bs.casso
#' @export 
#' @importFrom casso casso
bs.casso <- function(
	# casso-object
	fit,
	# How many bootstrapped datasets to generate
	bootstrap = 100,
	# RNG seed (integer) that can be set for exact reproducibility
	seed = NULL,
	# Level of verbosity (<1 supresses everything in R; parameter also passed to Fortran subroutine)
	verb = 0,
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
			ftemp <- casso::casso(x = xtemp, y = ytemp, k = fit@k, w = fit@w, family = fit@family, kmax = fit@kmax, print = verb, start = fit@start, verb = verb)		
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
	

#' @title Modified glm.control allowing 0 in maxit
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @param weights PARAM_DESCRIPTION, Default: rep.int(1, nobs)
#' @param start PARAM_DESCRIPTION, Default: NULL
#' @param etastart PARAM_DESCRIPTION, Default: NULL
#' @param mustart PARAM_DESCRIPTION, Default: NULL
#' @param offset PARAM_DESCRIPTION, Default: rep.int(0, nobs)
#' @param family PARAM_DESCRIPTION, Default: gaussian()
#' @param control PARAM_DESCRIPTION, Default: list()
#' @param intercept PARAM_DESCRIPTION, Default: TRUE
#' @param singular.ok PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{character(0)}}
#' @rdname casso:::.glm.fit.mod
#' @importFrom stats C_Cdqrls
.glm.control.mod <- function (epsilon = 1e-08, maxit = 25, trace = FALSE) 
{
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of 'epsilon' must be > 0")
    ## MODIFIED
    if (!is.numeric(maxit) || maxit < 0) 
        stop("maximum number of iterations must be >= 0")
    ## MODIFIED ENDS
    list(epsilon = epsilon, maxit = maxit, trace = trace)
}
	
#' Modified function from stats::glm.fit, allowing 0 maxit in control
.glm.fit.mod <- function (x, y, weights = rep.int(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep.int(0, nobs), family = gaussian(), 
    control = list(), intercept = TRUE, singular.ok = TRUE) 
{
    ## MODIFIED
    control <- do.call(".glm.control.mod", control)
    ## MODIFIED ENDS
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x)) 
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("invalid linear predictor values in empty model", 
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("invalid fitted means in empty model", 
                call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- sqrt((weights * mu.eta(eta)^2)/variance(mu))
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep_len(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                  nvars, paste(deparse(xnames), collapse = ", ")), 
                  domain = NA)
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1L) 
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("cannot find valid starting values: please specify some", 
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        ## MODIFIED
        iter <- 0L
        #for (iter in 1L:control$maxit) {
        ## MODIFIED ENDS
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (anyNA(varmu)) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning(gettextf("no observations informative at iteration %d", 
                  iter), domain = NA)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            ## MODIFIED
            fit <- .Call(stats:::C_Cdqrls, x[good, , drop = FALSE] * 
                w, z * w, min(1e-07, control$epsilon/1000), check = FALSE)
            ## MODIFIED ENDS    
                
            ## MODIFIED
            fit$coefficients <- start
            ## MODIFIED ENDS
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d", 
                  iter), domain = NA)
                ## MODIFIED
                #break
                ## MODIFIED ENDS
            }
            if (nobs < fit$rank) 
                stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation", 
                  "X matrix has rank %d, but only %d observations"), 
                  fit$rank, nobs), domain = NA)
            if (!singular.ok && fit$rank < nvars) 
                stop("singular fit encountered")
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance = ", dev, " Iterations - ", 
                  iter, "\n", sep = "")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance = ", dev, 
                    "\n", sep = "")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  ## MODIFIED
                  #start <- (start + coefold)/2
                  ## MODIFIED ENDS
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance = ", dev, 
                    "\n", sep = "")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                ## MODIFIED
                #break
                ## MODIFIED ENDS
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        ## MODIFIED
        #}
        ## MODIFIED ENDS
        if (!conv) 
            warning("glm.fit: algorithm did not converge", 
                call. = FALSE)
        if (boundary) 
            warning("glm.fit: algorithm stopped at boundary value", 
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("glm.fit: fitted probabilities numerically 0 or 1 occurred", 
                  call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps)) 
                warning("glm.fit: fitted rates numerically 0 occurred", 
                  call. = FALSE)
        }
        if (fit$rank < nvars) 
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY) 
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", 
            sum(good) - fit$rank))
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
        rank = rank, qr = if (!EMPTY) structure(fit[c("qr", 
            "rank", "qraux", "pivot", "tol")], 
            class = "qr"), family = family, linear.predictors = eta, 
        deviance = dev, aic = aic.model, null.deviance = nulldev, 
        iter = iter, weights = wt, prior.weights = weights, df.residual = resdf, 
        df.null = nulldf, y = y, converged = conv, boundary = boundary)
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
	bs	# Bootstrapped list from bs.casso
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
