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
	family <- fit@family
	# K-steps
	ks <- 1:nrow(fit@k)
	# Loop over the CV folds, make a list of predictions and the real values
	cvs <- lapply(1:fold, FUN=function(z){
		if(verb>=1) print(paste("CV fold", z))
		fittmp <- casso::casso(x=fit@x[cvsets$train[[z]],], y=fit@y[cvsets$train[[z]],], family=family, k=fit@k)
		# Perform predictions over all k-value fits
		pred <- lapply(fittmp@fits, FUN=function(f){
			predict(f, type="response", newdata=fit@x[cvsets$test[[z]],])
		})
		# Check predictions vs. real y-values
		list(
			# Predicted values
			pred = pred, 
			# True values
			true = fit@y[cvsets$test[[z]],]
		)
	})
	# Return cv results
	cvs
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
	verb = 3,
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
		ftemp <- casso(x = xtemp, y = ytemp, k = fit@k, w = fit@w, family = fit@family, print = verb, start = fit@start, verb = verb)
		
		ftemp@bperk # Return bootstrapped beta per ks
	})
	# Return bootstrapped beta per ks
	bperks
}
	

#' Modified glm.control allowing 0 in maxit
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