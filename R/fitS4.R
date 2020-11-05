####
#
# S4-class OO for saving models and relevant parameters
# DBDC-fit kit-structured L0-penalized model fits
#
####

#' S4-class for casso
#'
#'
#'
#' @export
setClass("casso", # abbreviation
	representation(
		# Set class representations for the S4-class @-slots
		## Results from Fortran
		bperk = "matrix",	# Rows: different kit k-values, columns: beta coefficients  in model fit
		fperk = "numeric",	# Target function values at each k-value
		kperk = "list",		# Which kits were picked at which k-step (named indices vector)
		cperk = "numeric",	# Cumulative kit costs per each k-value
		#cputime = "numeric",	# Total CPU time for computations
		rho = "numeric",	# Vector of utilized rho-values in the penalization
		start = "numeric",	# Method of generating starting values
		## Provided data, saved in R
		x = "matrix",		# Original data matrix for fitting
		y = "matrix",		# Response vector (two-column Surv-object for Cox modelling)
		k = "matrix",		# Kit indicator matrix
		w = "numeric",		# Kit value (cost ~ weights) vector, length ought to be same as nrow(k)
		## Additional
		family = "character",	# Utilized model family, in early versions only 'cox' supported (to be added: normal/gaussian, logistic, ...
		goodness = "numeric",	# Goodness metric at each k based on model fit; dependent on 'family': cox is concordance-index, ...
		fits = "list",		# Pre-fitted model objects, where each model fit only allows the optimal subset of features to be incorporated into the model; length should be equal to k sequence
		info = "character",	# Additional error messages, warnings or such reported e.g. during model fitting
		kmax = "integer",	# Number of maximum k tested
		metric = "character"	# Name of the goodness-of-fit metric used
	),	
	prototype(
		# Prototype base model object
		## Results from Fortran
		bperk = matrix(NA, nrow=0, ncol=0),
		fperk = NA_real_,
		kperk = list(),
		cperk = NA_real_,
		rho = NA_real_,
		start = NA_real_,
		## Provided data, saved in R
		x = matrix(NA, nrow=0, ncol=0),
		y = matrix(NA, nrow=0, ncol=0),
		k = matrix(NA, nrow=0, ncol=0),
		w = NA_real_,
		## Additional
		family = NA_character_,
		goodness = NA_real_,
		fits = list(),
		info = NA_character_,
		kmax = NA_integer_,
		metric = NA_character_
	),
	# Function for testing whether all S4-slots are legitimate for the S4 object 
	# TODO
	validity = function(
		object
	){
		# Return TRUE only if all slots fulfill the validity criteria 
		# (could be a list of if-statements that return FALSE prior to this if discrepancies are detected)
		TRUE
	}
)


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x Data matrix 'x'
#' @param y Response vector/two-column matrix 'y' (see: family); number of rows equal to nrow(x)
#' @param k Integer (0/1) kit indicator matrix; number of columns equal to ncol(x)
#' @param w Kit cost weight vector w of length nrow(k)
#' @param family Model family, should be one of: 'cox', 'mse'/'gaussian', or 'logistic, Default: 'gaussian'
#' @param print Level of verbosity in Fortran (may not be visible on all terminals); should be an integer between {range, range}, Default: 3
#' @param start Starting point generation method, see vignettes for details; should be an integer between {range,range}, Default: 2
#' @param verb Integer with additional integer values giving verbal feedback, Default: 1
#' @param kmax Maximum k step tested, by default all k are tested from k to maximum dimensionality
#' @return Fitted casso-object
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- casso(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit
#'  }
#' }
#' @seealso 
#'  \code{\link[survival]{coxph}},\code{\link[survival]{coxph.control}}
#'  \code{\link[stats]{glm}}
#'  \code{\link[casso]{character(0)}}
#' @rdname casso
#' @export 
#' @importFrom survival coxph coxph.control
#' @importFrom stats glm
casso <- function(
	# Data matrix x
	x, 
	# Response vector y (or 2-column survival response matrix)
	y, 
	# Kit indicator matrix k
	k, 
	# Kit cost weight vector w
	w, 
	# Model family (defines the likelihood function)
	family = "gaussian",
	## Tuning parameters
	print=3,# Level of verbosity (-1 for tidy output, 3 for debugging level verbosity)
	start=2,# Deterministic start point 
	verb=1, # Level of R verbosity (1 = standard, 2 = debug level, 3 = excessive debugging, 0<= none)
	kmax    # Maximum tested k-values
){
	# TODO: Sanity checks for input here
	x <- as.matrix(x)
	# ...
	if(verb>=2) print("Sanity checks ready")

	# TODO: Currently Cox assumes that event and time come in certain order in the 2-column y
	# Flip them depending on which way is expected
	if(any(class(y) %in% c("matrix", "array", "Surv"))){
		if(all(y[,1] %in% c(0,1))){

		}else if(all(y[,2] %in% c(0,1))){

		}
	}
	

	# If kit matrix is missing as input, assume that each variable is alone
	if(missing(k)){
		k <- matrix(0, nrow=ncol(x), ncol=ncol(x))
		diag(k) <- 1
	}
	## -> Intercept is an independent variable that is not subjected to penalization
	# If family is not Cox, (Intercept) requires its own row/column in K
	if(!family == "cox" & ncol(k) == ncol(x)){
		k <- rbind(0, cbind(0, k))
		k[1,1] <- 1
		if(!is.null(rownames(k)) & !is.null(colnames(k))) rownames(k)[1] <- colnames(k)[1] <- "(Intercept)"
	}
	# Checking k structure
	if(verb>=2){
		print("Input k:")
		print(k)
	}
	# If user has defined kmax, use it to pass to the Fortran function; otherwise kmax is the maximum number of kits in the input data
	if(missing(kmax)){
		kmax <- nrow(k)
	# Sanity checking for user provided parameter
	}else if(class(kmax) %in% "numeric"){
		kmax <- as.integer(kmax)
	}else if(!any(class(kmax) %in% c("integer", "numeric"))){
		stop("Provided kmax parameter ought to be of type 'integer' or 'numeric' cast to an integer")
	}
	if(verb>=2) print("Preprocessing k ready")

	# If kit weights are missing, assume them to be unit cost
	if(missing(w)){
		w <- rep(1, times=nrow(k))
	}
	# If cost for intercept has not been incorporated, add that as 0
	if(!family == "cox" & length(w) == ncol(x)){
		w <- c(0, w)
	}
	# Checking k structure
	if(verb>=2){
		print("Input w:")
		print(w)
	}
	if(verb>=2) print("Preprocessing w ready")

	# Call correct internal function based on the specified model family
	if(family=="cox"){
		# Call C function for Cox regression
		res <- .Call(c_casso_cox_f, 
			as.double(x), # Data matrix x
			as.double(y), # Response y
			as.integer(k), # Kit indicator matrix k 
			as.double(w), # Kit weights/costs
			as.integer(nrow(x)), # Number of samples (rows in x)
			as.integer(ncol(x)), # Number of variables (columns in x)
			as.integer(nrow(k)), # Number of kits
			as.integer(print), # Tuning parameter for verbosity
			as.integer(start), # Tuning parameter for starting values
			as.integer(kmax) # Tuning parameter for max k run
		)
		if(verb>=2){
			print(res)
		}
		# Beta per k steps
		bperk <- matrix(res[[1]], nrow = ncol(x), ncol = nrow(k))
		# Naming rows/cols/vector elements
		rownames(bperk) <- colnames(x)
	# Gaussian / normal distribution fit using mean-squared error	
	}else if(family %in% c("mse", "gaussian")){
		# Call C function for Mean-Squared Error regression
		res <- .Call(c_casso_mse_f, 
			as.double(x), # Data matrix x
			as.double(y), # Response y
			## Requires artificial addition of intercept kit?
			as.integer(k), # Kit indicator matrix k 
			as.double(w), # Kit weights/costs
			as.integer(nrow(x)), # Number of samples (rows in x)
			as.integer(ncol(x)), # Number of variables (columns in x)
			as.integer(nrow(k)), # Number of kits
			as.integer(print), # Tuning parameter for verbosity
			as.integer(start), # Tuning parameter for starting values
			as.integer(kmax) # Tuning parameter for max k run
		)
		# Beta per k steps
		# Add row for intercept
		bperk <- matrix(res[[1]], nrow = ncol(x)+1, ncol = nrow(k))
		# Naming rows/cols/vector elements
		rownames(bperk) <- c(colnames(x),"intercept")
		
	}else if(family == "logistic"){
		# Call C function for logistic regression
		res <- .Call(c_casso_logistic_f, 
			as.double(x), # Data matrix x
			as.integer(y), # Response y
			## Requires artificial addition of intercept kit?
			as.integer(k), # Kit indicator matrix k 
			as.double(w), # Kit weights/costs
			as.integer(nrow(x)), # Number of samples (rows in x)
			as.integer(ncol(x)), # Number of variables (columns in x)
			as.integer(nrow(k)), # Number of kits
			as.integer(print), # Tuning parameter for verbosity
			as.integer(start), # Tuning parameter for starting values
			as.integer(kmax) # Tuning parameter for max k run
		)
		# Beta per k steps
		# Add row for intercept
		bperk <- matrix(res[[1]], nrow = ncol(x)+1, ncol = nrow(k))
		# Naming rows/cols/vector elements
		rownames(bperk) <- c(colnames(x),"intercept")
	
	}

	# If verb is high enough, print midpoints
	if(verb>=3){
		print(res)
		print(dim(x))
		print(length(y))
		print(dim(k))
		print(length(w))
	}

	# Beta per k steps
	# Cox regression doesn't have intercept
	if(family == "cox"){
		bperk <- matrix(res[[1]], nrow = ncol(x), ncol = nrow(k))
		rownames(bperk) <- colnames(x)
	# All other model families have intercept
	}else{
		bperk <- matrix(res[[1]], nrow = ncol(x)+1, ncol = nrow(k))
		# Fortran subroutine returns intercept as the last element, swap that to front
		rownames(bperk) <- c(colnames(x), "(Intercept)")
		bperk <- bperk[c(nrow(bperk), 1:(nrow(bperk)-1)),]
	}

	# Target function values per k steps
	fperk <- as.numeric(res[[2]])
	# Naming rows/cols/vector elements
	names(fperk) <- colnames(bperk) <- paste("k_", 1:nrow(k), sep="")
	
	if(verb>=2){
		print("fperk")
		print(fperk)
	}
	# Transpose so that coefficients are columns and k-steps are rows in bperk
	bperk <- t(bperk)

	if(verb>=2){
		print("bperk")
		print(dim(bperk))
		print(bperk)
	}
	
	# Kits picked per each k-step
	kperk <- t(apply(bperk, MARGIN=1, FUN=function(z) as.integer(!z==0)))
	if(verb>=3){
		print("kperk1")
		print(dim(kperk))
		print(kperk)
	}
	# If dealing with non-Cox regression, omit (Intercept) from indicator matrix
	if(!family == "cox" & ncol(k) == ncol(x)){
		kperk <- kperk[,-1]
	}
	
	# Indices as a named vector
	kperk <- as.list(apply(kperk %*% t(k), MARGIN=1, FUN=function(z) { which(!z==0) }))
	if(verb>=3){
		print("kperk2")
		print(dim(kperk))
		print(kperk)
	}
	
	# Kit costs per each k-step
	cperk <- unlist(lapply(kperk, FUN=function(z) { sum(w[z], na.rm=TRUE) }))

	if(verb>=2){
		print("cperk")
		print(cperk)
	}

	# Return the freshly built S4 model object
	obj <- new("casso", 
		## Model fit results
		bperk = bperk[1:kmax,,drop=FALSE], # Beta coef per k-steps up to kmax
		fperk = fperk[1:kmax], # Target function values per k-steps up to kmax
		kperk = kperk[1:kmax], # Chosen kits per k-steps up to kmax
		cperk = cperk[1:kmax], # Total kit costs per each k-step up to kmax
		## Data slots
		x=as.matrix(x),	# Data matrix X
		y=as.matrix(y), # Response Y
		k=as.matrix(k), # Kit matrix K
		w=as.numeric(w),	# Vector of weights/costs for kits
		## Additional parameters
		family=as.character(family),	# Model family as a character string
		start = start,	# Method for generating starting points
		kmax = kmax	# Max run k-step
	)

	if(verb>=2){
		print("obj template created successfully")
	}
	
	# Fit lm/glm/coxph/... models per each estimated set of beta coefs (function call depends on 'family')
	try({
		obj@fits <- apply(obj@bperk, MARGIN=1, FUN=function(bs){
			if(family=="cox"){
				## Prefit a coxph-object
				survival::coxph(
					as.formula(paste("survival::Surv(time=obj@y[,1],event=obj@y[,2]) ~",paste(colnames(obj@x),collapse='+'))), # Formula for response 'y' modeled using data matrix 'x' 
					data=data.frame(obj@x), # Use data matrix 'x'
					init = bs, # Use model coefficients obtained using the DBDC optimization 
					control = survival::coxph.control(iter.max=0) # Prevent iterator from deviating from prior model parameters
				)
			}else if(family %in% c("mse", "gaussian")){
				## Prefit a linear glm-object with gaussian error; use heavily stabbed .glm.fit.mod allowing maxit = 0
				stats::glm(
					as.formula(paste("y ~",paste(colnames(obj@x),collapse='+')))
					, data = data.frame(obj@x), start = bs, family = gaussian(link="identity"), method = casso:::.glm.fit.mod
				)
			}else if(family=="logistic"){
				## Prefit a logistic glm-object with logistic link function; use heavily stabbed .glm.fit.mod allowing maxit = 0
				stats::glm(
					as.formula(paste("y ~",paste(colnames(obj@x),collapse='+')))
					, data = data.frame(obj@x),, start = bs, family = binomial(link="logit"), method = casso:::.glm.fit.mod
				)
			}
		})
	})

	if(verb>=2){
		print("fits-slot created successfully")
	}
		
	# Calculate/extract model goodness metric at each k
	try({
		# Cox regression
		if(family=="cox"){
			# Use c-index as the goodness measure
			obj@goodness <- unlist(lapply(obj@fits, FUN=function(z) { z$concordance["concordance"] }))
			obj@metric <- "cindex"
		}else if(family %in% c("mse", "gaussian")){
			# Use mean squared error as the goodness measure
			obj@goodness <- unlist(lapply(obj@fits, FUN=function(z) { mean((y - predict.glm(z, type="response"))^2) }))
			obj@metric <- "mse"
		}else if(family=="logistic"){
			# Use correct classification percent as the goodness measure
			obj@goodness <- unlist(lapply(obj@fits, FUN=function(z) { sum(as.numeric(y == (predict.glm(z, type="response")>0.5)))/length(y) }))
			obj@metric <- "accuracy"
		}
	})
	
	# Return the new model object
	obj
}
