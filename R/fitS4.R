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
		info = "character"	# Additional error messages, warnings or such reported e.g. during model fitting
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
		info = NA_character_
	),
	# Function for testing whether all S4-slots are legitimate for the S4 object 
	validity = function(
		object
	){
		# Return TRUE only if all slots fulfill the validity criteria (could be a list of if-statements that return FALSE prior to this if discrepancies are detected)
		TRUE
	}
)

#' Function for conveniently creating new model objects
#'
#' TODO: Extended explanation of the casso S4-class here.
#'
#' @param x Data matrix 'x'
#' @param y Response vector/two-column matrix 'y' (see: family); number of rows equal to nrow(x)
#' @param k Integer (0/1) kit indicator matrix; number of columns equal to ncol(x)
#' @param w Kit cost weight vector w of length nrow(k)
#' @param family Model family, should be one of: 'cox'
#' @param print Level of verbosity in Fortran (may not be visible on all terminals); should be an integer between {range, range}
#' @param start Starting point generation method, see vignettes for details; should be an integer between {range,range}
#'
#' @return casso S4-object fit using 
#'
#' @examples
#' data(ex)
#' fit <- casso(x=ex_X, y=ex_Y, k=ex_K, w=ex_c)
#' fit
#'
#' @export
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
	family = "cox",
	## Tuning parameters
	print=3,# Level of verbosity (-1 for tidy output, 3 for debugging level verbosity)
	start=2,#  
	verb=1  # Level of R verbosity (1 = standard, 2 = debug level, 0<= none)
){
	# Call correct internal function based on the specified model family
	if(family=="cox"){

		# Call C function
		res <- .Call(c_cassocox_f, 
			as.double(x), # Data matrix x
			as.double(y), # Response y
			as.integer(k), # Kit indicator matrix k 
			as.double(w), # Kit weights/costs
			as.integer(nrow(x)), # Number of samples (rows in x)
			as.integer(ncol(x)), # Number of variables (columns in x)
			as.integer(nrow(k)), # Number of kits
			as.integer(print), # Tuning parameter for verbosity
			as.integer(start) # Tuning parameter for starting values
		)
		if(verb>=2){
			print(res)
		}
	}
	# Beta per k steps
	bperk <- matrix(res[[1]], nrow = ncol(x), ncol = nrow(k))
	# Target function values per k steps
	fperk <- as.numeric(res[[2]])
	# Naming rows/cols/vector elements
	rownames(bperk) <- colnames(x)
	names(fperk) <- colnames(bperk) <- paste("k_", 1:nrow(k), sep="")
	# Transpose so that coefficients are columns and k-steps are rows in bperk
	bperk <- t(bperk)

	# Kits picked per each k-step
	kperk <- t(apply(bperk, MARGIN=1, FUN=function(z) as.integer(!z==0)))
	# Indices as a named vector
	kperk <- as.list(apply(kperk %*% t(k), MARGIN=1, FUN=function(z) { which(!z==0) }))
	
	# Kit costs per each k-step
	cperk <- unlist(lapply(kperk, FUN=function(z) { sum(w[z]) }))

	if(verb>=2){
		print(bperk)
		print(fperk)
		print(kperk)
	}

	# Return the freshly built S4 model object
	obj <- new("casso", 
		## Model fit results
		bperk = bperk, # Beta coef per k-steps
		fperk = fperk, # Target function values per k-steps
		kperk = kperk, # Chosen kits per k-steps
		cperk = cperk, # Total kit costs per each k-step
		start = start, # Method for generating starting points
		## Data slots
		x=as.matrix(x), 
		y=as.matrix(y), 
		k=as.matrix(k), 
		w=as.numeric(w), 
		## Additional
		family=as.character(family)
	)
	
	# Fit lm/glm/coxph/... models per each estimated set of beta coefs (function call depends on 'family')
	obj@fits <- apply(bperk, MARGIN=1, FUN=function(bs){
		if(family=="cox"){
			survival::coxph(
				survival::Surv(time=obj@y[,1], event=obj@y[,2]) ~ obj@x, # Formula for response 'y' modeled using data matrix 'x' 
				init = bs, # Use model coefficients obtained using the DBDC optimization 
				control = survival::coxph.control(iter.max=0) # Prevent iterator from deviating from prior model parameters
			)
		}
	})
	
	# Calculate/extract model goodness metric at each k
	if(family=="cox"){
		obj@goodness <- unlist(lapply(obj@fits, FUN=function(z) { z$concordance["concordance"] }))
	}
	
	# Return the new model object
	obj
}


#' Showing casso-objects
#'
#' @export
setMethod("show", "casso",
	function(object){
	cat("Cardinality-constrained Absolute Subset Selection Optimized model object\n")
	cat(paste("Model family: ", object@family,"\n", sep=""))
	cat(paste("k steps: ", nrow(object@k),"\n", sep=""))
	cat(paste("dim(x): [", paste(dim(object@x),collapse=","),"]\n", sep=""))
	cat(paste("dim(y): [", paste(dim(object@y),collapse=","),"]\n", sep=""))
})


