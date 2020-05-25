####
#
# S4-class OO for saving models and relevant parameters
# DBDC-fit kit-structured L0-penalized model fits
#
####

#' S4-class for basso
setClass("casso", # abbreviation
	representation(
		# Set class representations for the S4-class @-slots
		## Results from Fortran
		bperk = "matrix",	# Rows: different kit k-values, columns: beta coefficients  in model fit
		fperk = "numeric",	# Target function values at each k-value
		#cputime = "numeric",	# Total CPU time for computations
		kitorder = "integer",	# Integer order in which kits were selected
		kitnames = "character",	# Names of the corresponding kits (rownames of the provided kit indicator matrix)
		rho = "numeric",	# Vector of utilized rho-values in the penalization
		## Provided data, saved in R
		x = "matrix",		# Original data matrix for fitting
		y = "matrix",		# Response vector (two-column Surv-object for Cox modelling)
		k = "matrix",		# Kit indicator matrix
		w = "numeric",		# Kit value (cost ~ weights) vector, length ought to be same as nrow(k)
		## Additional
		family = "character",	# Utilized model family, in early versions only 'cox' supported (to be added: normal/gaussian, logistic, ...
		fits = "list",		# Pre-fitted model objects, where each model fit only allows the optimal subset of features to be incorporated into the model; length should be equal to k sequence
		info = "character"	# Additional error messages, warnings or such reported e.g. during model fitting
	),	
	prototype(
		# Prototype base model object
		## Results from Fortran
		bperk = matrix(NA, nrow=0, ncol=0),
		fperk = NA_real_,
		kitorder = NA_integer_,
		kitnames = NA_character_,
		rho = NA_real_,
		## Provided data, saved in R
		x = matrix(NA, nrow=0, ncol=0),
		y = matrix(NA, nrow=0, ncol=0),
		k = matrix(NA, nrow=0, ncol=0),
		w = NA_real_,
		## Additional
		family = NA_character_,
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
#'
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
	start=2 #  
){
	# Call correct internal function based on the specified model family
	if(family=="cox"){

		# Call C function
		res <- .Call(c_blassocox_f, 
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
		#print(fit)
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

	# Return the freshly built S4 model object
	new("casso", 
		## Model fit results
		bperk = bperk,
		fperk = fperk,
		## Data slots
		x=as.matrix(x), 
		y=as.matrix(y), 
		k=as.matrix(k), 
		w=as.numeric(w), 
		## Additional
		family=as.character(family)
	)
}





