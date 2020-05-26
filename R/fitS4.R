####
#
# S4-class OO for saving models and relevant parameters
# DBDC-fit kit-structured L0-penalized model fits
#
####

#' S4-class for casso
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
		kperk = list(),
		cperk = NA_real_,
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

	# Kits picked per each k-step
	kperk <- t(apply(bperk, MARGIN=1, FUN=function(z) as.integer(!z==0)))
	# Indices as a named vector
	kperk <- apply(kperk %*% t(k), MARGIN=1, FUN=function(z) { which(!z==0) })
	
	# Kit costs per each k-step
	cperk <- unlist(lapply(kperk, FUN=function(z) { sum(w[z]) }))

	# Return the freshly built S4 model object
	new("casso", 
		## Model fit results
		bperk = bperk, # Beta coef per k-steps
		fperk = fperk, # Target function values per k-steps
		kperk = kperk, # Chosen kits per k-steps
		cperk = cperk, # Total kit costs per each k-step
		## Data slots
		x=as.matrix(x), 
		y=as.matrix(y), 
		k=as.matrix(k), 
		w=as.numeric(w), 
		## Additional
		family=as.character(family)
	)
}


#' Showing casso-objects
#'
#' @export
setMethod("show", "casso",
	function(object){
	cat("Cardinality-constrained Absolute Subset Selection Optimizated model object\n")
	cat(paste("Model family: ", object@family,"\n", sep=""))
	cat(paste("k steps: ", nrow(object@k),"\n", sep=""))
	cat(paste("dim(x): [", paste(dim(object@x),collapse=","),"]\n", sep=""))
	cat(paste("dim(y): [", paste(dim(object@y),collapse=","),"]\n", sep=""))
})


#' Plotting casso-objects
#'
#' @export
setMethod("plot", "casso",
	function(
		object
	){
	legend <- TRUE
	mtexts <- TRUE
	par(las=2,  # All labels orthogonally to axes
		mar=c(7,4,1,4), # Inner margins
		oma=c(ifelse(mtexts, 2, 0), ifelse(mtexts, 2, 0), 0, ifelse(mtexts, 2, 0))) # Outer margins depend on additional labels with mtext
	y1 <- object@fperk
	y2 <- object@cperk
	x <- 1:nrow(object@k)
	# Two y-axes overlayed in a single graphics device, rigid first example
	plot.new()
	# First part (target function values)
	plot.window(xlim=c(1,length(x)), ylim=range(y1))
	axis(1, at=1:length(x), labels=x)
	axis(2, col.axis="red")
	box()
	points(1:length(x), y1, pch=16, col="red")
	points(1:length(x), y1, type="l", col="red")
	# Second part (cost accumulation)
	plot.window(xlim=c(1,length(x)), ylim=range(y2))
	axis(4, col.axis="blue")
	points(1:length(x), y2, pch=16, col="blue")
	points(1:length(x), y2, type="l", col="blue")
	if(legend){
		legend("top", col=c("red", "blue"), pch=16, lwd=1, legend=c("Target function value", "Cumulative kit cost"))
	}
	if(mtexts){
		mtext(side=1, text="K steps", las=0, outer=TRUE)
		mtext(side=2, text="Target function value", las=0, outer=TRUE)
		mtext(side=4, text="Cumulative kit costs", las=0, outer=TRUE)
	}
})



