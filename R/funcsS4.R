###
#
# Functions specific for the S4 class object 'casso'
# 'show', 'coef', etc base function replacements
#
###

###
#
# Replacing generics
#
###

#' Showing casso-objects
#'
#' @export
setMethod("show", "casso",
	function(object){
		cat("Cardinality-constrained Absolute Subset Selection Optimized model object\n")
		cat(paste("Model family: ", object@family,"\n", sep=""))
		cat(paste("k steps: ", object@kmax,"\n", sep=""))
		cat(paste("dim(x): [", paste(dim(object@x),collapse=","),"]\n", sep=""))
		cat(paste("dim(y): [", paste(dim(object@y),collapse=","),"]\n", sep=""))
	}
)

#' Taking coefficients of casso-objects
#'
#' @export
setMethod("coef", "casso",
	function(object, k){
		# Sanity checking for k-values
		if(missing(k)){
			stop("You need to provide parameter 'k' for obtaining coefficients at a certain k-step")
		}else if(k<1 | k>length(object@fits) | !is.numeric(k)){
			stop("Invalid k-value, should be an integer between {1,2, ..., kmax}")
		}
		# Returning the correct coef at k:th step
		stats::coef(object@fits[[k]])
	}
)

###
#
# casso-specific S4-functions
#
###


#' Return named vector of feature indices with a given k that are non-zero
#'
#' @export
setGeneric("feat",
	function(object, k){
		standardGeneric("feat")
	}
)
setMethod("feat", "casso",
	function(object, k){
		# Sanity checking for k-values
		if(missing(k)){
			stop("You need to provide parameter 'k' for obtaining coefficients at a certain k-step")
		}else if(k<1 | k>length(object@fits) | !is.numeric(k)){
			stop("Invalid k-value, should be an integer between {1,2, ..., kmax}")
		}
		# Returning the correct nonzero coefs at k:th step and name the indices according to data matrix column names
		tmp <- casso::coef(object@fits[[k]])
		# Return the non-zero coefs, named vector
		tmp[!tmp==0]
	}
)

#' Return named vector of kits indices with a given k that are non-zero
#'
#' @export
setGeneric("kits",
	function(object, k){
		standardGeneric("kits")
	}
)
setMethod("kits", "casso",
	function(object, k){
		# Sanity checking for k-values
		if(missing(k)){
			stop("You need to provide parameter 'k' for obtaining coefficients at a certain k-step")
		}else if(k<1 | k>length(object@fits) | !is.numeric(k)){
			stop("Invalid k-value, should be an integer between {1,2, ..., kmax}")
		}
		# Returning the correct kit indices while naming them
		kits <- object@kperk[[k]]
		# If data is without a kit structure, use names directly from variables (possibly including intercept)
		if(is.null(rownames(object@k)) & all(diag(object@k==1))){
			names(kits) <- colnames(object@bperk)[kits]
		# If data has a kit structure, use those for names
		}else{
			names(kits) <- rownames(object@k)[kits]
		}
		kits
	}
)
