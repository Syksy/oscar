###
#
# Functions specific for the S4 class object 'casso'
# 'show', 'coef', etc base function replacements
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

