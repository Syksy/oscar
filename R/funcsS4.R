###
#
# Functions specific for the S4 class object 'oscar'
# 'show', 'coef', etc base function replacements
#
###

###
#
# Replacing generics
#
###

#' Showing oscar-objects
#'
#' @rdname generics
#'
#' @export
setMethod("show", "oscar",
	function(object){
		cat("Optimal Subset Cardinality Regression (OSCAR) model object\n")
		cat(paste("Model family: ", object@family,"\n", sep=""))
		cat(paste("k steps: ", object@kmax,"\n", sep=""))
		cat(paste("dim(x): [", paste(dim(object@x),collapse=","),"]\n", sep=""))
		cat(paste("dim(y): [", paste(dim(object@y),collapse=","),"]\n", sep=""))
	}
)

#' Extract coefficients of oscar-objects
#'
#' @rdname generics
#'
#' @export
setMethod("coef", "oscar",
	function(object, k){
		# Sanity checking for k-values
		if(missing(k)){
			stop("You need to provide parameter 'k' for obtaining coefficients at a certain k-step")
		}else if(k<1 | k>object@kmax | !is.numeric(k)){
			stop("Invalid k-value, should be an integer between {1,2, ..., kmax}")
		}
		# Returning the correct coef at k:th step
		#stats::coef(object@fits[[k]])
		object@bperk[k,]
	}
)

#' Predicting based on oscar-objects
#'
#' @rdname generics
#'
#' @export
setMethod("predict", "oscar",
	# Showing possible options for types of responses
	function(object, k, type = c("response","link","nonzero","coefficients","label"), newdata = object@x){
		# Sanity checking for k-values
		if(missing(k)){
			stop("You need to provide parameter 'k' for obtaining coefficients at a certain k-step")
		}else if(k<1 | k>object@kmax | !is.numeric(k)){
			stop("Invalid k-value, should be an integer between {1,2, ..., kmax}")
		}else if(!object@family %in% c("cox", "mse", "gaussian", "logistic")){
			stop(paste("Invalid family-parameter in oscar-object:", object@family))
		}
		# Returning the correct coef at k:th step
		if(type[1] == "response"){
			if(object@family == "cox"){
				# Xb
				as.matrix(newdata) %*% t(object@bperk[k,,drop=FALSE])
			}else{
				# Need to add intercept to b0 + Xb in 
				cbind(1, as.matrix(newdata)) %*% t(object@bperk[k,,drop=FALSE])
			}
		# Non-zero coefficients
		}else if(type[1] == "nonzero"){
			object@bperk[k,which(!object@bperk[k,]==0)]
		# All coefficients
		}else if(type[1] == "coefficients"){
			object@bperk[k,]
		# Xb run through the link function
		}else if(type[1] == "link"){		
			if(object@family %in% "logistic"){
				Xb <- cbind(1, as.matrix(newdata)) %*% t(object@bperk[k,,drop=FALSE])
				1/(1+exp(-Xb))
			}else if(object@family %in% c("cox")){
				exp(as.matrix(newdata) %*% t(object@bperk[k,,drop=FALSE]))
			}
		# Class labels (binary for starters)
		}else if(type[1] == "label"){
			if(!object@family == "logistic"){
				stop("Parameter type == 'label' is intended for logistic or multinomial predictions")
			}
			Xb <- cbind(1, as.matrix(newdata)) %*% t(object@bperk[k,,drop=FALSE])
			as.integer(1/(1+exp(-Xb))>0.5)
		}
	}
)

#' Plot oscar-coefficients as a function of k and override default plot generic
#'
#' @rdname generics
#'
#' @export
setMethod("plot", "oscar",
	function(x, y, k=1:x@kmax, add=FALSE, intercept=FALSE, ...){
		# Should intercept be omitted from the plot
		if(!intercept & "(Intercept)" %in% colnames(x@bperk)){
			bperk <- x@bperk[,-which(colnames(x@bperk)=="(Intercept)")]
		}else{
			bperk <- x@bperk
		}
		# Not adding to an existing plot, making a fresh frame
		if(!add){
			plot.new() 
			plot.window(xlim=extendrange(k), ylim=extendrange(bperk))
			box()
			axis(1)
			axis(2)
			title(xlab="Cardinality 'k'", ylab="Beta coefficient")
			abline(h=0, col="grey", lwd=2)
		}
		# Plot each coefficient as its own line with a different colour
		ret <- lapply(1:ncol(bperk), FUN=function(i){
			vec <- bperk[,i]
			names(vec) <- paste("k_", k, sep="")
			points(x=k, y=vec, col=i, type="l")
			vec
		})
		names(ret) <- colnames(bperk)
		# Return a bperk coefficients list which is plotted
		invisible(ret)
	}
)


###
#
# oscar-specific S4-functions
#
###


#' Return named vector of feature indices with a given k that are non-zero
#'
#' @rdname generics
#'
#' @export
setGeneric("feat",
	function(object, k){
		standardGeneric("feat")
	}
)
setMethod("feat", "oscar",
	function(object, k){
		# Sanity checking for k-values
		if(missing(k)){
			stop("You need to provide parameter 'k' for obtaining coefficients at a certain k-step")
		}else if(k<1 | k>object@kmax | !is.numeric(k)){
			stop("Invalid k-value, should be an integer between {1,2, ..., kmax}")
		}
		# Returning the correct nonzero coefs at k:th step and name the indices according to data matrix column names
		tmp <- oscar::coef(object,k=k)
		# Return the non-zero coefs, named vector
		tmp[!tmp==0]
	}
)

#' Return named vector of kits indices with a given k that are non-zero
#'
#' @rdname generics
#'
#' @export
setGeneric("kits",
	function(object, k){
		standardGeneric("kits")
	}
)
setMethod("kits", "oscar",
	function(object, k){
		# Sanity checking for k-values
		if(missing(k)){
			stop("You need to provide parameter 'k' for obtaining coefficients at a certain k-step")
		}else if(k<1 | k>object@kmax | !is.numeric(k)){
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
