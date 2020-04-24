###
#
# First preliminary tests for R <-> C <-> Fortran interfaces
# Testing based on https://www.avrahamadler.com/2018/12/09/the-need-for-speed-part-1-building-an-r-package-with-fortran/
#
##

#' Wrap up matrix to Fortran
blassocox <- function(x, y, kits, costs, ...){
	if(!is.matrix(x)) stop("Function argument 'x' should be of 'matrix' class")
	if(!is.double(x)) { storage.mode(x) <- 'double' }	
	else{
		ncol <- ncol(x)
		nrow <- nrow(x)
		nkits <- nrow(kits)
		if(!is.integer(ncol)) { storage.mode(ncol) <- 'integer' }
		if(!is.integer(nrow)) { storage.mode(nrow) <- 'integer' }

		# Store betas per kits
		betakits <- matrix(
			.Call(c_blassocox_f, as.double(x), as.double(y), as.integer(kits), as.double(costs), as.integer(nrow), as.integer(ncol), as.integer(nkits)),
			nrow = ncol(x),
			ncol = nrow(kits)
		)
		colnames(betakits) <- colnames(x)
		rownames(betakits) <- rownames(kits)
		betakits
	}
}
