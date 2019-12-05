###
#
# First preliminary tests for R <-> C <-> Fortran interfaces
# Testing based on https://www.avrahamadler.com/2018/12/09/the-need-for-speed-part-1-building-an-r-package-with-fortran/
#
##

#' The first test blasso function
testblasso_f <- function(x){
	if(!is.double(x)) { storage.mode(x) <- 'double' }
	.Call(c_testblasso_f, x) # Calculate sum of vector values
}

#' Wrap up matrix to Fortran
blassomat_f <- function(x){
	if(!is.matrix(x)) stop("Function argument 'x' should be of 'matrix' class")
	else{
		ncol <- ncol(x)
		nrow <- nrow(x)
		# .Call(c_blassomat_f, as.double(x), as.integer(ncol), as.integer(nrow)
	}
}
