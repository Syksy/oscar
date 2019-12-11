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
	if(!is.double(x)) { storage.mode(x) <- 'double' }	
	else{
		ncol <- ncol(x)
		nrow <- nrow(x)
		if(!is.integer(ncol)) { storage.mode(ncol) <- 'integer' }
		if(!is.integer(nrow)) { storage.mode(nrow) <- 'integer' }
		print(paste("R: nrow:", nrow, "& ncol:", ncol))
		.Call(c_blassomat_f, as.double(x), as.integer(nrow), as.integer(ncol))
	}
}

