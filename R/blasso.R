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

