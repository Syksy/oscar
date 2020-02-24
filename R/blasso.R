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
		#print(paste("R: nrow:", nrow, "& ncol:", ncol))
		.Call(c_blassocox_f, as.double(x), as.double(y), as.integer(kits), as.double(costs), as.integer(nrow), as.integer(ncol), as.integer(nkits))
	}
}

#' R wrapper for Cox's DBDC
# Thus far:
#
#> library(blasso)
#> data(ex)
#> blasso::dbdc(
#x=     y=     k=     c=     start= print= ...=
#> blasso::dbdc(x=ex_X, y=ex_Y, c=ex_c, k=ex_K)
#Defining data constants...
#Preparing output variables...
#Reserving memory for output variables...
#Calling Fortran DBDC...
#
# ... crash.
dbdc <- function(
	x, 	# data matrix X
	y,	# two-column survival vector(s) Y (time, event)
	k, 	# kit matric K
	c, 	# cost vector c
	startPar = 1,	# Start control parameter
	printPar = 1,	# iprint control parameter
	... 	#
){
	nr <- nrow(x)
	nc <- ncol(x)
	nk <- nrow(k)
	if(!is.integer(nr)) { storage.mode(nr) <- 'integer' }
	if(!is.integer(nc)) { storage.mode(nc) <- 'integer' }
	if(!is.integer(nk)) { storage.mode(nk) <- 'integer' }
	# From the corresponding R <-> C wrapper:
	# ...
	# extern SEXP c_coxdbdc_loop_kits_f(SEXP nf, SEXP nr, SEXP nk, SEXP in_vX, SEXP in_vY, SEXP in_vK, SEXP in_vC, SEXP startPar, SEXP printPar){
	# ...
	# F77_CALL(coxdbdc_loop_kits_f)(REAL(beta_for_k), REAL(f_for_k), REAL(in_vX), INTEGER(in_vY), INTEGER(in_vK), REAL(in_vC), REAL(CPUtime), nft, nrecords, nkits, iprint, start);
	.Call(c_coxdbdc_loop_kits_f, as.integer(nc), as.integer(nr), as.integer(nk), as.double(x), as.integer(y), as.integer(k), as.double(c), as.integer(startPar), as.integer(printPar))
}

test <- function(
	x, 	# data matrix X
	y,	# two-column survival vector(s) Y (time, event)
	k, 	# kit matric K
	c, 	# cost vector c
	startPar = 1,
	printPar = 1,
	... 	#
){
	nr <- nrow(x)
	nc <- ncol(x)
	nk <- nrow(k)
	if(!is.integer(nr)) { storage.mode(nr) <- 'integer' }
	if(!is.integer(nc)) { storage.mode(nc) <- 'integer' }
	if(!is.integer(nk)) { storage.mode(nk) <- 'integer' }
	# extern SEXP c_coxdbdc_test_f(SEXP nf, SEXP nr, SEXP nk, SEXP in_vX, SEXP in_vY, SEXP in_vK, SEXP in_vC){
	.Call(c_coxdbdc_test_f, as.integer(nc), as.integer(nr), as.integer(nk), as.double(x), as.integer(y), as.integer(k), as.double(c), as.integer(startPar), as.integer(printPar))
}
