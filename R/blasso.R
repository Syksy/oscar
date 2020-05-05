###
#
# First preliminary tests for R <-> C <-> Fortran interfaces
# Testing based on https://www.avrahamadler.com/2018/12/09/the-need-for-speed-part-1-building-an-r-package-with-fortran/
#
##

#' Wrap up matrix to Fortran
blassocox <- function(x, y, kits, costs, print=3, start=2, ...){
	if(!is.matrix(x)) stop("Function argument 'x' should be of 'matrix' class")
	if(!is.double(x)) { storage.mode(x) <- 'double' }	
	# Sanity checks to prevent nonsensical data being sent to Fortran subroutine
	if(!ncol(kits)==ncol(x)) stop("'kits' matrix should have as many columns (variables) as the input data matrix 'x'")
	if(!all(kits) %in% c(0,1)) stop("'kits' should be an indicator matrix composing only of ones and zeroes")
	if(!nrow(kits)==length(costs)) stop("The number of rows in 'kits' (i.e. different kits) should be of same length as costs matrix (i.e. cost for each kits)")

	# Call C-interface and Fortran-subroutine	
	ncol <- ncol(x)
	nrow <- nrow(x)
	nkits <- nrow(kits)
	if(!is.integer(ncol)) { storage.mode(ncol) <- 'integer' }
	if(!is.integer(nrow)) { storage.mode(nrow) <- 'integer' }

	# Store betas per kits
	betakits <- matrix(
		.Call(c_blassocox_f, as.double(x), as.double(y), as.integer(kits), as.double(costs), as.integer(nrow), as.integer(ncol), as.integer(nkits),as.integer(print),as.integer(start)),
		nrow = ncol(x),
		ncol = nrow(kits)
	)
	rownames(betakits) <- colnames(x)
	colnames(betakits) <- paste("kstep_", 1:nrow(kits), sep="")
	betakits <- t(betakits)
	# Add complete set of kits in order as attributes to the beta/kits output
	kitorder <- apply(betakits, MARGIN=1, FUN=function(z) { which(!ex_K %*% z == 0) })
	names(kitorder) <- paste("kstep_",1:nrow(kits))
	attr(betakits, "kitorder") <- kitorder
	# Pick kit names and the order in which each individual kit was added
	attr(betakits, "kitnames") <- rownames(kits)[c(kitorder[[1]], unlist(lapply(2:length(kitorder), FUN=function(z) { base::setdiff(kitorder[[z]], kitorder[[z-1]]) })))]
	# Rows are the sequence in which kits are picked; columns are variables, which are picked in bundles according to kits
	betakits
}
