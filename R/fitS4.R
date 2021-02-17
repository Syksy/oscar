####
#
# S4-class OO for saving models and relevant parameters
# DBDC-fit kit-structured L0-penalized model fits
#
####

#' S4-class for oscar
#'
#'
#'
#' @export
setClass("oscar", # abbreviation
	representation(
		# Set class representations for the S4-class @-slots
		## Results from Fortran
		bperk = "matrix",	# Rows: different kit k-values, columns: beta coefficients  in model fit
		fperk = "numeric",	# Target function values at each k-value
		kperk = "list",		# Which kits were picked at which k-step (named indices vector)
		cperk = "numeric",	# Cumulative kit costs per each k-value
		#cputime = "numeric",	# Total CPU time for computations
		rho = "numeric",	# Vector of utilized rho-values in the penalization
		start = "numeric",	# Method of generating starting values
		## Provided data, saved in R
		x = "matrix",		# Original data matrix for fitting
		y = "matrix",		# Response vector (two-column Surv-object for Cox modelling)
		k = "matrix",		# Kit indicator matrix
		w = "numeric",		# Kit value (cost ~ weights) vector, length ought to be same as nrow(k)
		## Additional
		family = "character",	# Utilized model family, in early versions only 'cox' supported (to be added: normal/gaussian, logistic, ...
		goodness = "numeric",	# Goodness metric at each k based on model fit; dependent on 'family': cox is concordance-index, ...
		fits = "list",		# Pre-fitted model objects, where each model fit only allows the optimal subset of features to be incorporated into the model; length should be equal to k sequence
		info = "character",	# Additional error messages, warnings or such reported e.g. during model fitting
		kmax = "integer",	# Number of maximum k tested
		metric = "character",	# Name of the goodness-of-fit metric used
		AIC = "numeric"		# Akaike's Information Criterion for each member of @fits
	),	
	prototype(
		# Prototype base model object
		## Results from Fortran
		bperk = matrix(NA, nrow=0, ncol=0),
		fperk = NA_real_,
		kperk = list(),
		cperk = NA_real_,
		rho = NA_real_,
		start = NA_real_,
		## Provided data, saved in R
		x = matrix(NA, nrow=0, ncol=0),
		y = matrix(NA, nrow=0, ncol=0),
		k = matrix(NA, nrow=0, ncol=0),
		w = NA_real_,
		## Additional
		family = NA_character_,
		goodness = NA_real_,
		fits = list(),
		info = NA_character_,
		kmax = NA_integer_,
		metric = NA_character_,
		AIC = NA_real_
	),
	# Function for testing whether all S4-slots are legitimate for the S4 object 
	# TODO
	validity = function(
		object
	){
		# Return TRUE only if all slots fulfill the validity criteria 
		# (could be a list of if-statements that return FALSE prior to this if discrepancies are detected)
		TRUE
	}
)


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x Data matrix 'x'
#' @param y Response vector/two-column matrix 'y' (see: family); number of rows equal to nrow(x)
#' @param k Integer (0/1) kit indicator matrix; number of columns equal to ncol(x)
#' @param w Kit cost weight vector w of length nrow(k)
#' @param family Model family, should be one of: 'cox', 'mse'/'gaussian', or 'logistic, Default: 'gaussian'
#' @param print Level of verbosity in Fortran (may not be visible on all terminals); should be an integer between {range, range}, Default: 3
#' @param start Starting point generation method, see vignettes for details; should be an integer between {range,range}, Default: 2
#' @param verb Integer with additional integer values giving verbal feedback, Default: 1
#' @param kmax Maximum k step tested, by default all k are tested from k to maximum dimensionality
#' @param sanitize Whether input column names should be cleaned of potentially problematic symbols
#' @return Fitted oscar-object
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit
#'  }
#' }
#' @seealso 
#'  \code{\link[survival]{coxph}},\code{\link[survival]{coxph.control}}
#'  \code{\link[stats]{glm}}
#' @rdname oscar
#' @export 
#' @importFrom survival coxph coxph.control
#' @importFrom stats glm
oscar <- function(
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
	# Goodness metric used in model goodness-of-fit 
	metric,
	## Tuning parameters
	print=3,# Level of verbosity (-1 for tidy output, 3 for debugging level verbosity)
	start=2,# Strategy for choosing starting points at each n_k iteration
	verb=1, # Level of R verbosity (1 = standard, 2 = debug level, 3 = excessive debugging, 0<= none)
	kmax,    # Maximum tested k-values
	sanitize = TRUE,	# Whether column (i.e. variable) names are cleaned of potentially hazardous symbols
	in_mrounds = 5000, # The number of rounds in one main iteration 
	in_mit = 5000, # The number of main iteration 
	in_mrounds_esc = 5000, # The number of rounds in escape procedure
	in_b1, # Bundle B1, default min(n+5,1000) depends on the problem -> defined later
	in_b2 = 3, # Bundle B2
	in_b, # Bundle B in escape procedure, default 2n depends on the problem -> defined later
	in_m = 0.01, # Descent parameter
	in_m_clarke = 0.01, # Descent parameter in escape procedure
	in_c = 0.1, # Extra decrease parameter
	in_r_dec, # Decrease parameter, default depends on the problem -> defined later
	in_r_inc = 10^5, # Increase parameter
	in_eps1 = 5*10^(-5), # Enlargement parameter
	in_eps, # Stopping tolerance: Proximity measure, default depends on the problem -> defined later
	in_crit_tol =10^(-5) # Stopping tolerance: Criticality tolerance
){
	# TODO: Sanity checks for input here
	
	# Sanity checks: Are input matrices x and y available
	if(missing(x)){
		stop(paste("Input matrix x missing."))
	}
	if(missing(y)){
		stop(paste("Input matrix y missing."))
	}
	
	# Sanity checks: Are input x is a matrix. If not, try to convert into a matrix.
	if(!is.matrix(x)){
		try(x <- as.matrix(x))
	}

	
	# ...
	if(verb>=2) print("Sanity checks ready")

	# If no custom metric defined, using default goodness metrics:
	# mse/gaussian: mean-squared error
	# cox: c-index
	# logistic: (roc-)auc
	if(missing(metric)){
		if(family %in% c("mse","gaussian")){
			metric <- "mse"
		}else if(family == "cox"){
			metric <- "concordance"
		}else if(family == "logistic"){
			metric <- "auc"
		}else{
			stop(paste("Invalid parameter 'family':", family))
		}
	}else{
		# TODO; User provided custom metrics; 
		# check if legitimate (for logistic e.g. 'accuracy' and 'auc' ought to work)
	}

	# TODO: Currently Cox assumes that event and time come in certain order in the 2-column y
	# Flip them depending on which way is expected
	if(any(class(y) %in% c("matrix", "array", "Surv"))){
		if(all(y[,1] %in% c(0,1))){

		}else if(all(y[,2] %in% c(0,1))){

		}
	}
	
	# Sanity checks: Check that imput matrices x and y hav equal number of rows
	if(family=="cox"&&nrow(y)!=nrow(x)){
		stop(paste("Number of observations in the response matrix y (",nrow(y),") is not equal to number of observations in the predictor matrix x (",nrow(x),"). Please check both."))
	}
	
	###############################
	#### CHECKING kit matrix k ####
	# If kit matrix is missing as input, assume that each variable is alone
	if(missing(k)){
		k <- matrix(0, nrow=ncol(x), ncol=ncol(x))
		diag(k) <- 1
	}
	# Check that input k is a matrix
	if(!is.matrix(k)){
		try(k <- as.matrix(k))
	}
	
	# Check that kit matrix k and predictor matrix x have equal number of columns (features)
	if(family=="cox"&& ncol(k)!=ncol(x)){
		stop(paste("Number of columns in kit matrix k (",ncol(k),") should be equal to the amount of features (number of columns in x, ",ncol(x),"). Check that correct features are included."))
	}
	## -> Intercept is an independent variable that is not subjected to penalization
	# If family is not Cox, (Intercept) requires its own row/column in K
	if(!family == "cox" && ncol(k) == ncol(x)){
		k <- rbind(0, cbind(0, k))
		k[1,1] <- 1
		if(!is.null(rownames(k)) & !is.null(colnames(k))) rownames(k)[1] <- colnames(k)[1] <- "(Intercept)"
	}
	# Check that the kit matrix has only 0 or 1
	if(!all(k %in% c(0,1))){
		stop(paste("Values in the kit matrix k should be 0 or 1. Please check the kit matrix."))
	}
	# Check that there are no zero-columns in the kit matrix
	if(any(apply(k,MARGIN=2,FUN=sum)==0)){
		stop(paste("Zero column in the kit matrix k. Please check that each feature is exactly in one kit."))  ## For now only one kit per feature allowed.
	}
	# Check that each feature is only in one kit
	if(any(apply(k,MARGIN=2,FUN=sum)>1)){
		stop(paste("Some feature(s) are in multiple kits. Please check that each feature is exactly in one kit."))  ## For now only one kit per feature allowed.
	}
	# Check that each kit has at least one feature (not sure if necessary?))
	if(any(apply(k,MARGIN=1,FUN=sum)<1)){
		stop(paste("Some kit(s)) don't have any features. Please check that each kit has at least one feature."))
	}
	
	# Checking k structure
	if(verb>=2){
		print("Input k:")
		print(k)
	}

	########################
	#### CHECKING kmax  ####
	# If user has defined kmax, use it to pass to the Fortran function; otherwise kmax is the maximum number of kits in the input data
	if(missing(kmax)){
		kmax <- nrow(k)
	# Sanity checking for user provided parameter
	}else if(class(kmax) %in% "numeric"){
		kmax <- as.integer(kmax)
	}else if(!any(class(kmax) %in% c("integer", "numeric"))){
		stop("Provided kmax parameter ought to be of type 'integer' or 'numeric' cast to an integer")
	}
	if(verb>=2) print("Preprocessing k ready")


	############################
	#### CHECKING weights w ####
	# If kit weights are missing, assume them to be unit cost
	if(missing(w)){
		w <- rep(1, times=nrow(k))
	}
	
	# Check that w length is equal to number of rows in the kit matrix k (number of kits)
	if(family=="cox"&& length(w)!=nrow(k)){
		stop(paste("Number of kit costs (",length(w),") is not equal to number of kits (number of rows in the kit matrix k, ",nrow(k),"). Check that correct kits are included."))
	}
	
	# If cost for intercept has not been incorporated, add that as 0
	if(!family == "cox" & length(w) == ncol(x)){
		w <- c(0, w)
	}
	# Checking w structure
	if(verb>=2){
		print("Input w:")
		print(w)
	}
	if(verb>=2) print("Preprocessing w ready")
	
	#########################################################
	#### CHECKING tuning parameters and setting defaults ####
	if(in_mrounds <=0){ # The number of rounds in one main iteration 
		in_mrounds <- 5000
		warnings("Input in_mrounds should be >0. Default value 5000 is used instead.")
	}
	if(!is.integer(in_mrounds)){ warnings("Input in_mrounds should be an integer. The value is rounded to the nearest integer.")
	}
	if(in_mit<=0){ # The number of main iteration
		in_mit <- 5000
		warnings("Input in_mit should be >0. Default value 5000 is used instead.")
	}
	if(!is.integer(in_mit)){ warnings("Input in_mit should be an integer. The value is rounded to the nearest integer.")
	}
	if(in_mrounds_esc <=0){ # The number of rounds in escape procedure
		in_mrounds_esc <- 5000
		warnings("Input in_mrounds_esc should be >0. Default value 5000 is used instead.")
	}
	if(!is.integer(in_mrounds_esc)){ warnings("Input in_mrounds_esc should be an integer. The value is rounded to the nearest integer.")
	}
	if(missing(in_b1)){ # Bundle B1
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		in_b1 <- min(user_n+5,1000)	
	}
	if(in_b1 <=0){
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		in_b1 <- min(user_n+5,1000)
		warnings(paste("Input in_b1 should be >0. Default value ", in_b1, " is used instead",sep=""))
	}
	if(!is.integer(in_b1)){warnings("Input in_b1 should be an integer. The value is rounded to the nearest integer.")
	} 
	if(in_b2 <=0){ # Bundle B2
		in_b2 <- 3
		warnings("Input in_b2 should be >0. Default value 3 is used instead.")
	}
	if(!is.integer(in_b2)){ warnings("Input in_b2 should be an integer. The value is rounded to the nearest integer.")
	}
	if(missing(in_b)){ # Bundle B in escape procedure
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		in_b <- 2*user_n
	}
	if(in_b <=1){
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		in_b <- 2*user_n
		warnings(paste("Input in_b should be >1. Default value ", in_b," is used instead.",sep=""))
	}
	if(!is.integer(in_b)){ warnings("Input in_b should be an integer. The value is rounded to the nearest integer.")
	}
	if(in_m <=0 | in_m>=1){ # Descent parameter
		in_m <-0.2
		warnings("Input in_m should be in the interval (0,1). Default value 0.2 is used instead.")
	}
	if(in_m_clarke <=0 | in_m_clarke>=1){ # Descent parameter in escape procedure
		in_m_clarke <-0.01
		warnings("Input in_m_clarke should be in the interval (0,1). Default value 0.01 is used instead.")
	}
	if(in_c <=0 | in_c>=1){ # Extra decrease parameter
		in_c <-0.1
		warnings("Input in_c should be in the interval (0,1). Default value 0.1 is used instead.")
	}
	if(missing(in_r_dec)){ # Decrease parameter
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		if(user_n <10){in_r_dec <- 0.75
		}else if(user_n >=300){in_r_dec <- 0.99
		}else{in_r_dec <- trunc(user_n/(user_n+5)*100)/100  # Default is two first decimal from user_n/(user_n+5)
		}
	}
	if(in_r_dec <=0 | in_r_dec >=1){
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		if(user_n <10){in_r_dec <- 0.75
		}else if(user_n >=300){in_r_dec <- 0.99
		}else{in_r_dec <- trunc(user_n/(user_n+5)*100)/100  # Default is two first decimal from user_n/(user_n+5)
		}
		warnings(paste("Input in_r_dec should be in the interval (0,1). Default value ", in_r_dec," is used instead.",sep=""))
	}
	if(in_r_inc<=1){ # Increase parameter
		warnings("Input in_r_inc should be >1. Default value 10^7 is used instead.")
	}
	if(in_eps1<=0){ # Enlargement parameter
		warnings("Input in_eps1 should be >0. Default value 5*10^(-5) is used instead.")
	}
	if(missing(in_eps)){# Stopping tolerance: Proximity measure
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		if(user_n <=50){in_eps <- 10^(-6)
		}else{in_eps <- 10^(-5)
		}
	}
	if(in_eps<=0){
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		if(user_n <=50){in_eps <- 10^(-6)
		}else{in_eps <- 10^(-5)
		}
		warnings(paste("Input in_eps should be >0. Default value ", in_eps," is used instead.",sep=""))
	}
	#if(missing(in_crit_tol)){# Stopping tolerance: Criticality tolerance
	#	user_n <- ncol(x)
	#	if(family %in% c("mse","logistic","gaussian","normal")){
	#	user_n <- user_n +1}
	#	if(user_n <=200){in_crit_tol <- 10^(-5)
	#	}else{in_crit_tol <- 10^(-4)
	#	}
	#}
	if(in_crit_tol<=0){
		in_crit_tol <- 10^(-5)
		#user_n <- ncol(x)
		#if(family %in% c("mse","logistic","gaussian","normal")){
		#user_n <- user_n +1}
		#if(user_n <=200){in_crit_tol<- 10^(-5)
		#}else{in_crit_tol <- 10^(-4)
		#}
		warnings(paste("Input in_crit_tol should be >0. Default value ", in_crit_tol," is used instead.",sep=""))
	}	
	
	## Sanitize column names, replacing '+' with 'plus', '-' with 'minus', and ' ', '(' and ')' with '_'
	if(sanitize){
		colnames(x) <- gsub("\\ |\\(|)", "_", gsub("\\+", "plus", gsub("\\-", "minus", colnames(x))))
	}
	
	###############################
	#### Calling C and Fortran ####
	# Call correct internal function based on the specified model family
	if(family=="cox"){
		# Call C function for Cox regression
		res <- .Call(c_oscar_cox_f, 
			as.double(x), # Data matrix x
			as.double(y), # Response y
			as.integer(k), # Kit indicator matrix k 
			as.double(w), # Kit weights/costs
			as.integer(nrow(x)), # Number of samples (rows in x)
			as.integer(ncol(x)), # Number of variables (columns in x)
			as.integer(nrow(k)), # Number of kits
			as.integer(print), # Tuning parameter for verbosity
			as.integer(start), # Tuning parameter for starting values
			as.integer(kmax), # Tuning parameter for max k run
			as.integer(in_mrounds), # The number of rounds in one main iteration 
			as.integer(in_mit), # The number of main iteration 
			as.integer(in_mrounds_esc), # The number of rounds in escape procedure
			as.integer(in_b1), # Bundle B1
			as.integer(in_b2), # Bundle B2
			as.integer(in_b), # Bundle B in escape procedure
			as.double(in_m), # Descent parameter
			as.double(in_m_clarke), # Descent parameter in escape procedure
			as.double(in_c), # Extra decrease parameter
			as.double(in_r_dec), # Decrease parameter
			as.double(in_r_inc), # Increase parameter
			as.double(in_eps1), # Enlargement parameter
			as.double(in_eps), # Stopping tolerance: Proximity measure
			as.double(in_crit_tol) # Stopping tolerance: Criticality tolerance
		)
		if(verb>=2){
			print(res)
		}
		# Beta per k steps
		bperk <- matrix(res[[1]], nrow = ncol(x), ncol = nrow(k))
		# Naming rows/cols/vector elements
		rownames(bperk) <- colnames(x)
	# Gaussian / normal distribution fit using mean-squared error	
	}else if(family %in% c("mse", "gaussian")){
		# Call C function for Mean-Squared Error regression
		res <- .Call(c_oscar_mse_f, 
			as.double(x), # Data matrix x
			as.double(y), # Response y
			## Requires artificial addition of intercept kit?
			as.integer(k), # Kit indicator matrix k 
			as.double(w), # Kit weights/costs
			as.integer(nrow(x)), # Number of samples (rows in x)
			as.integer(ncol(x)), # Number of variables (columns in x)
			as.integer(nrow(k)), # Number of kits
			as.integer(print), # Tuning parameter for verbosity
			as.integer(start), # Tuning parameter for starting values
			as.integer(kmax), # Tuning parameter for max k run
			as.integer(in_mrounds), # The number of rounds in one main iteration 
			as.integer(in_mit), # The number of main iteration 
			as.integer(in_mrounds_esc), # The number of rounds in escape procedure
			as.integer(in_b1), # Bundle B1
			as.integer(in_b2), # Bundle B2
			as.integer(in_b), # Bundle B in escape procedure
			as.double(in_m), # Descent parameter
			as.double(in_m_clarke), # Descent parameter in escape procedure
			as.double(in_c), # Extra decrease parameter
			as.double(in_r_dec), # Decrease parameter
			as.double(in_r_inc), # Increase parameter
			as.double(in_eps1), # Enlargement parameter
			as.double(in_eps), # Stopping tolerance: Proximity measure
			as.double(in_crit_tol) # Stopping tolerance: Criticality tolerance
		)
		# Beta per k steps
		# Add row for intercept
		bperk <- matrix(res[[1]], nrow = ncol(x)+1, ncol = nrow(k))
		# Naming rows/cols/vector elements
		rownames(bperk) <- c(colnames(x),"intercept")
		
	}else if(family == "logistic"){
		# Call C function for logistic regression
		res <- .Call(c_oscar_logistic_f, 
			as.double(x), # Data matrix x
			as.integer(y), # Response y
			## Requires artificial addition of intercept kit?
			as.integer(k), # Kit indicator matrix k 
			as.double(w), # Kit weights/costs
			as.integer(nrow(x)), # Number of samples (rows in x)
			as.integer(ncol(x)), # Number of variables (columns in x)
			as.integer(nrow(k)), # Number of kits
			as.integer(print), # Tuning parameter for verbosity
			as.integer(start), # Tuning parameter for starting values
			as.integer(kmax), # Tuning parameter for max k run
			as.integer(in_mrounds), # The number of rounds in one main iteration 
			as.integer(in_mit), # The number of main iteration 
			as.integer(in_mrounds_esc), # The number of rounds in escape procedure
			as.integer(in_b1), # Bundle B1
			as.integer(in_b2), # Bundle B2
			as.integer(in_b), # Bundle B in escape procedure
			as.double(in_m), # Descent parameter
			as.double(in_m_clarke), # Descent parameter in escape procedure
			as.double(in_c), # Extra decrease parameter
			as.double(in_r_dec), # Decrease parameter
			as.double(in_r_inc), # Increase parameter
			as.double(in_eps1), # Enlargement parameter
			as.double(in_eps), # Stopping tolerance: Proximity measure
			as.double(in_crit_tol) # Stopping tolerance: Criticality tolerance
		)
		# Beta per k steps
		# Add row for intercept
		bperk <- matrix(res[[1]], nrow = ncol(x)+1, ncol = nrow(k))
		# Naming rows/cols/vector elements
		rownames(bperk) <- c(colnames(x),"intercept")
	
	}

	# If verb is high enough, print midpoints
	if(verb>=3){
		print(res)
		print(dim(x))
		print(length(y))
		print(dim(k))
		print(length(w))
	}

	# Beta per k steps
	# Cox regression doesn't have intercept
	if(family == "cox"){
		bperk <- matrix(res[[1]], nrow = ncol(x), ncol = nrow(k))
		rownames(bperk) <- colnames(x)
	# All other model families have intercept
	}else{
		bperk <- matrix(res[[1]], nrow = ncol(x)+1, ncol = nrow(k))
		# Fortran subroutine returns intercept as the last element, swap that to front
		rownames(bperk) <- c(colnames(x), "(Intercept)")
		bperk <- bperk[c(nrow(bperk), 1:(nrow(bperk)-1)),]
	}

	# Target function values per k steps
	fperk <- as.numeric(res[[2]])
	# Naming rows/cols/vector elements
	names(fperk) <- colnames(bperk) <- paste("k_", 1:nrow(k), sep="")
	
	if(verb>=2){
		print("fperk")
		print(fperk)
	}
	# Transpose so that coefficients are columns and k-steps are rows in bperk
	bperk <- t(bperk)

	if(verb>=2){
		print("bperk")
		print(dim(bperk))
		print(bperk)
	}
	
	# Kits picked per each k-step
	kperk <- t(apply(bperk, MARGIN=1, FUN=function(z) as.integer(!z==0)))
	if(verb>=3){
		print("kperk1")
		print(dim(kperk))
		print(kperk)
	}
	# If dealing with non-Cox regression, omit (Intercept) from indicator matrix
	if(!family == "cox" & ncol(k) == ncol(x)){
		kperk <- kperk[,-1]
	}
	
	# Indices as a named vector
	kperk <- as.list(apply(kperk %*% t(k), MARGIN=1, FUN=function(z) { which(!z==0) }))
	if(verb>=3){
		print("kperk2")
		print(dim(kperk))
		print(kperk)
	}
	
	# Kit costs per each k-step
	cperk <- unlist(lapply(kperk, FUN=function(z) { sum(w[z], na.rm=TRUE) }))

	if(verb>=2){
		print("cperk")
		print(cperk)
	}

	# Return the freshly built S4 model object
	obj <- new("oscar", 
		## Model fit results
		bperk = bperk[1:kmax,,drop=FALSE], # Beta coef per k-steps up to kmax
		fperk = fperk[1:kmax], # Target function values per k-steps up to kmax
		kperk = kperk[1:kmax], # Chosen kits per k-steps up to kmax
		cperk = cperk[1:kmax], # Total kit costs per each k-step up to kmax
		## Data slots
		x=as.matrix(x),	# Data matrix X
		y=as.matrix(y), # Response Y
		k=as.matrix(k), # Kit matrix K
		w=as.numeric(w),	# Vector of weights/costs for kits
		## Additional parameters
		family=as.character(family),	# Model family as a character string
		metric=as.character(metric),	# Model goodness metric
		start=start,	# Method for generating starting points
		kmax=kmax	# Max run k-step
	)

	if(verb>=2){
		print("obj template created successfully")
	}
	
	# Fit lm/glm/coxph/... models per each estimated set of beta coefs (function call depends on 'family')
	#try({
	#	# Model fits as a function of beta coefs
	#	obj@fits <- apply(obj@bperk, MARGIN=1, FUN=function(bs){
	#		# Debugging
	#		if(verb>=2) print("Performing obj@fits ...")
	#		if(family=="cox"){
	#			## Prefit a coxph-object
	#			survival::coxph(
	#				#as.formula(paste("survival::Surv(time=obj@y[,1],event=obj@y[,2]) ~",paste(colnames(obj@x),collapse='+'))), # Formula for response 'y' modeled using data matrix 'x' 
	#				as.formula("survival::Surv(time=obj@y[,1],event=obj@y[,2]) ~ ."), # Formula for response 'y' modeled using data matrix 'x' 
	#				data=data.frame(obj@x), # Use data matrix 'x'
	#				#data = data.frame(obj@x[,names(bs)]), # Use data matrix 'x'
	#				init = bs, # Use model coefficients obtained using the DBDC optimization 
	#				control = survival::coxph.control(iter.max=0) # Prevent iterator from deviating from prior model parameters
	#			)
	#		}else if(family %in% c("mse", "gaussian")){
	#			## Prefit a linear glm-object with gaussian error; use heavily stabbed .glm.fit.mod allowing maxit = 0
	#			stats::glm(
	#				#as.formula(paste("y ~",paste(colnames(obj@x),collapse='+')))
	#				as.formula("y ~ .")
	#				, data = data.frame(obj@x), start = bs, family = gaussian(link="identity"), method = oscar:::.glm.fit.mod
	#			)
	#		}else if(family=="logistic"){
	#			## Prefit a logistic glm-object with logistic link function; use heavily stabbed .glm.fit.mod allowing maxit = 0
	#			stats::glm(
	#				#as.formula(paste("y ~",paste(colnames(obj@x),collapse='+')))
	#				as.formula("y ~ .")
	#				, data = data.frame(obj@x), start = bs, family = binomial(link="logit"), method = oscar:::.glm.fit.mod
	#			)
	#			
	#			### Alternative function instead of glm (not tested here)
	#			#log.pred <- function(new.data){
	#			#	log.pred <-  bs%*%t(cbind(rep(1,nrow(new.data)),new.data))  ## NOTE! columnnames should be checked!
	#			#	log.pred <- exp(-log.pred)
	#			#	prob <- 1/(1+log.pred)
	#			#	pred <- lapply(prob,FUN=function(x){if(x>0.5){1}else{0}})  ## Cut-off 0.5 here
	#			#	return(pred)
	#			#}
	#
	#		}
	#	})
	#	# Extract corresponding model AICs as a function of k
	#	obj@AIC <- unlist(lapply(obj@fits, FUN=function(z) { stats::extractAIC(z)[2] }))
	#})
	#
	#if(verb>=2){
	#	print("fits-slot created successfully")
	#}
		
	# Calculate/extract model goodness metric at each k
	try({
		# Cox regression
		if(family=="cox"){
			# Use c-index as the goodness measure
			obj@goodness <- unlist(lapply(obj@fits, FUN=function(z) { z$concordance["concordance"] }))
		}else if(family %in% c("mse", "gaussian")){
			# Use mean squared error as the goodness measure
			obj@goodness <- unlist(lapply(obj@fits, FUN=function(z) { mean((y - predict.glm(z, type="response"))^2) }))
		}else if(family=="logistic"){
			# Use correct classification percent as the goodness measure
			# ROC-AUC
			if(metric=="auc"){
			
			# Accuracy
			}else if(metric=="accuracy"){
				#obj@goodness <- unlist(lapply(1:kmax, FUN=function(z) { sum(as.numeric(y == (predict.glm(z, type="response")>0.5)))/length(y) }))
			}
		}
	})
	
	# Return the new model object
	obj
}
