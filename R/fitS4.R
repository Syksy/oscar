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
		solver = "character",  	# Name of the solver used in the optimization (DBDC = 1 or LMBM = 2)
		in_selection = "integer", # Starting point selection method used in the optimization
		percentage = "numeric", # Percentage of used starting points (and potential predictors per k with in_selection 3 and 4)
		call = "call"		# Function call
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
		solver = NA_character_,
		in_selection ) = NA_integer_,
	 	percentage = NA_real_,
		call = call("oscar")
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


#' @title Main OSCAR fitting function
#' @description This function fits an OSCAR model object to the provided training data with the desired model family.
#' @param x Data matrix 'x'
#' @param y Response vector/two-column matrix 'y' (see: family); number of rows equal to nrow(x)
#' @param k Integer (0/1) kit indicator matrix; number of columns equal to ncol(x), Default: Unit diagonal indicator matrix
#' @param w Kit cost weight vector w of length nrow(k), Default: Equal cost for all variables
#' @param family Model family, should be one of: 'cox', 'mse'/'gaussian', or 'logistic, Default: 'cox'
#' @param metric Goodness metric, Default(s): Concordance index for Cox, MSE for Gaussian, and AUC for logistic regression
#' @param solver Solver used in the optimization, should be  1/'DBDC' or 2/'LMBM', Default: 1.
#' @param verb Level of verbosity in R, Default: 1
#' @param print Level of verbosity in Fortran (may not be visible on all terminals); should be an integer between {range, range}, Default: 3
#' @param kmax Maximum k step tested, by default all k are tested from k to maximum dimensionality, Default: ncol(x)
#' @param sanitize Whether input column names should be cleaned of potentially problematic symbols, Default: TRUE
#' @param control Tuning parameters for the optimizers, see function oscar.control(), Default: see ?oscar.control
#' @param ... Additional parameters
#'
#' @return Fitted oscar-object
#'
#' @details OSCAR utilizes the L0-pseudonorm, also known as the best subset selection, and makes use of a DC-formulation of the discrete feature selection task into a continuous one. Then an appropriate optimization algorithm is utilized to find optima at different cardinalities (k). The S4 model objects 'oscar' can then be passed on to various down-stream functions, such as oscar.pareto, oscar.cv, and oscar.bs, along with their supporting visualization functions.
#' @examples 
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit
#' }
#' @rdname oscar
#' @seealso \code{\link{oscar.cv}} \code{\link{oscar.bs}} \code{\link{oscar.pareto}} \code{\link{oscar.visu}} \code{\link{oscar.cv.visu}} \code{\link{oscar.bs.visu}} \code{\link{oscar.pareto.visu}} \code{\link{oscar.binplot}}
#' @export 
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
	# Solver used in the optimization (default = DBDC)
	solver = 1,
	# Tuning parameters in R
	verb=1, # Level of R verbosity (1 = standard, 2 = debug level, 3 = excessive debugging, 0<= none)
	print=3,# Level of Fortran verbosity (-1 for tidy output, 3 for debugging level verbosity)
	kmax,    # Maximum tested k-values
	sanitize = TRUE,	# Whether column (i.e. variable) names are cleaned of potentially hazardous symbols
	# Tuning parameter generation for the optimizers DBDC and LMBM; default values are generated with oscar.control-function, and can be replaced in the list
	percentage = 1, # Which percentage of possible starting points is used (1=100%)
	in_selection = 1, # Which starting point selection strategy is used (1 or 2)
	control,
	# Additional parameters, passed on to oscar.control(...)
	...
){
	# Save function call
	Call <- match.call()

	####
	#
	# R-side sanity checks for parameters and input
	#
	####

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
		stop("User provided custom goodness metrics not yet supported.")
	}

	# Flip Cox event/time columns to correct column order
	if(family == "cox" & inherits(y, c("matrix", "array", "Surv"))){
		# For Cox, we expect events to be second column
		if(all(y[,1] %in% c(0,1))){
			tmp <- y[,1]
			y[,1] <- y[,2]
			y[,2] <- tmp
		}else if(!all(y[,2] %in% c(0,1))){
			stop("Second column for y is expected to only have values '0' (no event) or '1' (observed event)")
		}
	}	
	# Check that input matrices x and y have equal number of rows
	if(family=="cox" && nrow(y)!=nrow(x)){
		stop(paste("Number of observations in the response matrix y (",nrow(y),") is not equal to number of observations in the predictor matrix x (",nrow(x),"). Please check both."))
	}
	# If kit matrix is missing as input, assume that each variable is alone
	if(missing(k)){
		k <- matrix(0, nrow=ncol(x), ncol=ncol(x))
		diag(k) <- 1
		rownames(k) <- colnames(x)
	}
	# Check that input k is a matrix
	if(!is.matrix(k)){
		try(k <- as.matrix(k))
		if(inherits(k, "try-error")) stop("Error casting kit matrix 'k' into a matrix; should consist only of indicator values 0 and 1")
	}
	# Check that kit matrix k and predictor matrix x have equal number of columns (features)
	if(family=="cox" && ncol(k)!=ncol(x)){
		stop(paste("Number of columns in kit matrix k (",ncol(k),") should be equal to the amount of features (number of columns in x, ",ncol(x),"). Check that correct features are included."))
	}
	# Check that the kit matrix has only 0 or 1
	if(!all(k %in% c(0,1))){
		stop(paste("Values in the kit matrix k should be 0 or 1. Please check the kit matrix."))
	}
	# Check that there are no zero-columns in the kit matrix
	if(any(apply(k,MARGIN=2,FUN=sum)==0)){
		stop(paste("Zero column in the kit matrix k. Please check that each feature is in atleast one kit."))  ## For now only one kit per feature allowed.
	}
	# Check that each kit has at least one feature (not sure if necessary?))
	if(any(apply(k,MARGIN=1,FUN=sum)<1)){
		stop(paste("Some kit(s) don't have any features. Please check that each kit has at least one feature."))
	}	
	# Check that no duplicate kits aka kits with exactly the same features
	if(nrow(unique(k))!=nrow(k)){
		stop(paste("Some kits have exactly the same features. Please check that each kit has a different combination of features."))
	}
	# Print output for k structure
	if(verb>=2){
		print("Input k:")
		print(k)
	}
	# If user has defined kmax, use it to pass to the Fortran function; otherwise kmax is the maximum number of kits in the input data
	if(missing(kmax)){
		kmax <- nrow(k)
	# Sanity checking for user provided parameter
	}else if(inherits(kmax, "numeric")){
		kmax <- as.integer(kmax)
	}else if(!inherits(kmax, c("integer", "numeric"))){
		stop("Provided kmax parameter ought to be of type 'integer' or 'numeric' cast to an integer")
	}
	# Check allowed kmax limits
	if(kmax<=0){
		kmax <- 1
	}else if(kmax > nrow(k)){
		kmax <- nrow(k)
	}
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
	# Check if solver is one of the four options 1/'DBDC' or 2/'LMBM'
	if(!(solver %in% c("dbdc", "DBDC", "lmbm","LMBM", 1, 2))){
		stop(paste("Solver should be either 1 (or 'DBDC') or 2 (or 'LMBM')."))
	}
	# If character given, change into numerical to be used when calling C.
	if(solver %in% c("DBDC", "dbdc")){
		solver <- 1
	}
	if(solver %in% c("LMBM", "lmbm")){
		solver <- 2
	}
	## Sanitize column names, replacing '+' with 'plus', '-' with 'minus', and ' ', '(' and ')' with '_'
	if(sanitize){
		colnames(x) <- gsub("\\ |\\(|)", "_", gsub("\\+", "plus", gsub("\\-", "minus", colnames(x))))
	}
	
	
	####
	#
	# Optimizer specific parameters; partially overlapping with Fortran checks
	#
	####
	# Extracted from oscar.control while passing user custom tuning parameters
	if(missing(control)){
		control <- oscar.control(x=x, family=family, ...)
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
			as.integer(control$start), # Tuning parameter for starting values
			as.integer(kmax), # Tuning parameter for max k run
			# DBDC tuning parameters
			as.integer(control$in_mrounds), # The number of rounds in one main iteration 
			as.integer(control$in_mit), # The number of main iteration 
			as.integer(control$in_mrounds_esc), # The number of rounds in escape procedure
			as.integer(control$in_b1), # Bundle B1
			as.integer(control$in_b2), # Bundle B2
			as.integer(control$in_b), # Bundle B in escape procedure
			as.double(control$in_m), # Descent parameter
			as.double(control$in_m_clarke), # Descent parameter in escape procedure
			as.double(control$in_c), # Extra decrease parameter
			as.double(control$in_r_dec), # Decrease parameter
			as.double(control$in_r_inc), # Increase parameter
			as.double(control$in_eps1), # Enlargement parameter
			as.double(control$in_eps), # Stopping tolerance: Proximity measure
			as.double(control$in_crit_tol), # Stopping tolerance: Criticality tolerance
			as.integer(sum(k)), #Number of ones in the kit matrix
			as.integer(solver), # The solver used in optimization
			# LMBM tuning parameters
			as.integer(control$na), # Size of the bundle
			as.integer(control$mcu), # Upper limit for maximum number of stored corrections
			as.integer(control$mcinit), # Initial maximum number of stored corrections
			as.double(control$tolf), # Tolerance for change of function values
			as.double(control$tolf2), # Second tolerance for change of function values.
			as.double(control$tolg), # Tolerance for the first termination criterion
			as.double(control$tolg2), # Tolerance for the second termination criterion
			as.double(control$eta), # Distance measure parameter
			as.double(control$epsL), # Line search parameter
			as.double(percentage), # Percentage of starting points used
			as.integer(in_selection) # Starting point selection method
		)
		if(verb>=2){
			print(res)
		}
		# Beta per k steps
		#bperk <- matrix(res[[1]], nrow = ncol(x), ncol = kmax)
		# Naming rows/cols/vector elements
		#rownames(bperk) <- colnames(x)
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
			as.integer(control$start), # Tuning parameter for starting values
			as.integer(kmax), # Tuning parameter for max k run
			# DBDC tuning parameters
			as.integer(control$in_mrounds), # The number of rounds in one main iteration 
			as.integer(control$in_mit), # The number of main iteration 
			as.integer(control$in_mrounds_esc), # The number of rounds in escape procedure
			as.integer(control$in_b1), # Bundle B1
			as.integer(control$in_b2), # Bundle B2
			as.integer(control$in_b), # Bundle B in escape procedure
			as.double(control$in_m), # Descent parameter
			as.double(control$in_m_clarke), # Descent parameter in escape procedure
			as.double(control$in_c), # Extra decrease parameter
			as.double(control$in_r_dec), # Decrease parameter
			as.double(control$in_r_inc), # Increase parameter
			as.double(control$in_eps1), # Enlargement parameter
			as.double(control$in_eps), # Stopping tolerance: Proximity measure
			as.double(control$in_crit_tol), # Stopping tolerance: Criticality tolerance
			as.integer(sum(k)), #Number of ones in the kit matrix
			as.integer(solver), # The solver used in optimization
			# LMBM tuning parameters
			as.integer(control$na), # Size of the bundle
			as.integer(control$mcu), # Upper limit for maximum number of stored corrections
			as.integer(control$mcinit), # Initial maximum number of stored corrections
			as.double(control$tolf), # Tolerance for change of function values
			as.double(control$tolf2), # Second tolerance for change of function values.
			as.double(control$tolg), # Tolerance for the first termination criterion
			as.double(control$tolg2), # Tolerance for the second termination criterion
			as.double(control$eta), # Distance measure parameter
			as.double(control$epsL) # Line search parameter
		)
		# Beta per k steps
		# Add row for intercept
		#bperk <- matrix(res[[1]], nrow = ncol(x)+1, ncol = nrow(k))
		# Naming rows/cols/vector elements
		#rownames(bperk) <- c(colnames(x),"intercept")
		
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
			as.integer(control$start), # Tuning parameter for starting values
			as.integer(kmax), # Tuning parameter for max k run
			# DBDC tuning parameters
			as.integer(control$in_mrounds), # The number of rounds in one main iteration 
			as.integer(control$in_mit), # The number of main iteration 
			as.integer(control$in_mrounds_esc), # The number of rounds in escape procedure
			as.integer(control$in_b1), # Bundle B1
			as.integer(control$in_b2), # Bundle B2
			as.integer(control$in_b), # Bundle B in escape procedure
			as.double(control$in_m), # Descent parameter
			as.double(control$in_m_clarke), # Descent parameter in escape procedure
			as.double(control$in_c), # Extra decrease parameter
			as.double(control$in_r_dec), # Decrease parameter
			as.double(control$in_r_inc), # Increase parameter
			as.double(control$in_eps1), # Enlargement parameter
			as.double(control$in_eps), # Stopping tolerance: Proximity measure
			as.double(control$in_crit_tol), # Stopping tolerance: Criticality tolerance
			as.integer(sum(k)), #Number of ones in the kit matrix
			as.integer(solver), # The solver used in optimization
			# LMBM tuning parameters
			as.integer(control$na), # Size of the bundle
			as.integer(control$mcu), # Upper limit for maximum number of stored corrections
			as.integer(control$mcinit), # Initial maximum number of stored corrections
			as.double(control$tolf), # Tolerance for change of function values
			as.double(control$tolf2), # Second tolerance for change of function values.
			as.double(control$tolg), # Tolerance for the first termination criterion
			as.double(control$tolg2), # Tolerance for the second termination criterion
			as.double(control$eta), # Distance measure parameter
			as.double(control$epsL) # Line search parameter
		)
		# Beta per k steps
		# Add row for intercept
		#bperk <- matrix(res[[1]], nrow = ncol(x)+1, ncol = nrow(k))
		# Naming rows/cols/vector elements
		#rownames(bperk) <- c(colnames(x),"intercept")
	
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
		bperk <- matrix(res[[1]], nrow = ncol(x), ncol = kmax)
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
	names(fperk) <- colnames(bperk) <- paste("k_", 1:kmax, sep="")
	
	if(verb>=2){
		print("fperk")
		print(fperk)
	}
	# Transpose so that coefficients are columns and k-steps are rows in bperk
	bperk <- t(bperk)
	# Print mid-point for beta per kit k matrix
	if(verb>=2){
		print("bperk")
		print(dim(bperk))
		print(bperk)
	}
	## Get kperk from Fortran
	kperk <- t(matrix(res[[3]], nrow = nrow(k), ncol = kmax))
	if(verb>=3){
		print("kperk1")
		print(dim(kperk))
		print(kperk)
	}
	
	kperk<- as.list(apply(kperk,MARGIN=1,FUN=function(z){z[which(!z==0)]}))
	kperk <- kperk[1:kmax]  # Take only until kmax already here, since for i>kmax, kperk[[i]] could be problematic
	# If dealing with non-Cox regression, omit intercept from kit names
	for(i in 1:length(kperk)){
		names(kperk[[i]])<-rownames(k)[kperk[[i]]]
	}
	
	# Kit costs per each k-step
	cperk <- unlist(lapply(kperk, FUN=function(z) { sum(w[z], na.rm=TRUE) }))

	if(verb>=2){
		print("cperk")
		print(cperk)
	}

	# Set solver as character into oscar object
	if(solver == 1)
	{
		solver <- "DBDC"
	}
	if(solver == 2)
	{
		solver <- "LMBM"
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
		start=control$start,	# Method for generating starting points
		kmax=kmax,	# Max run k-step
		solver = solver,# Solver used in optimization (DBDC or LMBM)
		in_selection = in_selection, # Used starting point selection method
		percentage = percentage, # Percentage of starting points
		call = Call	# Used function call
	)

	if(verb>=2){
		print("obj template created successfully")
	}
	
	# Calculate/extract model goodness metric at each k
	try({
		# Cox regression
		if(family=="cox"){
			# Use c-index as the goodness measure
			obj@goodness <- unlist(lapply(1:kmax, FUN=function(z) { survival::coxph(Surv(time=y[,1], event=y[,2]) ~ x %*% t(bperk[z,,drop=FALSE]))$concordance["concordance"] }))
		}else if(family %in% c("mse", "gaussian")){
			# Use mean squared error as the goodness measure
			obj@goodness <- unlist(lapply(1:kmax, FUN=function(z) { mean((y - (cbind(1, x) %*% t(bperk[z,,drop=FALSE])))^2) }))
		}else if(family=="logistic"){
			# Use correct classification percent as the goodness measure
			# ROC-AUC
			if(metric=="auc"){
				# TODO
				stop("ROC-AUC yet to be implemented.")
			# Accuracy
			}else if(metric=="accuracy"){
				#obj@goodness <- unlist(lapply(1:kmax, FUN=function(z) { sum(as.numeric(y == (predict.glm(z, type="response")>0.5)))/length(y) }))
				stop("TODO for the generic case")
			}else{
				stop(paste("Invalid goodness-of-fit metric:", metric))
			}
		}
	})
	
	# Return the new model object
	obj
}

#' @title Control OSCAR optimizer parameters
#'
#' @description Fine-tuning the parameters available for the DBDC and LMBM optimizers. See oscar documentation for the optimization algorithms for further details.
#'
#' @param x Input data matrix 'x'; will be used for calculating various control parameter defaults.
#' @param family Model family; should be one of 'cox', 'logistic', or 'gaussian'/'mse'
#' @param start Starting point generation method, see vignettes for details; should be an integer between {range,range}, Default: 2
#' @param in_mrounds DBDC: The maximum number of rounds in one main iteration, Default: 5000
#' @param in_mit DBDC: The maximum number of main iterations, Default: 5000
#' @param in_mrounds_esc DBDC: The maximum number of rounds in escape procedure, Default: 5000
#' @param in_b1 DBDC: The size of bundle B1, Default: min(n_feat+5,1000)
#' @param in_b2 DBDC: The size of bundle B2, Default: 3
#' @param in_b DBDC: Bundle B in escape procedure, Default: 2*n_feat
#' @param in_m DBDC: The descent parameter in main iteration, Default: 0.01
#' @param in_m_clarke DBDC: The descent parameter in escape procedure, Default: 0.01
#' @param in_c DBDC: The extra decrease parameter in main iteration, Default: 0.1
#' @param in_r_dec DBDC: The decrease parameter in main iteration, Default: 0.75, 0.99, or larger depending on n_obs (thresholds 10, 300, and above)
#' @param in_r_inc DBDC: The increase parameter in main iteration, Default: 10^5
#' @param in_eps1 DBDC: The enlargement parameter, Default: 5*10^(-5)
#' @param in_eps DBDC: The stopping tolerance (proximity measure), Default: 10^(-6) if number of features is <= 50, otherwise 10^(-5)
#' @param in_crit_tol DBDC: The stopping tolerance (criticality tolerance), Default: 10^(-5)
#' @param na LMBM: Size of the bundle, Default: 4
#' @param mcu LMBM: Upper limit for maximum number of stored corrections, Default: 7
#' @param mcinit LMBM: Initial maximum number of stored corrections, Default: 7
#' @param tolf LMBM: Tolerance for change of function values, Default: 10^(-5)
#' @param tolf2 LMBM: Second tolerance for change of function values, Default: 10^4
#' @param tolg LMBM: Tolerance for the first termination criterion, Default: 10^(-5)
#' @param tolg2 LMBM: Tolerance for the second termination criterion, Default: same as 'tolg'
#' @param eta LMBM: Distance measure parameter (>0), Default: 0.5
#' @param epsL LMBM: Line search parameter (0 < epsL < 0.25), Default: 0.125
#'
#' @return A list of sanity checked parameter values for the OSCAR optimizers.
#'
#' @details This function sanity checks and provides reasonable DBDC ('Double Bundle method for nonsmooth DC optimization' as described in Joki et al. (2018) <doi:10.1137/16M1115733>) and LMBM ('Limited Memory Bundle Method for large-scale nonsmooth optimization' as presented in Haarala et al. (2004) <doi:10.1080/10556780410001689225>) optimization tuning parameters. User may override custom values, though sanity checks will prevent unreasonable values and replace them. The returned list of parameters can be provided for the 'control' parameter when fitting oscar-objects.
#'
#' @examples 
#' if(interactive()){
#'   oscar.control() # Return a list of default parameters
#' }
#' @rdname oscar.control
#' @export 
oscar.control <- function(
	# Required input from main function to define control parameters
	x,
	family,
	# Default parameter construction part
	start=2,# Strategy for choosing starting points at each n_k iteration
	# Tuning parameters for DBDC
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
	in_crit_tol = 10^(-5), # Stopping tolerance: Criticality tolerance
	# Tuning parameters for LMBM
	na = 4, # Size of the bundle in LMBM, na >= 2
	mcu = 7, # Upper limit for maximum number of stored corrections in LMBM, mcu >= 3
	mcinit = 7, # Initial maximum number of stored corrections in LMBM, mcu >= mcinit >= 3 
	tolf = 10^(-5), # Tolerance for change of function values in LMBM
	tolf2 = 10^4, # Second tolerance for change of function values.
           #                          - If tolf2 < 0 the the parameter and the corresponding termination criterion will be ignored. 
           #                          - If tolf2 = 0 the default value 1.0E+4 will be used
	tolg = 10^(-5), # Tolerance for the first termination criterion in LMBM
	tolg2 = tolg, # Tolerance for the second termination criterion in LMBM (default = tolg)
	eta = 0.5, # Distance measure parameter in LMBM, eta > 0
	epsL = 0.125 # Line search parameter in LMBM, 0 < epsL < 0.25
){
	# Required data input
	if(missing(x) || missing(family)){
		stop("In order to generate feasible default parameters, user has to provide data ('x') and the model family ('family') as parameters to this function")
	}
	#### CHECKING tuning parameters and setting defaults ####
	if(in_mrounds <= 0){ # The number of rounds in one main iteration 
		in_mrounds <- 5000
		warnings("Input in_mrounds should be >0. Default value 5000 is used instead.")
	}
	if(!is.integer(in_mrounds)){ warnings("Input in_mrounds should be an integer. The value is rounded to the nearest integer.")
		in_mrounds <- round(in_mrounds)
	}
	if(in_mit <= 0){ # The number of main iteration
		in_mit <- 5000
		warnings(paste("Input in_mit should be >0. Default value", in_mit, "is used instead."))
	}
	if(!is.integer(in_mit)){ 
		in_mit <- round(in_mit)
		warnings("Input in_mit should be an integer. The value is rounded to the nearest integer.")
	}
	if(in_mrounds_esc <= 0){ # The number of rounds in escape procedure
		in_mrounds_esc <- 5000
		warnings(paste("Input in_mrounds_esc should be >0. Default value", in_mrounds_esc, "is used instead."))
	}
	if(!is.integer(in_mrounds_esc)){ warnings("Input in_mrounds_esc should be an integer. The value is rounded to the nearest integer.")
		in_mrounds_esc <- round(in_mrounds_esc)
	}
	if(missing(in_b1)){ # Bundle B1
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal"))
		{
			user_n <- user_n +1
		}
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
		in_b1 <- round(in_b1)
	} 
	if(in_b2 <=0){ # Bundle B2
		in_b2 <- 3
		warnings("Input in_b2 should be >0. Default value 3 is used instead.")
	}
	if(!is.integer(in_b2)){ warnings("Input in_b2 should be an integer. The value is rounded to the nearest integer.")
		in_b2 <- round(in_b2)
	}
	if(missing(in_b) || in_b <=1){ # Bundle B in escape procedure
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal"))
		{
			user_n <- user_n +1
		}
		in_b <- 2*user_n
		warnings(paste("Input in_b should be >1. Default value ", in_b," is used instead.",sep=""))
	}
	if(!is.integer(in_b)){ 
		in_b <- round(in_b)
		warnings("Input in_b should be an integer. The value is rounded to the nearest integer.")
	}
	if(in_m <=0 | in_m>=1){ # Descent parameter
		in_m <- 0.2
		warnings(paste("Input in_m should be in the interval (0,1). Default value", in_m, "is used instead."))
	}
	if(in_m_clarke <=0 | in_m_clarke>=1){ # Descent parameter in escape procedure
		in_m_clarke <- 0.01
		warnings(paste("Input in_m_clarke should be in the interval (0,1). Default value", in_m_clarke, "is used instead."))
	}
	if(in_c <=0 | in_c>=1){ # Extra decrease parameter
		in_c <- 0.1
		warnings(paste("Input in_c should be in the interval (0,1). Default value", in_c, "is used instead."))
	}
	if(missing(in_r_dec)){ # Decrease parameter
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal"))
		{
			user_n <- user_n +1
		}
		if(user_n <10){
			in_r_dec <- 0.75
		}else if(user_n >=300){
			in_r_dec <- 0.99
		}else{
			in_r_dec <- trunc(user_n/(user_n+5)*100)/100  # Default is two first decimal from user_n/(user_n+5)
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
		in_r_inc <- 10^5
		warnings("Input in_r_inc should be >1. Default value 10^5 is used instead.")
	}
	if(in_eps1<=0){ # Enlargement parameter
		in_eps1 <- 5*10^(-5)
		warnings("Input in_eps1 should be >0. Default value 5*10^(-5) is used instead.")
	}
	if(missing(in_eps) || in_eps<=0){# Stopping tolerance: Proximity measure
		user_n <- ncol(x)
		if(family %in% c("mse","logistic","gaussian","normal")){
		user_n <- user_n +1}
		if(user_n <=50){
			in_eps <- 10^(-6)
		}else{
			in_eps <- 10^(-5)
		}
		warnings(paste("Input in_eps should be >0. Default value ", in_eps," is used instead.",sep=""))
	}
	if(in_crit_tol<=0){
		in_crit_tol <- 10^(-5)
		warnings(paste("Input in_crit_tol should be >0. Default value ", in_crit_tol," is used instead.",sep=""))
	}	
	
	##############################################
	#### CHECKING tuning parameters for LMBM  ####
	if(na <2){ # Size of the bundle
		na <- as.integer(4)
		warnings(paste("Input na should be >=2. Default value", na, "is used instead."))
	}
	if(!is.integer(na)){ 
		na <- as.integer(round(na))
		warnings("Input na should be an integer. The value is rounded to the nearest integer >=2.")
	}
	if(mcu <3){ # Upper limit for maximum number of stored correctionse
		mcu <- as.integer(7)
		warnings("Input mcu should be >=3. Default value 7 is used instead.")
	}
	if(!is.integer(mcu)){ warnings("Input mcu should be an integer. The value is rounded to the nearest integer >=2.")
		mcu <- as.integer(round(mcu))
	}
	if(mcinit <3 || mcinit >mcu){ # Initial maximum number of stored corrections
		if(mcu <7){
			mcinit <- mcu
		}else{
			mcinit <- 7
		}
		
		warnings(paste("Input mcinit should be <=mcu and >=3. Default value", mcinit, "used if mcu >=7, otherwise mcinit=mcu."))
	}
	if(!is.integer(mcinit)){ 
		mcinit <- round(mcinit)
		warnings("Input mcinit should be an integer. The value is rounded to the nearest integer >=2.")
	}
	if(eta <=0){ # Distance measure parameter
		eta <- 0.5
		warnings(paste("Input eta should be >0. Default value", eta, "is used instead."))
	}
	if(epsL <=0 || epsL>=0.25){ # Line search parameter
		epsL <- 0.125
		warnings(paste("Input epsl should be >0 and <0.25. Default value", epsL, "is used instead."))
	}
	# Return all tuning parameters
	list(
		start = start,	# Strategy for choosing starting points at each n_k iteration
		# Tuning parameters for DBDC
		in_mrounds = in_mrounds,	# The number of rounds in one main iteration 
		in_mit = in_mit,	# The number of main iteration 
		in_mrounds_esc = in_mrounds_esc,	# The number of rounds in escape procedure
		in_b1 = in_b1,	# Bundle B1, default min(n+5,1000) depends on the problem -> defined later
		in_b2 = in_b2,	# Bundle B2
		in_b = in_b,	# Bundle B in escape procedure, default 2n depends on the problem -> defined later
		in_m = in_m,	# Descent parameter
		in_m_clarke = in_m_clarke,	# Descent parameter in escape procedure
		in_c = in_c,	# Extra decrease parameter
		in_r_dec = in_r_dec,	# Decrease parameter, default depends on the problem -> defined later
		in_r_inc = in_r_inc,	# Increase parameter
		in_eps1 = in_eps1,	# Enlargement parameter
		in_eps = in_eps,	# Stopping tolerance: Proximity measure, default depends on the problem -> defined later
		in_crit_tol = in_crit_tol,	# Stopping tolerance: Criticality tolerance
		# Tuning parameters for LMBM
		na = na,	# Size of the bundle in LMBM, na >= 2
		mcu = mcu,	# Upper limit for maximum number of stored corrections in LMBM, mcu >= 3
		mcinit = mcinit,	# Initial maximum number of stored corrections in LMBM, mcu >= mcinit >= 3 
		tolf = tolf,	# Tolerance for change of function values in LMBM
		tolf2 = tolf2,	# Second tolerance for change of function values.
		tolg = tolg,	# Tolerance for the first termination criterion in LMBM
		tolg2 = tolg2,	# Tolerance for the second termination criterion in LMBM (default = tolg)
		eta = eta,	# Distance measure parameter in LMBM, eta > 0
		epsL = epsL	# Line search parameter in LMBM, 0 < epsL < 0.25 (default = 1.0E-4.),
	)
}
