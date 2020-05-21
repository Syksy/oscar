####
#
# S4-class OO for saving models and relevant parameters
# DBDC-fit kit-structured L0-penalized (Cox) model fits
#
####

#' S4-class for className
setClass("className",
	representation(
		# Set class representations for the S4-class @-slots
		## Results from Fortran
		betaperkits = "matrix",	# Rows: different kit k-values, columns: beta coefficients  in model fit
		#cputime = "numeric",	# Total CPU time for computations
		kitorder = "integer",	# Integer order in which kits were selected
		kitnames = "character",	# Names of the corresponding kits (rownames of the provided kit indicator matrix)
		fkits = "numeric",	# Target function values at each k-value
		rho = "numeric",	# Vector of utilized rho-values in the penalization
		## Provided data, saved in R
		x = "matrix",		# Original data matrix for fitting
		#y = "Surv",	# Response vector (two-column Surv-object for Cox modelling)
		k = "matrix",		# 
		w = "numeric",		# Kit value (cost ~ weights) vector, length ought to be same as nrow(k)
		## Extra
		info = "character"	# Additional error messages, warnings or such reported e.g. during model fitting
	),	
	prototype(
		# Prototype base model object
		## Results from Fortran
		betaperkits = matrix(NA, nrow=0, ncol=0),
		kitorder = NA_integer_,
		kitnames = NA_character_,
		fkits = NA_real_,
		rho = NA_real_,
		## Provided data, saved in R
		x = matrix(NA, nrow=0, ncol=0),
		#y = NA,
		k = matrix(NA, nrow=0, ncol=0),
		w = NA_real_,
		## Extra
		info = NA_character_
	)
)

# Function that checks whether the provided model object is valid
setValidity("className", .check_validity)

# Test whether all S4-slots are legitimate		
.check_validity <- function(
	object
){
	# Return TRUE only if all slots fulfill the validity criteria
	TRUE
}
		
		