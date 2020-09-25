####
#
# Visualization functions for blasso-objects and other relevant objects and variables
#
####

#' Target function value and total kit cost as a function of number of kits included
#'
#' @param y
#' @param cols
#' @param legend
#' @param mtexts
#' @param ...
#'
#' @examples
#' data(ex)
#' fit <- casso(x=ex_X, y=ex_Y, k=ex_K, w=ex_c)
#' visu(fit, y=c("goodness", "cost") # Goodness-of-fit vs. cost of kits
#' visu(fit, y=c("target", "cost") # Target function value vs. cost of kits
#'
#' @export
visu <- function(
	object,	# casso-object (with corresponding slots available)
	## Options for plotting on the y-axes:
	# target: Target objective function value at each k step
	# cost: model kit cost at each k step
	# goodness: model goodness-of-fit measure at each k step
	# cv: cross-validated model generalization goodness-measure at each k step
	## Notice only 1st and 2nd element of the vector is used; if vector is of length 1, only first y-axis is used
	y = c("target", "cost", "goodness", "cv"), 
	# Associated y-axis colors
	cols = c("red", "blue"),
	legend = "top", # Legend on top, FALSE/NA omits legend, otherwise it's used for placing the legend
	mtexts = TRUE, # Outer margin texts
	...
){
	if(!class(object) %in% "casso") stop("'object' should be of class 'casso'")
	par(las=2,  # All labels orthogonally to axes
		mar=c(7,4,1,4), # Inner margins
		oma=c(ifelse(mtexts, 2, 0), ifelse(mtexts, 2, 0), 0, ifelse(mtexts, 2, 0))) # Outer margins depend on additional labels with mtext
	x <- 1:nrow(object@k)
	# Maximum of two y-axes overlayed in a single graphics device
	plot.new()

	leg <- c()
	# First y-axis options
	if(y[1]=="target"){
		y1 <- object@fperk	
		leg <- c(leg, "Target function value")
	}else if(y[1]=="cost"){
		y1 <- object@cperk
		leg <- c(leg, "Total kit cost")		
	}else if(y[1]=="goodness"){
		y1 <- object@goodness
		leg <- c(leg, "Model goodness-of-fit")		
	}else if(y[1]=="cv"){
		# TODO
		leg <- c(leg, "Cross-validated goodness-of-fit")		
	}else{
		stop(paste("Invalid y[1] parameter (", y[1],"), should be one of: 'target', 'cost', 'goodness', 'cv'", sep=""))
	}
	plot.window(xlim=c(1,length(x)), ylim=range(y1))
	axis(1, at=1:length(x), labels=x)
	axis(2, col.axis=cols[1])
	box()
	points(1:length(x), y1, pch=16, col=cols[1])
	points(1:length(x), y1, type="l", col=cols[1])
	# If length(y)>1, then plot a second y-axis
	if(length(y)>1){
		if(y[2]=="target"){
			y2 <- object@fperk	
			leg <- c(leg, "Target function value")
		}else if(y[2]=="cost"){
			y2 <- object@cperk
			leg <- c(leg, "Total kit cost")		
		}else if(y[2]=="goodness"){
			y2 <- object@goodness
			leg <- c(leg, "Model goodness-of-fit")		
		}else if(y[2]=="cv"){
			# TODO
			leg <- c(leg, "Cross-validated goodness-of-fit")		
		}else{
			stop(paste("Invalid y[2] parameter (", y[2],"), should be one of: 'target', 'cost', 'goodness', 'cv'", sep=""))
		}
		plot.window(xlim=c(1,length(x)), ylim=range(y2))
		axis(4, col.axis=cols[2])
		box()
		points(1:length(x), y2, pch=16, col=cols[2])
		points(1:length(x), y2, type="l", col=cols[2])
	}
	# Plot legend according to desired y
	if(!is.na(legend) & !legend==""){ # If model is included (NA or FALSE omits it)
		legend(legend, 
			col=cols[1:2], 
			pch=16, lwd=1, 
			legend=leg
		)
	}
	if(mtexts){
		mtext(side=1, text="K steps", las=0, outer=TRUE)
		mtext(side=2, text=leg[1], las=0, outer=TRUE)
		if(length(y)>1) mtext(side=4, text=leg[2], las=0, outer=TRUE)
	}
}


# Visualize bootstrap coefficient paths for a casso model object
bs.visu <- function(
	bs	# Bootstrapped list from bs.casso
){
	# Omit entries with try-errors
	bs <- bs[-which(unlist(lapply(bs, FUN=class))=="try-error")]
	
	
}

