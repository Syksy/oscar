####
#
# Visualization functions for blasso-objects and other relevant objects and variables
#
####

#' Target function value and total kit cost as a function of number of kits included
#'
#' @export
funcvalcost <- function(
	blassocoxfit, # Fit object for blasso Cox model
	x = "kits", # Which x-axis to use; number of "kits" or ...
	legend = TRUE, # Whether to add the default legend
	mtexts = TRUE, # Whether default label texts should be added to bottom, left and right margins
	... # Additional parameters
){
	par(las=2,  # All labels orthogonally to axes
		mar=c(7,4,1,4), # Inner margins
		oma=c(ifelse(mtexts, 2, 0), ifelse(mtexts, 2, 0), 0, ifelse(mtexts, 2, 0))) # Outer margins depend on additional labels with mtext
	y1 <- attr(blassocoxfit, "fkits")
	y2 <- attr(blassocoxfit, "costkits")
	if(x=="kits"){
		x = paste(paste("Kit #",attr(blassocoxfit, "kitindices"), sep=""), paste("(",attr(blassocoxfit, "kitnames"),")",sep=""), sep="\n")	
	}else{
		stop("TODO")
	}
	# Two y-axes overlayed in a single graphics device, rigid first example
	plot.new()
	# First part (target function values)
	plot.window(xlim=c(1,length(x)), ylim=range(y1))
	axis(1, at=1:length(x), labels=x)
	axis(2, col.axis="red")
	box()
	points(1:length(x), y1, pch=16, col="red")
	points(1:length(x), y1, type="l", col="red")
	# Second part (cost accumulation)
	plot.window(xlim=c(1,length(x)), ylim=range(y2))
	axis(4, col.axis="blue")
	points(1:length(x), y2, pch=16, col="blue")
	points(1:length(x), y2, type="l", col="blue")
	if(legend){
		legend("top", col=c("red", "blue"), pch=16, lwd=1, legend=c("Target function value", "Cumulative kit cost"))
	}
	if(mtexts){
		mtext(side=1, text="Kits added to the model in optimal order", las=0, outer=TRUE)
		mtext(side=2, text="Target function value", las=0, outer=TRUE)
		mtext(side=4, text="Cumulative kit costs", las=0, outer=TRUE)
	}
}

