####
#
# Visualization functions for blasso-objects and other relevant objects and variables
#
####

#' Target function value and total kit cost as a function of number of kits included
#'
#' @export
visu <- function(
	object,	# casso-object (with corresponding slots available)
	legend = TRUE, # Legend on top
	mtexts = TRUE, # Outer margin texts
	gtexts = TRUE, # Goodness measure texts
	...
){
	if(!class(object) %in% "casso") stop("'object' should be of class 'casso'")
	par(las=2,  # All labels orthogonally to axes
		mar=c(7,4,1,4), # Inner margins
		oma=c(ifelse(mtexts, 2, 0), ifelse(mtexts, 2, 0), 0, ifelse(mtexts, 2, 0))) # Outer margins depend on additional labels with mtext
	y1 <- object@fperk
	y2 <- object@cperk
	y3 <- object@goodness
	x <- 1:nrow(object@k)
	# Two y-axes overlayed in a single graphics device, rigid first example
	plot.new()
	# Target function values
	plot.window(xlim=c(1,length(x)), ylim=range(y1))
	axis(1, at=1:length(x), labels=x)
	axis(2, col.axis="red")
	box()
	points(1:length(x), y1, pch=16, col="red")
	points(1:length(x), y1, type="l", col="red")
	# Model 
	plot.window(xlim=c(1,length(x)), ylim=range(y3))
	points(1:length(x), y3, pch=16, col="green")	
	points(1:length(x), y3, type="l", col="green")	
	text(x=1:length(x), y=y3, labels=y3, pos=3, col="green") # Write goodness of measure values as text
	# Model kit costs at each k
	plot.window(xlim=c(1,length(x)), ylim=range(y2))
	axis(4, col.axis="blue")
	points(1:length(x), y2, pch=16, col="blue")
	points(1:length(x), y2, type="l", col="blue")
	if(legend){
		legend("top", col=c("red", "green", "blue"), pch=16, lwd=1, legend=c("Target function value", "Model goodness-of-fit", "Cumulative kit cost"))
	}
	if(mtexts){
		mtext(side=1, text="K steps", las=0, outer=TRUE)
		mtext(side=2, text="Target function value", las=0, outer=TRUE)
		mtext(side=4, text="Cumulative kit costs", las=0, outer=TRUE)
	}
}





