####
#
# Visualization functions for blasso-objects and other relevant objects and variables
#
####

#' @title Target function value and total kit cost as a function of number of kits included
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION, Default: c("target", "cost", "goodness", "cv")
#' @param cols PARAM_DESCRIPTION, Default: c("red", "blue")
#' @param legend PARAM_DESCRIPTION, Default: 'top'
#' @param mtexts PARAM_DESCRIPTION, Default: TRUE
#' @param add PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname visu
#' @export
visu <- function(
	object,	# oscar-object (with corresponding slots available)
	## Options for plotting on the y-axes:
	# target: Target objective function value at each k step
	# cost: model kit cost at each k step
	# goodness: model goodness-of-fit measure at each k step
	# cv: cross-validated model generalization goodness-measure at each k step
	## Notice only 1st and 2nd element of the vector is used; if vector is of length 1, only first y-axis is used
	y = c("target", "cost", "goodness", "cv", "AIC"), 
	# Associated y-axis colors
	cols = c("red", "blue"),
	legend = "top", # Legend on top, FALSE/NA omits legend, otherwise it's used for placing the legend
	mtexts = TRUE, # Outer margin texts
	add = FALSE, # Should plot be added into an existing frame / plot
	main = ""
){
	if(!class(object) %in% "oscar") stop("'object' should be of class 'oscar'")
	par(las=2,  # All labels orthogonally to axes
		mar=c(7,4,1,4), # Inner margins
		oma=c(ifelse(mtexts, 2, 0), ifelse(mtexts, 2, 0), 0, ifelse(mtexts, 2, 0))) # Outer margins depend on additional labels with mtext
	#x <- 1:nrow(object@k)
	# Use kmax to truncate k-path
	x <- 1:object@kmax
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
		leg <- c(leg, paste("Model goodness-of-fit (", object@metric, ")"))		
	}else if(y[1]=="cv"){
		# TODO
		leg <- c(leg, "Cross-validated goodness-of-fit")		
	}else if(y[1]=="AIC"){
		y1 <- object@AIC
		leg <- c(leg, "AIC")
	}else{
		stop(paste("Invalid y[1] parameter (", y[1],"), should be one of: 'target', 'cost', 'goodness', 'cv', 'AIC'", sep=""))
	}
	if(!add){
		plot.window(xlim=c(1,length(x)), ylim=range(y1))
		axis(1, at=1:length(x), labels=x)
		axis(2, col.axis=cols[1])
		box(); title(main=main)
	}
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
			leg <- c(leg, paste("Model goodness-of-fit", object@metric, ")"))		
		}else if(y[2]=="cv"){
			# TODO
			leg <- c(leg, paste("Cross-validated goodness-of-fit", object@metric, ")"))		
		}else if(y[2]=="aic"){
			y2 <- object@aic
			leg <- c(leg, "AIC")
		}else{
			stop(paste("Invalid y[2] parameter (", y[2],"), should be one of: 'target', 'cost', 'goodness', 'cv'", sep=""))
		}
		if(!add){
			plot.window(xlim=c(1,length(x)), ylim=range(y2))
			axis(4, col.axis=cols[2])
			box()
			points(1:length(x), y2, pch=16, col=cols[2])
			points(1:length(x), y2, type="l", col=cols[2])
		}
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

#' @title Visualize bootstrapping of a fit oscar object
#' @description FUNCTION_DESCRIPTION
#' @param bs PARAM_DESCRIPTION
#' @param intercept PARAM_DESCRIPTION, Default: FALSE
#' @param add PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname bs
#' @export
bs.visu <- function(
	bs, # Bootstrap array as produced by bs.oscar
	intercept = FALSE, # Whether intercept coefficient ought to be plotted as well
	add = FALSE # Should plot be added into an existing frame / plot
){
	# Sanity checking
	if(!length(dim(bs))==3) stop("Parameter 'bs' ought to be a 3-dim array as produced by bs.oscar")
	# Remove intercept if needed
	if(!intercept & "(Intercept)" %in% dimnames(bs)[[2]]){
		bs <- bs[,-which("(Intercept)" == dimnames(bs)[[2]]),]
	}
	
	# Plot bootstrapped runs
	# Plot new or add
	if(!add){
		# If not adding to an existing plot, setup the graphics device
		plot.new()
		plot.window(xlim=range(1:dim(bs)[1]), ylim=extendrange(bs))
		box(); axis(1); axis(2)
		abline(h=0, lwd=2, col="grey")
		title(xlab="Cardinality 'k'", ylab="Beta coefficients", main="Bootstrapped coefficients")
	}
	# Iterate over coefficients
	for(i in 1:dim(bs)[2]){
		# Iterate over bootstraps
		for(b in 1:dim(bs)[3]){
			points(x=1:dim(bs)[1], y=bs[,i,b], col=i, lwd=2, type="l")
		}
	}
}

#' @title Visualize cross-validation as a function of k
#' @description FUNCTION_DESCRIPTION
#' @param cvs PARAM_DESCRIPTION
#' @param add PARAM_DESCRIPTION
#' @param main PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname cv
#' @export
cv.visu <- function(
	cvs, # Matrix produced by cv.oscar; rows are cv-folds, cols are k-values
	add = FALSE, # Should plot be added into an existing frame / plot
	main = "OSCAR cross-validation", # Main title
	xlab = "Cardinality 'k'",
	ylab = "CV performance",
	...
){
	# Compute statistics for the CV-curve
	means <- apply(cvs, MARGIN=2, FUN=mean)
	sds <- apply(cvs, MARGIN=2, FUN=sd)
	# x-coordinates
	x <- 1:ncol(cvs)
	# Plot new or add
	if(!add){
		# If not adding to an existing plot, setup the graphics device
		plot.new()
		plot.window(xlim=range(x), ylim=extendrange(c(means-sds, means+sds)))
		box(); axis(1); axis(2)
		title(main=main, xlab=xlab, ylab=ylab)
	}
	# Plotting
	points(x, means, type="l")
	# Means as a function of k
	points(x, means, pch=16, col="red")
	# Standard errors as a function of k
	arrows(x0=x, y0=means-sds, x1=x, y1=means+sds, code=3, angle=90, length=0.1)
}

#' @title Bootstrap visualization with boxplot, percentage of new additions
#'
#' @export
bs.boxplot <- function(
	bs, # Bootstrap array as produced by bs.oscar
	...
){
	nbootstrap <- dim(bs)[3]
	nft <- dim(bs)[2]
	nkits <- dim(bs)[1]
	howoften.new <- matrix(rep(0,nkits*nft),nrow=nkits)
	colnames(howoften.new) <- colnames(bs[,,1])
	for(k in 1:nkits){  #Kits
	  for(i in 1:nft){ #Measures
	    lkm <-0
	    valittu <- 0
	    for(j in 1:nbootstrap){  #Bootsrap iteration
	      if(bs[k,i,j]!=0){
		valittu <- valittu +1
		if(k>1&&bs[k-1,i,j]==0){
		  lkm <- lkm +1
		}
	      }
	    }
	    if(k==1){  ## If only 1 kit, no previous to compare, so overall percent is used
		howoften.new[k,i] <- valittu/nbootstrap
	    }else{   ## Compared to previous, how often each feature is chosen as the new variable
		howoften.new[k,i] <- lkm/nbootstrap}
	  }
	}
	barplot(t(howoften.new)[,nkits:1,drop=FALSE],horiz=TRUE,col=rainbow(38))

}

#' @title Visualize binary indicator matrix optionally coupled with cross-validation performance
#'
#' @export
binplot <- function(
	fit, # Model fit object
	cv, # Cross-validation performance if calculated
	kmax, # Maximum k to draw to; if missing, using kmax from fit@kmax
	collines = TRUE, # Draw vertical lines to bottom part
	rowlines = TRUE, # Draw horizontal lines to highlight variables
	cex.axis = 0.6, 
	heights=c(0.2, 0.8),
	... # Additional parameters passed on to hamlet::hmap
){
	if(!class(fit) %in% c("oscar")){
		stop("'fit' should be a fit oscar object")
	}
	if(missing(kmax)){
		kmax <- fit@kmax
	}
	
	if(!missing(cv)){
		# Extra CV annotation on top
		par(mar=c(0,4,1,1), las=1)
		# Set up paneling
		layout(matrix(c(1,2), nrow=2), heights=heights)
		oscar::cv.visu(cv[,1:kmax])
		axis(1, at=1:kmax)
	}

	# Bottom binarized indicator panel
	par(mar=c(4,4,1,1))
	bfit <- binarize(fit)[,1:kmax]
	plot.new()
	plot.window(xlim=c(0.2,0.8), ylim=c(0.2,0.8))
	
	h <- hamlet::hmap(oscar::binarize(fit), Colv=NA, Rowv=NA, toplim=0.8, bottomlim=0.2, col=c("lightgrey", "darkgrey"), nbins=2, namerows=FALSE, namecols=FALSE, add=TRUE)
	#title(ylab="Non-zero coefficients", xlab="Cardinality 'k'")
	title(xlab="Cardinality 'k'")
	axis(1, at=h$coltext$xseq[1:kmax], labels=1:kmax)
	axis(2, at=h$rowtext$yseq, labels=h$rowtext$rownam, cex.axis=cex.axis, las=1)
	# Line annotations
	if(collines) abline(v=h$coltext$xseq[1:kmax], col="grey")
	if(rowlines) abline(h=h$rowtext$yseq, col="grey")
	# Legend
	legend("topright", col=c("lightgrey", "darkgrey"), pch=15, legend=c("Excluded", "Included"), bg="white")
	box()
}

#' @title Bootstrap + cross-validation heatmap plot
#'
#' @export
bs.plot <- function(
	fit, # Model fit object
	#cv, # Cross-validation performance if calculated
	bs, # Bootstrapped data estimates if calculated
	kmax, # Maximum k to draw to; if missing, using kmax from fit@kmax
	cex.axis = 0.6, 	
	...
){
	if(!class(fit) %in% c("oscar")){
		stop("'fit' should be a fit oscar object")
	}
	if(missing(kmax)){
		kmax <- fit@kmax
	}
	par(xpd=NA)
	# Cross-validation
	#if(!missing(cv)){
	#	# Extra CV annotation on top
	#	par(mar=c(0,4,1,1), las=1)
	#	# Set up paneling
	#	layout(matrix(c(1,2), nrow=2), heights=heights)
	#	oscar::cv.visu(cv[,1:kmax])
	#	axis(1, at=1:kmax)
	#}
	# Bootstrapping
	if(!missing(bs)){
		bmat <- Reduce("+", lapply(1:dim(bs)[3], FUN=function(z) { !bs[,,z]==0 }))/dim(bs)[3]		
	}else{
		bmat <- oscar::binarize(fit)
	}
	# Order in which variables were first picked in the original fit
	varorder <- unique(unlist(apply(fit@bperk, MARGIN=1, FUN=function(z) which(!z==0))))
	bmat <- t(bmat[,varorder])

	h <- hamlet::hmap(bmat, Colv=NA, Rowv=NA, toplim=0.8, bottomlim=0.2, col=colorRampPalette(c("orange", "red", "black", "blue", "cyan"))(100), nbins=100, namerows=FALSE, namecols=FALSE)
	title(xlab="Cardinality 'k'")
	axis(1, at=h$coltext$xseq[1:kmax], labels=1:kmax)
	axis(2, at=h$rowtext$yseq, labels=h$rowtext$rownam, cex.axis=cex.axis, las=1)
	hamlet::hmap.key(h)
}
