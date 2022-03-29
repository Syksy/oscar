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
#' @param main Main title
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname oscar.visu
#' @export
oscar.visu <- function(
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

#' @title Visualize pareto-optimal front in respect to cross-validation and variable/kit costs
#' @description FUNCTION_DESCRIPTION
#' @param object PARAM_DESCRIPTION
#' @param cv PARAM_DESCRIPTION
#' @param k PARAM_DESCRIPTION
#' @param plot PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#'
#' @rdname pareto
#' @export
pareto <- function(
	object, # OSCAR model fit
	cv, # OSCAR cross-validation object
	k = 1:object@kmax, # k cardinalities to plot the cost (x-axis) versus cv metric (y-axis)
	plot = TRUE # Should a pareto front be plotted in addition to returning those points
){
	

	invisible(paretopoints)
}

#' @title Visualize bootstrapping of a fit oscar object
#'
#' @description This function visualizes bootstrapped model coefficients over multiple bootstrap runs as lines in a graph
#'
#' @param bs Bootstrapped 3-dimensional array for an oscar object as produced by oscar.bs
#' @param intercept Whether model intercept should be plotted also as a coefficient, Default: FALSE
#' @param add Should plot be added on top of an existing plot device
#'
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname oscar.bs.visu
#' @export
oscar.bs.visu <- function(
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
#'
#' @description This function plots the model performance as a function of cardinality for k-fold cross-validation. Performance metric depends on user choice and model family (i.e. lower MSE is good, higher C-index is good).
#'
#' @param cvs Matrix produced by oscar.cv; rows are cv-folds, cols are k-values
#' @param add Should plot be added on top of an existing plot device
#' @param main Main title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param ... Additional parameters passed on top the CV points
#'
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname oscar.cv.visu
#' @export
oscar.cv.visu <- function(
	cvs, # Matrix produced by oscar.cv; rows are cv-folds, cols are k-values
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
	points(x, means, type="l", ...)
	# Means as a function of k
	points(x, means, pch=16, col="red", ...)
	# Standard errors as a function of k
	arrows(x0=x, y0=means-sds, x1=x, y1=means+sds, code=3, angle=90, length=0.1)
}

#' @title Bootstrap visualization with boxplot, percentage of new additions
#'
#' @description This function plots as barplots as a function of k-cardinality in what proporties certain coefficients were chosen as non-zero over the bootstrap runs.
#' 
#' @param bs Bootstrapped 3-dimensional array for an oscar object as produced by oscar.bs
#' @param ... Additional parameters passed on to barplot
#'
#' @rdname oscar.bs.boxplot
#' @export
oscar.bs.boxplot <- function(
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
	barplot(t(howoften.new)[,nkits:1,drop=FALSE],horiz=TRUE,col=rainbow(38), ...)

}

#' @title Visualize binary indicator matrix optionally coupled with cross-validation performance
#'
#' @description TODO
#'
#' @param fit TODO
#' @param cv TODO
#' @param kmax TODO
#' @param collines TODO
#' @param rowlines TODO
#' @param cex.axis TODO
#' @param heights TODO
#' @param ... Additional parameters passed on to hamlet::hmap
#'
#' @importFrom hamlet hmap hmap.key
#' @export
oscar.binplot <- function(
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
#' @description TODO
#'
#' @param fit TODO
#' @param bs TODO
#' @param kmax TODO
#' @param cex.axis TODO
#' @param palet TODO
#' @param nbins TODO
#' @param Colv TODO
#' @param Rowv TODO
#' @param ... TODO
#'
#' @importFrom hamlet hmap hmap.key
#' @export
oscar.bs.plot <- function(
	fit, # Model fit object
	bs, # Bootstrapped data estimates if calculated
	kmax, # Maximum k to draw to; if missing, using kmax from fit@kmax
	cex.axis = 0.6, # Axis cex, for scaling
	palet = colorRampPalette(c("orange", "red", "black", "blue", "cyan"))(dim(bs)[3]), # color palette used in the heatmap, by default the number of bootstrapped datasets as separate colors
	nbins = 100, # Number of color bins
	Colv=NA, # Sorting of columns
	Rowv=NA, # Sorting of rows
	...
){
	if(!class(fit) %in% c("oscar")){
		stop("'fit' should be a fit oscar object")
	}
	if(missing(kmax)){
		kmax <- fit@kmax
	}
	par(mar=c(4,4,0,0))
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

	h <- hamlet::hmap(bmat, Colv=Colv, Rowv=Rowv, xlim=c(0.25, 1), ylim=c(0, 0.75), col=palet, nbins=nbins, namerows=FALSE, namecols=FALSE)
	title(xlab="Cardinality 'k'")
	axis(1, at=h$coltext$xseq[1:kmax], labels=1:kmax)
	axis(2, at=h$rowtext$yseq, labels=h$rowtext$rownam, cex.axis=cex.axis, las=1)
	hamlet::hmap.key(h)
}
					
#' @title Visualize cost and C-index with paretofront without oscar object
#'
#' @description TODO
#'
#' @param costs Vector including costs values
#' @param cindexes Vector including C-index values
#' @param nzeros Vector including the number of non-zeros of the model in each point
#' @param col Vector including desired colors for each number of non-zeros, Default : 'black'
#' @param ncolors Number of required colors, Default: max(nzeros)-min(nzeros)+1
#' @param mar.1st Marginals for the main plot, Default: c(4,4,2,0)
#' @param mar.2nd Marginals for the legend, Default: c(2.5,0,2,2)
#' @param monochrome If the plot should be with colors (F) or monochrome (T), Default: F
#'
#' @importFrom rPref psel
#'					
#' @rdname plot.smoothing.spline
#' @export
plot.pareto.costCI <- function(costs, cindexes, nzeros, col,ncolors,ylab="cindex",xlab="cost",pch=19,bg=col,cex=0.2,
                               ylim = c(min(cindexes),max(cindexes)), xlim=c(min(costs),max(costs)),
                               las=1, title="", width.of.layout =  c(2.2,0.8),height.of.layout= c(1,1),
                               mar.1st=c(4,4,2,0),mar.2nd =c(2.5,0,2,2),monochrome=F, cex.pareto=0.5, col.paretoline="grey"){
  if(monochrome==T){
     # If plotting with only one color, plot without legend
    if(missing(col)){  #Default color to black
      col <- "black"
    }
  }else{ # If plotting with colors, leaving area for legend
     # Also deal with color variables if missing
    layout(matrix(1:2,ncol=2), width = width.of.layout,height = height.of.layout)
  
  
  if(missing(ncolors)){
    ncolors <- max(nzeros)-min(nzeros)+1
  }
  if(missing(col)){
    col <- c()
    if(min(nzeros)==0){
      for(i in nzeros){
        col <- c(col,rainbow(ncolors)[i+1])
      }
    }else{
      for(i in nzeros){
        col <- c(col,rainbow(ncolors)[i])
      }
    }
  }
  
  }
  
  par(mar=mar.1st)
  
  if(monochrome==T){
    plot(costs,cindexes,ylab=ylab,xlab=xlab,col=col,cex=cex,pch=pch,bg=bg,las=las,ylim=ylim, xlim=xlim)
  }else{
    plot(costs,cindexes,col=col,ylab=ylab,xlab=xlab,cex=cex,pch=pch,las=las,ylim=ylim, xlim=xlim)
  }
  ## TDL: Packages shouldn't be loaded, but instead called with ::
  #library(rPref) ## !!!
  pareto.points <- rPref::psel(data.frame("cost"=costs,"CI"=cindexes),pref=rPref::low(cost)*rPref::high(CI))
  points(pareto.points[,1],pareto.points[,2],type='l',col=col.paretoline,lwd=1.5)
  
  if(monochrome==T){
    points(pareto.points[,1],pareto.points[,2],pch=pch,bg=bg,cex=cex.pareto,col=col)
    title(title)
  }else{ #Plot with colors and color legend
    points(pareto.points[,1],pareto.points[,2],pch=pch,cex=cex.pareto,col=col[as.numeric(rownames(pareto.points))])
    title(title)
    par(mar=mar.2nd)
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Cardinality',cex.main=1)
    rect(rep(0.2,ncolors),seq(0,1-1/ncolors,by=1/ncolors),rep(1.2,ncolors),seq(1/ncolors,1,by=1/ncolors),col=rainbow(ncolors))
    # Plot only every second tick label to avoid crowded look 
    tick.y <-seq(1/(2*ncolors),1-1/(2*ncolors),2/(ncolors))
    #browser()
    text(x=1.5, y =tick.y, labels = seq(min(nzeros),max(nzeros),1)[ceiling(tick.y*(ncolors))])
  }
  
}

#' @title Visualize smoothing spline with the first and second derivativests
#' @description FUNCTION_DESCRIPTION
#' @param x Values for x-axis
#' @param y Values for y-axis
#' @param ylab Title for y-axis
#' @param xlab Title for x-axis
#' @param a.title Main title for the scatter plot with spline
#' @param df Degrees of freedom, Default: 'ceiling(length(unique(x)))/2)'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#'
#' @rdname plot.smoothing.spline
#' @export
plot.smoothing.spline <- function(
	x, # X-values
	y, # Y-values
	ylab="y",  # Title of y-axis
	xlab="x", # Title of x-axis
	a.title="", # Title of the scatter plot with spline
	df=ceiling(length(unique(x)))/2)
{ # Degrees of freedom
  spline<-smooth.spline(x=x,y=y,df=df)  
  d2.spline <- predict(spline,seq(min(x),max(x),by=(max(x)-min(x))/100),deriv=2)
  d1.spline <- predict(spline,seq(min(x),max(x),by=(max(x)-min(x))/100),deriv=1)
  par(mfrow=c(1,3))
  plot(x,y,col=rainbow(length(x)),pch=19,cex=0.2,las=1,ylab=ylab,xlab=xlab)
  axis(side=2,col="red",las=1)
  title(a.title)
  lines(predict(spline.costCI,seq(100,720,by=10)),lwd=1)
  
  plot(d2.spline$x,d2.spline$y,type='l',lwd=2,xlab=xlab,ylab="D2(spline)",las=1)
  abline(h=0)
  title("b) Spline second derivative")
  plot(d1.spline$x,d1.spline$y,type='l',lwd=2,xlab=xlab,ylab="D(spline)",las=1)
  abline(h=0)
  title("c) Spline first derivative")
  
}

