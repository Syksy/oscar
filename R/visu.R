####
#
# Visualization functions for blasso-objects and other relevant objects and variables
#
####

#' @title Target function value and total kit cost as a function of number of kits included
#'
#' @description Plot oscar S4-object goodness-of-fit, kit costs, and similar performance metrics.
#'
#' @param fit Fitted oscar S4-class object
#' @param y Plotted y-axes supporting two simultaneous axes with different scales, Default: c("target", "cost", "goodness", "cv")
#' @param cols Colours for drawn lines, Default: c("red", "blue")
#' @param legend Location of legend or omission of legend with NA, Default: 'top'
#' @param mtexts Outer margin texts
#' @param add Should plot be added into an existing frame / plot
#' @param main Main title
#'
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   oscar.visu(fit, y=c("target", "cost"))
#'  }
#' }
#' @rdname oscar.visu
#' @export
oscar.visu <- function(
	fit,	# oscar-object (with corresponding slots available)
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
	if(!inherits(fit, "oscar")) stop("'fit' should be of class 'oscar'")
	opar <- par(no.readonly = TRUE)

	par(las=2,  # All labels orthogonally to axes
		mar=c(7,4,1,4), # Inner margins
		oma=c(ifelse(mtexts, 2, 0), ifelse(mtexts, 2, 0), 0, ifelse(mtexts, 2, 0))) # Outer margins depend on additional labels with mtext
	# Use kmax to truncate k-path
	x <- 1:fit@kmax
	# Maximum of two y-axes overlayed in a single graphics device
	plot.new()

	leg <- c()
	# First y-axis options
	if(y[1]=="target"){
		y1 <- fit@fperk	
		leg <- c(leg, "Target function value")
	}else if(y[1]=="cost"){
		y1 <- fit@cperk
		leg <- c(leg, "Total kit cost")		
	}else if(y[1]=="goodness"){
		y1 <- fit@goodness
		leg <- c(leg, paste("Model goodness-of-fit (", fit@metric, ")"))		
	}else if(y[1]=="cv"){
		# TODO
		leg <- c(leg, "Cross-validated goodness-of-fit")		
	}else if(y[1]=="AIC"){
		y1 <- fit@AIC
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
			y2 <- fit@fperk	
			leg <- c(leg, "Target function value")
		}else if(y[2]=="cost"){
			y2 <- fit@cperk
			leg <- c(leg, "Total kit cost")		
		}else if(y[2]=="goodness"){
			y2 <- fit@goodness
			leg <- c(leg, paste("Model goodness-of-fit", fit@metric, ")"))		
		}else if(y[2]=="cv"){
			# TODO
			leg <- c(leg, paste("Cross-validated goodness-of-fit", fit@metric, ")"))		
		}else if(y[2]=="aic"){
			y2 <- fit@aic
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
	par(opar)
}

#' @title Visualize bootstrapping of a fit oscar object
#'
#' @description This function visualizes bootstrapped model coefficients over multiple bootstrap runs as lines in a graph
#'
#' @param bs Bootstrapped 3-dimensional array for an oscar object as produced by oscar.bs
#' @param intercept Whether model intercept should be plotted also as a coefficient, Default: FALSE
#' @param add Should plot be added on top of an existing plot device
#'
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_bs <- oscar.bs(fit, bootstrap = 20, seed = 123)
#'   oscar.bs.visu(fit_bs)
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
	if(!length(dim(bs))==3) stop("Parameter 'bs' ought to be a 3-dim array as produced by oscar.bs")
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
#' @param cv Matrix produced by oscar.cv; rows are cv-folds, cols are k-values
#' @param add Should plot be added on top of an existing plot device
#' @param main Main title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param ... Additional parameters passed on top the CV points
#'
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_cv <- oscar.cv(fit, fold = 10, seed = 123)
#'   oscar.cv.visu(fit_cv)
#'  }
#' }
#' @rdname oscar.cv.visu
#' @export
oscar.cv.visu <- function(
	cv, # Matrix produced by oscar.cv; rows are cv-folds, cols are k-values
	add = FALSE, # Should plot be added into an existing frame / plot
	main = "OSCAR cross-validation", # Main title
	xlab = "Cardinality 'k'",
	ylab = "CV performance",
	...
){
	# Compute statistics for the CV-curve
	means <- apply(cv, MARGIN=2, FUN=mean)
	sds <- apply(cv, MARGIN=2, FUN=sd)
	# x-coordinates
	x <- 1:ncol(cv)
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
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_bs <- oscar.bs(fit, bootstrap = 20, seed = 123)
#'   oscar.bs.boxplot(fit_bs)
#'  }
#' }
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

#' @title Visualize binary indicator matrix optionally coupled with cross-validation performance for oscar models
#'
#' @description This visualization function makes use of the sparsified beta-coefficient matrix form as a function of cardinality. Optionally, user may showcase cross-validation performance alongside at the same cardinality values.
#'
#' @param fit Fitted oscar S4-class object
#' @param cv Matrix produced by oscar.cv; rows are cv-folds, cols are k-values
#' @param kmax Maximum cardinality 'k'
#' @param collines Should vertical lines be drawn to bottom part
#' @param rowlines Should horizontal lines be drawn to highlight variables
#' @param cex.axis Axis magnification
#' @param heights Paneling proportions as a numeric vector of length 2
#' @param ... Additional parameters passed on to hamlet::hmap
#'
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_cv <- oscar.cv(fit, fold = 10, seed = 123)
#'   oscar.binplot(fit=fit, cv=fit_cv)
#'  }
#' }
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
	if(!inherits(fit, "oscar")){
		stop("'fit' should be an oscar S4-object")
	}
	if(missing(kmax)){
		kmax <- fit@kmax
	}
	
	if(!missing(cv)){
		opar <- par(no.readonly = TRUE)
		# Extra CV annotation on top
		par(mar=c(0,4,1,1), las=1)
		# Set up paneling
		layout(matrix(c(1,2), nrow=2), heights=heights)
		oscar::oscar.cv.visu(cv[,1:kmax])
		axis(1, at=1:kmax)
		par(opar)
	}

	opar <- par(no.readonly = TRUE)
	# Bottom binarized indicator panel
	par(mar=c(4,4,1,1))
	bfit <- oscar::oscar.binarize(fit)[,1:kmax]
	plot.new()
	plot.window(xlim=c(0.2,0.8), ylim=c(0.2,0.8))
	par(opar)
	
	h <- hamlet::hmap(oscar::oscar.binarize(fit), Colv=NA, Rowv=NA, toplim=0.8, bottomlim=0.2, col=c("lightgrey", "darkgrey"), nbins=2, namerows=FALSE, namecols=FALSE, add=TRUE, ...)
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

#' @title Bootstrap heatmap plot for oscar models
#'
#' @description This function neatly plots a colourized proportion of variables chosen as a function of cardinalities over a multitude of bootstrap runs. This helps model diagnostics in assesssing variable importance.
#'
#' @param fit Fitted oscar S4-class object
#' @param bs Bootstrapped 3-dimensional array for an oscar object as produced by oscar.bs
#' @param kmax Maximum cardinality 'k'
#' @param cex.axis Axis magnification
#' @param palet Colour palette
#' @param nbins Number of bins (typically ought to be same as number of colours in the palette)
#' @param Colv Column re-ordering indices or a readily built dendrogram
#' @param Rowv Row re-ordering indices or a readily built dendrogram
#' @param ... Additional parameters passed on to the hamlet::hmap function
#'
#' @details Further heatmap parameters available from ?hmap
#'
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_bs <- oscar.bs(fit, bootstrap = 20, seed = 123)
#'   oscar.bs.plot(fit, fit_bs)
#'  }
#' }
#' @importFrom hamlet hmap hmap.key
#' @export
oscar.bs.plot <- function(
	fit, # Model fit object
	bs, # Bootstrapped data estimates if calculated
	kmax, # Maximum k to draw to; if missing, using kmax from fit@kmax
	cex.axis = 0.6, # Axis cex, for scaling
	palet = colorRampPalette(c("orange", "red", "black", "blue", "cyan"))(dim(bs)[3]), # color palette used in the heatmap, by default the number of bootstrapped datasets as separate colors
	nbins = dim(bs)[3], # Number of color bins
	Colv=NA, # Sorting of columns
	Rowv=NA, # Sorting of rows
	...
){
	if(!inherits(fit, "oscar")){ 
		stop("'fit' should be an oscar S4-object")
	}
	if(missing(kmax)){
		kmax <- fit@kmax
	}
	opar <- par(no.readonly = TRUE)
	par(mar=c(4,4,0,0))
	# Bootstrapping
	if(!missing(bs)){
		bmat <- Reduce("+", lapply(1:dim(bs)[3], FUN=function(z) { !bs[,,z]==0 }))/dim(bs)[3]		
	}else{
		bmat <- oscar::oscar.binarize(fit)
	}
	# Order in which variables were first picked in the original fit
	varorder <- unique(unlist(apply(fit@bperk, MARGIN=1, FUN=function(z) which(!z==0))))
	bmat <- t(bmat[,varorder])

	h <- hamlet::hmap(bmat, Colv=Colv, Rowv=Rowv, xlim=c(0.25, 1), ylim=c(0, 0.75), col=palet, nbins=nbins, namerows=FALSE, namecols=FALSE, ...)
	title(xlab="Cardinality 'k'")
	axis(1, at=h$coltext$xseq[1:kmax], labels=1:kmax)
	axis(2, at=h$rowtext$yseq, labels=h$rowtext$rownam, cex.axis=cex.axis, las=1)
	hamlet::hmap.key(h)
	par(opar)
}

#' Visualize oscar model pareto front
#'
#' Visualization function for showing the pareto front for cardinality 'k' and model goodness metric, either from goodness-of-fit or from cross-validation
#'
#' @param fit Fit oscar S4-object
#' @param cv A cross-validation matrix as produced by oscar.cv; if CV is not provided, then goodness-of-fit from fit object itself is used rather than cross-validation generalization metric
#' @param xval The x-axis to construct pareto front based on; by default 'cost' vector for features/kits, can also be 'cardinality'/'k'
#' @param weak If weak pareto-optimality is allowed; by default FALSE.
#' @param summarize Function that summarizes over cross-validation folds; by default, this is the mean over the k-folds.
#' @param add If the fit should be added on top of an existing plot; in that case leaving out labels etc. By default new plot is called.
#' @param ... Additional parameters provided for the plotting functions
#'
#' @examples 
#' \dontrun{
#' if(interactive()){
#'   data(ex)
#'   fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
#'   fit_cv <- oscar.cv(fit, fold = 10, seed = 123)
#'   par(mfrow=c(1,2))
#'   oscar.pareto.visu(fit=fit) # Model goodness-of-fit
#'   oscar.pareto.visu(fit=fit, cv=fit_cv) # Model cross-validation performance
#'  }
#' }
#' @export
oscar.pareto.visu <- function(
	fit,
	cv,
	xval = "cost",
	weak = FALSE,
	summarize = mean,
	add = FALSE,
	...
){
	# Use oscar.pareto-function to obtain the pareto front, and then focus on plotting
	if(missing(cv)){
		paretos <- oscar::oscar.pareto(fit = fit, xval = xval, weak = weak)
	}else{
		paretos <- oscar::oscar.pareto(fit = fit, cv = cv, xval = xval, weak = weak, summarize = summarize)
	}
	# If we don't add to an existing plot, an informative canvas is created
	if(!add){
		plot.new()
		plot.window(xlim=range(paretos[,1]), ylim=range(paretos[,2]))
		box(); axis(1); axis(2)
		if(xval == "cost"){
			title(xlab="Variable/kit cost")
		}else if(xval %in% c("cardinality", "k")){
			title(xlab="Cardinality 'k'")
		}
		if(!missing(cv)){
			title(ylab=paste("CV performance (", fit@metric, ")", sep=""))
		}else{
			title(ylab=paste("Model performance (", fit@metric, ")", sep=""))
		}
		title(main="Pareto front for oscar fit")
	}
	# Actual pareto front, with model fits in black and front with a red line
	points(x=paretos[,1], y=paretos[,2], pch=16)
	points(x=paretos[which(paretos[,3]),1], y=paretos[which(paretos[,3]),2], type="l", col="red", lwd=2, ...)
}

