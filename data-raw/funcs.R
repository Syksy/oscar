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
  lines(predict(spline,seq(100,720,by=10)),lwd=1)
  
  plot(d2.spline$x,d2.spline$y,type='l',lwd=2,xlab=xlab,ylab="D2(spline)",las=1)
  abline(h=0)
  title("b) Spline second derivative")
  plot(d1.spline$x,d1.spline$y,type='l',lwd=2,xlab=xlab,ylab="D(spline)",las=1)
  abline(h=0)
  title("c) Spline first derivative")
  
}

