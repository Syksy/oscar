% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitS4.R
\name{oscar}
\alias{oscar}
\title{Main OSCAR fitting function}
\usage{
oscar(
  x,
  y,
  k,
  w,
  family = "cox",
  metric,
  solver = 1,
  verb = 1,
  print = 3,
  kmax,
  sanitize = TRUE,
  percentage = 1,
  in_selection = 1,
  storeX = TRUE,
  storeY = TRUE,
  control,
  ...
)
}
\arguments{
\item{x}{Data matrix 'x'}

\item{y}{Response vector/two-column matrix 'y' (see: family); number of rows equal to nrow(x)}

\item{k}{Integer (0/1) kit indicator matrix; number of columns equal to ncol(x), Default: Unit diagonal indicator matrix}

\item{w}{Kit cost weight vector w of length nrow(k), Default: Equal cost for all variables}

\item{family}{Model family, should be one of: 'cox', 'mse'/'gaussian', or 'logistic, Default: 'cox'}

\item{metric}{Goodness metric, Default(s): Concordance index for Cox, MSE for Gaussian, and AUC for logistic regression}

\item{solver}{Solver used in the optimization, should be  1/'DBDC' or 2/'LMBM', Default: 1.}

\item{verb}{Level of verbosity in R, Default: 1}

\item{print}{Level of verbosity in Fortran (may not be visible on all terminals); should be an integer between {range, range}, Default: 3}

\item{kmax}{Maximum k step tested, by default all k are tested from k to maximum dimensionality, Default: ncol(x)}

\item{sanitize}{Whether input column names should be cleaned of potentially problematic symbols, Default: TRUE}

\item{percentage}{Percentage of possible starting points used within range [0,1], Default: 1}

\item{in_selection}{Which starting point selection strategy is used (1, 2 or 3), Default: 1}

\item{storeX}{If data matrix X should be saved in the model object; turning this off might would help with memory, Default: TRUE}

\item{storeY}{If data response Y should be saved in the model object; turning this off might would help with memory, Default: TRUE}

\item{control}{Tuning parameters for the optimizers, see function oscar.control(), Default: see ?oscar.control}

\item{...}{Additional parameters}
}
\value{
Fitted oscar-object
}
\description{
This function fits an OSCAR model object to the provided training data with the desired model family.
}
\details{
OSCAR utilizes the L0-pseudonorm, also known as the best subset selection, and makes use of a DC-formulation of the discrete feature selection task into a continuous one. Then an appropriate optimization algorithm is utilized to find optima at different cardinalities (k). The S4 model objects 'oscar' can then be passed on to various down-stream functions, such as oscar.pareto, oscar.cv, and oscar.bs, along with their supporting visualization functions.
}
\examples{
if(interactive()){
  data(ex)
  fit <- oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family='cox')
  fit
}
}
\seealso{
\code{\link{oscar.cv}} \code{\link{oscar.bs}} \code{\link{oscar.pareto}} \code{\link{oscar.visu}} \code{\link{oscar.cv.visu}} \code{\link{oscar.bs.visu}} \code{\link{oscar.pareto.visu}} \code{\link{oscar.binplot}}
}
