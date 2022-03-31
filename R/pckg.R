###
#
# Package information for roxygen and external functionality via C and Fortran
#
###

#' oscar: Optimal Subset Cardinality Regression
#'
#' OSCAR models utilize the L0-pseudonorm to select an optimal subset of features that generalizes linear regression models to a variety of families. Currently supported models include conventional Gaussian regression (family="mse" or family="gaussian"), Binomial/Logistic regression (family="logistic"), and Cox proportional hazards modeling (family="cox").
#' 
#' @section oscar functionality, list of functions
#'
#' @docType package
#' @name oscar
#' @importFrom grDevices colorRampPalette extendrange rainbow
#' @importFrom graphics abline arrows axis barplot box layout legend lines mtext par plot.new
#' @importFrom methods new
#' @importFrom stats sd smooth.spline
#' @importfrom utils str
#' @useDynLib oscar, .registration=TRUE
NULL
