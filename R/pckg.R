###
#
# Package information for roxygen and external functionality via C and Fortran
#
###

#' oscar: Optimal Subset Cardinality Regression
#'
#' OSCAR models utilize the L0-pseudonorm to select an optimal subset of features that generalizes linear regression models to a variety of families. Currently supported models include conventional Gaussian regression (family="mse" or family="gaussian"), Binomial/Logistic regression (family="logistic"), and Cox proportional hazards modeling (family="cox").
#' 
#' @references Halkola AS, Joki K, Mirtti T, Mäkelä MM, Aittokallio T, Laajala TD (2023) OSCAR: Optimal subset cardinality regression using the L0-pseudonorm with applications to prognostic modelling of prostate cancer. PLoS Comput Biol 19(3): e1010333. \doi{10.1371/journal.pcbi.1010333}
#'
#' @docType package
#' @name oscar-package
#' @importFrom grDevices colorRampPalette extendrange rainbow
#' @importFrom graphics abline arrows axis barplot box layout legend lines mtext par plot.new plot.window points rect text title
#' @importFrom methods new show
#' @importFrom stats sd smooth.spline predict.glm
#' @importFrom utils str
#' @importFrom pROC auc roc
#' @useDynLib oscar, .registration=TRUE
NULL
