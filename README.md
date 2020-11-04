# casso
Cardinality-constrained Absolute Subset Selection Optimized (CASSO) regression models

Checking in command terminal using (no ERRORs, WARNINGs, or NOTEs except CRAN first submission NOTE):
* R CMD build casso
* R CMD check casso_x.y.z.tar.gz --as-cran

Current R-package folder structure:

Folders

* ./man/  : All R-code manuals (\*.Rd) go here
* ./src/  : All source code (\*.c, \*.f95) goes here; e.g. Fortran implementation that is passed through C, including makefile
* ./R/    : All R-code (\*.R) goes here
* ./data/ : TO DO; example dataset goes here as a \*.RData
* ./inst/ : TO DO
* ...

Root files
* .gitignore  : Ignored files such as pre-compiled files and/or files specific e.g. for R studio
* DESCRIPTION : Mandatory R-package description file
* NAMESPACE   : Mandatory R-package export for which functions are available in the package's namespace
* README.md   : This file
* ...

Version plan:
* 0.4 - Implement cross-validation (cv) and/or bootstrap (bs) for model diagnostics
* 0.5 - User-friendliness in R-package; vignette, helper tools, casso S4-class applicable functions
* ...
* 0.x - Roxygenize R package documentation
* 0.x - Use omp_lib for parallelization in Fortran (or a similar approach)
* 0.x - Testing for final modeling procedure (Ideally R unit testing folder)
* 0.x - Variants allowed for goodness measures for the three model families
* ...

Compiling the code and running examples:
* > R CMD build casso
* > R CMD INSTALL casso_VERSION.tar.gz
* > Alternatively from inside an on-going R-session run: > install.packages("casso_VERSION.tar.gz", source=TRUE)
* > R -> Start R session inside command prompt

\> library(casso) # Load blasso-package  

\> data(ex) # Load example dataset (consists of ex_X, ex_Y, ex_K and ex_c) for Cox regression  

\> fit <- casso::casso(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family="cox") # Test run, notice this uses all the data! Smaller test would be feasible  

\> fit # Show model results and other attributes  

\> casso::visu(fit, y=c("target", "goodness")) # Visualize fit as a function of allowed kits  

\> \# Example of cross-validation based k-selection  

\> cv <- casso::cv.casso(fit, fold=5, seed=123) # 5-fold cross-validation, with fixed seed  

\> casso::cv.visu(cv) # Visualize the cross-validation curve (highest point in c-index in optimal generalizable k-value)  

\> \# Naive example using events in logistic regression (not good modelling! placeholder)  

\> fit2 <- casso::casso(x=ex_X, y=ex_Y[,2], k=ex_K, w=ex_c, family="logistic")  

\> casso::visu(fit2, y=c("target", "goodness"))  

\> \# Example swiss fertility data  

\> data(swiss)  

\> fit_mse <- casso::casso(x=swiss[,-1], y=swiss[,1], family="gaussian")  

\> fit_mse@fits  

...
