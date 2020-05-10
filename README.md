# blasso
Budget-Optimized LASSO private repo (R-package, using Fortran through C)

Checking in command terminal using (no ERRORs, WARNINGs, or NOTEs except CRAN first submission NOTE):
* R CMD build blasso
* R CMD check blasso_x.y.z.tar.gz --as-cran

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
* 0.x - Roxygenize R package documentation
* 0.x - Testing for final modeling procedure (Ideally R unit testing folder)
* 0.x - Include a small R test data (currently ex_X and ex_Y are already in ePCR, maybe include only ex_K and ex_c in blasso)
* 0.2 - Basic diagnostics and results interpretation of the modeling results; visualizations, fitting of coxph-models in R based on the results from Fortran, prediction, model coefficient extraction etc
* 0.3 - Additional DBDC fine-tuning; starting-point robustness, speed, 
* 0.4 - An S4-class for blasso-objects instead of saving helper variables as attributes
* 0.5 - User-friendliness in R-package; vignette, helper tools
* 0.6 - Comfortable to other similar packages, e.g. predict.blasso-function to with type = "coefficients" and with adjustable "k" parameter, type = "response" for particular Cox-model predictions etc
* ...

Compiling the code:
* R CMD build blasso
* R CMD INSTALL blasso_VERSION.tar.gz
* Alternatively from inside an on-going R-session run: > install.packages("blasso_VERSION.tar.gz", source=TRUE)
* R -> Start R session inside command prompt
* > library(blasso) # Load blasso-package
* > data(ex) # Load example dataset (consists of ex_X, ex_Y, ex_K and ex_c)
* > bc <- blassocox::blassocox(x=ex_X, y=ex_Y, kits=ex_K, costs=ex_c) # Test run, notice this uses all the data! Smaller test would be feasible
* > bc # Show model results and other attributes
* ...

Working with git:
* git status
* git pull
* git checkout branchname
* git add .
* git commit -m "Did X changes"
* git push
