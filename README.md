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

TO DO:
* Include Kaisa's preliminary Fortran implementation
* Extend current naive example
* Roxygenize R package documentation
* Make sure Fortran makefile is correct
* Testing for final modeling procedure (Ideally R unit testing folder)
* Extending R functionality to actually work
* Include a small R test data
* Many many things...
* ...

Compiling the code:
* R CMD build blasso
* R CMD INSTALL blasso_VERSION.tar.gz
* Alternatively from inside an on-going R-session run: > install.packages("blasso_VERSION.tar.gz", source=TRUE)
* R -> Start R session inside command prompt
* > library(blasso) # Load blasso-package
* > data(ex) # Load example dataset (consists of ex_X, ex_Y, ex_K and ex_c)
* > blassocox(x=ex_X, y=ex_Y, kits=ex_K, costs=ex_c) # Test run, notice this uses all the data! Smaller test would be feasible
* ...

Working with git:
* git status
* git pull
* git checkout branchname
* git add .
* git commit -m "Did X changes"
* git push
