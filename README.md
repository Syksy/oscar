# blasso
Budget-Optimized LASSO private repo (R-package, using Fortran through C)

Structure:

Folders
./man/ : All R-code manuals (*.Rd) go here
./src/ : All source code (*.c, *.f95) go here; e.g. Fortran implementation that is passed through C, including makefile
./R/   : All R-code goes here

Root files
.gitignore  : Ignored files such as pre-compiled files and/or files specific e.g. for R studio
DESCRIPTION : Mandatory R-package description file
NAMESPACE   : Mandatory R-package export for which functions are available in the package's namespace
README.md   : This file


TO DO:
- Include Kaisa's preliminary Fortran implementation
- Extend current naive example
- Make sure Fortran makefile is correct
- Testing
- Extending R functionality to actually work
- Include a small R test data
- Many many things...
