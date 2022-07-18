
# OSCAR models with R

Optimal Subset CArdinality Regression (OSCAR) models with R

The R extension package `oscar` offers L0-pseudonorm based cardinality constrained regression models, with a user-friendly R interface tied to Fortran implemented optimizers via C. The R package comes with a suite of useful model building tools, such as visualizations, diagnotics, cross-validation, and bootstrapping functions. 

# Citation

Anni S. Halkola, Kaisa Joki, Tuomas Mirtti, Marko M. Mäkelä, Tero Aittokallio, Teemu D. Laajala. OSCAR: Optimal subset cardinality regression using the L0-pseudonorm with applications to prognostic modelling of prostate cancer. bioRxiv 2022.06.29.498064; doi: https://doi.org/10.1101/2022.06.29.498064

<a href="https://www.biorxiv.org/content/10.1101/2022.06.29.498064v1" class="uri">Direct bioarXiv link to the preprint</a>

# CRAN

The `oscar` package is available from the Central R Archive Network
(CRAN), and can be installed with:

    install.packages("oscar")

Canonical direct CRAN link:

<a href="https://CRAN.R-project.org/package=oscar" class="uri">https://CRAN.R-project.org/package=oscar</a>

# GitHub installation

While the CRAN installation may be most convenient for user
installation, the latest modifications and package release come through
GitHub. Therefore there may be a slight lag in delivering the latest
release to CRAN, and the user may wish to install the latest
GitHub-version by:

    install.packages("devtools")
    devtools::install_github("Syksy/oscar")

# Brief examples

Below are brief examples on `oscar` usage. More elaborate examples can
be found from within the package in the examples vignette.


    # Load the oscar-package

    library(oscar)

    # Load example dataset (consists of ex_X, ex_Y, ex_K and ex_c) for
    # Cox regression

    data(ex)

    # Test run, notice this uses all the data! Smaller test would be
    # feasible

    fit <- oscar::oscar(x = ex_X, y = ex_Y, k = ex_K, w = ex_c, family = "cox")

    # Show model results and other attributes

    fit

    # Visualize fit as a function of allowed kits

    oscar::oscar.visu(fit, y = c("target", "goodness"))

    # Example of cross-validation based k-selection

    cv <- oscar::oscar.cv(fit, fold = 5, seed = 123)  # 5-fold cross-validation, with fixed seed  

    # Visualize the cross-validation curve (highest point in c-index in
    # optimal generalizable k-value)

    oscar::oscar.cv.visu(cv)

    # Naive example using events in logistic regression (placeholder PoC,
    # end-point ought to be time-dependent rather)

    fit2 <- oscar::oscar(x = ex_X, y = ex_Y[, 2], k = ex_K, w = ex_c, family = "logistic")

    oscar::oscar.visu(fit2, y = c("target", "goodness"))

    # Example swiss fertility data

    data(swiss)

    fit_mse <- oscar::oscar(x = swiss[, -1], y = swiss[, 1], family = "gaussian")

    fit_mse
