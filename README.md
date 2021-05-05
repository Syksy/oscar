---
output: github_document
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# OSCAR models with R
Optimal Subset CArdinality Regression (OSCAR) models with R

## Version plan / feature milestones:
* ...
* 0.x - Make sure there are no WARNINGs or NOTEs in 'R CMD check pckg.tar.gz --as-cran'
* 0.x - Roxygenize R package documentation
* 0.x - Use omp_lib for parallelization in Fortran (or a similar approach)
* 0.x - Testing for final modeling procedure (Ideally R unit testing folder)
* 0.x - Variants allowed for goodness measures for the three model families
* ...

```{r example, warning = FALSE, message = FALSE}

# devtools::install_github("Syksy/oscar")

library(oscar) # Load the oscar-package  

data(ex) # Load example dataset (consists of ex_X, ex_Y, ex_K and ex_c) for Cox regression  

fit <- oscar::oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family="cox") # Test run, notice this uses all the data! Smaller test would be feasible  

fit # Show model results and other attributes  

oscar::visu(fit, y=c("target", "goodness")) # Visualize fit as a function of allowed kits  

# Example of cross-validation based k-selection  

cv <- oscar::cv.oscar(fit, fold=5, seed=123) # 5-fold cross-validation, with fixed seed  

oscar::cv.visu(cv) # Visualize the cross-validation curve (highest point in c-index in optimal generalizable k-value)  

# Naive example using events in logistic regression (not good modelling! placeholder)  

fit2 <- oscar::oscar(x=ex_X, y=ex_Y[,2], k=ex_K, w=ex_c, family="logistic")  

oscar::visu(fit2, y=c("target", "goodness"))  

# Example swiss fertility data  

data(swiss)  

fit_mse <- oscar::oscar(x=swiss[,-1], y=swiss[,1], family="gaussian")  

fit_mse@fits  

```