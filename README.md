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

```{r example, warning = FALSE, message = FALSE}

# devtools::install_github("Syksy/oscar")

library(oscar) # Load the oscar-package  

data(ex) # Load example dataset (consists of ex_X, ex_Y, ex_K and ex_c) for Cox regression  

fit <- oscar::oscar(x=ex_X, y=ex_Y, k=ex_K, w=ex_c, family="cox") # Test run, notice this uses all the data! Smaller test would be feasible  

fit # Show model results and other attributes  

oscar::oscar.visu(fit, y=c("target", "goodness")) # Visualize fit as a function of allowed kits  

# Example of cross-validation based k-selection  

cv <- oscar::oscar.cv(fit, fold=5, seed=123) # 5-fold cross-validation, with fixed seed  

oscar::oscar.cv.visu(cv) # Visualize the cross-validation curve (highest point in c-index in optimal generalizable k-value)  

# Naive example using events in logistic regression (placeholder PoC, end-point ought to be time-dependent rather)  

fit2 <- oscar::oscar(x=ex_X, y=ex_Y[,2], k=ex_K, w=ex_c, family="logistic")  

oscar::oscar.visu(fit2, y=c("target", "goodness"))  

# Example swiss fertility data  

data(swiss)  

fit_mse <- oscar::oscar(x=swiss[,-1], y=swiss[,1], family="gaussian")  

fit_mse

```
