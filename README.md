
<!-- README.md is generated from README.Rmd. Please edit that file -->
PLmixed
=======

The purpose of `PLmixed` is to extend the capabilities of `lme4` to allow factor structures (i.e., factor loadings and discrimination parameters) to be freely estimated. Thus, factor analysis and item response theory models with multiple hierarchical levels and/or crossed random effects can be estimated using code that requires little more input than that required by `lme4`. All of the strengths of `lme4`, including the ability to add (possibly random) covariates and an arbitrary number of crossed random effects, are encompassed within `PLmixed`. In fact, `PLmixed` uses `lme4` and `optim` to estimate the model using nested maximizations. Details of this approach can be found in Jeon and Rabe-Hesketh (2012). A manuscript documenting the use of `PLmixed` is currently in preparation.

Installation
------------

`PLmixed` can be installed from CRAN with:

``` r
install.packages("PLmixed")
```

<!-- You can install the developmental version of PLmixed from github with: -->
<!-- # ```{r gh-installation, eval = FALSE} -->
<!-- # install.packages("devtools") -->
<!-- devtools::install_github("njrockwood/PLmixed") -->
<!-- ``` -->
