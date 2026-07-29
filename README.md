
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GEMini (Generating Equivalent Models)

<!-- badges: start -->

<!-- badges: end -->

GEMini presents a series of functions to help researchers generate and
evaluate equivalent models in structural equation modelling (SEM). By
providing a *lavaan* model specification (external data not required),
primary function `equiv_gen()` generates a list of equivalent models
according to rules laid out by Lee and Hershberger (1990). Supplemental
functions `param_summary()` and `equiv_trim()` allow users to summarise
the list of equivalent models generated and trim the list of models
respectively (see REDACTED, in prep. for details).

## Installation

You can install the development version of GEMini from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("SarcylY/GEMini")
```

or

``` r
# install.packages("devtools")
devtools::install_github("SarcylY/GEMini")
```

## References

\[REDACTED\]. (2026). *GEMini: An R Package for Generating Equivalent
Models in Structural Equation Modelling* \[Manuscript in preparation\].
