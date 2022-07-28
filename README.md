
<!-- README.md is generated from README.Rmd. Please edit that file -->

# malt

<!-- badges: start -->
<!-- badges: end -->

This package implements the sampling algorithm: Metropolis Adjusted
Langevin Trajectories (Riou-Durand and Vogrinc 2022). Details available
on [arXiv](https://arxiv.org/abs/2202.13230). The MALT algorithm is a
robust extension of Hamiltonian Monte Carlo, for which robustness of
tuning is enabled by a positive choice of damping parameter (a.k.a
friction).

## Installation

You can install the development version of malt from
[R-universe](https://r-universe.dev/) with:

``` r
# Enable repository from lrioudurand
options(repos = c(
  lrioudurand = 'https://lrioudurand.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
# Download and install malt in R
install.packages('malt')
# Browse the malt manual pages
help(package = 'malt')
```
