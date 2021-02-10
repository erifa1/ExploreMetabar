
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ExploreMetabar
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4317188.svg)](https://doi.org/10.5281/zenodo.4317188)

<!-- badges: start -->

<!-- [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) -->

<!-- badges: end -->

ExploreMetabar is a shiny application used to explore metabarcoding data
(16S, ITS) with interactive plots, integrated statistical test, and
differential analysis.

## Installation

R3.6.3 is required. For Windows, [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and [git](https://git-scm.com/download/win) are required.

You can install the development version of ExploreMetabar from this
repository with:

``` r
install.packages("devtools")
devtools::install_git("https://forgemia.inra.fr/umrf/exploremetabar")
```

## To run Shiny app

``` r
library(ExploreMetabar)
ExploreMetabar::run_app()
```

## On ShinyApps.io

To test ExploreMetabar:
https://erifa1.shinyapps.io/exploremetabar/

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->
