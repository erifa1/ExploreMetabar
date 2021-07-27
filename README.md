
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

On ubuntu 18.04 those libraries are needed:
```{r}
remotes::system_requirements(path = '.', os="ubuntu", os_release = "18.04")
 [1] "apt-get install -y software-properties-common" "add-apt-repository -y ppa:cran/libgit2"       
 [3] "apt-get update"                                "apt-get install -y make"                      
 [5] "apt-get install -y zlib1g-dev"                 "apt-get install -y libicu-dev"                
 [7] "apt-get install -y git"                        "apt-get install -y pandoc"                    
 [9] "apt-get install -y libxml2-dev"                "apt-get install -y libcurl4-openssl-dev"      
[11] "apt-get install -y libssl-dev"                 "apt-get install -y libgit2-dev"               
[13] "apt-get install -y libglpk-dev"                "apt-get install -y libgmp3-dev"               
[15] "apt-get install -y libjpeg-dev"                "apt-get install -y libpng-dev"
```


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
https://erifa1.shinyapps.io/exploremetabar/ (old version)

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->
