
# ExploreMetabar

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4317188.svg)](https://doi.org/10.5281/zenodo.4317188)


ExploreMetabar is a shiny application used to explore metabarcoding data
(16S, ITS) with interactive plots, integrated statistical tests, and
differential analysis.

## Online instance 

ExploreMetabar is hosted on [migale bioinformatics facility](https://migale.inrae.fr/) and available at this link:
 
**[https://shiny.migale.inrae.fr/app/exploremetabar](https://shiny.migale.inrae.fr/app/exploremetabar)**


## Installation

R3.6.3 or upper is required.


* **Linux (recommended)**

On ubuntu 18.04 those libraries are needed:

```bash
apt-get update && apt-get install -y  git-core libcurl4-openssl-dev libgit2-dev libglpk-dev libgmp-dev libicu-dev libpng-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev libtiff-dev libjpeg-dev libbz2-dev libgmp3-dev software-properties-common
```

* Windows

[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and [git](https://git-scm.com/download/win) are required.

You can install the development version of ExploreMetabar from this
repository with:

``` r
install.packages("devtools")
devtools::install_git("https://forgemia.inra.fr/umrf/exploremetabar")
```

### To run Shiny app

``` r
library(ExploreMetabar)
ExploreMetabar::run_app()
```

## Docker

To install ExploreMetabar via docker environment:

```bash
sudo docker pull erifa1/exploremetabar:v1.002
sudo docker run -it -p 3838:3838 erifa1/exploremetabar:v1.002
```

## Citation 

Etienne RIFA, & Sebastien Theil. (2020). ExploreMetabar: v1.0.0, https://forgemia.inra.fr/umrf/exploremetabar. Zenodo. https://doi.org/10.5281/zenodo.4317188






