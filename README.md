
# ExploreMetabar

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5245195.svg)](https://doi.org/10.5281/zenodo.5245195)


ExploreMetabar is a shiny application used to explore metabarcoding data
(16S, ITS) with interactive plots, integrated statistical tests, and
differential analysis.

## Online instances

ExploreMetabar is hosted thanks to [SK8 INRAE](https://sk8.inrae.fr/), [UCA mesocentre](https://mesocentre.uca.fr/) and [migale bioinformatics facility](https://migale.inrae.fr/):


**[https://explore-metabar.sk8.inrae.fr](https://explore-metabar.sk8.inrae.fr)**

**[https://shiny.mesocentre.uca.fr/app/exploremetabar](https://shiny.mesocentre.uca.fr/app/exploremetabar)**

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
sudo docker pull erifa1/exploremetabar:latest
sudo docker run -it -p 3838:3838 erifa1/exploremetabar:latest
```

## Citation

Etienne RIFA, & Sebastien Theil. (2020). ExploreMetabar: v1.0.1, https://forgemia.inra.fr/umrf/exploremetabar. Zenodo. https://doi.org/10.5281/zenodo.5245195
