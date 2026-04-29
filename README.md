
# ExploreMetabar

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14610142.svg)](https://doi.org/10.5281/zenodo.4317187)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R ≥ 4.4.3](https://img.shields.io/badge/R-%E2%89%A54.4.3-blue.svg)](https://cran.r-project.org/)

Interactive Shiny application for exploring metabarcoding (16S, ITS) data.

## Online instance

ExploreMetabar is hosted thanks to [SK8 INRAE](https://sk8.inrae.fr/):

**[https://explore-metabar.sk8.inrae.fr](https://explore-metabar.sk8.inrae.fr)**

## Overview

ExploreMetabar takes a [phyloseq](https://joey711.github.io/phyloseq/) object as input and guides you through a complete microbiome analysis workflow: filter and normalize your data, then explore it through 10 independent analysis modules covering community composition, diversity, differential abundance, and multivariate statistics. All results are downloadable as tables, plots, or R objects.

## Analysis modules

| Module | Key methods |
|---|---|
| **Data Loading** | Phyloseq filtering, normalization (TSS, CLR, VST, Hellinger), metadata editing, color scheme upload |
| **Composition** | Stacked bar plots by taxonomic rank, top-N taxa, group merging |
| **Alpha Diversity** | Observed, Chao1, ACE, Shannon, Simpson, InvSimpson · ANOVA/Tukey · linear regression |
| **Beta Diversity** | PCoA, NMDS, dbRDA · Bray-Curtis, Jaccard, UniFrac · PERMANOVA, pairwise Adonis, betadisper, envfit |
| **Taxa Boxplots** | Kruskal-Wallis, pairwise Wilcoxon, Pearson/Spearman/Kendall · FDR correction |
| **Heatmap** | Hierarchical clustering · TSS, CLR, VST, Hellinger, log10 normalization · metadata annotations |
| **Differential Analysis** | DESeq2, metagenomeSeq, MetaCoder heat trees · consensus table across methods |
| **ASV Venn** | Venn/nVennR diagrams · alluvial shared-taxon charts |
| **Clustering** | Ward/complete/average linkage · silhouette k-estimation · indicator species (indicspecies) |
| **MixOmics sPLS-DA** | 3-stage workflow: Initial → Tune (cross-validation) → Final · individuals, loadings, biplot |

## Input data

ExploreMetabar expects an **RData file containing a phyloseq object named `data`**. The object must include:

- OTU/ASV abundance table (required)
- Sample metadata (required)
- Taxonomy table (required)
- Phylogenetic tree (optional — enables UniFrac distances)
- Reference sequences (optional — enables FASTA exports)

## Installation

R 4.4.3 or higher is required.

**Linux (recommended)**

On Ubuntu 20.04, install system dependencies first:

```bash
apt-get update && apt-get install -y git-core libcurl4-openssl-dev libgit2-dev \
  libglpk-dev libgmp-dev libicu-dev libpng-dev libssl-dev libxml2-dev make \
  pandoc pandoc-citeproc zlib1g-dev libtiff-dev libjpeg-dev libbz2-dev \
  libgmp3-dev software-properties-common
```

**Windows**

[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and [git](https://git-scm.com/download/win) are required.

**Install from repository**

```r
install.packages("renv")
options(renv.config.gitlab.host = "https://forge.inrae.fr")

renv::install("gitlab::umrf/exploremetabar@master")
```

**Run the app**

```r
library(ExploreMetabar)
ExploreMetabar::run_app()
```

## Docker (older versions)

```bash
sudo docker pull erifa1/exploremetabar:latest
sudo docker run -it -p 3838:3838 erifa1/exploremetabar:latest
```

## Citation

Etienne RIFA, & Sebastien Theil. (2025). ExploreMetabar: v3.0, https://forge.inrae.fr/umrf/exploremetabar. Zenodo. https://doi.org/10.5281/zenodo.4317187

## License

MIT © Etienne Rifa, Sebastien Theil
