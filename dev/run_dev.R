  # Set options here
  options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

  # Detach all loaded packages and clean your environment
  golem::detach_all_attached()
  # rm(list=ls(all.names = TRUE))

  # Document and reload your package
  golem::document_and_reload()

  # Run the application
  devtools::load_all(".")
  ExploreMetabar::run_app()

# one line test
options(golem.app.prod = FALSE);golem::detach_all_attached();golem::document_and_reload();devtools::load_all(".");ExploreMetabar::run_app()
  # # Deploy
  # options(repos = BiocManager::repositories()); getOption("repos")
  # rsconnect::deployApp("~/Repository/LRF/shiny_app/ExploreMetabar/")
  # rsconnect::appDependencies("~/Repository/LRF/shiny_app/ExploreMetabar/")


  # package     version       source
  # 1              AlgDesign       1.2.0         CRAN
  # 2          AnnotationDbi      1.48.0 Bioconductor
  # 3                     BH    1.72.0-3         CRAN
  # 4                Biobase      2.46.0 Bioconductor
  # 5           BiocGenerics      0.32.0 Bioconductor
  # 6           BiocParallel      1.20.1 Bioconductor
  # 7             Biostrings      2.54.0 Bioconductor
  # 8                    DBI       1.1.0         CRAN
  # 9                 DESeq2     1.27.31       github
  # 10                    DT        0.13         CRAN
  # 11          DelayedArray      0.12.2 Bioconductor
  # 12                    GA         3.2         CRAN
  # 13          GenomeInfoDb      1.22.1 Bioconductor
  # 14      GenomeInfoDbData       1.2.2 Bioconductor
  # 15         GenomicRanges      1.38.0 Bioconductor
  # 16                   IHW      1.14.0 Bioconductor
  # 17               IRanges      2.20.2 Bioconductor
  # 18            KernSmooth     2.23-17         CRAN
  # 19                  MASS    7.3-51.6         CRAN
  # 20                Matrix      1.2-18         CRAN
  # 21             PLNmodels 0.10.5-9000       github
  # 22                    R6       2.4.1         CRAN
  # 23          RColorBrewer       1.1-2         CRAN
  # 24                 RCurl    1.98-1.1         CRAN
  # 25               RSQLite       2.2.0         CRAN
  # 26                  Rcpp     1.0.4.6         CRAN
  # 27         RcppArmadillo 0.9.860.2.0         CRAN
  # 28              Rhdf5lib       1.8.0 Bioconductor
  # 29                 Rtsne        0.15         CRAN
  # 30             S4Vectors      0.24.3 Bioconductor
  # 31  SummarizedExperiment      1.16.1 Bioconductor
  # 32           VennDiagram      1.6.20         CRAN
  # 33             WikidataR       1.4.0         CRAN
  # 34             WikipediR       1.5.0         CRAN
  # 35                Wrench       1.4.0 Bioconductor
  # 36                   XML    3.99-0.3         CRAN
  # 37               XVector      0.26.0 Bioconductor
  # 38                  ade4      1.7-15         CRAN
  # 39             agricolae       1.3-2         CRAN
  # 40              annotate      1.64.0 Bioconductor
  # 41                   ape         5.3         CRAN
  # 42               askpass         1.1         CRAN
  # 43            assertthat       0.2.1         CRAN
  # 44               attempt       0.3.0         CRAN
  # 45             backports       1.1.6         CRAN
  # 46             base64enc       0.1-3         CRAN
  # 47            biomformat      1.14.0 Bioconductor
  # 48                   bit    1.1-15.2         CRAN
  # 49                 bit64       0.9-7         CRAN
  # 50                bitops       1.0-6         CRAN
  # 51                  blob       1.2.1         CRAN
  # 52                  bold       0.9.0         CRAN
  # 53                  brew       1.0-6         CRAN
  # 54                 broom       0.5.5         CRAN
  # 55               caTools      1.18.0         CRAN
  # 56                 callr       3.4.3         CRAN
  # 57                 class      7.3-17         CRAN
  # 58              classInt       0.4-3         CRAN
  # 59                   cli       2.0.2         CRAN
  # 60                 clipr       0.7.0         CRAN
  # 61               cluster       2.1.0         CRAN
  # 62             codetools      0.2-16         CRAN
  # 63            colorspace       1.4-1         CRAN
  # 64              combinat       0.0-8         CRAN
  # 65            commonmark         1.7         CRAN
  # 66                config         0.3         CRAN
  # 67              corrplot        0.84         CRAN
  # 68               cowplot       1.0.0         CRAN
  # 69                crayon       1.3.4         CRAN
  # 70             crosstalk     1.1.0.1         CRAN
  # 71                  crul       0.9.0         CRAN
  # 72                  curl         4.3         CRAN
  # 73            data.table      1.12.8         CRAN
  # 74                  desc       1.2.0         CRAN
  # 75                digest      0.6.25         CRAN
  # 76           dockerfiler       0.1.3         CRAN
  # 77                 dplyr       0.8.5         CRAN
  # 78                 e1071       1.7-3         CRAN
  # 79              ellipsis       0.3.0         CRAN
  # 80              evaluate        0.14         CRAN
  # 81                 fansi       0.4.1         CRAN
  # 82                farver       2.0.3         CRAN
  # 83               fastmap       1.0.1         CRAN
  # 84             fastmatch       1.1-0         CRAN
  # 85               fdrtool      1.2.15         CRAN
  # 86               forcats       0.5.0         CRAN
  # 87               foreach       1.5.0         CRAN
  # 88               formatR         1.7         CRAN
  # 89                    fs       1.4.1         CRAN
  # 90         futile.logger       1.4.3         CRAN
  # 91        futile.options       1.0.1         CRAN
  # 92                 gdata      2.18.0         CRAN
  # 93               gdtools       0.2.2         CRAN
  # 94            genefilter      1.68.0 Bioconductor
  # 95           geneplotter      1.64.0 Bioconductor
  # 96              generics       0.0.2         CRAN
  # 97             ggfittext       0.8.1         CRAN
  # 98               ggplot2       3.3.0         CRAN
  # 99               ggrepel       0.8.2         CRAN
  # 100                   gh       1.1.0         CRAN
  # 101                git2r      0.26.1         CRAN
  # 102           glassoFast         1.0         CRAN
  # 103               glmnet       3.0-2         CRAN
  # 104                 glue  1.4.0.9000       github
  # 105                golem       0.2.1         CRAN
  # 106               gplots       3.0.3         CRAN
  # 107            gridExtra         2.3         CRAN
  # 108               gtable       0.3.0         CRAN
  # 109               gtools       3.8.2         CRAN
  # 110                haven       2.2.0         CRAN
  # 111                 here         0.1         CRAN
  # 112               hexbin      1.28.1         CRAN
  # 113                highr         0.8         CRAN
  # 114                  hms       0.5.3         CRAN
  # 115               hoardr       0.5.2         CRAN
  # 116            htmltools       0.4.0         CRAN
  # 117          htmlwidgets       1.5.1         CRAN
  # 118             httpcode       0.3.0         CRAN
  # 119               httpuv       1.5.2         CRAN
  # 120                 httr       1.4.1         CRAN
  # 121               igraph       1.2.5         CRAN
  # 122                  ini       0.3.1         CRAN
  # 123              isoband       0.2.1         CRAN
  # 124            iterators      1.0.12         CRAN
  # 125             jsonlite       1.6.1         CRAN
  # 126                 klaR      0.6-15         CRAN
  # 127                knitr        1.28         CRAN
  # 128             labeling         0.3         CRAN
  # 129             labelled       2.2.2         CRAN
  # 130             lambda.r       1.2.4         CRAN
  # 131                later       1.0.0         CRAN
  # 132              lattice     0.20-41         CRAN
  # 133             lazyeval       0.2.2         CRAN
  # 134            lifecycle       0.2.0         CRAN
  # 135                limma      3.42.2 Bioconductor
  # 136               locfit     1.5-9.4         CRAN
  # 137           lpsymphony      1.14.0 Bioconductor
  # 138             magrittr         1.5         CRAN
  # 139             markdown         1.1         CRAN
  # 140          matrixStats      0.56.0         CRAN
  # 141              memoise       1.1.0         CRAN
  # 142            metacoder  0.3.3.9002       github
  # 143        metagenomeSeq      1.28.2 Bioconductor
  # 144                 mgcv      1.8-31         CRAN
  # 145           microbiome       1.8.0 Bioconductor
  # 146                 mime         0.9         CRAN
  # 147               miniUI     0.1.1.1         CRAN
  # 148             multtest      2.42.0 Bioconductor
  # 149              munsell       0.5.0         CRAN
  # 150              natserv       0.3.0         CRAN
  # 151                 nlme     3.1-147         CRAN
  # 152               nloptr     1.2.2.1         CRAN
  # 153              openssl       1.4.1         CRAN
  # 154              permute       0.9-5         CRAN
  # 155             phangorn       2.5.5         CRAN
  # 156             phyloseq      1.30.0 Bioconductor
  # 157            phylotate         1.3         CRAN
  # 158               pillar       1.4.3         CRAN
  # 159               pixmap      0.4-11         CRAN
  # 160             pkgbuild       1.0.6         CRAN
  # 161            pkgconfig       2.0.3         CRAN
  # 162              pkgload       1.0.2         CRAN
  # 163                plogr       0.2.0         CRAN
  # 164               plotly     4.9.2.1         CRAN
  # 165                 plyr       1.8.6         CRAN
  # 166               praise       1.0.0         CRAN
  # 167          prettyunits       1.1.1         CRAN
  # 168             processx       3.4.2         CRAN
  # 169             progress       1.2.2         CRAN
  # 170             promises       1.1.0         CRAN
  # 171                   ps       1.3.2         CRAN
  # 172                purrr       0.3.3         CRAN
  # 173             quadprog       1.5-8         CRAN
  # 174            questionr       0.7.0         CRAN
  # 175             rappdirs       0.3.1         CRAN
  # 176                readr       1.3.1         CRAN
  # 177             rematch2       2.1.1         CRAN
  # 178              remotes       2.1.1         CRAN
  # 179              rentrez       1.2.2         CRAN
  # 180              reshape       0.8.8         CRAN
  # 181             reshape2       1.4.4         CRAN
  # 182                rhdf5      2.30.1 Bioconductor
  # 183                ritis       0.8.0         CRAN
  # 184                rlang       0.4.5         CRAN
  # 185            rmarkdown         2.1         CRAN
  # 186                 rncl       0.8.4         CRAN
  # 187                 rotl      3.0.10         CRAN
  # 188             roxygen2       7.1.0         CRAN
  # 189            rprojroot       1.3-2         CRAN
  # 190             rredlist       0.6.0         CRAN
  # 191           rstudioapi        0.11         CRAN
  # 192                rvest       0.3.5         CRAN
  # 193               scales       1.1.0         CRAN
  # 194            segmented       1.1-0         CRAN
  # 195              selectr       0.4-2         CRAN
  # 196               seqinr       3.6-1         CRAN
  # 197               shades       1.4.0         CRAN
  # 198                shape       1.4.4         CRAN
  # 199                shiny     1.4.0.2         CRAN
  # 200       shinydashboard       0.7.1         CRAN
  # 201                 slam      0.1-47         CRAN
  # 202                 snow       0.4-3         CRAN
  # 203              solrium       1.1.4         CRAN
  # 204          sourcetools       0.1.7         CRAN
  # 205                   sp       1.4-1         CRAN
  # 206              stringi       1.4.6         CRAN
  # 207              stringr       1.4.0         CRAN
  # 208             survival      3.1-12         CRAN
  # 209              svglite       1.2.3         CRAN
  # 210                  sys         3.3         CRAN
  # 211          systemfonts       0.2.0         CRAN
  # 212                 taxa       0.3.3         CRAN
  # 213               taxize      0.9.94         CRAN
  # 214             testthat       2.3.2         CRAN
  # 215               tibble       3.0.0         CRAN
  # 216                tidyr       1.0.2         CRAN
  # 217           tidyselect       1.0.0         CRAN
  # 218              tinytex        0.21         CRAN
  # 219               traits       0.4.2         CRAN
  # 220            triebeard       0.3.0         CRAN
  # 221             urltools       1.7.3         CRAN
  # 222              usethis       1.6.0         CRAN
  # 223                 utf8       1.1.4         CRAN
  # 224                vctrs       0.3.1         CRAN
  # 225                vegan       2.5-6         CRAN
  # 226          viridisLite       0.3.0         CRAN
  # 227              whisker         0.4         CRAN
  # 228             wikitaxa       0.3.0         CRAN
  # 229                withr       2.2.0         CRAN
  # 230               worrms       0.4.0         CRAN
  # 231                 xfun        0.13         CRAN
  # 232                 xml2       1.3.1         CRAN
  # 233               xtable       1.8-4         CRAN
  # 234                 yaml       2.2.1         CRAN
  # 235             zlibbioc      1.32.0 Bioconductor
  # 236                  zoo       1.8-7         CRAN
