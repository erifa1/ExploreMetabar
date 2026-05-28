# Test helpers shared across all module tests.
#
# Loads one of the bundled phyloseq fixtures from inst/data_test/ and returns
# a shiny::reactiveValues object that mirrors the shape `app_server.R` builds.
# Module servers can be tested with:
#
#   shiny::testServer(
#     mod_xxx_server,
#     args = list(r = mock_r()),
#     expr = { ... }
#   )

load_test_phyloseq <- function(fixture = "phy_test_numeric") {
  path <- system.file("data_test", paste0(fixture, ".rdata"),
                      package = "ExploreMetabar")
  if (!nzchar(path) || !file.exists(path)) {
    # Fallback for in-source testing (devtools::test() before install)
    path <- file.path("..", "..", "inst", "data_test",
                      paste0(fixture, ".rdata"))
  }
  if (!file.exists(path)) {
    # Try the .Rdata variant
    path <- sub("\\.rdata$", ".Rdata", path)
  }
  env <- new.env()
  load(path, envir = env)
  env$data
}


mock_r <- function(fixture = "phy_test_numeric") {
  phy <- load_test_phyloseq(fixture)
  # Force a *plain* data.frame — as.data.frame() on a sample_data object can
  # return an S4 sample_data, which then breaks tidyr::unite() downstream.
  sdat <- data.frame(phyloseq::sample_data(phy), check.names = FALSE,
                     stringsAsFactors = FALSE)
  num_vars <- colnames(sdat)[vapply(sdat, is.numeric, logical(1L))]
  fact_vars <- colnames(sdat)[!vapply(sdat, is.numeric, logical(1L))]
  ranks <- phyloseq::rank_names(phy)

  # Pre-build deterministic color maps so plotting modules don't choke.
  factor_colors <- lapply(stats::setNames(fact_vars, fact_vars), function(v) {
    mods <- unique(as.character(sdat[[v]]))
    stats::setNames(
      grDevices::hcl.colors(length(mods), palette = "Set2"),
      mods
    )
  })
  taxa_colors <- lapply(stats::setNames(ranks, ranks), function(rk) {
    taxa <- unique(as.character(phyloseq::tax_table(phy)[, rk]))
    taxa <- taxa[!is.na(taxa)]
    stats::setNames(
      grDevices::hcl.colors(length(taxa), palette = "Set3"),
      taxa
    )
  })
  numeric_palettes <- stats::setNames(
    as.list(rep("viridis::viridis", length(num_vars))),
    num_vars
  )

  shiny::reactiveValues(
    phyloseq_filtered      = shiny::reactive(phy),
    phyloseq_filtered_norm = shiny::reactive(phy),
    phyloseq_data          = shiny::reactive(phy),
    sdat                   = shiny::reactive(sdat),
    var_list               = shiny::reactive(colnames(sdat)),
    factor_list            = shiny::reactive(fact_vars),
    data_ready             = shiny::reactive(TRUE),
    rank_glom              = shiny::reactive("ASV"),
    norm_method            = shiny::reactive(0L),
    fdata                  = phy,
    factor_colors          = shiny::reactive(factor_colors),
    numeric_palettes       = shiny::reactive(numeric_palettes),
    taxa_colors            = shiny::reactive(taxa_colors),
    tabs                   = shiny::reactiveValues(tabselected = "data_loading")
  )
}
