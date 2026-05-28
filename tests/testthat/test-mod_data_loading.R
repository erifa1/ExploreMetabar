context("mod_data_loading: upload + filter + normalize pipeline")

# mod_data_loading is the most stateful module — it owns reactiveValues
# (r_values$phyobj_*), embeds datamods::filter_data_server, handles file
# uploads, and writes the consumer reactives onto the shared `r`.
# These tests stay at smoke-test level: they confirm the module mounts
# without erroring and exposes the expected interface.

test_that("mod_data_loading_ui returns a tag list", {
  ui <- mod_data_loading_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_data_loading_server: phyloseq_data loads the default fixture when no upload", {
  # The fixture is bundled with the installed package; if it's not installed
  # yet (e.g. running tests in source tree), skip rather than fail.
  fixture_path <- system.file("data_test", "robjects.Rdata",
                              package = "ExploreMetabar")
  skip_if(!nzchar(fixture_path) || !file.exists(fixture_path),
          "robjects.Rdata not installed; skipping")

  r <- mock_r()
  shiny::testServer(
    mod_data_loading_server,
    args = list(r = r),
    expr = {
      phy <- phyloseq_data()
      expect_s4_class(phy, "phyloseq")
    }
  )
})

test_that("mod_data_loading_server: sdat_initial returns a data.frame with sample.id", {
  fixture_path <- system.file("data_test", "robjects.Rdata",
                              package = "ExploreMetabar")
  skip_if(!nzchar(fixture_path) || !file.exists(fixture_path),
          "robjects.Rdata not installed; skipping")

  shiny::testServer(
    mod_data_loading_server,
    args = list(r = mock_r()),
    expr = {
      phyloseq_data()  # trigger initial load
      sdat <- sdat_initial()
      expect_s3_class(sdat, "data.frame")
      expect_true("sample.id" %in% colnames(sdat))
    }
  )
})
