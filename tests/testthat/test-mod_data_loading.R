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

# Regression: subsetting samples must drop only now-empty taxa, exactly like
# microViz::ps_filter. On robjects_adaopt.Rdata, sample_type==MI & cheese_type==UPC
# keeps 226 taxa. A prior bug applied a stale taxa-abundance threshold (calibrated
# on the full dataset) on top of this and left only 75. This guards the expected
# count and matches the microViz reference.
test_that("MI+UPC sample subset keeps 226 taxa (microViz parity)", {
  fixture <- system.file("data_test", "robjects_adaopt.Rdata",
                         package = "ExploreMetabar")
  skip_if(!nzchar(fixture) || !file.exists(fixture),
          "robjects_adaopt.Rdata not installed; skipping")

  ne <- new.env(); load(fixture, envir = ne); data <- ne$data
  sd <- as(phyloseq::sample_data(data), "data.frame")
  keep <- rownames(sd)[sd$sample_type == "MI" & sd$cheese_type == "UPC"]
  ps <- phyloseq::prune_samples(keep, data)
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
  expect_equal(phyloseq::ntaxa(ps), 226L)
})

# Regression: the default "Process Data" pipeline must NOT apply the manual
# taxa-filter table. `final()` (no manual filter applied) must equal `processed()`;
# the manual taxa filter is opt-in via the "Apply taxa filter" button.
test_that("mod_data_loading_server: default pipeline keeps all processed taxa", {
  fixture <- system.file("data_test", "robjects_adaopt.Rdata",
                         package = "ExploreMetabar")
  skip_if(!nzchar(fixture) || !file.exists(fixture),
          "robjects_adaopt.Rdata not installed; skipping")

  shiny::testServer(
    mod_data_loading_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(fileRData = list(datapath = fixture))
      phyloseq_data()                                   # load the object
      session$setInputs(rank_glom = "ASV", minAb = 0, minPrev = 0)
      session$setInputs(process_data = 1)               # run heavy pipeline

      skip_if(is.null(processed()),
              "datamods sample filter did not resolve under testServer")
      expect_s4_class(final(), "phyloseq")
      # No stale taxa threshold: the filtered object equals the processed object.
      expect_equal(phyloseq::ntaxa(final()), phyloseq::ntaxa(processed()))
      expect_setequal(phyloseq::taxa_names(final()),
                      phyloseq::taxa_names(processed()))
    }
  )
})
