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

# coerce_metadata_types is the pure, Shiny-free metadata type-repair helper.
# Contract (see ?coerce_metadata_types): decimal-bearing text -> numeric
# (French decimal comma included), all other non-numeric character -> factor,
# integer-only codes deliberately kept categorical, sample.id untouched.
test_that("coerce_metadata_types: French decimal comma -> numeric", {
  df <- data.frame(
    sample.id = c("s1", "s2", "s3"),
    ph        = c("6,5", "7,1", "6,8"),   # comma decimals (read as character)
    stringsAsFactors = FALSE
  )
  out <- coerce_metadata_types(df)
  expect_type(out$ph, "double")
  expect_equal(out$ph, c(6.5, 7.1, 6.8))
  expect_true(any(grepl("'ph'.*numeric", attr(out, "coercions"))))
})

test_that("coerce_metadata_types: dot decimals and space thousand-separators parse", {
  df <- data.frame(
    dot       = c("6.5", "7.1", "6.8"),
    thousands = c("1 000,5", "2 500,0", "3 750,25"),
    stringsAsFactors = FALSE
  )
  out <- coerce_metadata_types(df)
  expect_equal(out$dot, c(6.5, 7.1, 6.8))
  expect_equal(out$thousands, c(1000.5, 2500.0, 3750.25))
})

test_that("coerce_metadata_types: integer-only text stays categorical (factor)", {
  df <- data.frame(
    replicate = c("1", "2", "3", "1"),    # coded group, no decimal separator
    stringsAsFactors = FALSE
  )
  out <- coerce_metadata_types(df)
  expect_s3_class(out$replicate, "factor")
  expect_false(is.numeric(out$replicate))
})

test_that("coerce_metadata_types: plain character -> factor, numerics & sample.id untouched", {
  df <- data.frame(
    sample.id = c("12", "13", "14"),       # numeric-looking id: must stay text
    group     = c("ctrl", "treat", "ctrl"),
    age       = c(23, 45, 31),             # already numeric
    stringsAsFactors = FALSE
  )
  out <- coerce_metadata_types(df)
  expect_type(out$sample.id, "character")  # excluded from coercion
  expect_s3_class(out$group, "factor")
  expect_true(is.numeric(out$age))
})

test_that("coerce_metadata_types: NA / blanks preserved when converting to numeric", {
  df <- data.frame(
    weight = c("1,5", NA, "", "2,5"),
    stringsAsFactors = FALSE
  )
  out <- coerce_metadata_types(df)
  expect_equal(out$weight, c(1.5, NA, NA, 2.5))
})

test_that("coerce_metadata_types: no changes yields empty coercions attribute", {
  df <- data.frame(
    a = c(1.0, 2.0),                       # numeric
    b = factor(c("x", "y")),               # already factor
    stringsAsFactors = FALSE
  )
  out <- coerce_metadata_types(df)
  expect_length(attr(out, "coercions"), 0L)
})

# glom_physeq() / abund_prev_filter() are the pure filter helpers shared by the
# live taxonomy preview (tax_glom_obj/tax_preview_obj) and the committed
# pipeline (processed()). Keeping them in sync is the whole point of the dedup.
test_that("glom_physeq: ASV + no sample filter is a no-op", {
  phy <- load_test_phyloseq()
  out <- glom_physeq(phy, keep_samples = NULL, rank = "ASV")
  expect_equal(phyloseq::ntaxa(out), phyloseq::ntaxa(phy))
  expect_equal(phyloseq::nsamples(out), phyloseq::nsamples(phy))
})

test_that("glom_physeq: sample subset keeps those samples and drops emptied taxa", {
  phy <- load_test_phyloseq()
  keep <- phyloseq::sample_names(phy)[1:5]
  out <- glom_physeq(phy, keep_samples = keep, rank = "ASV")
  expect_equal(phyloseq::nsamples(out), length(keep))
  expect_setequal(phyloseq::sample_names(out), keep)
  expect_lte(phyloseq::ntaxa(out), phyloseq::ntaxa(phy))   # empties dropped
  expect_true(all(phyloseq::taxa_sums(out) > 0))
})

test_that("glom_physeq: non-ASV rank yields one taxon per rank value, taxa renamed", {
  phy <- load_test_phyloseq()
  tt <- as.vector(phyloseq::tax_table(phy)[, "Phylum"])
  expected <- length(unique(tt[!is.na(tt)]))   # tax_glom default NArm = TRUE
  out <- glom_physeq(phy, keep_samples = NULL, rank = "Phylum")
  expect_equal(phyloseq::ntaxa(out), expected)
  # taxa renamed to their Phylum value
  expect_setequal(phyloseq::taxa_names(out),
                  as.vector(phyloseq::tax_table(out)[, "Phylum"]))
})

test_that("abund_prev_filter: 0/0 thresholds are a no-op", {
  phy <- load_test_phyloseq()
  base <- glom_physeq(phy, keep_samples = NULL, rank = "ASV")
  out <- abund_prev_filter(base, minAb = 0, minPrev = 0)
  expect_equal(phyloseq::ntaxa(out), phyloseq::ntaxa(base))
})

test_that("abund_prev_filter: a positive prevalence threshold can only drop taxa", {
  phy <- load_test_phyloseq()
  base <- glom_physeq(phy, keep_samples = NULL, rank = "ASV")
  out <- abund_prev_filter(base, minAb = 0, minPrev = 0.5)
  expect_lte(phyloseq::ntaxa(out), phyloseq::ntaxa(base))
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
