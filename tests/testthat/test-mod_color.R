context("color helpers and mod_color contract")

# ---- pure helpers ----------------------------------------------------------

test_that("qualitative_palette_pool is non-empty and well-shaped", {
  pool <- qualitative_palette_pool()
  expect_s3_class(pool, "data.frame")
  expect_true(all(c("package", "palette", "length") %in% colnames(pool)))
  expect_gt(nrow(pool), 0)
})


test_that("assign_qualitative_colors is deterministic across calls", {
  mods <- list(group = c("a", "b", "c"), batch = c("x", "y"))
  a <- assign_qualitative_colors(mods)
  b <- assign_qualitative_colors(mods)
  expect_identical(a, b)
  expect_named(a, c("group", "batch"))
  expect_named(a$group, c("a", "b", "c"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}", a$group)))
})


test_that("assign_qualitative_colors uses distinct palettes per variable", {
  mods <- list(g1 = letters[1:3], g2 = letters[1:3])
  res <- assign_qualitative_colors(mods)
  # the two variables share modality names but should not share colors
  expect_false(identical(unname(res$g1), unname(res$g2)))
})


test_that("assign_qualitative_colors handles NA modalities", {
  mods <- list(g = c("a", NA, "b"))
  res <- assign_qualitative_colors(mods)
  expect_true("NA" %in% names(res$g))
  expect_equal(unname(res$g[["NA"]]), "#000000")
})


test_that("assign_qualitative_colors falls back deterministically for very high cardinality", {
  # cardinality larger than any palette in the pool forces the
  # high_cardinality_palette fallback
  big <- list(g = as.character(seq_len(500L)))
  a <- assign_qualitative_colors(big)
  b <- assign_qualitative_colors(big)
  expect_identical(a, b)
  expect_length(a$g, 500L)
  expect_true(all(grepl("^#", a$g)))
})


test_that("high_cardinality_palette stays qualitative (regression: no Spectral red ramp)", {
  hex <- high_cardinality_palette(100L)
  expect_length(hex, 100L)
  # at this cardinality we expect well over 50 distinct colors —
  # a sequential ramp (the previous Spectral fallback) gives barely
  # any visually distinct neighbors at the start of the vector
  expect_gt(length(unique(hex)), 80L)
  # adjacent colors should not be quasi-identical (sequential ramp
  # would have neighbors within a few RGB units of each other)
  rgb_mat <- t(grDevices::col2rgb(hex))
  step_dist <- sqrt(rowSums((rgb_mat[-1, ] - rgb_mat[-nrow(rgb_mat), ])^2))
  expect_gt(median(step_dist), 50)
})


test_that("assign_numeric_palettes returns one palette id per variable, cycling pool", {
  res <- assign_numeric_palettes(c("age", "ph", "salinity"))
  expect_named(res, c("age", "ph", "salinity"))
  expect_true(all(vapply(res, is.character, logical(1L))))
  expect_true(all(grepl("::", unlist(res, use.names = FALSE))))
})


test_that("assign_numeric_palettes returns an empty list when given no vars", {
  expect_identical(assign_numeric_palettes(character(0)), list())
})


test_that("read_color_csv returns NULL on missing columns", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  utils::write.csv(data.frame(foo = 1, bar = 2), tmp, row.names = FALSE)
  expect_null(read_color_csv(tmp))
})


test_that("apply_color_csv routes __numeric__ rows to numerics, others to colors", {
  mods <- list(group = c("a", "b"))
  df <- data.frame(
    variable = c("group", "group", "__numeric__"),
    modality = c("a",     "b",      "age"),
    color    = c("#111111", "#222222", "viridis::viridis"),
    stringsAsFactors = FALSE
  )
  res <- apply_color_csv(df, mods, accept_numeric = TRUE)
  expect_named(res$colors, "group")
  expect_named(res$numerics, "age")
  expect_equal(unname(res$colors$group["a"]), "#111111")
  expect_equal(res$numerics$age, "viridis::viridis")
  expect_equal(nrow(res$dropped), 0L)
})


test_that("apply_color_csv drops rows with unknown variable / invalid color / unknown palette", {
  mods <- list(known = c("a", "b"))
  df <- data.frame(
    variable = c("unknown",  "known",       "known",  "known",   "__numeric__"),
    modality = c("a",        "a",           "ghost",  "b",       "x"),
    color    = c("#111111",  "not_a_color", "#cccccc", "#222222", "fake::pal"),
    stringsAsFactors = FALSE
  )
  res <- apply_color_csv(df, mods, accept_numeric = TRUE)
  expect_equal(nrow(res$dropped), 4L)
  expect_true("unknown variable / rank"      %in% res$dropped$reason)
  expect_true("invalid color"                %in% res$dropped$reason)
  expect_true("unknown modality"             %in% res$dropped$reason)
  expect_true("unknown continuous palette"   %in% res$dropped$reason)
  # only the valid (known, b, #222222) row survived
  expect_equal(unname(res$colors$known["b"]), "#222222")
  expect_false("a" %in% names(res$colors$known))
})


test_that("apply_color_csv with accept_numeric = FALSE rejects __numeric__ rows", {
  mods <- list(Phylum = c("Firmicutes"))
  df <- data.frame(
    variable = c("Phylum", "__numeric__"),
    modality = c("Firmicutes", "age"),
    color    = c("#aabbcc", "viridis::viridis"),
    stringsAsFactors = FALSE
  )
  res <- apply_color_csv(df, mods, accept_numeric = FALSE)
  expect_length(res$numerics, 0L)
  expect_true(any(grepl("numeric row not accepted", res$dropped$reason)))
})


test_that("fill_missing_colors fills only the absent modalities and reports them", {
  loaded <- list(g = c(a = "#111111"))
  full   <- list(g = c("a", "b", "c"))
  res <- fill_missing_colors(loaded, full)
  expect_named(res$g, c("a", "b", "c"))
  expect_equal(unname(res$g["a"]), "#111111")
  missing <- attr(res, "missing")
  expect_setequal(missing$g, c("b", "c"))
})


test_that("build_factor_color_df produces long-format suitable for re-import", {
  fc <- list(group = c(a = "#111111", b = "#222222"))
  np <- list(age = "viridis::viridis")
  df <- build_factor_color_df(fc, np)
  expect_named(df, c("variable", "modality", "color"))
  expect_equal(nrow(df), 3L)
  expect_true("__numeric__" %in% df$variable)
  # round-trip: re-parsing the data.frame against the same ground-truth
  # restores both maps
  res <- apply_color_csv(df, list(group = c("a", "b")), accept_numeric = TRUE)
  expect_equal(res$colors$group, fc$group)
  expect_equal(res$numerics$age, np$age)
})


test_that("build_taxa_color_df handles an empty taxa_colors list", {
  df <- build_taxa_color_df(list())
  expect_equal(nrow(df), 0L)
  expect_named(df, c("variable", "modality", "color"))
})


# ---- mod_color server contract --------------------------------------------

test_that("mod_color_server writes factor_colors, numeric_palettes and taxa_colors onto r", {
  skip_if_not_installed("shiny")
  fixture_path <- system.file("data_test", "phy_test_numeric.rdata",
                              package = "ExploreMetabar")
  if (!nzchar(fixture_path) || !file.exists(fixture_path)) {
    skip("phy_test_numeric.rdata not installed")
  }
  env <- new.env()
  load(fixture_path, envir = env)
  phy <- env$data

  shiny::testServer(
    mod_color_server,
    args = list(r = shiny::reactiveValues(
      phyloseq_filtered = shiny::reactive(phy),
      sdat = shiny::reactive(as.data.frame(phyloseq::sample_data(phy))),
      var_list = shiny::reactive(
        colnames(as.data.frame(phyloseq::sample_data(phy)))
      ),
      factor_list = shiny::reactive({
        sdat <- as.data.frame(phyloseq::sample_data(phy))
        colnames(sdat)[!vapply(sdat, is.numeric, logical(1L))]
      })
    )),
    expr = {
      session$flushReact()
      sdat <- as.data.frame(phyloseq::sample_data(phy))
      num_vars <- colnames(sdat)[vapply(sdat, is.numeric, logical(1L))]
      fact_vars <- colnames(sdat)[!vapply(sdat, is.numeric, logical(1L))]

      fc <- r$factor_colors()
      np <- r$numeric_palettes()
      tc <- r$taxa_colors()

      # numerics never appear in factor_colors
      expect_true(all(names(fc) %in% fact_vars))
      expect_false(any(names(fc) %in% num_vars))

      # every numeric variable gets a continuous palette
      expect_setequal(names(np), num_vars)
      expect_true(all(grepl("::", unlist(np, use.names = FALSE))))

      # taxa colors: one entry per rank
      expect_setequal(names(tc), phyloseq::rank_names(phy))

      # determinism — recomputing returns identical output
      session$flushReact()
      expect_identical(fc, r$factor_colors())
      expect_identical(np, r$numeric_palettes())
      expect_identical(tc, r$taxa_colors())
    }
  )
})
