context("mod_beta: beta diversity / ordinations")

test_that("mod_beta_ui returns a tag list", {
  ui <- mod_beta_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_beta_server: get_species_table returns a numeric matrix sized to the phyloseq", {
  shiny::testServer(
    mod_beta_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(beta_factor = "modalite")
      m <- get_species_table()
      expect_true(is.matrix(m) || is.data.frame(m))
      expect_equal(nrow(m), phyloseq::nsamples(r$phyloseq_filtered()))
    }
  )
})

test_that("mod_beta_server: ord() is gated on launch_beta button", {
  shiny::testServer(
    mod_beta_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        beta_factor = "modalite",
        ordination  = "PCOA",
        metrics     = "bray"
      )
      expect_error(ord(), class = "shiny.silent.error")
    }
  )
})

test_that("mod_beta_server: ord() returns a PCoA result when launched", {
  shiny::testServer(
    mod_beta_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        beta_factor  = "modalite",
        ordination   = "PCOA",
        metrics      = "bray",
        launch_beta  = 1
      )
      res <- ord()
      # vegan::capscale returns a 'capscale' object (subclass of 'cca')
      expect_true(inherits(res, c("capscale", "cca", "rda", "metaMDS")))
    }
  )
})

test_that("mod_beta_server: PERMANOVA reactives compute aligned results when launched", {
  shiny::testServer(
    mod_beta_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        beta_factor = "modalite",
        ordination  = "PCOA",
        metrics     = "bray",
        launch_beta = 1
      )

      # adonis2: a data.frame with the partitioned-variance columns. The
      # fixture has no NA, so every sample is retained (none dropped).
      ad <- get_adonis_res()
      expect_s3_class(ad, "data.frame")
      expect_true("Df" %in% names(ad))
      expect_gte(nrow(ad), 3L)                       # Model / Residual / Total

      # pairwise adonis: one row per pair of the 3 'modalite' levels.
      pw <- get_pairwise_res()
      expect_s3_class(pw, "data.frame")
      expect_true(all(c("pairs", "p.adjusted") %in% names(pw)))
      expect_equal(nrow(pw), choose(nlevels(local_metadata()[["modalite"]]), 2))

      # dispersion: a betadisper object whose groups align to the dist samples.
      bd <- get_dispersion_res()
      expect_s3_class(bd, "betadisper")
      expect_equal(length(bd$distances), phyloseq::nsamples(r$phyloseq_filtered()))
    }
  )
})

test_that("mod_beta_server: pairwise PERMANOVA rejects a numeric factor", {
  shiny::testServer(
    mod_beta_server,
    args = list(r = mock_r()),
    expr = {
      # 'pH' is a numeric metadata column in the fixture; pairwise comparisons
      # and dispersion only make sense for categorical groups.
      session$setInputs(
        beta_factor = "pH",
        ordination  = "PCOA",
        metrics     = "bray",
        launch_beta = 1
      )
      expect_true(isNumFactor())
      expect_error(get_pairwise_res(), class = "shiny.silent.error")
      expect_error(get_dispersion_res(), class = "shiny.silent.error")
    }
  )
})

test_that("mod_beta_server: get_constr_formula defaults to 'spe ~ 1' with no terms", {
  shiny::testServer(
    mod_beta_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(param_mode = "picker", constr_picker = character(0))
      expect_equal(get_constr_formula(), "spe ~ 1")
    }
  )
})

test_that("mod_beta_server: get_constr_formula assembles 'spe ~ a + b' from picker selection", {
  shiny::testServer(
    mod_beta_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        param_mode    = "picker",
        constr_picker = c("modalite", "Affinage")
      )
      expect_equal(get_constr_formula(), "spe ~ modalite + Affinage")
    }
  )
})
