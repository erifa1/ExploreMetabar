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
