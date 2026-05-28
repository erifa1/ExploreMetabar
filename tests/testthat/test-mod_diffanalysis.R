context("mod_diffanalysis: differential abundance")

test_that("mod_diffanalysis_ui returns a tag list", {
  ui <- mod_diffanalysis_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_diffanalysis_server: local_physeq subsets to two conditions", {
  shiny::testServer(
    mod_diffanalysis_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        diff_factor = "modalite",
        Cond1       = "CL",
        Cond2       = "La"
      )
      phy <- local_physeq()
      expect_s4_class(phy, "phyloseq")
      # only samples from the two requested modalities should remain
      mods <- unique(as.character(
        phyloseq::sample_data(phy)[["modalite"]]
      ))
      expect_setequal(mods, c("CL", "La"))
    }
  )
})

test_that("mod_diffanalysis_server: deseqDA is gated on launch_diff", {
  shiny::testServer(
    mod_diffanalysis_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        diff_factor = "modalite",
        Cond1       = "CL",
        Cond2       = "La"
      )
      expect_error(deseqDA(), class = "shiny.silent.error")
    }
  )
})
