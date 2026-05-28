context("mod_alpha: alpha diversity")

test_that("mod_alpha_ui returns a tag list", {
  ui <- mod_alpha_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_alpha_server: alpha1() is gated on launch_alpha button", {
  shiny::testServer(
    mod_alpha_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(Fact1 = "modalite", checkbox1 = TRUE, metrics = "Shannon")
      expect_error(alpha1(), class = "shiny.silent.error")
    }
  )
})

test_that("mod_alpha_server: alpha1() returns a richness table when launched", {
  shiny::testServer(
    mod_alpha_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        Fact1        = "modalite",
        checkbox1    = TRUE,
        metrics      = "Shannon",
        launch_alpha = 1
      )
      LL <- alpha1()
      expect_named(LL, c("alphatab", "data"))
      expect_s3_class(LL$alphatab, "data.frame")
      expect_true(all(c("Observed", "Shannon", "Chao1") %in% colnames(LL$alphatab)))
      # one row per sample
      expect_equal(nrow(LL$alphatab), phyloseq::nsamples(r$phyloseq_filtered()))
    }
  )
})

test_that("mod_alpha_server: alphagrp_table summarises by categorical factor", {
  shiny::testServer(
    mod_alpha_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        Fact1        = "modalite",
        checkbox1    = TRUE,
        metrics      = "Shannon",
        launch_alpha = 1
      )
      tab <- alphagrp_table()
      expect_s3_class(tab, "data.frame")
      # categorical factor → grouped summary: one row per modalite level (CL, La, SP)
      expect_equal(nrow(tab), 3L)
      expect_true(any(grepl("^mean_", colnames(tab))))
    }
  )
})

test_that("mod_alpha_server: reacalpha returns ANOVA + Tukey for a categorical factor", {
  shiny::testServer(
    mod_alpha_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        Fact1        = "modalite",
        checkbox1    = TRUE,
        metrics      = "Shannon",
        launch_alpha = 1
      )
      LL <- reacalpha()
      expect_named(LL, c("form1", "aov1", "groups1"))
      expect_true(grepl("Shannon", as.character(LL$form1)))
      expect_s3_class(LL$groups1, "data.frame")
    }
  )
})
