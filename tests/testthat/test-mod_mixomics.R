context("mod_mixomics: sPLS-DA")

test_that("mod_mixomics_ui returns a tag list", {
  ui <- mod_mixomics_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_mixomics_server: x() returns an OTU matrix oriented sample x taxa", {
  shiny::testServer(
    mod_mixomics_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(factor_spls_da = "modalite")
      mat <- x()
      expect_equal(nrow(mat),
                   phyloseq::nsamples(r$phyloseq_filtered()))
      expect_equal(ncol(mat),
                   phyloseq::ntaxa(r$phyloseq_filtered()))
    }
  )
})

test_that("mod_mixomics_server: y() returns the requested factor as a vector", {
  shiny::testServer(
    mod_mixomics_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(factor_spls_da = "modalite")
      fact <- y()
      expect_length(fact, phyloseq::nsamples(r$phyloseq_filtered()))
      expect_setequal(unique(as.character(fact)), c("CL", "La", "SP"))
    }
  )
})

test_that("mod_mixomics_server: splsda_initial is gated on launch_initial", {
  shiny::testServer(
    mod_mixomics_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        factor_spls_da   = "modalite",
        nb_comp_initial  = 2
      )
      expect_error(splsda_initial(), class = "shiny.silent.error")
    }
  )
})

test_that("mod_mixomics_server: list_colors rejects numeric factors", {
  shiny::testServer(
    mod_mixomics_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(factor_spls_da = "pH")
      expect_error(list_colors(), class = "shiny.silent.error")
    }
  )
})

test_that("tune script template ships and fully substitutes its tokens", {
  tmpl_path <- app_sys("templates", "tune_splsda_template.R")
  expect_true(file.exists(tmpl_path))

  tmpl <- readLines(tmpl_path)
  expect_true(any(grepl("{{FACTOR}}",   tmpl, fixed = TRUE)))
  expect_true(any(grepl("{{DISTANCE}}", tmpl, fixed = TRUE)))
  expect_true(any(grepl("{{MEASURE}}",  tmpl, fixed = TRUE)))

  rendered <- tmpl
  rendered <- gsub("{{FACTOR}}",   "modalite", rendered, fixed = TRUE)
  rendered <- gsub("{{DISTANCE}}", "max.dist", rendered, fixed = TRUE)
  rendered <- gsub("{{MEASURE}}",  "BER",      rendered, fixed = TRUE)
  expect_false(any(grepl("{{", rendered, fixed = TRUE)))
})
