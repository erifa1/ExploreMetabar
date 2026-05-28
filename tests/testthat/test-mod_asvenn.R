context("mod_asvenn: Venn / alluvial")

test_that("mod_asvenn_ui returns a tag list", {
  ui <- mod_asvenn_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_asvenn_server: resVenn() is gated on go1", {
  shiny::testServer(
    mod_asvenn_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(Fact1 = "modalite", lvls1 = c("CL", "La", "SP"), minAb = 0)
      expect_error(resVenn(), class = "shiny.silent.error")
    }
  )
})

test_that("mod_asvenn_server: resVenn() builds the TF list and v.table when go1 fires", {
  shiny::testServer(
    mod_asvenn_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        Fact1 = "modalite",
        lvls1 = c("CL", "La", "SP"),
        minAb = 0,
        go1   = 1
      )
      res <- resVenn()
      expect_named(res, c("TF", "v.table"))
      expect_named(res$TF, c("CL", "La", "SP"))
      expect_s3_class(res$v.table, "data.frame")
      # one column per modality, plus 'taxa' and 'taxo'
      expect_true(all(c("taxa", "taxo", "CL", "La", "SP") %in% colnames(res$v.table)))
    }
  )
})

test_that("mod_asvenn_server: getPal returns one color per selected level", {
  shiny::testServer(
    mod_asvenn_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        Fact1 = "modalite",
        lvls1 = c("CL", "La"),
        minAb = 0,
        go1   = 1
      )
      pal <- getPal()
      expect_true(length(pal) >= 2L)
    }
  )
})
