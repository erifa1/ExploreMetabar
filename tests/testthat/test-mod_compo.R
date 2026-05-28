context("mod_compo: composition bar plots")

test_that("mod_compo_ui returns a tag list", {
  ui <- mod_compo_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_compo_server starts without errors and exposes empty state", {
  shiny::testServer(
    mod_compo_server,
    args = list(r = mock_r()),
    expr = {
      # The internal compo() reactive is gated on input$go1 — without a
      # button click it must not have computed anything yet.
      expect_error(compo(), class = "shiny.silent.error")
    }
  )
})

test_that("mod_compo_server computes both raw and relative plots when go1 fires", {
  shiny::testServer(
    mod_compo_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        Ord1        = "modalite",
        RankCompo   = "Genus",
        topTax      = 5,
        radio1      = "1",
        autoorder1  = TRUE,
        go1         = 1
      )
      LL <- compo()
      expect_type(LL, "list")
      expect_named(LL, c("p1", "p2"))
      # bars_fun returns a ggplot/plotly object; assert both are non-null
      expect_false(is.null(LL$p1))
      expect_false(is.null(LL$p2))
    }
  )
})

test_that("mod_compo_server merges samples when radio1 = 3", {
  shiny::testServer(
    mod_compo_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        Ord1        = "modalite",
        RankCompo   = "Genus",
        topTax      = 5,
        radio1      = "3",        # merge samples
        autoorder1  = TRUE,
        go1         = 1
      )
      LL <- compo()
      expect_type(LL, "list")
      expect_length(LL, 2L)
    }
  )
})
