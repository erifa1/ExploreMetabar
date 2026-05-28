context("mod_taxaboxplot: per-taxon boxplots and stats")

test_that("mod_taxaboxplot_ui returns a tag list", {
  ui <- mod_taxaboxplot_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_taxaboxplot_server: get_pval_table is gated on go1", {
  shiny::testServer(
    mod_taxaboxplot_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(boxplot_fact1 = "modalite", order1 = TRUE)
      expect_error(get_pval_table(), class = "shiny.silent.error")
    }
  )
})

test_that("mod_taxaboxplot_server: Kruskal-Wallis p-value table for categorical factor", {
  shiny::testServer(
    mod_taxaboxplot_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        boxplot_fact1 = "modalite",
        order1        = TRUE,
        go1           = 1
      )
      tab <- get_pval_table()
      expect_s3_class(tab, "data.frame")
      expect_true(all(c("taxa", "p.value", "p.adj") %in% colnames(tab)))
      # one row per taxon (joined with tax_table)
      expect_equal(
        nrow(tab),
        phyloseq::ntaxa(r$phyloseq_filtered_norm())
      )
    }
  )
})

test_that("mod_taxaboxplot_server: correlation table for numeric factor", {
  shiny::testServer(
    mod_taxaboxplot_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        boxplot_fact1 = "pH",
        order1        = FALSE,
        cor_test      = "spearman",
        go1           = 1
      )
      tab <- get_pval_table()
      expect_s3_class(tab, "data.frame")
      expect_true("cor.coef" %in% colnames(tab))
    }
  )
})
