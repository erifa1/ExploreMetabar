context("shared reactive factories and helpers")

# ---- round_df (pure helper) ------------------------------------------------

test_that("round_df rounds numeric columns to N digits", {
  df <- data.frame(a = c(1.23456, 7.89012), b = c("x", "y"),
                   stringsAsFactors = FALSE)
  out <- round_df(df, digits = 2)
  expect_equal(out$a, c(1.23, 7.89))
  expect_equal(out$b, c("x", "y"))
})

test_that("round_df leaves non-data.frame inputs untouched", {
  expect_identical(round_df(1:3), 1:3)
  expect_identical(round_df("a"), "a")
  expect_null(round_df(NULL))
})

test_that("round_df preserves column order and non-numeric types", {
  df <- data.frame(x = 1.111, y = 2L, z = TRUE, s = "a",
                   stringsAsFactors = FALSE)
  out <- round_df(df, digits = 1)
  expect_equal(colnames(out), c("x", "y", "z", "s"))
  expect_type(out$z, "logical")
  expect_type(out$s, "character")
})


# ---- reactive factories ----------------------------------------------------
# Each factory returns a reactive, so calling it requires an active reactive
# context. We drive them through testServer() using a tiny throwaway module.
# The factories are bound as locals inside the moduleServer function, which
# puts them in scope of the testServer expr block.

factories_module <- function(r) {
  function(id) {
    shiny::moduleServer(id, function(input, output, session) {
      factor_input <- shiny::reactive(input$facts)
      get_meta_col <<- make_get_meta_col(factor_input, r)
      local_metadata <<- make_local_metadata(factor_input, get_meta_col, r)
      local_metadata_all <<- make_local_metadata(
        factor_input, get_meta_col, r, keep_all_cols = TRUE
      )
      local_physeq <<- make_local_physeq(local_metadata, r)
      is_num <<- make_is_num_factor(get_meta_col, local_metadata_all)
    })
  }
}


test_that("make_get_meta_col returns a single factor name unchanged", {
  shiny::testServer(
    app = factories_module(mock_r()),
    expr = {
      session$setInputs(facts = "modalite")
      expect_equal(get_meta_col(), "modalite")
    }
  )
})

test_that("make_get_meta_col concatenates multiple factors with '_'", {
  shiny::testServer(
    app = factories_module(mock_r()),
    expr = {
      session$setInputs(facts = c("modalite", "Affinage"))
      expect_equal(get_meta_col(), "modalite_Affinage")
    }
  )
})

test_that("make_get_meta_col rejects multiple numeric factors", {
  shiny::testServer(
    app = factories_module(mock_r()),
    expr = {
      session$setInputs(facts = c("pH", "Salt.."))
      expect_error(get_meta_col(), class = "shiny.silent.error")
    }
  )
})

test_that("make_local_metadata unites multi-factor columns and casts to factor", {
  shiny::testServer(
    app = factories_module(mock_r()),
    expr = {
      session$setInputs(facts = c("modalite", "Affinage"))
      md <- local_metadata()
      expect_true("modalite_Affinage" %in% colnames(md))
      expect_s3_class(md$modalite_Affinage, "factor")
      expect_setequal(colnames(md), c("sample.id", "modalite_Affinage"))
    }
  )
})

test_that("make_local_metadata with keep_all_cols keeps every metadata column", {
  shiny::testServer(
    app = factories_module(mock_r()),
    expr = {
      session$setInputs(facts = "modalite")
      md_all <- local_metadata_all()
      expect_true(ncol(md_all) > 2)
      expect_true("pH" %in% colnames(md_all))
    }
  )
})

test_that("make_local_physeq swaps sample_data on the phyloseq object", {
  shiny::testServer(
    app = factories_module(mock_r()),
    expr = {
      session$setInputs(facts = "modalite")
      phy <- local_physeq()
      expect_s4_class(phy, "phyloseq")
      md <- as.data.frame(phyloseq::sample_data(phy))
      expect_true("modalite" %in% colnames(md))
    }
  )
})

test_that("make_is_num_factor distinguishes numeric and categorical metadata", {
  shiny::testServer(
    app = factories_module(mock_r()),
    expr = {
      session$setInputs(facts = "modalite")
      expect_false(is_num())
    }
  )
  shiny::testServer(
    app = factories_module(mock_r()),
    expr = {
      session$setInputs(facts = "pH")
      expect_true(is_num())
    }
  )
})
