context("mod_heatmap: heatmap + feature selection")

test_that("mod_heatmap_ui returns a tag list", {
  ui <- mod_heatmap_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_heatmap_server: agglom_data() returns a phyloseq at ASV rank", {
  shiny::testServer(
    mod_heatmap_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(rank = "ASV")
      d <- agglom_data()
      expect_s4_class(d, "phyloseq")
      expect_equal(phyloseq::ntaxa(d),
                   phyloseq::ntaxa(r$phyloseq_filtered()))
    }
  )
})

test_that("mod_heatmap_server: agglom_data() agglomerates at a higher rank", {
  shiny::testServer(
    mod_heatmap_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(rank = "Phylum")
      d <- agglom_data()
      expect_s4_class(d, "phyloseq")
      # phylum-level agglomeration must reduce or preserve taxa count
      expect_lte(phyloseq::ntaxa(d),
                 phyloseq::ntaxa(r$phyloseq_filtered()))
    }
  )
})

test_that("mod_heatmap_server: selected_data() returns full agglom when select_features is FALSE", {
  shiny::testServer(
    mod_heatmap_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        rank             = "ASV",
        select_features  = FALSE,
        norm             = 0L
      )
      sd <- selected_data()
      expect_s4_class(sd, "phyloseq")
      expect_equal(phyloseq::ntaxa(sd),
                   phyloseq::ntaxa(r$phyloseq_filtered()))
    }
  )
})

test_that("mod_heatmap_server: agglom_normalized_data() applies TSS (norm=1)", {
  shiny::testServer(
    mod_heatmap_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        rank            = "ASV",
        select_features = FALSE,
        norm            = 1L     # total-sum scaling
      )
      d <- agglom_normalized_data()
      expect_s4_class(d, "phyloseq")
      # TSS rescales each sample to sum to 1
      sample_sums <- phyloseq::sample_sums(d)
      expect_true(all(abs(sample_sums - 1) < 1e-8))
    }
  )
})
