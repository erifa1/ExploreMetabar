context("mod_cluster: hierarchical clustering")

test_that("mod_cluster_ui returns a tag list", {
  ui <- mod_cluster_ui("test")
  expect_s3_class(ui, c("shiny.tag", "shiny.tag.list"), exact = FALSE)
})

test_that("mod_cluster_server: compute.dist returns a dist object", {
  shiny::testServer(
    mod_cluster_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        dist.meth        = "bray",
        hclust.meth      = "ward.D2",
        clstr_rank_glom  = "ASV",
        clust_fact1      = "modalite"
      )
      d <- compute.dist()
      expect_s3_class(d, "dist")
      expect_equal(attr(d, "Size"),
                   phyloseq::nsamples(r$phyloseq_filtered_norm()))
    }
  )
})

test_that("mod_cluster_server: compute.clust returns an hclust", {
  shiny::testServer(
    mod_cluster_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        dist.meth   = "bray",
        hclust.meth = "ward.D2",
        clust_fact1 = "modalite"
      )
      hc <- compute.clust()
      expect_s3_class(hc, "hclust")
    }
  )
})

test_that("mod_cluster_server: compute.k is gated on launch_clust", {
  shiny::testServer(
    mod_cluster_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        dist.meth   = "bray",
        hclust.meth = "ward.D2",
        clust_fact1 = "modalite",
        k.meth      = "silhouette"
      )
      expect_error(compute.k(), class = "shiny.silent.error")
    }
  )
})

test_that("mod_cluster_server: compute.k returns a positive integer when launched", {
  shiny::testServer(
    mod_cluster_server,
    args = list(r = mock_r()),
    expr = {
      session$setInputs(
        dist.meth     = "bray",
        hclust.meth   = "ward.D2",
        clust_fact1   = "modalite",
        k.meth        = "silhouette",
        launch_clust  = 1
      )
      k <- compute.k()
      expect_true(is.numeric(k))
      expect_gte(k, 1L)
    }
  )
})
