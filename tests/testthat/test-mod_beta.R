# tests/testthat/test-mod_beta.R
#
# Tests for mod_beta — validates the beta diversity analysis pipeline.
#
# Since mod_beta_server uses the legacy callModule() pattern (not moduleServer),
# testServer() cannot access internal reactives directly. We test the
# underlying computation pipeline that the module executes.
#
# NOTE: Some environments (notably R 4.5.x + vegan 2.7.x + glibc 2.39)
# experience glibc malloc crashes when running multiple vegan/phyloseq
# operations in a single testthat session. This is a C-level memory issue
# in compiled dependencies, not a bug in ExploreMetabar code.
# If you see "malloc(): invalid next size" errors, run tests individually:
#   Rscript -e 'source("R/pairwise_adonis.R"); testthat::test_file("tests/testthat/test-mod_beta.R")'
#
# Usage:
#   devtools::test(filter = "mod_beta")
#   testthat::test_file("tests/testthat/test-mod_beta.R")

library(testthat)
library(phyloseq)


# ── Load test data ──
test_physeq <- local({
  path <- system.file("data_test/robjects_600.Rdata", package = "ExploreMetabar")
  if (path == "") {
    path <- file.path("../../inst/data_test/robjects_600.Rdata")
  }
  env <- new.env()
  load(path, envir = env)
  env$data
})

test_sdat <- local({
  sd <- as.data.frame(sample_data(test_physeq))
  sd$sample.id <- rownames(sd)
  sd
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 1: veganifyOTU helper function
# ══════════════════════════════════════════════════════════════════════════════
test_that("veganifyOTU returns a samples x taxa matrix", {
  mat <- veganifyOTU(test_physeq)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), nsamples(test_physeq))
  expect_equal(ncol(mat), ntaxa(test_physeq))
  expect_equal(sort(rownames(mat)), sort(sample_names(test_physeq)))
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 2: pairwise.adonis standalone
# ══════════════════════════════════════════════════════════════════════════════
test_that("pairwise.adonis returns correct structure", {
  mat <- veganifyOTU(test_physeq)
  d <- vegan::vegdist(mat, method = "bray")

  res <- pairwise.adonis(d, factors = test_sdat$SampleType, p.adjust.m = "fdr")

  expect_s3_class(res, "data.frame")
  expect_true("pairs" %in% colnames(res))
  expect_true("p.adjusted" %in% colnames(res))
  expect_true("sig" %in% colnames(res))
  # 4 levels of SampleType -> C(4,2) = 6 pairwise comparisons
  expect_equal(nrow(res), 6)
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 3: Hellinger transformation
# ══════════════════════════════════════════════════════════════════════════════
test_that("Hellinger transformation produces valid output", {
  mat <- veganifyOTU(test_physeq)
  spe <- vegan::decostand(mat, method = "hell")

  expect_true(is.matrix(spe))
  expect_equal(dim(spe), dim(mat))
  expect_true(all(spe >= 0))
  expect_true(all(spe <= 1))
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 4: Bray-Curtis distance
# ══════════════════════════════════════════════════════════════════════════════
test_that("Bray-Curtis distance matrix is computed correctly", {
  dist_mat <- phyloseq::distance(test_physeq, method = "bray")

  expect_s3_class(dist_mat, "dist")
  expect_equal(attr(dist_mat, "Size"), nsamples(test_physeq))
  expect_true(all(as.vector(dist_mat) >= 0))
  expect_true(all(as.vector(dist_mat) <= 1))
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 5: PCOA ordination (capscale)
# ══════════════════════════════════════════════════════════════════════════════
test_that("PCOA ordination via capscale runs correctly", {
  mat <- veganifyOTU(test_physeq)
  spe <- vegan::decostand(mat, method = "hell")

  ord_res <- vegan::capscale(spe ~ 1, spe, distance = "bray")
  expect_true(inherits(ord_res, "capscale"))

  sites <- vegan::scores(ord_res, display = "sites", correlation = TRUE)
  expect_equal(nrow(sites), nsamples(test_physeq))

  axes <- colnames(sites)
  expect_true(length(axes) >= 2)
  expect_true(all(grepl("MDS", axes)))

  eig <- vegan::eigenvals(ord_res)
  expect_true(length(eig) > 0)
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 6: NMDS ordination
# ══════════════════════════════════════════════════════════════════════════════
test_that("NMDS ordination runs correctly", {
  mat <- veganifyOTU(test_physeq)
  spe <- vegan::decostand(mat, method = "hell")

  nmds_res <- vegan::metaMDS(spe, k = 2, distance = "bray",
                              trace = FALSE, autotransform = FALSE)
  expect_true(inherits(nmds_res, "metaMDS"))

  sites <- vegan::scores(nmds_res, display = "sites")
  expect_equal(nrow(sites), nsamples(test_physeq))
  expect_true(nmds_res$stress >= 0)
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 7: RDA constrained ordination
# ══════════════════════════════════════════════════════════════════════════════
test_that("RDA constrained ordination runs correctly", {
  mat <- veganifyOTU(test_physeq)
  spe <- vegan::decostand(mat, method = "hell")

  rda_res <- vegan::rda(spe ~ SampleType, data = test_sdat)
  expect_true(inherits(rda_res, "rda"))

  sites <- vegan::scores(rda_res, display = "sites", correlation = TRUE)
  expect_equal(nrow(sites), nsamples(test_physeq))

  species <- vegan::scores(rda_res, display = "species", correlation = TRUE)
  expect_equal(nrow(species), ntaxa(test_physeq))
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 8: Betadisper (dispersion analysis)
# ══════════════════════════════════════════════════════════════════════════════
test_that("Betadisper dispersion analysis runs correctly", {
  dist_mat <- phyloseq::distance(test_physeq, method = "bray")

  disp <- vegan::betadisper(dist_mat, test_sdat$SampleType)
  expect_s3_class(disp, "betadisper")

  disp_anova <- anova(disp)
  expect_s3_class(disp_anova, "data.frame")

  disp_tukey <- TukeyHSD(disp)
  expect_true("group" %in% names(disp_tukey))
  expect_true(all(disp$distances >= 0))
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 9: PERMANOVA (adonis2)
# ══════════════════════════════════════════════════════════════════════════════
test_that("PERMANOVA adonis2 runs correctly", {
  dist_mat <- phyloseq::distance(test_physeq, method = "bray")
  sdat_test <- test_sdat
  sdat_test$Depth <- sample_sums(test_physeq)

  adonis_res <- vegan::adonis2(
    dist_mat ~ Depth + SampleType,
    data = sdat_test,
    permutations = 99
  )

  expect_s3_class(adonis_res, "data.frame")
  expect_true("Pr(>F)" %in% colnames(adonis_res))
  expect_true("R2" %in% colnames(adonis_res))
  expect_true(adonis_res["SampleType", "R2"] > 0)
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 10: Envfit
# ══════════════════════════════════════════════════════════════════════════════
test_that("Envfit runs correctly on PCOA ordination", {
  mat <- veganifyOTU(test_physeq)
  spe <- vegan::decostand(mat, method = "hell")
  ord_res <- vegan::capscale(spe ~ 1, spe, distance = "bray")

  en <- vegan::envfit(ord_res, test_sdat[, "SampleType", drop = FALSE], na.rm = TRUE)
  expect_true(!is.null(en$factors) || !is.null(en$vectors))
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 11: Unifrac distance (requires phylogenetic tree)
# ══════════════════════════════════════════════════════════════════════════════
test_that("Unifrac distance works when tree is present", {
  skip_if(is.null(phy_tree(test_physeq, errorIfNULL = FALSE)),
          "No phylogenetic tree in test data")

  dist_uf <- phyloseq::distance(test_physeq, method = "unifrac")
  expect_s3_class(dist_uf, "dist")
  expect_equal(attr(dist_uf, "Size"), nsamples(test_physeq))
})


# ══════════════════════════════════════════════════════════════════════════════
# Test 12: Full PCOA pipeline produces a valid ggplot
# ══════════════════════════════════════════════════════════════════════════════
test_that("Full PCOA pipeline produces a valid ggplot", {
  meta_col <- "SampleType"
  spe <- vegan::decostand(veganifyOTU(test_physeq), method = "hell")
  ord_res <- vegan::capscale(spe ~ 1, spe, distance = "bray")

  sites_df <- tibble::as_tibble(
    vegan::scores(ord_res, display = "sites", correlation = TRUE),
    rownames = "sample.id"
  )
  sites_df <- dplyr::inner_join(
    sites_df,
    test_sdat[, c("sample.id", meta_col)],
    by = "sample.id"
  )

  axes <- colnames(vegan::scores(ord_res, display = "sites"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = sites_df,
      ggplot2::aes(
        x = .data[[axes[1]]],
        y = .data[[axes[2]]],
        fill = .data[[meta_col]]
      ),
      shape = 23, size = 5
    ) +
    ggplot2::stat_ellipse(
      data = sites_df,
      ggplot2::aes(
        x = .data[[axes[1]]],
        y = .data[[axes[2]]],
        group = .data[[meta_col]],
        color = .data[[meta_col]]
      )
    )

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(sites_df), nsamples(test_physeq))
  expect_true(meta_col %in% colnames(sites_df))
})
