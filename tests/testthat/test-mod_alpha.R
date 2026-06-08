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

# ---- pure helpers ------------------------------------------------------------

test_that("box_whisker_stats matches Tukey quartiles and clamps whiskers", {
  v <- c(1:10, 100)  # 100 is an outlier beyond the upper fence
  s <- box_whisker_stats(v)
  expect_named(s, c("q1", "median", "q3", "lowerfence", "upperfence"))
  q <- stats::quantile(v, c(.25, .5, .75), type = 7, names = FALSE)
  expect_equal(c(s$q1, s$median, s$q3), q)
  # whiskers stop at the most extreme value within 1.5*IQR (not at the outlier)
  expect_equal(s$lowerfence, 1)
  expect_equal(s$upperfence, 10)
  expect_lt(s$upperfence, max(v))
})

test_that("box_whisker_stats handles a single-value group and empty input", {
  s <- box_whisker_stats(c(5, 5, 5))
  expect_equal(unname(unlist(s)), rep(5, 5))
  empty <- box_whisker_stats(numeric(0))
  expect_true(all(is.na(unlist(empty))))
})

test_that("p_to_stars maps p-values to the expected thresholds", {
  expect_equal(p_to_stars(c(0.0005, 0.005, 0.04, 0.2)),
               c("***", "**", "*", "ns"))
})

test_that("make_sig_brackets builds one bracket per significant pair", {
  lvls <- c("A", "B", "C")
  sig <- data.frame(
    comparison = c("B-A", "C-A"),
    `p adj`    = c(0.0004, 0.03),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  br <- make_sig_brackets(lvls, sig, ymax = 10, yrange = 5)
  expect_length(br$shapes, 2)
  expect_length(br$annotations, 2)
  # stars come from p adj
  expect_equal(vapply(br$annotations, `[[`, "", "text"), c("***", "*"))
  # brackets are stacked at increasing heights, ytop sits above the top one
  heights <- vapply(br$annotations, `[[`, numeric(1), "y")
  expect_true(all(diff(heights) > 0))
  expect_gt(br$ytop, max(heights))
  # B-A spans category indices 0..1 (0-based), C-A spans 0..2
  expect_match(br$shapes[[1]]$path, "^M 0,.* L 0,.* L 1,.* L 1,")
  expect_match(br$shapes[[2]]$path, "^M 0,.* L 0,.* L 2,.* L 2,")
})

test_that("make_sig_brackets resolves level names containing a hyphen", {
  lvls <- c("ctrl-1", "ctrl-2")
  sig <- data.frame(
    comparison = "ctrl-2-ctrl-1",  # naive strsplit on '-' would break this
    `p adj`    = 0.01,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  br <- make_sig_brackets(lvls, sig, ymax = 1, yrange = 1)
  expect_length(br$shapes, 1)
  expect_match(br$shapes[[1]]$path, "^M 0,.* L 0,.* L 1,.* L 1,")
})

test_that("make_sig_brackets returns empty geometry when nothing is significant", {
  br <- make_sig_brackets(c("A", "B"),
                          sig_df = data.frame(comparison = character(0),
                                              `p adj` = numeric(0),
                                              check.names = FALSE),
                          ymax = 10, yrange = 5)
  expect_length(br$shapes, 0)
  expect_length(br$annotations, 0)
  expect_equal(br$ytop, 10)
})
