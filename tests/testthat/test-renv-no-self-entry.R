# Guard: the deploy image installs ExploreMetabar from a source clone at the
# release tag (see Dockerfile in the SK8 deploy repo). An `ExploreMetabar`
# self-entry in renv.lock is a leftover of the old install path; `renv::restore`
# would then try to clone+install the package itself, masking the source build
# and pinning a stale SHA. `renv::snapshot()` normally omits the project package,
# so a self-entry only reappears if added by hand -- this test catches that.
#
# renv.lock is .Rbuildignored, so it is absent under R CMD check / installed
# tests; the guard runs only in a source checkout (where snapshots happen).

test_that("renv.lock has no ExploreMetabar self-entry", {
  lock <- testthat::test_path("..", "..", "renv.lock")
  skip_if_not(file.exists(lock), "renv.lock not present (not a source checkout)")

  lines <- readLines(lock, warn = FALSE)
  expect_false(
    any(grepl('^\\s*"ExploreMetabar"\\s*:', lines)),
    info = paste(
      "renv.lock contains an ExploreMetabar self-entry. Remove it:",
      "the SK8 image installs the package from a source clone at the tag,",
      "not via renv.lock. A self-entry pins a stale SHA and breaks that flow."
    )
  )
})
