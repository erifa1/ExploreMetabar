# Guard: vegan 2.7-3 causes glibc heap corruption (`malloc(): invalid next
# size`) inside adonis2 / pairwiseAdonis permutation C code under R 4.6 — it
# crashes the whole R process during PERMANOVA in the Beta diversity tab. The
# fix is to pin vegan to 2.7-1 (validated in the sister project micropipeline).
# This test fails loudly if a future `renv::snapshot()` / dependency refresh
# re-bumps the lockfile to the known-bad 2.7-3. The matching upper bound lives
# in DESCRIPTION (`vegan (<= 2.7-1)`).
#
# renv.lock is .Rbuildignored, so it is absent under R CMD check / installed
# tests; the guard runs only in a source checkout (where snapshots happen).

test_that("renv.lock does not pin the heap-corrupting vegan 2.7-3", {
  lock <- testthat::test_path("..", "..", "renv.lock")
  skip_if_not(file.exists(lock), "renv.lock not present (not a source checkout)")

  lines <- readLines(lock, warn = FALSE)
  vi <- grep('^\\s*"vegan"\\s*:', lines)
  skip_if(length(vi) == 0, "no vegan entry in renv.lock")

  ver_line <- grep('"Version"', lines[vi[1]:(vi[1] + 6L)], value = TRUE)[1]
  ver <- sub('.*"Version"\\s*:\\s*"([^"]+)".*', "\\1", ver_line)

  expect_false(
    identical(ver, "2.7-3"),
    info = paste(
      "renv.lock pins vegan 2.7-3, which causes a `malloc(): invalid next size`",
      "heap-corruption crash in adonis2/pairwiseAdonis under R 4.6. Pin 2.7-1:",
      "renv::install('vegan@2.7-1'); renv::snapshot()."
    )
  )
})
