# ExploreMetabar 2.1.0

## Bug fixes
* Fixed `browser()` calls left in production code (`mod_beta.R`, `mod_source_tracker.R`).
* Fixed `refseq(table, ...)` bug in `mod_diffanalysis.R` — was calling base R `table` instead of the phyloseq object.
* Fixed filtered ASV table download exporting raw (unfiltered) data instead of filtered data in `mod_data_loading.R`.
* Fixed unreachable `setProgress()` call after `return()` in `mod_alpha.R`.

## Internal changes
* Version number is now dynamically read from DESCRIPTION via `packageVersion()` — no more hardcoded version strings.
* Synchronized version across DESCRIPTION, golem-config.yml, and UI header.
* Replaced hardcoded `/tmp/` paths with `tempdir()` for portability and multi-user safety.
* Removed dead/commented-out code blocks (unused module references, debug print statements).

# ExploreMetabar 0.0.0.9000

* Initial golem scaffold.
* Added a `NEWS.md` file to track changes to the package.
