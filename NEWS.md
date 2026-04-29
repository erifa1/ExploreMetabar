# ExploreMetabar 3.0

## UI overhaul
* Full migration to **bslib Bootstrap 5**: all modules now use `layout_sidebar`, `navset_card_underline`, and `accordion` panels for a consistent, modern layout.

## New features
* sPLS-DA workflow (MixOmics module) restructured into three explicit stages: Initial exploration → Parameter tuning (ncomp + keepX) → Final model. Tuning results auto-fill the final stage inputs.

## Bug fixes
* Fixed NMDS screeplot crash in `mod_beta.R` (`return(p)` called outside branch where `p` was defined).
* Fixed unifrac validation logic in `mod_beta.R` (wrong boolean operator blocked all metrics when no tree was present).
* Fixed `req()` calls in `mod_alpha.R` passing reactive function objects instead of reactive values, bypassing NULL guards.
* Fixed `updateSelectInput()` missing `session` argument in `mod_cluster.R`.

## Internal changes
* Switched all module roxygen tags from `@export` to `@noRd` per golem convention (module functions are not part of the public API).

* Added missing Bioconductor packages (`Biobase`, `Biostrings`, `DESeq2`, `metagenomeSeq`, `metacoder`, `microbiome`, `mixOmics`, `reshape2`) to DESCRIPTION Imports.

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
