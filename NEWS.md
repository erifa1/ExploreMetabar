# ExploreMetabar 3.1.0

## New features
* **Color management module (`mod_color`)**: a dedicated "Colors" panel centralizes every plot color. Assign per-modality colors for categorical metadata, continuous palettes for numeric metadata, and per-taxon colors for each taxonomic rank. Schemes can be exported and re-imported as long-format CSV.
* **Beta diversity**: ordination plots now overlay the top-N taxa contributions as labelled arrows.
* **Branding via `_brand.yml`**: theme colors and typography are driven by a single `inst/_brand.yml` (bslib `brand` integration). The navbar now carries the INRAE logo and external links (SK8, source code, issues) pushed to the right edge.

## Internal changes
* Color logic extracted into pure, Shiny-free helpers (`color_helpers.R`) consumed by `mod_color`; plotting modules read the shared `r$factor_colors()`, `r$numeric_palettes()`, and `r$taxa_colors()` reactives.
* Fixed numeric / high-cardinality variable handling in color assignment.
* Added a `testthat` suite covering the alpha, beta, cluster, color, composition, data-loading, differential-analysis, heatmap, mixomics, taxa-boxplot, and shared-utils modules.

## Cleanup
* Removed obsolete scratch files (`inst/old_files/`, `inst/test/`), a duplicate test fixture in `data-raw/`, and the unused `ranks_ref` dataset (dropped `LazyData`).

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
