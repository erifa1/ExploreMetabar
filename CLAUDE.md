# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

ExploreMetabar is a production-grade **Shiny web application** built as an R package using the [Golem](https://golemverse.org/) framework. It provides interactive exploration of metabarcoding (16S, ITS) data via a phyloseq-based workflow.

## Development Commands

```r
# Load all and launch the app locally
pkgload::load_all()
ExploreMetabar::run_app()

# Regenerate documentation (NAMESPACE, man/)
devtools::document(roclets = c('rd', 'collate', 'namespace'))

# Run tests
devtools::test()

# Run a single test file
devtools::test(filter = "mod_beta")

# Full package check
devtools::check()
```

Dev scripts in `dev/` orchestrate these tasks: `dev/02_dev.R` for adding modules/dependencies, `dev/03_deploy.R` for pre-deployment checks, `dev/run_dev.R` for quick local runs.

## Architecture

### Module system

The app is a navbar with 10 independent Shiny modules, all sharing a single `r = reactiveValues()` object passed from `app_server.R`:

| Module | File | Purpose |
|--------|------|---------|
| Data Loading | `mod_data_loading.R` | Upload RData, inspect/filter phyloseq objects |
| Composition | `mod_compo.R` | Stacked bar plots by taxa/group |
| Alpha diversity | `mod_alpha.R` | Observed, Shannon, Chao1, etc. + ANOVA/Tukey |
| Beta diversity | `mod_beta.R` | PCoA, NMDS, adonis2 PERMANOVA |
| Taxa boxplots | `mod_taxaboxplot.R` | Abundance boxplots + statistical tests |
| Heatmap | `mod_heatmap.R` | Heatmap + dendextend clustering |
| Diff analysis | `mod_diffanalysis.R` | DESeq2, metagenomeSeq, metacoder |
| ASV Venn | `mod_asvenn.R` | Venn diagrams for ASV overlap |
| Clustering | `mod_cluster.R` | Hierarchical clustering + dendrograms |
| MixOmics | `mod_mixomics.R` | sPLS-DA multivariate analysis |

### Shared reactive state

The `r` object is the backbone of the application:
- `r$phyloseq_filtered()` — the active phyloseq object after user filtering
- `r$sdat()` — sample metadata data frame
- `r$var_list()` — character vector of available metadata variable names
- `r$data_ready()` — gating reactive; `TRUE` once the user has loaded and processed data in the Input Data tab
- `r$tabs$tabselected` — tracks the active nav panel (set via an `observe` on `input$tabs`)
- `r$fdata` — initialized to `NULL` in `reactiveValues`; populated by `mod_data_loading`

### Navigation invariant

`app_server.R` installs an `observeEvent(input$tabs, …)` that bounces the user back to the Input Data tab via `bslib::nav_select(id = "tabs", selected = "data_loading")` whenever they try to switch tabs before `r$data_ready()` is `TRUE`. New modules must not bypass this — assume `r$phyloseq_filtered()` is available by the time the module's reactives run, but still `req()` it to be safe.

### Global server setup

`app_server.R` also:
- bumps `shiny.maxRequestSize` to 30 MB for RData uploads
- enables `thematic::thematic_shiny()` so ggplot/base plots inherit bslib theme colors automatically — do **not** hard-code palette/theme inside modules; rely on thematic + `_brand.yml`

### Reactive factories (`mod_utils_shared.R`)

Four factories create reusable reactives for module-level data access. Each takes the previously built reactives as arguments (composition, not a single `r` handle):

- `make_get_meta_col(factor_input, r)` — reactive returning the metadata column name for the selected factor(s); collapses multiple factors with `_` and rejects multi-numeric selections via `validate()`
- `make_is_num_factor(get_meta_col, local_metadata)` — reactive returning `TRUE` when the chosen factor column is numeric
- `make_local_metadata(factor_input, get_meta_col, r, keep_all_cols = FALSE)` — reactive returning the module's local metadata data frame, with multi-factor columns united via `tidyr::unite` and converted to factor. `mod_beta` passes `keep_all_cols = TRUE`
- `make_local_physeq(local_metadata, r)` — reactive returning a phyloseq object whose `sample_data` has been swapped for the local metadata

Each module composes these in its server function and uses the resulting reactives rather than touching `r` directly for reads.

Also in `mod_utils_shared.R`: `round_df(x, digits = 4)` — used at the `DT::renderDataTable` rendering layer to round numeric columns for display without mutating upstream reactive data that other reactives consume at full precision.

### Performance pattern

Heavy computations (ordinations, differential tests, permutation tests) are gated behind `eventReactive()` tied to an action button, never triggered on every input change. This is intentional — do not convert these to reactive expressions.

## Key Files

- `app_ui.R` — navbar scaffold with Bootstrap 5 theme (`_brand.yml` colors: primary `#00a3a6`, secondary `#423089`)
- `app_server.R` — wires all 10 module servers together with shared `r`
- `phyloseq_extended_graphical_methods.R` — domain-specific plotting utilities (stacked bars, rarefaction curves)
- `bars_fun.R` — bar plot generation helpers
- `color_functions.R` — custom palette management
- `pairwise_adonis.R` — pairwise PERMANOVA utility used by `mod_beta`

## Visual Testing

Do not attempt to launch or visually test the app. Only verify syntax and logic correctness (R CMD check, `devtools::check()`, `devtools::test()`). The user handles all final UI/UX testing.

## Testing

The test suite is currently minimal: `tests/testthat/` contains only the default `test-golem-recommended.R` stub from `golem::use_recommended_tests()`. There are no module-level tests and no snapshot directory yet. When adding tests for a module, follow the standard `testServer()` pattern and load fixtures from `inst/data_test/`.

Test datasets for local development:
- `inst/data_test/phy_test_numeric.rdata` — small object (24 KB)
- `inst/data_test/robjects_600.Rdata` — medium (55 KB)
- `inst/data_test/robjects.Rdata` — full (269 KB)

## Dependencies and reproducibility

R 4.6.0, Bioconductor 3.23. The project uses `renv` with `renv.lock` as the **single source of truth** for the dependency tree. To bootstrap:

```r
renv::restore()
```

Only three packages are pulled from GitHub (`DESCRIPTION` `Remotes:`): `vmikk/metagMisc`, `mikemc/speedyseq`, `vqf/nVennR`. Everything else resolves from CRAN or Bioc 3.23 — do not add new GitHub remotes unless the package is unavailable from a canonical source.

When you add or update a dependency in `DESCRIPTION`:

```r
pak::local_install_deps(dependencies = TRUE)   # install the new package(s)
renv::snapshot()                                # rewrite renv.lock
file.copy("renv.lock", "../explore-metabar/renv.lock", overwrite = TRUE)
```

Never add a CRAN-only package that duplicates functionality already provided by the Bioconductor ecosystem (e.g., prefer `phyloseq` transforms over reimplementing them).

## Deployment (sister repo)

Deployment lives in a **separate repository** at `../explore-metabar/` — it holds the `Dockerfile` and `.gitlab-ci-sk8.yml` for the SK8 (INRAE) Shiny hosting platform. That Dockerfile builds on `rocker/r-ver:4.6.0`, installs `renv 1.2.2`, then runs `renv::restore()` against a copy of this repo's `renv.lock`. The two `renv.lock` files **must stay byte-identical** — sync after every `renv::snapshot()` (see command above). The local `Dockerfile` in this repo is obsolete and not used for SK8 deployment.

## Golem Conventions

- New modules are scaffolded with `golem::add_module("name")`, which creates `mod_name.R` and a test stub.
- New dependencies are added with `usethis::use_package("pkg")` or `golem::add_js_file()` / `golem::add_css_file()` for assets.
- `inst/golem-config.yml` controls the app name and production flag (`golem.app.prod`).
- The `NAMESPACE` and `man/` files are auto-generated by roxygen2 — never edit them by hand.
