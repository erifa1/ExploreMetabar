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

An 11th module, `mod_color.R`, is mounted inside `mod_data_loading_ui` (in the "Colors" nav panel). It owns all color maps and writes the `r$factor_colors()` / `r$numeric_palettes()` / `r$taxa_colors()` reactives consumed by every plotting module. See "Color contract" below.

### Shared reactive state

The `r` object is the backbone of the application:
- `r$phyloseq_filtered()` — the active phyloseq object after user filtering
- `r$sdat()` — sample metadata data frame
- `r$var_list()` — character vector of available metadata variable names
- `r$data_ready()` — gating reactive; `TRUE` once the user has loaded and processed data in the Input Data tab
- `r$tabs$tabselected` — tracks the active nav panel (set via an `observe` on `input$tabs`)
- `r$fdata` — initialized to `NULL` in `reactiveValues`; populated by `mod_data_loading`
- `r$factor_list()` — character vector of metadata variables whose column is **not** numeric (i.e. factors / characters). Populated by `mod_data_loading`.
- `r$factor_colors()` — named list keyed by factor variable; each entry is a named character vector `c(<modality> = "#hex", ...)`. **Factors only** — numeric variables are not in this list.
- `r$numeric_palettes()` — named list keyed by numeric metadata variable; each entry is a `"package::palette"` string accepted by `paletteer::paletteer_c`. Consumers that color by a numeric variable must read from here, not `r$factor_colors()`.
- `r$taxa_colors()` — named list keyed by taxonomic rank; each entry is a named character vector `c(<taxon> = "#hex", ...)`.

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
- `color_helpers.R` — pure (Shiny-free) palette helpers used by `mod_color`
- `mod_color.R` — server-only module that owns the color reactives and the "Colors" UI panel inside `mod_data_loading_ui`
- `pairwise_adonis.R` — pairwise PERMANOVA utility used by `mod_beta`

### Color contract

`mod_color_server("color_ui_1", r = r)` is mounted as a nested module inside `mod_data_loading_server` and writes three reactives onto the shared `r`:

- `r$factor_colors()` — categorical metadata variables only (named hex vector per variable).
- `r$numeric_palettes()` — numeric metadata variables (`"package::palette"` string per variable). **Never** appears in `r$factor_colors()`.
- `r$taxa_colors()` — taxonomic ranks (named hex vector per rank).

When a consumer module colors a plot by a metadata variable that may be numeric, it must branch on `var %in% names(r$numeric_palettes())` (or equivalent) and wire `paletteer::scale_*_paletteer_c()` in the numeric path instead of `scale_*_manual(values = r$factor_colors()[[var]])`. Plots that are inherently categorical (Venn, stacked bars, sPLS-DA, pairwise differential contrasts) `validate(need(...))` upfront.

CSV import/export uses long-format with columns `variable, modality, color`. Numeric assignments are rows with `variable = "__numeric__"` where `modality` is the variable name and `color` is a `package::palette` string.

## Visual Testing

Do not attempt to launch or visually test the app. Only verify syntax and logic correctness (R CMD check, `devtools::check()`, `devtools::test()`). The user handles all final UI/UX testing.

## Testing

`tests/testthat/` covers each module (`test-mod_*.R`), the shared reactive factories (`test-mod_utils_shared.R`), and a deploy guard (`test-renv-no-self-entry.R`), alongside the golem `test-golem-recommended.R` app-launch stub; shared mocks live in `helper-mock-r.R`. When adding tests for a module, follow the standard `testServer()` pattern and load fixtures from `inst/data_test/`.

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
```

Do **not** hand-copy `renv.lock` into the deploy repo — `make release` (see *Releasing*) syncs it as part of every release. `renv.lock` must **not** contain an `ExploreMetabar` self-entry: `renv::snapshot()` omits the project package by design, and `tests/testthat/test-renv-no-self-entry.R` fails if one reappears (the deploy installs the package from a source clone at the release tag, not from the lockfile).

Never add a CRAN-only package that duplicates functionality already provided by the Bioconductor ecosystem (e.g., prefer `phyloseq` transforms over reimplementing them).

## Releasing a new version

Releases are driven by one command, run from this repo's `master` with a clean tree:

```bash
make release VERSION=3.2.0
```

`make release` ([Makefile](Makefile)):
1. bumps `Version:` in `DESCRIPTION` and `golem_version:` in `inst/golem-config.yml` (idempotent),
2. commits and creates the annotated tag `v3.2.0`, then pushes `master` + the tag to forge,
3. runs the `sync-deploy` target: copies `renv.lock` and pins `ARG EM_REF=v3.2.0` in the deploy repo's `Dockerfile`, then commits + pushes the deploy repo — which fires the SK8 build and auto-deploy.

Add a `NEWS.md` section for the version first (the target warns if it is missing). The UI header and startup log read the version from `DESCRIPTION` via `packageVersion()`, so nothing else needs editing. `DEPLOY_REPO` defaults to `$HOME/git-repo/explore-metabar` and is overridable (`make release VERSION=… DEPLOY_REPO=…`).

## Deployment (sister repo)

Deployment lives in a **separate repository** at `$HOME/git-repo/explore-metabar` (`sk8/sk8-apps/tbi/explore-metabar` on forge) — it holds the `Dockerfile` and `.gitlab-ci-sk8.yml` for the SK8 (INRAE) Shiny hosting platform. SK8 has **no inbound trigger**: it rebuilds whenever a commit/tag is pushed to that repo (which `make release` does) and then auto-deploys (`MISE_A_JOUR_AUTOMATIQUE: "true"`).

The deploy `Dockerfile` builds on `rocker/r-ver:4.6.0`, installs `renv 1.2.2`, then:
1. **deps** — `renv::restore()` from the synced `renv.lock` (dependency tree only, **no** `ExploreMetabar` self-entry);
2. **app** — `git clone`s this (public) repo at the tag in `ARG EM_REF` and `R CMD INSTALL`s from source. The tag — not the lockfile — is the single source of truth for the deployed version.

Never hand-edit the deploy repo's `renv.lock` or `EM_REF`; `make release` writes both. SK8's CI also reads `renv.lock` for the R version and dependency cache. The local `Dockerfile` in this repo is obsolete and not used for SK8 deployment.

## Golem Conventions

- New modules are scaffolded with `golem::add_module("name")`, which creates `mod_name.R` and a test stub.
- New dependencies are added with `usethis::use_package("pkg")` or `golem::add_js_file()` / `golem::add_css_file()` for assets.
- `inst/golem-config.yml` controls the app name and production flag (`golem.app.prod`).
- The `NAMESPACE` and `man/` files are auto-generated by roxygen2 — never edit them by hand.
