#' @import shiny
#' @import bslib
#' @importFrom bsicons bs_icon
#' @importFrom base64enc dataURI
#' @importFrom brand.yml read_brand_yml

# `brand.yml` is only a Suggests of bslib, but bslib delegates all `_brand.yml`
# parsing to it, so it is a hard runtime dependency of this app. The importFrom
# above declares that (and keeps R CMD check's dependency note quiet).

app_ui <- function() {
  # Resolve these at runtime, NOT at package top level. A top-level
  # `system.file(package = "ExploreMetabar")` is evaluated while the lazy-load
  # database is built during install, when it points at the staging directory
  # (.../00LOCK-ExploreMetabar/00new/...). That absolute path then gets baked
  # into the installed package and trips R's staged-install guard with
  # "ERROR: hard-coded installation path". Resolving inside the function keeps
  # installation safe and still points at the real install location at runtime.
  SK8img <- base64enc::dataURI(file=system.file(file.path('app/www', 'SK8.png'), package='ExploreMetabar'), mime="image/png")
  # White INRAE logo for the (teal) navbar; mime must be set so the SVG renders.
  INRAEimg <- base64enc::dataURI(file=system.file(file.path('app/www', 'inrae-logo-white.svg'), package='ExploreMetabar'), mime="image/svg+xml")

  # INRAE branding lives in inst/_brand.yml so it ships with the installed
  # package; bslib reads colors + typography from it (single source of truth).
  brand_path <- system.file("_brand.yml", package = "ExploreMetabar")

  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),

    # App-wide "computing" cue: Shiny's built-in busy-indicator page pulse
    # (a thin animated bar that sweeps across the top while the session is
    # recalculating). Spinners are disabled so we keep only the top sweep;
    # its look is themed in inst/app/www/style.css via the --shiny-pulse-*
    # custom properties. Shiny renders it at z-index 9999, above the navbar.
    useBusyIndicators(spinners = FALSE, pulse = TRUE),

    # Bootstrap 5 theme with INRAE styling
    page_navbar(
      id = "tabs",
      title = span(
        "Explore Metabar ",
        tags$small(
          as.character(utils::packageVersion("ExploreMetabar")),
          style = "opacity: 0.7;"
        )
      ),
      # flatly base, with colors + typography driven by inst/_brand.yml
      # (single source of truth; falls back to plain flatly if file missing).
      theme = bs_theme(
        version = 5,
        preset = "flatly",
        brand = if (nzchar(brand_path)) brand_path else FALSE
      ),
      fillable = TRUE,
      window_title = NA,

      # Navigation items for each module
      nav_panel(
        title = span(bs_icon("1-circle"), " Input Data"),
        value = "data_loading",
        mod_data_loading_ui("data_loading_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("1-circle"), " Community Composition"),
        value = "tab_compo",
        mod_compo_ui("compo_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("2-circle"), " Alpha diversity"),
        value = "tab_alpha",
        mod_alpha_ui("alpha_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("3-circle"), " Beta diversity"),
        value = "tab_beta",
        mod_beta_ui("beta_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("4-circle"), " Boxplot/Tests"),
        value = "tab_boxplot",
        mod_taxaboxplot_ui("taxaboxplot_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("5-circle"), " Heatmap"),
        value = "heatmap",
        mod_heatmap_ui("heatmap_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("6-circle"), " Differential Analysis"),
        value = "tab_diff",
        mod_diffanalysis_ui("diffanalysis_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("7-circle"), " ASVenn"),
        value = "tab_asvenn",
        mod_asvenn_ui("asvenn_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("8-circle"), " Cluster Analysis"),
        value = "cluster",
        mod_cluster_ui("cluster_ui_1")
      ),

      nav_panel(
        title = span(bs_icon("9-circle"), " MixOmics - SPLS-DA"),
        value = "tab_mixomics",
        mod_mixomics_ui("mixomics_1")
      ),

      # External links + INRAE branding, pushed to the right edge
      nav_spacer(),
      nav_item(tags$a(
        tags$img(src = SK8img, alt = "SK8", height = "32px"),
        href = "https://sk8.inrae.fr/", target = "_blank",
        class = "nav-link py-1", title = "SK8 hosting platform (home)"
      )),
      nav_item(tags$a(
        bs_icon("git"), " Source",
        href = "https://forge.inrae.fr/umrf/exploremetabar", target = "_blank",
        class = "nav-link", title = "GitLab repository"
      )),
      nav_item(tags$a(
        bs_icon("bug"), " Issues",
        href = "https://forge.inrae.fr/umrf/ExploreMetabar/-/issues", target = "_blank",
        class = "nav-link", title = "Report an issue on GitLab"
      )),
      nav_item(tags$a(
        tags$img(src = INRAEimg, alt = "INRAE", height = "32px"),
        href = "https://www.inrae.fr", target = "_blank",
        class = "nav-link py-1", title = "INRAE"
      ))
    )
  )
}

#' @import shiny
golem_add_external_resources <- function(){
  addResourcePath(
    'www', system.file('app/www', package = 'ExploreMetabar')
  )

  tags$head(
    golem::activate_js(),
    golem::favicon(),
    # App-wide non-blocking "computing" bar (keyed on Shiny's shiny-busy class).
    tags$link(rel = "stylesheet", type = "text/css", href = "www/style.css")
  )
}