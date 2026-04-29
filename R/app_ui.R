#' @import shiny
#' @import bslib
#' @importFrom bsicons bs_icon
#' @importFrom base64enc dataURI


SK8img <- base64enc::dataURI(file=system.file(file.path('app/www', 'SK8.png'), package='ExploreMetabar'))
UCAimg <- base64enc::dataURI(file=system.file(file.path('app/www', 'uca2.png'), package='ExploreMetabar'))
MIGimg <- base64enc::dataURI(file=system.file(file.path('app/www', 'migale2.png'), package='ExploreMetabar'))

description <- read.dcf(system.file("DESCRIPTION", package = "ExploreMetabar"))
package_name <- description[1, "Package"]
package_version <- description[1, "Version"]


app_ui <- function() {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),

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
      theme = bs_theme(
        version = 5,
        preset = "flatly",
        primary = "#00a3a6",
        secondary = "#423089",
        success = "#9dc544",
        info = "#9ed6e3",
        warning = "#ed6e6c",
        danger = "#ed6e6c"
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

      # Navbar menu for external links
      nav_menu(
        title = "Links",
        nav_item(tags$a("SK8", href = "https://sk8.inrae.fr/", target = "_blank", class = "nav-link")),
        nav_item(tags$a(bs_icon("code"), " Source code", href = "https://forge.inrae.fr/umrf/exploremetabar", target = "_blank", class = "nav-link")),
        nav_item(tags$a(bs_icon("bug"), " Issues", href = "https://forge.inrae.fr/umrf/ExploreMetabar/-/issues", target = "_blank", class = "nav-link"))
      ),

      # Dark mode toggle
      nav_spacer(),
      nav_item(input_dark_mode())
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
    golem::favicon()
  )
}