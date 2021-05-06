#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import shinyjs

app_ui <- function() {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # shinyjs::useShinyjs(), # used for onclick function from shinyjs package.
    # List the first level UI elements here
    dashboardPage(skin = "red",
                  dashboardHeader(
                                  title = "Explore Metabar",
                                  tags$li(class="dropdown",tags$a(icon("gitlab"),href="https://forgemia.inra.fr/umrf/exploremetabar")), #, target="_blank"
                                  tags$li(class="dropdown",tags$a(icon("clinic-medical"),href="https://forgemia.inra.fr/umrf/exploremetabar/-/issues"))
                                  ),

                  dashboardSidebar(
                    sidebarMenu(
                      id="tabs",
                      style = "position: fixed; overflow: visible",
                      menuItem("Input Data", tabName= 'data_loading', icon=icon("diagnoses")),
                      menuItem("Community Composition", tabName = "tab_compo", icon = icon("chart-pie")),
                      menuItem("Alpha diversity", tabName = "tab_alpha", icon = icon("chart-bar")),
                      menuItem("Beta diversity ", tabName = "tab_beta", icon = icon("chart-bar")),
                      menuItem("Boxplot/Tests", tabName = "tab_boxplot", icon = icon("microscope")),
                      menuItem('Heatmap', tabName = 'heatmap', icon = icon('chart-bar')),
                      menuItem("Differential Analysis", tabName = "tab_diff", icon = icon("microscope")),
                      menuItem("ASVenn", tabName = "tab_asvenn", icon = icon("microscope")),
                      menuItem("SourceTracker", tabName = "source_tracker", icon = icon("sourcetree"))
                      # menuItem("PLN network", tabName = "tab_networkpln", icon = icon("project-diagram")),
                      # menuItem("DiffExplore", tabName = "tab_diffexplore", icon = icon("leaf"))
                    )
                  ),

                  dashboardBody(

                    tabItems(
                      tabItem(tabName = 'data_loading',
                              mod_data_loading_ui("data_loading_ui_1")
                      ),
                      tabItem(tabName = "tab_compo",
                              mod_compo_ui("compo_ui_1")
                      ),
                      tabItem(tabName = "tab_alpha",
                              mod_alpha_ui("alpha_ui_1")
                      ),
                      tabItem(tabName = "tab_beta",
                              mod_beta_ui("beta_ui_1")
                      ),
                      tabItem(tabName = "tab_boxplot",
                              mod_taxaboxplot_ui("taxaboxplot_ui_1")
                      ),
                      tabItem(tabName = "tab_diff",
                              mod_diffanalysis_ui("diffanalysis_ui_1")
                      ),
                      tabItem(tabName = "tab_asvenn",
                              mod_asvenn_ui("asvenn_ui_1")
                      ),
                      tabItem(tabName = "source_tracker",
                              mod_source_tracker_ui("source_tracker_ui_1")
                      ),
                      tabItem(tabName = "heatmap",
                              mod_heatmap_ui("heatmap_ui_1")
                      )
                      # tabItem(tabName = "tab_networkpln",
                      #         mod_networkpln_ui("networkpln_ui_1")
                      # ),
                      # tabItem(tabName = "tab_diffexplore",
                      #         mod_diffexplore_ui("diffexplore_ui_1")
                      # )
                    )
                  )

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
    # Add here all the external resources
    # If you have a custom.css in the inst/app/www
    # Or for example, you can add shinyalert::useShinyalert() here
    #tags$link(rel="stylesheet", type="text/css", href="www/custom.css")
  )
}
