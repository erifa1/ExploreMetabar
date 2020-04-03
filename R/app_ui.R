#' @import shiny
#' @import shinydashboard
app_ui <- function() {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here 
    dashboardPage(skin = "red",
                  dashboardHeader(title = "Explore Metabar"),
                  
                  dashboardSidebar(
                    sidebarMenu(
                      # menuItem("Transform phyloseq object", tabName = "Transform", icon = icon("leaf")),
                      menuItem("Metadatas/Subset", tabName = "input_select", icon = icon("leaf")),
                      menuItem("ASVtable", tabName = "export_asvtable", icon = icon("leaf")),
                      menuItem("Compo", tabName = "tab_compo", icon = icon("leaf")),
                      menuItem("Alpha", tabName = "tab_alpha", icon = icon("leaf")),
                      menuItem("Boxplot/Tests", tabName = "tab_boxplot", icon = icon("leaf")),
                      menuItem("Differential Analysis", tabName = "tab_diff", icon = icon("leaf")),
                      menuItem("DiffExplore", tabName = "tab_diffexplore", icon = icon("leaf"))
                    )
                  ),
                  
                  dashboardBody(
                    tabItems(
                      tabItem(tabName = "input_select",
                              mod_metadata_subset_ui("metadata_subset_ui_1")
                      ),
                      tabItem(tabName = "export_asvtable",
                              mod_export_asvtaxtable_ui("export_asvtaxtable_ui_1")
                      ),
                      tabItem(tabName = "tab_compo",
                              mod_compo_ui("compo_ui_1")
                      ),
                      tabItem(tabName = "tab_alpha",
                              mod_alpha_ui("alpha_ui_1")
                      ),
                      tabItem(tabName = "tab_boxplot",
                              mod_taxaboxplot_ui("taxaboxplot_ui_1")
                      ),
                      tabItem(tabName = "tab_diff",
                                 mod_diffanalysis_ui("diffanalysis_ui_1")
                      ),
                      tabItem(tabName = "tab_diffexplore",
                              mod_diffexplore_ui("diffexplore_ui_1")
                      )
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
