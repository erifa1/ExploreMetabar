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
                      tags$div(
                        title = "RData where 'data' is a phyloseq object.",
                        fileInput("fileRData",
                                  label = "RData : ",
                                  placeholder = "data.RData")
                      ),
                      selectInput(
                        "RankGlom",
                        label = "Select rank to glom : ",
                        choices = "",
                        selected = 1
                      ),
                      radioButtons(
                        "NORM",
                        label = "Normalization : ",
                        inline = TRUE,
                        choices = list(
                          `No Norm` = "Raw",
                          `CLR` = "CLR norm.",
                          `TSS` = "TSS norm.",
                          `VST` = "VST norm."
                        ), selected = "CLR norm."
                      ),
                      # menuItem("Transform phyloseq object", tabName = "Transform", icon = icon("leaf")),
                      menuItem("Metadatas/Subset", tabName = "Metadatas", icon = icon("leaf"))   
                    )
                  ),
                  
                  dashboardBody(
                    
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
