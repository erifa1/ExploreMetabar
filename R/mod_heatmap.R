#' heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' @import phyloseq
#' @import ggplot2
mod_heatmap_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        infoBox("",
                HTML(
                  paste(
                    h4("New Heatmap module."),
                    h5("Ecology organized heatmap. Row an columns are organized using ordination methods (NMDS, PCA) and computed on ecological distances (bray, unifrac, jaccard)."),
                    h5("Please cite: "), 
                    h5(tags$a(href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-45")) )
                ),
                icon = icon("info-circle"), fill=TRUE, width = 12
        )
      ),
      fluidRow(
        box(title = "Settings", width = 12, status = "warning", solidHeader = TRUE,
          selectInput(
            ns("src_fact1"),
            label = "Sample label",
            choices = ""
          ),
          selectInput(
            ns("ord.method"),
            label = "Ordiantion method",
            choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
            selected = 'NMDS'
          ),
          selectInput(
            ns("dist"),
            label = "Ecological distance method ",
            choices = c("bray", "jaccard", "dpcoa", "unifrac", "wunifrac"),
            selected = 'NMDS'
          )
        )
      ),
      fluidRow(
        box(title = "Heatmap", width=12, height = 12, status = 'primary', solidHeader = FALSE,
          shinycustomloader::withLoader(
            plotOutput(ns('heatmap_plot')),
            type = "html", loader = "loader1"
          )
        )
      )
    )
  )
}
    
#' heatmap Server Function
#'
#' @noRd 
mod_heatmap_server <- function(input, output, session, r){
  ns <- session$ns
  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "src_fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })
  
  get_heatmap <- reactive({
    req(input$src_fact1, r$phyloseq_filtered())
    phyloseq::plot_heatmap(r$phyloseq_filtered(), method = input$ord.method, distance = input$dist, sample.label = input$src_fact1)
  })
  
  output$heatmap_plot <- renderPlot(
    withProgress(message = 'Computing heatmap...',{
      get_heatmap()
    })
    
  )
}
    
## To be copied in the UI
# mod_heatmap_ui("heatmap_ui_1")
    
## To be copied in the server
# callModule(mod_heatmap_server, "heatmap_ui_1")
 
