#' asvenn UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_asvenn_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      infoBox("",
              "Select conditions to highlight shared taxa",
              icon = icon("info-circle"), fill=TRUE, width = 10),
      box(
        selectInput(
          ns("Fact1"),
          label = "Select factor to test: ",
          choices = ""
        ),
        uiOutput(ns("lvls1")),
        numericInput(ns("minAb"), "Minimum raw abundance:", 1, min = 1, max = NA),
        actionButton(ns("go1"), "Run/Update ASVenn", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(
        plotOutput(ns("venn1"), height = "800px"),
        title = "Venn Diagram:", width = 12, status = "primary", solidHeader = TRUE, height = 900
      ),
      box(
        downloadButton(outputId = ns("otable_download"), label = "Download Table"),
        DT::dataTableOutput(ns("tabvenn1")),
        title = "Venn table:", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE
      )
    )
  )
}

#' asvenn Server Function
#'
#' @noRd
#'
#' @importFrom futile.logger flog.threshold
#' @importFrom grid grid.draw
#' @importFrom grDevices rainbow
#' @importFrom VennDiagram calculate.overlap
#' @importFrom VennDiagram venn.diagram
#' @importFrom ranomaly ASVenn_fun
#'

mod_asvenn_server <- function(input, output, session, r=r){
  ns <- session$ns

  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "Fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })

# output$lvls1 <- reactive({
#   level1 <- na.omit(levels(as.factor(sample_data(data)[,input$Fact1]@.Data[[1]])) )
#   level1
# })

  output$lvls1 = renderUI({
    req(input$Fact1, r$phyloseq_filtered())
    level1 <- na.omit(levels(as.factor(sample_data(r$phyloseq_filtered())[,input$Fact1]@.Data[[1]])) )
    checkboxGroupInput(ns("lvls1"), label = "Select up to 7 levels :",
                       choices = level1, inline = TRUE, selected = level1[1:3])

  })

  resVenn <- eventReactive(input$go1, {   #TF
    req(r$phyloseq_filtered())
    resVenn = ASVenn_fun(data = r$phyloseq_filtered(), output = NULL, rank = "ASV", column1 = input$Fact1, lvls = input$lvls1, shared = TRUE)

    resVenn
  })

  output$venn1 <- renderPlot({
    invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
    if(length(input$lvls1) >= 2 & length(input$lvls1) <= 7){
        resVenn()$venn_plot
    }else{showNotification("Choose 2 to 7 levels...", type="error", duration = 5)
          return(NULL)
          }
  })

  output$tabvenn1 <-DT::renderDataTable({
    resVenn()$TABf #tabvenn1()
  }, filter="top", options = list(scrollX = TRUE))

  output$otable_download <- downloadHandler(
    filename = "venn_table.csv",
    content = function(file) {
      req(tabvenn1())
      write.table(tabvenn1(), file, sep="\t", row.names=FALSE)
    }
  )



}

## To be copied in the UI
# mod_asvenn_ui("asvenn_ui_1")

## To be copied in the server
# callModule(mod_asvenn_server, "asvenn_ui_1")
