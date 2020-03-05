# Module UI
  
#' @title   mod_diffanalysis_ui and mod_diffanalysis_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_diffanalysis
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList 
mod_diffanalysis_ui <- function(id){
  ns <- NS(id)
  tagList(
    selectInput(
      ns("Fact1"),
      label = "Select factor to test: ",
      choices = ""
    ),
    verbatimTextOutput(ns("print1")),
    
    selectInput(
      ns("Cond1"),
      label = "Select Condition1 to test: ",
      choices = ""
    )
  )
}
    
# Module Server
    
#' @rdname mod_diffanalysis
#' @export
#' @keywords internal
    
mod_diffanalysis_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  observe({
    updateSelectInput(session, "Fact1",
                      choices = r$data16S()@sam_data@names)
    
    # updateSelectInput(session, "Cond1",
    #                   choices = Cond1())
  })
  
  output$print1 <- renderPrint({
    if(is.null(input$Fact1)){return(NULL)}
    print(input$Fact1)
  })
  
  # output$Cond1 <- renderPrint({
  #   if(is.null(input$Fact1)){return(NULL)}
  #   cond1 = levels(r$data16S()@sam_data[,input$Fact1])
  #   print(cond1)
  #   cond1
  # })
  
  output$Cond1 <- reactive({
    if(is.null(input$Fact1)){return(NULL)}
      cond1 = levels(r$data16S()@sam_data[,input$Fact1])
      print(input$Fact1)
      print(cond1)
      cond1
  })
  
}
    
## To be copied in the UI
# mod_diffanalysis_ui("diffanalysis_ui_1")
    
## To be copied in the server
# callModule(mod_diffanalysis_server, "diffanalysis_ui_1")
 
