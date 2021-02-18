#' source_tracker UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_source_tracker_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      box(
        fluidRow(
          selectInput(
            ns("src_fact1"),
            label = "Select factor to test: ",
            choices = ""
          )
        ),
        fluidRow(
          checkboxGroupInput(
            ns('sources_box'), label = 'Choose your sources', choices=''
          )
          radioButtons(
            ns('sink_radio'), label='Choose your sink', choices = ''
          ),
        )
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ) 
    )
  )
}
    
#' source_tracker Server Function
#'
#' @noRd 
mod_source_tracker_server <- function(input, output, session, r = r){
  ns <- session$ns

  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "src_fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })
  
  
  observe({
    req(r$sdat(), input$src_fact1)
  })
  
  observe({
    req(r$sdat(), input$src_fact1, input$sources_box)
    tmp <- r$sdat()
    uniq.name <- unique(tmp[,input$src_fact1])
    indices <- 1:length(uniq.name)
    print(uniq.name)
    print(indices)
    choices.list <- as.list(setNames(indices, uniq.name))
    print(choices.list)
    updateRadioButtons(session, 'sink_radio',
                      choices = choices.list)
  })
  
}
    
## To be copied in the UI
# mod_source_tracker_ui("source_tracker_ui_1")
    
## To be copied in the server
# callModule(mod_source_tracker_server, "source_tracker_ui_1")
 
