#' source_tracker UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session,r Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' 
mod_source_tracker_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      box(
        fluidPage(
          fluidRow(
            selectInput(
              ns("src_fact1"),
              label = "Select factor to test: ",
              choices = ""
            )
          ),
          fluidRow(
            uiOutput(ns("sources_box")),
            uiOutput(ns("sink_radio")),
          ),
          fluidRow(
            actionButton(ns('src_go'), label='Launch', style = "material-circle", color = "primary", icon = icon('play'))
          )
        ),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      )
    )
  )
}
    
#' source_tracker Server Function
#'
#' @noRd 
#' 
#' 
mod_source_tracker_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  r_values <- reactiveValues(factor_list=NULL)
  
  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "src_fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })
  
  output$sources_box <- renderUI({
    req(r$sdat(), input$src_fact1)
    levels <- na.omit(levels(r$sdat()[,input$src_fact1]))
    checkboxGroupButtons(ns('sources_box'), "Choose your sources:", choices=levels, justified = TRUE, checkIcon = list(yes = icon("ok", lib = "glyphicon")))
  })
  
  output$sink_radio <- renderUI({
    req(input$sources_box)
    levels <- na.omit(levels(r$sdat()[,input$src_fact1]))
    ch <- dplyr::setdiff(levels, input$sources_box )
    radioGroupButtons(ns('sink_radio'), "Choose your sink", choices=ch, justified = TRUE, checkIcon = list(yes = icon("ok", lib = "glyphicon")))
  })
  
  
  launch_sourceTracker <- reactive({
    # phy_obj <- r$phyloseq_filtered()
    ff <- glue::glue("tt <- r$sdat()${input$src_fact1} %in% c(input$sources_box, input$sink_radio )")
    print(ff)
    eval(parse(text=ff))
  })
  
  srcTracker_reactive <- eventReactive(input$src_go, {
    cat(file=stderr(), 'launch sourceTracker...')
    launch_sourceTracker()
  })
  
}
    
## To be copied in the UI
# mod_source_tracker_ui("source_tracker_ui_1")
    
## To be copied in the server
# callModule(mod_source_tracker_server, "source_tracker_ui_1")
