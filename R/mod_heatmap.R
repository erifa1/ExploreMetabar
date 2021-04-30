#' heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_heatmap_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' heatmap Server Function
#'
#' @noRd 
mod_heatmap_server <- function(input, output, session, r){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_heatmap_ui("heatmap_ui_1")
    
## To be copied in the server
# callModule(mod_heatmap_server, "heatmap_ui_1")
 
