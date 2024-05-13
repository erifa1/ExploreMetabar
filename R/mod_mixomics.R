#' mixomics UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_mixomics_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' mixomics Server Functions
#'
#' @noRd 
mod_mixomics_server <- function(input, output, session, r){
    ns <- session$ns
}
    
## To be copied in the UI
# mod_mixomics_ui("mixomics_1")
    
## To be copied in the server
# mod_mixomics_server("mixomics_1")
