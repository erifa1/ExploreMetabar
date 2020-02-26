# Module UI
  
#' @title   mod_metadata_subset_ui and mod_metadata_subset_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_metadata_subset
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList 
mod_metadata_subset_ui <- function(id){
  ns <- NS(id)
  tagList(
  
  )
}
    
# Module Server
    
#' @rdname mod_metadata_subset
#' @export
#' @keywords internal
    
mod_metadata_subset_server <- function(input, output, session){
  ns <- session$ns
}
    
## To be copied in the UI
# mod_metadata_subset_ui("metadata_subset_ui_1")
    
## To be copied in the server
# callModule(mod_metadata_subset_server, "metadata_subset_ui_1")
 
