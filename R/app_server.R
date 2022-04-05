#' @import shiny
#' @import rhdf5
options(shiny.maxRequestSize=30*1024^2)

app_server <- function(input, output,session) {

  r <- reactiveValues(
    tabs = reactiveValues()
  )

  observe({
    r$tabs$tabselected <- input$tabs
  })


  # List the first level callModules here
  callModule(mod_data_loading_server, "data_loading_ui_1", session=session, r=r)
  callModule(mod_compo_server, "compo_ui_1", session = session, r = r)
  callModule(mod_alpha_server, "alpha_ui_1", session = session, r = r)
  callModule(mod_beta_server, "beta_ui_1", session = session, r = r)
  callModule(mod_taxaboxplot_server, "taxaboxplot_ui_1", session = session, r = r)
  callModule(mod_diffanalysis_server, "diffanalysis_ui_1", session = session, r = r)
  callModule(mod_asvenn_server, "asvenn_ui_1", session = session, r = r)
  callModule(mod_source_tracker_server, "source_tracker_ui_1", session = session, r = r)
  callModule(mod_heatmap_server, "heatmap_ui_1", session=session,r=r)
  
  # callModule(mod_diffexplore_server, "diffexplore_ui_1", session = session, r = r)
  # callModule(mod_networkpln_server, "networkpln_ui_1", session = session, r = r)

}
