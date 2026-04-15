#' @import shiny
#' @import rhdf5
#' @importFrom thematic thematic_shiny


app_server <- function(input, output,session) {
  options(shiny.maxRequestSize=30*1024^2)

  # Enable thematic for automatic plot theming
  thematic::thematic_shiny()

  r <- reactiveValues(
    tabs = reactiveValues(),
    fdata = NULL
  )

  observe({
    r$tabs$tabselected <- input$tabs
  })


  # List the first level modules here
  mod_data_loading_server("data_loading_ui_1", r = r)
  mod_compo_server("compo_ui_1", r = r)
  mod_alpha_server("alpha_ui_1", r = r)
  mod_beta_server("beta_ui_1", r = r)
  mod_taxaboxplot_server("taxaboxplot_ui_1", r = r)
  mod_diffanalysis_server("diffanalysis_ui_1", r = r)
  mod_asvenn_server("asvenn_ui_1", r = r)
  mod_heatmap_server("heatmap_ui_1", r = r)
  mod_cluster_server("cluster_ui_1", r = r)
  mod_mixomics_server("mixomics_1", r = r)

}