# Module UI

#' @title   mod_compo_ui and mod_compo_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_compo
#'
#' @keywords internal
#' @export
#' @importFrom shiny NS tagList
#' @import plotly
mod_compo_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(

      infoBox("",
              "Use phyloseq object without taxa merging step.",
              icon = icon("info-circle"), fill=TRUE, width = 10),

      box(
        selectInput(
          ns("RankCompo"),
          label = "Select rank to plot: ",
          choices = ""
        ),

        selectInput(
          ns("Ord1"),
          label = "Select variable to order/split samples (X axis): ",
          choices = ""
        ),

        selectInput(
          ns("Fact1"),
          label = "Select variable for changing X axis tick labels (for unsplitted plots):",
          choices = ""
        ),
        numericInput(ns("topTax"), "Number of top taxa to plot:", 10, min = 1, max = NA),
        checkboxInput(ns("checkbox1"), label = "Splitted plots", value = FALSE),
        actionButton(ns("go1"), "Run Composition Plot", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(plotlyOutput(ns("compo2")),
          title = "Relative abundance:", width = 12, status = "primary", solidHeader = TRUE),
      box(plotlyOutput(ns("compo1")),
          title = "Raw abundance:", width = 12, status = "primary", solidHeader = TRUE),
      box(verbatimTextOutput(ns("totalsum1")),
          title = "Total sum per samples:", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE)
  )
  )
}

# Module Server

#' @rdname mod_compo
#' @export
#' @keywords internal
#' @import plotly
#' @importFrom microbiome aggregate_top_taxa
#' @importFrom reshape2 melt
#' @importFrom ranomaly bars_fun

mod_compo_server <- function(input, output, session, r = r){
  ns <- session$ns

  observe({
    updateSelectInput(session, "RankCompo",
                      choices = rank_names(r$data16S()),
                      selected = "Genus")
    updateSelectInput(session, "Fact1",
                      choices = r$data16S()@sam_data@names)
    updateSelectInput(session, "Ord1",
                      choices = r$data16S()@sam_data@names)
  })

  compo <- eventReactive(input$go1, {
    LL=list()

    # print("compo")
    # print(r$rowselect())
    # print(r$data16S())
    # print(r$asvselect())
    withProgress({
      Fdata <- prune_samples(sample_names(r$data16S())[r$rowselect()], r$data16S())
      Fdata <- prune_taxa(taxa_sums(Fdata) > 0, Fdata)
      if( r$RankGlom() == "ASV"){
        Fdata <- prune_taxa(r$asvselect(), Fdata)
      }

      LL$p1 = bars_fun(data = Fdata, top = input$topTax, Ord1 = input$Ord1, Fact1 = input$Fact1, rank=input$RankCompo, relative = FALSE, outfile=NULL, split = input$checkbox1) %>% config(toImageButtonOptions = list(format = "svg"))
      

      LL$p2 = bars_fun(data = Fdata, top = input$topTax, Ord1 = input$Ord1, Fact1 = input$Fact1, rank=input$RankCompo, relative = TRUE, outfile=NULL, split = input$checkbox1) %>% config(toImageButtonOptions = list(format = "svg"))

      LL
    }, message="Processing, please wait...")

  })

  output$compo1 <- renderPlotly({
    LL <- compo()
    LL$p1
  })

  output$compo2 <- renderPlotly({
    LL <- compo()
    LL$p2
  })

  output$totalsum1 <- renderPrint({
    Fdata <- prune_samples(sample_names(r$data16S())[r$rowselect()], r$data16S())
    Fdata <- prune_taxa(taxa_sums(Fdata) > 0, Fdata)
    if( r$RankGlom() == "ASV"){
      Fdata <- prune_taxa(r$asvselect(), Fdata)
    }

    print(sample_sums(Fdata))
  })


}


## To be copied in the UI
# mod_compo_ui("compo_ui_1")

## To be copied in the server
# callModule(mod_compo_server, "compo_ui_1")
