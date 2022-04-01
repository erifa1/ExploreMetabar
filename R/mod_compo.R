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
#' @importFrom plotly plotlyOutput
mod_compo_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      infoBox("",
              "Use phyloseq object without taxa merging step.",
              icon = icon("info-circle"), fill=TRUE, width = 10),

      box(
        # radioButtons(
        #   ns("compo_norm_bool"),
        #   label = "Use normalized data",
        #   inline = TRUE,
        #   choices = list(
        #     "Raw" = 0 ,
        #     "Normalized" = 1
        #   ), selected = 0
        # ),

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
        numericInput(ns("topTax"), "Number of top taxa to plot:", 10, min = 1, max = NA),
        radioButtons(ns("radio1"), label = ("Plot display:"), choices = list("Default" = 1, "Splitted groups" = 2, "Merge samples" = 3),
        selected = 1),

        actionButton(ns("go1"), "Run Composition Plot", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
        ),

      box(
        downloadButton(outputId = ns("DLcompo2"), label = "Download plot"),
        plotlyOutput(ns("compo2")),
        title = "Relative abundance:", width = 12, status = "primary", solidHeader = TRUE),
      # box(plotlyOutput(ns("compo3")),
      #     title = "VST Normalized abundance:", width = 12, status = "primary", solidHeader = TRUE),
      box(
        downloadButton(outputId = ns("DLcompo1"), label = "Download plot"),
        plotlyOutput(ns("compo1")),
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
#' @importFrom plotly renderPlotly
#' @importFrom reshape2 melt

mod_compo_server <- function(input, output, session, r = r){
  ns <- session$ns
  observeEvent(r$tabs$tabselected, {
    if(r$tabs$tabselected=='tab_compo' && !isTruthy(r$phyloseq_filtered())){
      shinyalert(title = "Oops", text="Phyloseq object not present. Return to input data and validate all steps.", type='error')
    }
  })

  observeEvent(input$go1,{
    if(input$RankCompo==''){
      shinyalert(title = "Oops", text="You must provide a rank to plot.", type='error')
    }
  })

  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "RankCompo",
                      choices = phyloseq::rank_names(r$phyloseq_filtered()),
                      selected = r$rank_glom())
    updateSelectInput(session, "Ord1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })

  compo <- eventReactive(input$go1, {
    cat(file=stderr(),'Creating plots...',"\n")
    req(input$topTax, input$Ord1, input$RankCompo, r$phyloseq_filtered(), r$phyloseq_filtered_norm) #input$compo_norm_bool,
    LL=list()
    # if(input$compo_norm_bool==0){
      Fdata <- r$phyloseq_filtered()
      print(Fdata)
    # }
    # if(input$compo_norm_bool==1){
      # Fdatanorm <- Fdata
      # otable <- Fdatanorm@otu_table@.Data+1
      # otableVST <- DESeq2::varianceStabilizingTransformation(otable, fitType='local')
      # Fdatanorm@otu_table@.Data <- otableVST
      # print(Fdatanorm)
    # }

    withProgress({
      if(input$radio1 == 3){
        cat(file=stderr(),'Merged...',"\n")
        Fdata <- phyloseq::merge_samples(Fdata, group=input$Ord1, fun=mean)
        sample_data(Fdata)[[input$Ord1]] <- sample_names(Fdata)
        split1 = FALSE

      }else{
        if(input$radio1 == 1 | input$radio1 == 3){split1 = FALSE}else{split1 = TRUE}
        cat(file=stderr(),'Std...',"\n")
      }

      LL$p1 = bars_fun(Fdata, rank=input$RankCompo, top = input$topTax, Ord1 = input$Ord1, relative = FALSE, outfile = NULL, split = split1, autoorder = FALSE, verbose = FALSE, split_sid_order = FALSE, ylab = "Raw abundance")
      LL$p2 = bars_fun(Fdata, rank=input$RankCompo, top = input$topTax, Ord1 = input$Ord1, relative = TRUE, outfile = NULL, split = split1, autoorder = FALSE, verbose = FALSE, split_sid_order = FALSE, ylab = "Relative abundance")

      LL

    }, message="Processing, please wait...")

  })



  output$compo1 <- renderPlotly({
    LL <- compo()
    LL$p1 %>% config(toImageButtonOptions = list(format = "svg"))
  })

  output$compo2 <- renderPlotly({
    LL <- compo()
    LL$p2 %>% config(toImageButtonOptions = list(format = "svg"))
  })

  # output$compo3 <- renderPlotly({
  #   LL <- compo()
  #   LL$p3
  # })

  output$totalsum1 <- renderPrint({
      Fdata <- r$phyloseq_filtered()
    print(sample_sums(Fdata))
  })

  output$DLcompo2 <- downloadHandler(
    filename = "plot_compo_relative_abundance.html",
    content = function(file) {
      req(compo())
      LL <- compo()
      plot1 <- LL$p2 %>% config(toImageButtonOptions = list(format = "svg"))
      saveWidget(plot1, file= file)
    }
  )

  output$DLcompo1 <- downloadHandler(
    filename = "plot_compo_raw_abundance.html",
    content = function(file) {
      req(compo())
      LL <- compo()
      plot1 <- LL$p1 %>% config(toImageButtonOptions = list(format = "svg"))
      saveWidget(plot1, file= file)
    }
  )


}


## To be copied in the UI
# mod_compo_ui("compo_ui_1")

## To be copied in the server
# callModule(mod_compo_server, "compo_ui_1")
