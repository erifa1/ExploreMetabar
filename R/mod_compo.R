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
#' @import bslib
#' @import bsicons

mod_compo_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Info card
    card(
      full_screen = TRUE,
      card_header(class = "bg-info"),
      "Use the phyloseq object without performing the taxa merging step."
    ),

    br(),

    # Settings card
    card(
      full_screen = TRUE,
      card_header(bs_icon("gear"), " Settings"),
      selectInput(
        ns("RankCompo"),
        label = "Select rank to plot: ",
        choices = ""
      ),
      shinyWidgets::pickerInput(
        ns("Ord1"),
        label = "Select one or more categorial variable to order/split samples (X axis): ",
        choices = "",
        multiple = TRUE
      ),
      numericInput(ns("topTax"), "Number of top taxa to plot:", 10, min = 1, max = NA),
      radioButtons(ns("radio1"), label = ("Plot display:"), choices = list("Default" = 1, "Splitted groups" = 2, "Merge samples" = 3),
      selected = 1, inline = TRUE),
      checkboxInput(ns("autoorder1"), "Autoorder samples", value = TRUE),
      actionButton(ns("go1"), "Run Composition Plot", icon = bs_icon("play"))
    ),

    br(),

    # Relative abundance plot
    card(
      full_screen = TRUE,
      card_header(bs_icon("pie-chart"), " Relative abundance"),
      downloadButton(outputId = ns("DLcompo2"), label = "Download plot"),
      plotlyOutput(ns("compo2"), height = "600px")
    ),

    br(),

    # Raw abundance plot
    card(
      full_screen = TRUE,
      card_header(bs_icon("bar-chart"), " Raw abundance"),
      downloadButton(outputId = ns("DLcompo1"), label = "Download plot"),
      plotlyOutput(ns("compo1"), height = "600px")
    ),

    br(),

    # Total sum info
    card(
      full_screen = TRUE,
      card_header(bs_icon("info-circle"), " Total sum per samples"),
      accordion(
        open = FALSE,
        accordion_panel(
          "Details",
          verbatimTextOutput(ns("totalsum1"))
        )
      )
    )
  )
}

# Module Server

#' @rdname mod_compo
#' @export
#' @keywords internal
#' @importFrom plotly renderPlotly
#' @importFrom tidyr pivot_longer
#' @import futile.logger
#' @importFrom futile.logger flog.info

mod_compo_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns

  observeEvent(input$go1,{
    if(input$RankCompo==''){
      shinyalert(title = "Oops", text="You must provide a rank to plot.", type='error')
      req(FALSE)
    }
  })

  observe({
    req(r$phyloseq_filtered())
    ranks1 <- phyloseq::rank_names(r$phyloseq_filtered())
    updateSelectInput(session, "RankCompo",
                      choices = ranks1,
                      selected = ranks1[length(ranks1)])
    shinyWidgets::updatePickerInput(session, "Ord1",
                      choices = r$var_list(),selected = r$var_list()[1])
  })

  observeEvent(input$Ord1, {
    if ("sample.id" %in% input$Ord1) {
      updateRadioButtons(session, "radio1",
        choices = list("Default" = 1),
        selected = 1, inline = TRUE)
    } else {
      updateRadioButtons(session, "radio1",
        choices = list("Default" = 1, "Splitted groups" = 2, "Merge samples" = 3),
        selected = input$radio1, inline = TRUE)
    }
  })

  factor_input <- reactive({ input$Ord1 })
  get_meta_col  <- make_get_meta_col(factor_input, r)
  local_metadata <- make_local_metadata(factor_input, get_meta_col, r)
  local_physeq  <- make_local_physeq(local_metadata, r)


  compo <- eventReactive(input$go1, {
    flog.info('compo - Creating plots...')
    req(input$topTax, get_meta_col(), input$RankCompo, local_physeq())

    LL=list()
    Fdata <- local_physeq()

    withProgress({
      if(input$radio1 == 3){
        flog.info('compo - Merged...')
        Fdata <- phyloseq::merge_samples(Fdata, group=get_meta_col(), fun=mean)
        sample_data(Fdata)[[get_meta_col()]] <- sample_names(Fdata)
        split1 = FALSE
      } else{
        if(input$radio1 == 1 | input$radio1 == 3){split1 = FALSE}else{split1 = TRUE}
        flog.info('compo - Std...')
      }

      LL$p1 = bars_fun(Fdata, rank=input$RankCompo, top = input$topTax, Ord1 = get_meta_col(), relative = FALSE, outfile = NULL, split = split1, autoorder = input$autoorder1, verbose = FALSE, split_sid_order = FALSE, ylab = "Raw abundance", pal = c(r$factor_colors()[[get_meta_col()]], r$taxa_colors()[[input$RankCompo]], 'Other' = 'grey'))
      LL$p2 = bars_fun(Fdata, rank=input$RankCompo, top = input$topTax, Ord1 = get_meta_col(), relative = TRUE, outfile = NULL, split = split1, autoorder = input$autoorder1, verbose = FALSE, split_sid_order = FALSE, ylab = "Relative abundance", pal = c(r$factor_colors()[[get_meta_col()]], r$taxa_colors()[[input$RankCompo]], 'Other' = 'grey'))

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
  })
}

## To be copied in the UI
# mod_compo_ui("compo_ui_1")

## To be copied in the server
# mod_compo_server("compo_ui_1", r = r)