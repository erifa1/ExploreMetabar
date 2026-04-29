# Module UI

#' @title   mod_alpha_ui and mod_alpha_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_alpha
#'
#' @keywords internal
#' @export
#' @importFrom shiny NS tagList
#' @importFrom plotly plotlyOutput
#' @importFrom shinyalert shinyalert useShinyalert
#' @import bslib
#' @import bsicons
mod_alpha_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "Settings",
      open = "desktop",
      width = "350px",
      htmltools::p(
        "Use the phyloseq object without performing the taxa merging step.",
        style = "font-size: 0.9em; color: grey;"
      ),
      tags$hr(),
      uiOutput(ns('ui_alpha_factor')),
      checkboxInput(ns("checkbox1"), label = "Automatic order factor", value = TRUE),
      actionButton(ns("launch_alpha"), "Run Alpha Diversity", icon = bs_icon("play-fill"), class = "btn-primary w-100 btn-lg")
    ),

    navset_card_underline(
      title = "Alpha Indexes",
      full_screen = TRUE,
      nav_panel(
        title = "Alpha indexes table",
        icon = bs_icon("table"),
        DT::dataTableOutput(ns("alphaout")),
        downloadButton(outputId = ns("alpha_download"), label = "Download Table")
      ),
      nav_panel(
        title = "Alpha indexes by group",
        icon = bs_icon("people"),
        DT::dataTableOutput(ns("alphagrp")),
        downloadButton(outputId = ns("alphagrp_download"), label = "Download group Table")
      ),
      nav_panel(
        title = "Boxplot",
        icon = bs_icon("bar-chart"),
        radioButtons(ns("metrics"), "Choose one index:", inline = TRUE,
                     choices = list("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"),
                     selected = c("Shannon")
        ),
        plotly::plotlyOutput(ns("boxplot"))
      ),
      nav_panel(
        title = "Statistics",
        icon = bs_icon("calculator"),
        uiOutput(ns('anovaBox'))
      )
    )
  )
}

# Module Server

#' @rdname mod_alpha
#' @export
#' @keywords internal
#' @import phyloseq
#' @importFrom DT renderDataTable
#' @importFrom plotly renderPlotly config layout
#' @importFrom agricolae HSD.test
#' @importFrom gtools mixedsort
#' @import futile.logger
#' @importFrom futile.logger flog.info

mod_alpha_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns

  factor_input <- reactive({ input$Fact1 })
  get_meta_col  <- make_get_meta_col(factor_input, r)
  local_metadata <- make_local_metadata(factor_input, get_meta_col, r)
  isNumFactor   <- make_is_num_factor(get_meta_col, local_metadata)
  local_physeq  <- make_local_physeq(local_metadata, r)


  observe({
    req(r$phyloseq_filtered(), r$var_list())
    shinyWidgets::updatePickerInput(session, "Fact1",
                      choices = r$var_list())
  })

  output$ui_alpha_factor <- renderUI({
        shinyWidgets::pickerInput(
          ns("Fact1"),
          label = "Select one or more factor to test (when multiple selection, qualitative factors are concatenated): ",
          choices = r$var_list(),
          selected = r$var_list()[1],
          multiple = FALSE
        )
  })


  alpha1 <- eventReactive(input$launch_alpha, {
    withProgress(message = 'Computing alpha diversity tables', min=0, max=10, value = 0,{
      flog.info('computing alpha1...')
      req(local_physeq())

      data <- local_physeq()
      setProgress(value = 5, detail = 'estimate richness')
      alphatab <- estimate_richness(data, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                                                       "InvSimpson") )
      row.names(alphatab) = sample_names(data)

      LL=list()
      LL$alphatab = as.data.frame(alphatab)
      LL$data = data
      flog.info('computing alpha1 done.')
      setProgress(value = 10, detail = 'done')
      return(LL)
    })
  })


  output$alphaout <- DT::renderDataTable({
    LL = alpha1()
    round_df(LL$alphatab, 4)
  }, options = list(
    pageLength = 5,
    scrollX = TRUE
  ))


  alphagrp_table <- eventReactive(input$launch_alpha, {
    req(alpha1(), local_metadata(), get_meta_col())
    withProgress(message = 'Group table', min=0, max=10, value = 0,{
    alpha.table <- alpha1()$alphatab

    alpha.table =  tibble::rownames_to_column(alpha.table)
    alpha.table <- dplyr::left_join(local_metadata(), alpha.table, by = c( "sample.id" = "rowname"))

    alpha.table[,'rowname'] <- NULL

    if(!isNumFactor()){
      alpha.table <- alpha.table %>%
        group_by_at(get_meta_col()) %>%
        summarise(
          tibble(
            across(where(is.numeric), ~mean(.x), .names = "mean_{.col}"),
            across(where(is.numeric), ~median(.x), .names = "median_{.col}")
          )
        )
    }
    setProgress(value = 10, detail = 'done')
    return(alpha.table)
    })
  })

  output$alphagrp <- DT::renderDataTable({
    round_df(as.data.frame(alphagrp_table()), 4)
  }, options = list(
    pageLength = 5,
    scrollX = TRUE
  ))

  output$alphagrp_download <- downloadHandler(
    filename = "alphagrp_index.csv",
    content = function(file) {
      write.table(alphagrp_table(), file, sep="\t", col.names=NA)}
  )

  output$alpha_download <- downloadHandler(
    filename = "alpha_index.csv",
    content = function(file) {
      LL = alpha1()
      write.table(LL$alphatab, file, sep="\t", col.names=NA)}
  )


  boxtab <- eventReactive(input$launch_alpha, {
    req(r$sdat(), input$Fact1)
    withProgress(message = 'Boxplot table', min=0, max=10, value = 0,{
    flog.info('boxtab function...')
    LL = alpha1()
    alphatab =  tibble::rownames_to_column(LL$alphatab)

    boxtab <- dplyr::left_join(r$sdat(), alphatab, by = c('sample.id' = "rowname"))

    if(! is.numeric(boxtab[, input$Fact1])){
      if(input$checkbox1){
        boxtab[[input$Fact1]] <- factor(boxtab[[input$Fact1]], levels = gtools::mixedsort(levels(as.factor(boxtab[[input$Fact1]]))))
      }
    }

    if( !any(names(boxtab)=="sample.id") ) {
      boxtab <- dplyr::rename(boxtab, sample.id = rowname)
    }

    boxtab$Depth <- sample_sums(local_physeq())
    setProgress(value = 10, detail = 'done')
    flog.info('boxtab done.')
    return(boxtab)
    })
  }
)


  get_box_plot <- reactive({
    req(boxtab(), get_meta_col())
    flog.info('renderPloty...')
    withProgress(message = 'Rendering plot...', min=0, max=10, value = 0,{
      dt <- boxtab()
      if(is.numeric(dt[,get_meta_col()])){
        p <- plot_ly(dt, x = as.formula(glue("~{get_meta_col()}")), y = as.formula(glue("~{input$metrics}")),
                     color = as.formula(glue("~{get_meta_col()}")), type = 'scatter')
      } else{
        p <- plot_ly(dt, x = as.formula(glue("~{get_meta_col()}")), y = as.formula(glue("~{input$metrics}")),
                     color = as.formula(glue("~{get_meta_col()}")), type = 'box', colors = r$factor_colors()[[input$Fact1]])
      }
      p %>% layout(title=input$metrics, yaxis = list(title = glue('{input$metrics}')), barmode = 'stack') %>%
        config(toImageButtonOptions = list(format = "svg"))
    })
  })


  output$boxplot <- renderPlotly({
    get_box_plot()
  })

  output$alphalrm <- renderPrint(
    print(alpha_lm())
  )

  alpha_lm <- reactive({
    dt <- boxtab()
    form1 = glue::glue("{input$metrics} ~ Depth + {get_meta_col()}")
    fit = lm(as.formula(form1), data=dt)
    return(summary(fit))
  })

  output$anovaBox <- renderUI({
    if(isNumFactor()){
      tagList(
        h3("Linear regression model"),
        verbatimTextOutput(ns("alphalrm"))
      )
    } else {
      tagList(
        h3("ANOVA Results"),
        DT::dataTableOutput(ns("testalpha")),
        tags$hr(),
        h3("Post-hoc Tukey HSD Test"),
        downloadButton(outputId = ns("boxtab_download"), label = "Download Tukey Table"),
        DT::dataTableOutput(ns("boxstats"))
      )
    }
  })


  reacalpha <- reactive({
    req(input$metrics, get_meta_col(), boxtab())

    flog.info('Alpha tests...')
    withProgress(message = 'Statistics...', min=0, max=10, value = 0,{

    form1 = glue::glue("{input$metrics} ~ Depth + {get_meta_col()}")
    anova_res1 <- aov( as.formula(form1), boxtab())

    tukey_hsd <- TukeyHSD(anova_res1, get_meta_col())

    LL = list()
    LL$form1 = form1
    LL$aov1 = summary(anova_res1)
    LL$groups1 <- tukey_hsd[[get_meta_col()]]

    LL$groups1 <- LL$groups1 %>% as.data.frame() %>% rownames_to_column('comparison')

    setProgress(value = 10, detail = 'done')

    })
    flog.info('Done...')
    return(LL)
 })


 output$testalpha <- DT::renderDataTable({
   req(input$metrics)
   tt <- reacalpha()
   anova_df <- as.data.frame(tt$aov1[[1]])
   anova_df <- tibble::rownames_to_column(anova_df, "Term")
   datatable(round_df(anova_df, 4), options = list(
     pageLength = 10,
     scrollX = TRUE,
     dom = 't'
   ))
 }, server = FALSE)


 output$boxstats <- DT::renderDataTable({
   req(reacalpha)
   LL = reacalpha()
   round_df(as.data.frame(LL$groups1), 4)
 }, filter="top", options = list(
   pageLength = 5,
   scrollX = TRUE
 ))


 output$boxtab_download <- downloadHandler(
   filename = "alpha_boxplot_stats.csv",
   content = function(file) {
     req(reacalpha)
     LL = reacalpha()
     write.table(LL$groups1, file, sep="\t", col.names=NA)}
 )
  })
}

## To be copied in the UI
# mod_alpha_ui("alpha_ui_1")

## To be copied in the server
# mod_alpha_server("alpha_ui_1", r = r)