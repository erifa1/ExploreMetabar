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
mod_alpha_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      infoBox("",
              "Use the phyloseq object without performing the taxa merging step.",
              icon = icon("info-circle"), fill=TRUE, width = 10),

      box(
        uiOutput(ns('ui_alpha_factor')),
        checkboxInput(ns("checkbox1"), label = "Automatic order factor", value = TRUE),

        actionButton(ns("launch_alpha"), "Run Alpha Diversity", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),

      box(
        DT::dataTableOutput(ns("alphaout")),
        downloadButton(outputId = ns("alpha_download"), label = "Download Table", icon = icon("download"), class = "butt",
                       style="background-color: #3b9ef5"),
        width=12, status = "primary", solidHeader = TRUE, title = "Alpha indexes table", collapsible = TRUE, collapsed = FALSE
      ),
      box(
        DT::dataTableOutput(ns("alphagrp")),
        downloadButton(outputId = ns("alphagrp_download"), label = "Download group Table", icon = icon("download"), class = "butt",
                       style="background-color: #3b9ef5"),
        width=12, status = "primary", solidHeader = TRUE, title = "Alpha indexes by group", collapsible = TRUE, collapsed = FALSE
      ),
      box(
        radioButtons(ns("metrics"), "Choose one index:", inline = TRUE,
                     choices =
                       list("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                            "InvSimpson"),
                     selected = c("Shannon")
        ),
        plotly::plotlyOutput(ns("boxplot")),
        width=12, status = "primary", solidHeader = TRUE, title = "Boxplot"
      ),
      uiOutput(ns('anovaBox'))
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

mod_alpha_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns

  factor_input <- reactive({ input$Fact1 })
  get_meta_col  <- make_get_meta_col(factor_input, r)
  local_metadata <- make_local_metadata(factor_input, get_meta_col, r)
  isNumFactor   <- make_is_num_factor(get_meta_col, local_metadata)
  local_physeq  <- make_local_physeq(local_metadata, r)
  

  observeEvent(r$tabs$tabselected, {
    flog.info(paste0('tab - ', r$tabs$tabselected))
    if(r$tabs$tabselected!='data_loading' && !isTruthy(r$phyloseq_filtered())){
      shinyalert::shinyalert(title = "Oops", text="Phyloseq object not present. Return to input data and validate all steps.", type='error')
      req(FALSE)
    }
  })
  
  
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
    LL$alphatab
  }, filter="top",options = list(pageLength = 5, scrollX = TRUE))

  
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
            across(where(is.numeric), ~round(mean(.x),2), .names = "mean_{.col}"),
            across(where(is.numeric), ~round(median(.x),2), .names = "median_{.col}")
          )
        )
    }
    setProgress(value = 10, detail = 'done')
    return(alpha.table)
    })
  })

  output$alphagrp <- DT::renderDataTable({
    alphagrp_table()
  }, filter="top",options = list(pageLength = 5, scrollX = TRUE))

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
      box(
        h3("Linear regression model"),
        box(verbatimTextOutput(ns("alphalrm")), width=12, status = "primary"),
        width=12, status = "primary", solidHeader = TRUE, title = "Statistics and tests", collapsible = TRUE
      )
    }
    else{
      box(
            h3("ANOVA results"),
            box(verbatimTextOutput(ns("testalpha")), width=12, status = "primary"),

            h3("TukeyHSD test results"),
            downloadButton(outputId = ns("boxtab_download"), label = "Download Table", icon = icon("download")),
            DT::dataTableOutput(ns("boxstats")),

            width=12, status = "primary", solidHeader = TRUE, title = "Statistics and tests", collapsible = TRUE
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


 output$testalpha <- renderPrint({
  req(input$metrics)
   tt <- reacalpha()
   print(tt$form1)
   print(tt$aov1)
 })


 output$boxstats <- DT::renderDataTable({
   req(reacalpha)
   LL = reacalpha()
   LL$groups1
 }, filter="top", options = list(pageLength = 5, scrollX = TRUE))


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
