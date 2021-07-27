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
#' @import shinyalert
#' @import bslib
mod_alpha_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      theme = bslib::bs_theme(version=4),
      useShinyalert(),
      infoBox("",
              "Use phyloseq object without taxa merging step.",
              icon = icon("info-circle"), fill=TRUE, width = 10),

      box(
        selectInput(
          ns("Fact1"),
          label = "Select factor to test: ",
          choices = ""
        ),
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
        plotly::plotlyOutput(ns("plot2")),
        width=12, status = "primary", solidHeader = TRUE, title = "Boxplot"
      ),
      box(
        h3("ANOVA results"),
        box(verbatimTextOutput(ns("testalpha")), width=12, status = "primary"),

        h3("TukeyHSD test results"),
        downloadButton(outputId = ns("boxtab_download"), label = "Download Table", icon = icon("download")),
        DT::dataTableOutput(ns("boxstats")),

        width=12, status = "primary", solidHeader = TRUE, title = "Statistics and tests", collapsible = TRUE
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

mod_alpha_server <- function(input, output, session, r = r){
  ns <- session$ns

  observeEvent(r$tabs$tabselected, {
    if(r$tabs$tabselected=='tab_alpha' && !isTruthy(r$phyloseq_filtered())){
      shinyalert(title = "Oops", text="Phyloseq object not present. Return to input data and validate all steps.", type='error')
    }
  })


  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "Fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })




  alpha1 <- eventReactive(input$launch_alpha,{
    flog.info('computing alpha1...')
    req(r$phyloseq_filtered())

    data <- r$phyloseq_filtered()

    alphatab <- estimate_richness(data, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                                                     "InvSimpson") )
    row.names(alphatab) = sample_names(data)

    LL=list()
    LL$alphatab = as.data.frame(alphatab)
    LL$data = data
    flog.info('computing alpha1 done.')
    return(LL)
  })


  output$alphaout <- DT::renderDataTable({
    LL = alpha1()
    LL$alphatab
  }, filter="top",options = list(pageLength = 5, scrollX = TRUE))

  alphagrp_table <- reactive({
    alpha.table <- alpha1()$alphatab
    metadata = tibble::rownames_to_column(r$sdat())
    metadata <- select(metadata, rowname, input$Fact1)
    alpha.table =  tibble::rownames_to_column(alpha.table)
    alpha.table <- dplyr::left_join(metadata, alpha.table, by = "rowname")

    alpha.table[,'rowname'] <- NULL

    alpha.table <- alpha.table %>%
      group_by_at(input$Fact1) %>%
      summarise(
        tibble(
          across(where(is.numeric), ~round(mean(.x),2), .names = "mean_{.col}"),
          across(where(is.numeric), ~round(median(.x),2), .names = "median_{.col}")
        )
      )
    return(alpha.table)
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


  boxtab <- eventReactive(input$launch_alpha,{
    req(r$sdat(), input$checkbox1, input$Fact1, r$phyloseq_filtered())
    flog.info('boxtab function')
    LL = alpha1()

    metadata = tibble::rownames_to_column(r$sdat())
    alphatab =  tibble::rownames_to_column(LL$alphatab)


    boxtab <- dplyr::left_join(metadata, alphatab, by = "rowname")

    if(input$checkbox1){
      print("ORDER factor")
      fun = glue::glue( "boxtab${input$Fact1} <- forcats::fct_relevel(boxtab[[input$Fact1]])")
      eval(parse(text=fun))
    }

    if( !any(names(boxtab)=="sample.id") ) {
      print("change rowname to sample.id")
      dplyr::rename(boxtab, sample.id = rowname)
    }

    boxtab$Depth <- sample_sums(r$phyloseq_filtered())
    
    boxtab
  }
)


  output$plot2 <- renderPlotly({
   plot_ly(boxtab(), x = as.formula(glue("~{input$Fact1}")), y = as.formula(glue("~{input$metrics}")),
           color = as.formula(glue("~{input$Fact1}")), type = 'box') %>% #, name = ~variable, color = ~variable) %>% #, color = ~variable
     layout(title=input$metrics, yaxis = list(title = glue('{input$metrics}')), xaxis = list(title = 'Samples'), barmode = 'stack') %>%
    config(toImageButtonOptions = list(format = "svg"))
 })


  reacalpha <- eventReactive(input$launch_alpha,{
    req(input$metrics, input$Fact1)
    flog.info('reacalpha')
    
    anova_data = boxtab()

    form1 = glue::glue("{input$metrics} ~ Depth + {input$Fact1}")
    anova_res1 <- aov( as.formula(form1), anova_data)

    fun <- glue::glue("tukey_hsd <- TukeyHSD(anova_res1, \"{input$Fact1}\")")
    eval(parse(text=fun))

    LL = list()
    LL$form1 = form1
    LL$aov1 = summary(anova_res1)

    fun <- glue::glue("LL$groups1 <- tukey_hsd${input$Fact1}")
    eval(parse(text=fun))

    return(LL)
 })

 output$testalpha <- renderPrint({
   req(reacalpha)
   LL = reacalpha()
   cat("ANOVA\n##########\n")
   print(LL$form1)
   print(LL$aov1)
 })

 output$boxstats <- DT::renderDataTable({
   req(reacalpha)
   LL = reacalpha()
   LL$groups1
 }, filter="top",options = list(pageLength = 5, scrollX = TRUE))

 output$boxtab_download <- downloadHandler(
   filename = "alpha_boxplot_stats.csv",
   content = function(file) {
     req(reacalpha)
     LL = reacalpha()
     write.table(LL$groups1, file, sep="\t", col.names=NA)}
 )
}

## To be copied in the UI
# mod_alpha_ui("alpha_ui_1")

## To be copied in the server
# callModule(mod_alpha_server, "alpha_ui_1")
