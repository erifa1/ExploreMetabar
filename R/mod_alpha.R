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
              "Use phyloseq object without taxa merging step.",
              icon = icon("info-circle"), fill=TRUE, width = 10),

      box(
        selectInput(
          ns("Fact1"),
          label = "Select factor to test: ",
          choices = ""
        ),
        selectInput(
          ns("Fact2"),
          label = "Select second factor to combine: ",
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

        h3("TukeyHSD adjusted p-value heatmap"),
        plotOutput(ns("pvalue_matrix")),

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
    print(r$tabs$tabselected)
    if(r$tabs$tabselected!='data_loading' && !isTruthy(r$phyloseq_filtered())){
      shinyalert::shinyalert(title = "Oops", text="Phyloseq object not present. Return to input data and validate all steps.", type='error')
    }
  })


  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "Fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })

  observe({
    req(r$phyloseq_filtered(), input$Fact1)
    updateSelectInput(session, "Fact2",
                      choices = c('none',dplyr::setdiff(r$phyloseq_filtered()@sam_data@names, input$Fact1)))
  })

  alpha1 <- eventReactive(input$launch_alpha,{
    withProgress(message = 'Computing alpha diversity tables', min=0, max=10, value = 0,{
      flog.info('computing alpha1...')
      req(r$phyloseq_filtered())

      data <- r$phyloseq_filtered()
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

  alphagrp_table <- reactive({
    withProgress(message = 'Group table', min=0, max=10, value = 0,{
    alpha.table <- alpha1()$alphatab
    metadata <- as(r$sdat(), "data.frame")
    if(input$Fact2 != 'none' && ! is.numeric(metadata[, input$Fact1])){
      meta.col <- paste0(input$Fact1, '_', input$Fact2)
      metadata <- tidyr::unite(metadata, !!meta.col, input$Fact1, input$Fact2, na.rm=TRUE)
      metadata[, meta.col] <- as.factor(metadata[, meta.col])
      metadata <- select(metadata, rowname, meta.col)
    }
    else{
      meta.col <- input$Fact1
      metadata <- select(metadata, rowname, input$Fact1)
    }

    alpha.table =  tibble::rownames_to_column(alpha.table)
    alpha.table <- dplyr::left_join(as(r$sdat(), "data.frame"), alpha.table, by = c( "sample.id" = "rowname"))

    alpha.table[,'rowname'] <- NULL

    if(! is.numeric(alpha.table[, input$Fact1])){
      alpha.table <- alpha.table %>%
        group_by_at(input$Fact1) %>%
        summarise(
          tibble(
            across(where(is.numeric), ~round(mean(.x),2), .names = "mean_{.col}"),
            across(where(is.numeric), ~round(median(.x),2), .names = "median_{.col}")
          )
        )
    }
    return(alpha.table)
    setProgress(value = 10, detail = 'done')
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


  boxtab <- eventReactive(input$launch_alpha,{
    req(r$sdat(), input$Fact1, r$phyloseq_filtered())
    withProgress(message = 'Boxplot table', min=0, max=10, value = 0,{
    flog.info('boxtab function')
    LL = alpha1()

    alphatab =  tibble::rownames_to_column(LL$alphatab)
    metadata <- as(r$sdat(), "data.frame")
    if(input$Fact2 != 'none' && ! is.numeric(metadata[, input$Fact1])){
      meta.col <- paste0(input$Fact1, '_', input$Fact2)
      metadata <- tidyr::unite(metadata, !!meta.col, input$Fact1, input$Fact2, na.rm=TRUE)
      metadata[, meta.col] <- as.factor(metadata[, meta.col])
    } else{
      meta.col <- input$Fact1
    }
    boxtab <- dplyr::left_join(metadata, alphatab, by = c('sample.id' = "rowname"))

    if(! is.numeric(boxtab[, input$Fact1])){
      if(input$checkbox1){
        print("ORDER factor")
        fun = glue::glue( "boxtab${input$Fact1} = factor( boxtab${input$Fact1}, levels = gtools::mixedsort(levels(as.factor(boxtab${input$Fact1}))) ) ")

        eval(parse(text=fun))
      }
    }


    if( !any(names(boxtab)=="sample.id") ) {
      print("change rowname to sample.id")
      dplyr::rename(boxtab, sample.id = rowname)
    }

    boxtab$Depth <- sample_sums(r$phyloseq_filtered())
    setProgress(value = 10, detail = 'done')
    return(boxtab)
    })
  }
)


  output$plot2 <- renderPlotly({
    withProgress(message = 'Rendering plot...', min=0, max=10, value = 0,{
    dt <- boxtab()
    # browser()
    if(input$Fact2 != 'none'){
      meta.col <- paste0(input$Fact1, '_', input$Fact2)
    } else{
      meta.col <- input$Fact1
    }

    if(is.numeric(dt[,meta.col])){
      p <- plot_ly(dt, x = as.formula(glue("~{meta.col}")), y = as.formula(glue("~{input$metrics}")),
                   color = as.formula(glue("~{meta.col}")), type = 'scatter')
    } else{
      p <- plot_ly(dt, x = as.formula(glue("~{meta.col}")), y = as.formula(glue("~{input$metrics}")),
                   color = as.formula(glue("~{meta.col}")), type = 'box')
    }
     p %>% layout(title=input$metrics, yaxis = list(title = glue('{input$metrics}')), barmode = 'stack') %>%
    config(toImageButtonOptions = list(format = "svg"))

  })
 })


  reacalpha <- reactive({
    req(input$metrics, input$Fact1)
    metadata = tibble::rownames_to_column(r$sdat())
    if(input$Fact2 != 'none'){
      meta.col <- paste0(input$Fact1, '_', input$Fact2)
      metadata <- tidyr::unite(metadata, !!meta.col, input$Fact1, input$Fact2, na.rm=TRUE)
      metadata[, meta.col] <- as.factor(metadata[, meta.col])
    } else{
      meta.col <- input$Fact1
    }
    cat(file=stderr(),'Alpha tests...',"\n")
    withProgress(message = 'Statistics...', min=0, max=10, value = 0,{
    anova_data = boxtab()

    form1 = glue::glue("{input$metrics} ~ Depth + {meta.col}")
    anova_res1 <- aov( as.formula(form1), anova_data)

    fun <- glue::glue("tukey_hsd <- TukeyHSD(anova_res1, \"{meta.col}\")")
    eval(parse(text=fun))

    LL = list()
    LL$form1 = form1
    LL$aov1 = summary(anova_res1)
    fun <- glue::glue("LL$groups1 <- tukey_hsd${meta.col}")
    eval(parse(text=fun))

    setProgress(value = 10, detail = 'done')

    })
    cat(file=stderr(),'Done...',"\n")

    return(LL)
 })

 alpha_test <- reactive({
  t <- reacalpha()
  return(t)
 })

 output$testalpha <- renderPrint({
  req(input$metrics)
   tt <- alpha_test()
   print(tt$form1)
   print(tt$aov1)
 })


 output$pvalue_matrix <- renderPlot({
   req(reacalpha)
   LL = reacalpha()
   mat <- as.data.frame(LL$groups1) %>% rownames_to_column()
   mat <- separate(mat, rowname, c('var1', 'var2'), sep='-')
   mat <- dplyr::select(mat, var1, var2, `p adj`) %>% dplyr::filter(!is.na('p adj')) # %>% tidyr::pivot_wider(names_from = var1, values_from = 'p adj')
   mat <- mat %>% mutate_if(is.numeric, round, digits=3)

   # mat <- mat %>% mutate_if(is.numeric, round, digits=3) %>% column_to_rownames('var2') %>% as.matrix()
   # DT::datatable(mat, options = list(pageLength = 20, scrollX = TRUE)) %>% DT::formatStyle(colnames(mat), backgroundColor = styleInterval(c(0,0.01,0.05,1), c("white","greenyellow", "lightgreen","orange","red")))
   ggplot(mat, aes(x=var1, y=var2)) +
     geom_tile(aes(fill=`p adj`)) +
     geom_text(aes(label = `p adj`)) +
     scale_fill_gradientn(colors = c('green', 'red'), breaks = c(0, 0.05), limits = c(0, 0.05)) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
}

## To be copied in the UI
# mod_alpha_ui("alpha_ui_1")

## To be copied in the server
# callModule(mod_alpha_server, "alpha_ui_1")
