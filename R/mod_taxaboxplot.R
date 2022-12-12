# Module UI

#' @title   mod_taxaboxplot_ui and mod_taxaboxplot_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_taxaboxplot
#'
#' @keywords internal
#' @export
#' @importFrom plotly plotlyOutput
#' @importFrom shiny NS tagList
mod_taxaboxplot_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidPage(

      infoBox("Reminder :",
              "You can select specific sample in Metadatas/Subset module, and agglomerate to specific rank in ASVtable module",
              icon = icon("info-circle"), fill=TRUE, width = 10),

      box(
        selectInput(
          ns("boxplot_fact1"),
          label = "Select factor to test: ",
          choices = ""
        ),
        uiOutput(ns('ui_radio_tests')),
        actionButton(ns("go1"), "Run Test/Correlation", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),

      box(
      h2(icon("diagnoses"),"Click on feature below to generate plot:"),
      DT::dataTableOutput(ns("pvalout1")),
      title = "Features:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(
          checkboxInput(ns("order1"), label = "Automatic order factor", value = TRUE),
          plotlyOutput(ns("boxplot1")), #, height=500
          title = "Boxplot:", width = 12, status = "primary", solidHeader = TRUE
          ),
      uiOutput(ns('ui_pair_test'))
    )

  )
}

# Module Server

#' @rdname mod_taxaboxplot
#' @export
#' @keywords internal
#' @importFrom plotly plot_ly renderPlotly
#' @importFrom DT datatable
#' @importFrom DT formatStyle
#' @importFrom DT formatRound
#' @importFrom DT styleInterval
#' @import formulaic

mod_taxaboxplot_server <- function(input, output, session, r = r){
  ns <- session$ns

  observe({
    req(r$phyloseq_filtered(), r$var_list())
    updateSelectInput(session, "boxplot_fact1",
                      choices = r$var_list())

  })
  
  
  isNumFactor <- reactive({
    req(input$boxplot_fact1, r$sdat())
    metadata <- r$sdat()
    if(is.numeric(metadata[, input$boxplot_fact1])){
      return(TRUE)
    } else{
      return(FALSE)
    }
  })
  
  
  output$ui_radio_tests <- renderUI({
    if(isNumFactor()){
      radioButtons(ns('cor_test'),
                   'Select correlation test:',
                   choices = c('pearson', 'spearman', 'kendall'),
                   selected = 'pearson')
    }
  })
  
  
  output$ui_pair_test <- renderUI({
    if(! isNumFactor()){
      box(DT::dataTableOutput(ns("wilcoxDT")),
          title = "Results of pairwise wilcox test:", width = 12, status = "primary", solidHeader = TRUE)
      # box(verbatimTextOutput(ns("wilcoxprint")),
      #     title = "Raw Results of pairwise wilcox test:", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE)
    }
  })


  get_pval_table <- eventReactive(input$go1, {
    if(isNumFactor()){
      get_corr_pval_table()
    } else{
      get_kruskal_pval_table()
    }
  })
  
  
  get_merged_table <- reactive({
    otable <- otu_table(r$phyloseq_filtered_norm()) %>% t() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      rownames_to_column('sample.id')
    metadata <- r$sdat()
    metadata <- metadata[, input$boxplot_fact1, drop=FALSE] %>% rownames_to_column('sample.id')
    mtable <- left_join(otable, metadata, by='sample.id')
    return(mtable)
  })
  
  
  get_kruskal_pval_table <- reactive({
    mtable <- get_merged_table()
    results <- tibble('taxa' = as.character(),
                      'p.value' = as.numeric())
    for(taxa in taxa_names(r$phyloseq_filtered_norm())){
      res = kruskal.test(mtable[,taxa], mtable[,input$boxplot_fact1])
      results <- results %>% add_row('taxa' = taxa, 'p.value' = res$p.value)
    }
    return(results)
  })
  
  
  get_corr_pval_table <- reactive({
    mtable <- get_merged_table()
    results <- tibble('taxa' = as.character(),
                      'p.value' = as.numeric(),
                      'cor.coef' = as.numeric())
    for(taxa in taxa_names(r$phyloseq_filtered_norm())){
      res <- stats::cor.test(mtable[,taxa], mtable[,input$boxplot_fact1], method = input$cor_test)
      results <- results %>% add_row('taxa' = taxa, 'p.value' = res$p.value, cor.coef = res$estimate)
    }
    return(results)
  })


  output$pvalout1 <- DT::renderDataTable({
    req(get_pval_table())
    get_pval_table() %>% datatable(selection = "single", filter="top") %>%
      formatStyle(
        'p.value',
        backgroundColor=styleInterval(c(0,0.01,0.05,1), c("white","greenyellow", "lightgreen","yellow","red")))
  })

  ordertable1 <- reactive({
    mtable <- get_merged_table()
    if(input$order1){
      fun = glue::glue( "mtable${input$boxplot_fact1} = factor( mtable${input$boxplot_fact1}, levels = gtools::mixedsort(unique(mtable${input$boxplot_fact1})) ) " )
      eval(parse(text=fun))
    }
    return(mtable)
  })


  output$boxplot1 <- renderPlotly({
    if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
    
    if(isNumFactor()){
      mtable <- get_merged_table()
      ptype <- 'scatter'
      select1  <- get_corr_pval_table()[input$pvalout1_row_last_clicked,'taxa'] %>% pull
      p <- ggplotly(ggplot2::ggplot(data = mtable, aes_string(x = formulaic::add.backtick(select1), y = input$boxplot_fact1)) + 
                 geom_point() + 
                 geom_smooth(method = 'lm', se = T, na.rm = T, show.legend = T))
                   
    } else{
      mtable <- ordertable1()
      select1  <- get_kruskal_pval_table()[input$pvalout1_row_last_clicked,'taxa'] %>% pull
      p <- plot_ly(mtable, x = as.formula(glue("~ {input$boxplot_fact1}")), y = as.formula(paste0("~", formulaic::add.backtick(select1))),
                   color = as.formula(glue("~{input$boxplot_fact1}")), type = 'box')
    }
    
    return(p)
  })

  
  get_pairwise_test <- reactive({
    if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
    mtable <- get_merged_table()
    select1  <- get_kruskal_pval_table()[input$pvalout1_row_last_clicked,'taxa'] %>% pull
    res = pairwise.wilcox.test(mtable[,select1], mtable[,input$boxplot_fact1], p.adjust.method = 'none')
    return(res)
  })


  output$wilcoxDT <- DT::renderDataTable({
    req(get_pairwise_test())
    LL = get_pairwise_test()
    wtab = as.data.frame(LL$p.value)

    wtab %>%
      tibble::rownames_to_column() %>%
      reshape2::melt(value.name = "pvalue") %>%
      na.omit() %>%
      rename(Condition1 = rowname)%>%
      rename(Condition2 = variable) %>%
      datatable() %>%
      formatStyle("pvalue",
        backgroundColor = styleInterval(c(0,0.05), c("white","greenyellow", "white"))
    )
  })

}

## To be copied in the UI
# mod_taxaboxplot_ui("taxaboxplot_ui_1")

## To be copied in the server
# callModule(mod_taxaboxplot_server, "taxaboxplot_ui_1")
