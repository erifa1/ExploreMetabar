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
#' @noRd
#' @importFrom plotly plotlyOutput
#' @importFrom shiny NS tagList
#' @importFrom bslib layout_sidebar sidebar accordion accordion_panel navset_card_underline nav_panel
#' @importFrom bsicons bs_icon
mod_taxaboxplot_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "Settings",
      open = "desktop",
      width = "350px",
      htmltools::p(
        "Kruskal\u2013Wallis test on each taxon for all factors. P-values adjusted using FDR. For numerical factors, samples with zero abundance are omitted.",
        style = "font-size: 0.9em; color: grey;"
      ),
      tags$hr(),
      accordion(
        id = ns("config_accordion"),
        open = "Test settings",
        multiple = TRUE,
        accordion_panel(
          "Test settings",
          icon = bs_icon("sliders"),
          uiOutput(ns('ui_picker')),
          uiOutput(ns('ui_radio_tests'))
        )
      ),
      tags$hr(),
      actionButton(
        ns("go1"),
        "Run Test / Correlation",
        icon = bs_icon("play-fill"),
        class = "btn-primary w-100 btn-lg"
      )
    ),
    navset_card_underline(
      title = "Results",
      full_screen = TRUE,
      nav_panel(
        "Features",
        icon = bs_icon("table"),
        htmltools::p("Click on a feature below to generate the plot:"),
        DT::dataTableOutput(ns("pvalout1"))
      ),
      nav_panel(
        "Plot",
        icon = bs_icon("bar-chart-line"),
        checkboxInput(ns("order1"), label = "Automatic order factor", value = TRUE),
        plotlyOutput(ns("boxplot1"))
      ),
      nav_panel(
        "Pairwise Wilcox",
        icon = bs_icon("calculator"),
        uiOutput(ns('ui_pair_test'))
      )
    )
  )
}

# Module Server

#' @rdname mod_taxaboxplot
#' @noRd
#' @keywords internal
#' @importFrom plotly plot_ly renderPlotly
#' @importFrom DT datatable
#' @importFrom DT formatStyle
#' @importFrom DT formatRound
#' @importFrom DT styleInterval
#' @importFrom rlang .data

mod_taxaboxplot_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns

  
  output$ui_picker <- renderUI({
    req(r$var_list())
    shinyWidgets::pickerInput(ns("boxplot_fact1"),
                              label = "Select factor to test: ",
                              choices = r$var_list(),
                              selected = r$var_list()[1],
                              multiple = FALSE,
                              options = pickerOptions(
                                actionsBox = TRUE,
                                liveSearch = TRUE,
                                showContent = FALSE
                              ),
                              choicesOpt = list(
                                content = unlist(lapply(
                                  X = r$var_list(),
                                  FUN = function(x) {
                                    htmltools::doRenderTags(
                                      tags$div(
                                        splitLayout(cellWidths = 200,
                                                    tags$div(
                                                      style = htmltools::css(fontWeight = "bold"),
                                                      x
                                                    ),
                                                    tags$div(
                                                      style = htmltools::css(color = 'grey'),
                                                      class(r$sdat()[,x])
                                                    ),
                                                    tags$div(
                                                      style = htmltools::css(color = 'grey'),
                                                      paste0(sum(is.na(r$sdat()[,x])), '/', nrow(r$sdat()), ' NAs')
                                                    )
                                        )
                                      )
                                    )
                                  }
                                ))
                              )
    )
  })

  
  
  factor_input   <- reactive({ input$boxplot_fact1 })
  get_meta_col   <- make_get_meta_col(factor_input, r)
  local_metadata <- make_local_metadata(factor_input, get_meta_col, r)
  isNumFactor    <- make_is_num_factor(get_meta_col, local_metadata)

  local_physeq <- reactive({
    phy <- r$phyloseq_filtered_norm()
    phyloseq::sample_data(phy) <- phyloseq::sample_data(local_metadata())
    return(phy)
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
      DT::dataTableOutput(ns("wilcoxDT"))
    }
  })


  get_pval_table <- eventReactive(input$go1, {
    if(isNumFactor()){
      res <- get_corr_pval_table()
    } else{
      res <- get_kruskal_pval_table()
    }
    t_table <- as.data.frame(tax_table(r$phyloseq_filtered_norm())) %>% rownames_to_column()
    res <- left_join(res, t_table, by = c('taxa' = 'rowname')) %>% 
      mutate(p.adj = p.adjust(p.value, method = "fdr"), .after = p.value) %>% 
      arrange(p.value)
    return(res)
  })
  
  
  get_merged_table <- eventReactive(input$go1, {
    withProgress({
      otable <- otu_table(local_physeq()) %>% t() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        rownames_to_column('sample.id')
      incProgress(amount = 0.5)
      metadata <- local_metadata()
      incProgress(amount = 0.1)
      if(! any(colnames(metadata) %in% "sample.id")){
        metadata <- metadata[, get_meta_col(), drop=FALSE] %>% rownames_to_column('sample.id')
      }
      incProgress(amount = 0.1)
      mtable <- left_join(otable, metadata, by='sample.id')
      incProgress(amount = 0.3)
    }, message = 'Merging table...', min = 0, max = 1)
    return(mtable)
  })
  
  
  get_kruskal_pval_table <- reactive({
    mtable <- get_merged_table()
    withProgress({
      results <- tibble('taxa' = as.character(),
                        'p.value' = as.numeric())
      for(taxa in taxa_names(local_physeq())){
        res = kruskal.test(mtable[,taxa], mtable[,get_meta_col()])
        results <- results %>% add_row('taxa' = taxa, 'p.value' = res$p.value)
        incProgress(amount = 1/length(taxa_names(local_physeq())))
      }
    }, message = 'Computing kruskal.test...', min = 0, max = 1)
    return(results)
  })
  
  
  get_corr_pval_table <- reactive({
    mtable <- get_merged_table()
    results <- tibble('taxa' = as.character(),
                      'available.data' = as.numeric(),
                      'p.value' = as.numeric(),
                      'cor.coef' = as.numeric())
    withProgress({
      for(taxa in taxa_names(local_physeq())){
        tmp <- mtable[,c(taxa, get_meta_col())]
        tmp <- na.omit(tmp)
        tmp <- tmp[tmp[,taxa] > 0,]
        if(nrow(tmp) > 5){
          res <- stats::cor.test(tmp[, taxa], tmp[, get_meta_col()], method = input$cor_test)
          results <- results %>% add_row('taxa' = taxa, 'available.data' = nrow(tmp), 'p.value' = res$p.value, cor.coef = res$estimate)
        } else {
          results <- results %>% add_row('taxa' = taxa, 'available.data' = nrow(tmp), 'p.value' = NA, cor.coef = NA)
        }
        incProgress(amount = 1/length(taxa_names(local_physeq())))
      }
    }, message = 'Computing correlations...', min = 0, max = 1)
    
    return(results)
  })


  output$pvalout1 <- DT::renderDataTable({
    req(get_pval_table())
    get_pval_table() %>% datatable(selection = "single", filter="top", options = list(scrollX = TRUE)) %>%
      formatStyle(
        'p.value',
        backgroundColor=styleInterval(c(0,0.01,0.05,0.1,1), 
        c("white", "greenyellow", "lightgreen","yellow","orange","red") )) %>%
      formatStyle(
        'p.adj',
        backgroundColor=styleInterval(c(0,0.01,0.05,0.1,1), 
        c("white", "greenyellow", "lightgreen","yellow","orange","red") ))
  })

  ordertable1 <- reactive({
    mtable <- get_merged_table()
    if(input$order1){
      mtable[[get_meta_col()]] <- factor(mtable[[get_meta_col()]], levels = gtools::mixedsort(unique(mtable[[get_meta_col()]])))
    }
    return(mtable)
  })


  output$boxplot1 <- renderPlotly({
    if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
    
    
    if(isNumFactor()){
      mtable <- get_merged_table()
      ptype <- 'scatter'
      select1  <- get_pval_table()[input$pvalout1_row_last_clicked,'taxa'] %>% pull
      p <- ggplotly(ggplot2::ggplot(data = mtable, aes(x = .data[[select1]], y = .data[[get_meta_col()]])) + 
                 geom_point() + 
                 geom_smooth(method = 'lm', se = T, na.rm = T, show.legend = T))
                   
    } else{
      mtable <- ordertable1()
      select1  <- get_pval_table()[input$pvalout1_row_last_clicked,'taxa'] %>% pull
      p <- plot_ly(mtable, x = as.formula(glue("~ `{get_meta_col()}`")), y = as.formula(paste0("~ `", select1, "`")),
                   color = as.formula(glue("~ `{get_meta_col()}`")), type = 'box', colors=r$factor_colors()[[input$boxplot_fact1]])
    }
    
    return(p)
  })

  
  get_pairwise_test <- reactive({
    if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
    mtable <- get_merged_table()
    select1  <- get_kruskal_pval_table()[input$pvalout1_row_last_clicked,'taxa'] %>% pull
    res = pairwise.wilcox.test(mtable[,select1], mtable[,get_meta_col()], p.adjust.method = 'none')
    return(res)
  })


  output$wilcoxDT <- DT::renderDataTable({
    req(get_pairwise_test())
    LL = get_pairwise_test()
    wtab = as.data.frame(LL$p.value)

    wtab %>%
      tibble::rownames_to_column("Condition1") %>%
      tidyr::pivot_longer(-Condition1, names_to = "Condition2", values_to = "pvalue") %>%
      na.omit() %>%
      datatable() %>%
      formatStyle("pvalue",
        backgroundColor = styleInterval(c(0,0.05), c("white","greenyellow", "white"))
    )
  })

  })
}

## To be copied in the UI
# mod_taxaboxplot_ui("taxaboxplot_ui_1")

## To be copied in the server
# mod_taxaboxplot_server("taxaboxplot_ui_1", r = r)
