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
        actionButton(ns("go1"), "Run Test/Boxplot", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),

      box(
      h2(icon("diagnoses"),"Click on feature below to generate boxplot:"),
      DT::dataTableOutput(ns("pvalout1")),
      title = "Features:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(
          checkboxInput(ns("order1"), label = "Automatic order factor", value = TRUE),
          plotlyOutput(ns("boxplot1")), #, height=500
          title = "Boxplot:", width = 12, status = "primary", solidHeader = TRUE
          ),
      box(DT::dataTableOutput(ns("wilcoxDT")),
          title = "Results of pairwise wilcox test:", width = 12, status = "primary", solidHeader = TRUE),
      box(verbatimTextOutput(ns("wilcoxprint")),
          title = "Raw Results of pairwise wilcox test:", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE)
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

mod_taxaboxplot_server <- function(input, output, session, r = r){
  ns <- session$ns

  observe({
    req(r$phyloseq_filtered(), r$sdat())
    metadata <- as(r$sdat(), "data.frame")
    updateSelectInput(session, "boxplot_fact1",
                      choices = colnames(metadata))

  })
  
  isNumFactor <- reactive({
    req(input$boxplot_fact1, r$sdat())
    metadata <- as(r$sdat(), "data.frame")
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

  
  LjoinGlom <- reactive({
    req(r$phyloseq_filtered_norm(), r$rank_glom(), r$sdat())
    withProgress({
      # browser()
      Fdata <- r$phyloseq_filtered_norm() #r$dat()

      #If taxa names begin with a number
      if(any(grepl("^[0-9].*$", taxa_names(Fdata)))) {
        taxa_names(Fdata) <- paste("ASV_", taxa_names(Fdata), sep="")
      }

      # print("BP sdata")
      stable <- Fdata %>%
        sample_data() %>%
        as.matrix() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        tibble::rownames_to_column()

      # print("BP otable")
      otable <- Fdata %>%
        otu_table() %>%
        # as.matrix() %>%
        t() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        tibble::rownames_to_column()

      # print("BP otable ok")
      # print(r$rank_glom())
      if(r$rank_glom() != "ASV"){
        lvls <- names(otable)[-1]

      }
      else{
        lvls <- names(otable)
      }


      joinGlom <- dplyr::left_join(stable, otable, by = "rowname")
      if( !any(names(joinGlom)=="sample.id") ) { print("change rowname to sample.id"); dplyr::rename(joinGlom, sample.id = rowname) }

      LL <- list()
      LL$joinGlom <- joinGlom
      LL$lvls <- lvls
      LL
    }, message = "Construct table...")
  })


  listBP <- reactive({
    req(input$boxplot_fact1, LjoinGlom())
    withProgress({
      # browser()
      LL = LjoinGlom()
      joinGlom <- LL$joinGlom
      lvls <- LL$lvls

      # print(length(lvls))
      stock=NULL
      print("loop")
      stock=NULL; pval1=NULL; taxa1=NULL
      for(i in lvls[-1]){
        # print(i)
        if(mean(joinGlom[,i]) == 0){next}
        res = kruskal.test(joinGlom[,i], joinGlom[,input$boxplot_fact1])
        pval1 = c(pval1, res$p.value)
        taxa1 = c(taxa1, i)
        if(res$p.value < 0.05){stock = c(stock, i)}
      }
      # print("cbind")
      # print(length(taxa1))
      # print(length(pval1))
      respval <- cbind.data.frame(Taxa = taxa1, kruskal.pvalue = pval1)

      print(head(as.data.frame(respval)))

      LL = list()
      LL$joinGlom = joinGlom
      LL$pval = respval
      LL
    }, message="Kruskall test...")
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
    metadata <- as(r$sdat(), "data.frame")
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
      browser()
      fit <- lm(as.formula(paste0(input$boxplot_fact1, '~', formulaic::add.backtick(select1))), data = mtable)
      
                   
    } else{
      mtable <- ordertable1()
      select1  <- get_kruskal_pval_table()[input$pvalout1_row_last_clicked,'taxa'] %>% pull
      ptype <- 'box'
      # p <- plot_ly(mtable, x = as.formula(glue("~ {input$boxplot_fact1}")), y = as.formula(paste0("~", formulaic::add.backtick(select1))),
      #         color = as.formula(glue("~{input$boxplot_fact1}")), type = 'box') %>% 
      #   layout(title=select1, yaxis = list(title = glue('{input$NORM} abundance')), xaxis = list(title = 'Samples'), barmode = 'stack')
    }
    p <- plot_ly(mtable, x = as.formula(glue("~ {input$boxplot_fact1}")), y = as.formula(paste0("~", formulaic::add.backtick(select1))),
                 color = as.formula(glue("~{input$boxplot_fact1}")), type = ptype)
    return(p)
  })

  # statsBP1 <- reactive({
  #   if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
  #   LL = listBP()
  #   stab <- LL$pval
  #   joinGlom <- LL$joinGlom
  #   select1  <- stab[input$pvalout1_row_last_clicked,1]
  #   fun = glue( "res = pairwise.wilcox.test(joinGlom[,'{select1}'], joinGlom[,'{input$boxplot_fact1}'], p.adjust.method = 'none')" )
  #   eval(parse(text=fun))
  #   LL$res = res
  #   LL$select1 = select1
  #   LL
  # })

  # output$wilcoxprint <- renderPrint({
  #   if(! isNumFactor()){
  #     LL = statsBP1()
  #     print(LL$select1)
  #     print(LL$res)
  #   }
  # })

  
# output$wilcoxDT <- DT::renderDataTable({
#   if( isNumFactor()){
#     
#   } else{
#     LL = statsBP1()
#     wtab = as.data.frame(LL$res$p.value)
#   
#     wtab %>%
#       tibble::rownames_to_column() %>%
#       reshape2::melt(value.name = "pvalue") %>%
#       na.omit() %>%
#       rename(Condition1 = rowname)%>%
#       rename(Condition2 = variable) %>%
#       datatable() %>%
#       formatStyle("pvalue",
#         backgroundColor = styleInterval(c(0,0.05), c("white","greenyellow", "white"))
#     )
#   }
# })

# 
#   output$statsBP1 <- reactive({
#     req(input$pvalout1_row_last_clicked, input$boxplot_fact1)
#     # if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
#     joinGlom <- LL$joinGlom
#     select1  <- stab[input$pvalout1_row_last_clicked,1]
#     tab1  <- joinGlom[,c(input$boxplot_fact1, select1)]
# 
#     tt = tab1 %>%
#       group_by(SampleType) %>%
#       group_map(~ summary(.x))
#   })
# 

}

## To be copied in the UI
# mod_taxaboxplot_ui("taxaboxplot_ui_1")

## To be copied in the server
# callModule(mod_taxaboxplot_server, "taxaboxplot_ui_1")
