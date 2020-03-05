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
#' @import plotly
#' @importFrom shiny NS tagList 
mod_taxaboxplot_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    fluidPage(
      h1("Boxplots"),
      h3("You can select specific sample in Metadatas/Subset tab, and agglomerate to specific rank in ASVtable tab."),
      selectInput(
        ns("Fact1"),
        label = "Select factor to test: ",
        choices = ""
      ),
      h4("Clic on feature to generate boxplot:"),
      fluidRow(box(dataTableOutput(ns("pvalout1"))
                   , width = 8)),
      box(plotlyOutput(ns("boxplot1")), width=10, height=500),
      # verbatimTextOutput(ns("sids2")),
      box(verbatimTextOutput(ns("wilcoxprint")), width=10)
    )
    
  )
}
    
# Module Server
    
#' @rdname mod_taxaboxplot
#' @export
#' @keywords internal
#' @import plotly
#' @importFrom DT datatable
    
mod_taxaboxplot_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  observe({
    updateSelectInput(session, "Fact1",
                      choices = r$data16S()@sam_data@names)
    
  })
  
  LjoinGlom <- reactive({
    withProgress({
      print("melting table")
      # data.melt <- psmelt(Fdata)
      Fdata <- r$dat() #subglom()
      print("BP sdata")
      stable <- Fdata %>%
        sample_data() %>%
        as.matrix() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        tibble::rownames_to_column()
      
      print("BP otable")
      otable <- Fdata %>%
        otu_table() %>%
        # as.matrix() %>%
        t() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        tibble::rownames_to_column()
      
      print("BP otable ok")
      print(r$RankGlom())
      if(r$RankGlom() != "ASV"){
        # lvls <- paste(substr(tax_table(Fdata)[,r$RankGlom()],1,20), "_", names(otable)[-1],sep="")
        # names(otable)[-1] <- lvls
        lvls <- names(otable)[-1]
        
      }else(lvls <- names(otable))
      
      
      joinGlom <-
        dplyr::left_join(stable, otable, by = "rowname") %>%
        dplyr::rename(sample.id = rowname)
      
      LL=list()
      LL$joinGlom <- joinGlom
      LL$lvls <- lvls
      LL
    }, message = "Construct table...")
  })
  
  
  listBP <- reactive({
    withProgress({
      LL = LjoinGlom()
      joinGlom <- LL$joinGlom
      lvls <- LL$lvls
      
      print(length(lvls))
      stock=NULL
      print("loop")
      stock=NULL; pval1=NULL; taxa1=NULL
      for(i in lvls[-1]){
        # print(i)
        if(mean(joinGlom[,i]) == 0){next}
        res = kruskal.test(joinGlom[,i], joinGlom[,input$Fact1])
        pval1 = c(pval1, res$p.value)
        taxa1 = c(taxa1, i)
        if(res$p.value < 0.05){stock = c(stock, i)}
      }
      print("cbind")
      print(length(taxa1))
      print(length(pval1))
      respval <- cbind.data.frame(Taxa = taxa1, kruskal.pvalue = pval1)
      
      print(head(as.data.frame(respval)))
      
      LL = list()
      LL$joinGlom = joinGlom
      LL$pval = respval
      LL
    }, message="Kruskall test...")
  })
  
  output$pvalout1 <- renderDataTable({
    LL = listBP()
    # print(head(as.data.frame(LL$pval)))
    # print(str(as.data.frame(LL$pval)))
    datatable(as.data.frame(LL$pval), selection = "single", filter="top")
  })
  
  output$sids2 <- reactive({
    LL = listBP()
    stab <- LL$pval
    select1  <- stab[input$pvalout1_row_last_clicked,1]
    return(select1)
  })
  
  output$boxplot1 <- renderPlotly({
    if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
    print("Boxplot")
    LL = listBP()
    stab <- LL$pval
    joinGlom <- LL$joinGlom
    select1  <- stab[input$pvalout1_row_last_clicked,1]
    print(head(joinGlom))
    print(str(joinGlom))
    plot_ly(joinGlom, x = as.formula(glue("~{input$Fact1}")), y = as.formula(glue("~{select1}")),
            color = as.formula(glue("~{input$Fact1}")), type = 'box') %>% #, name = ~variable, color = ~variable) %>% #, color = ~variable
      layout(title=select1, yaxis = list(title = glue('{input$NORM} abundance')), xaxis = list(title = 'Samples'), barmode = 'stack')
  })

statsBP1 <- reactive({  
  if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
  LL = listBP()
  stab <- LL$pval
  joinGlom <- LL$joinGlom
  select1  <- stab[input$pvalout1_row_last_clicked,1]
  fun = glue( "res = pairwise.wilcox.test(joinGlom[,'{select1}'], joinGlom[,'{input$Fact1}'], p.adjust.method = 'none')" )
  eval(parse(text=fun))
  LL$res = res
  LL$select1 = select1
  LL
})

output$wilcoxprint <- renderPrint({
  LL = statsBP1()
  print(LL$select1)
  print(LL$res)
  })
  
  
  output$statsBP1 <- reactive({
    if(is.null(input$pvalout1_row_last_clicked)){return(NULL)}
    joinGlom <- LL$joinGlom
    select1  <- stab[input$pvalout1_row_last_clicked,1]
    tab1  <- joinGlom[,c(input$Fact1, select1)]
    
    tt = tab1 %>%
      group_by(SampleType) %>%
      group_map(~ summary(.x))
  })
  
  
}
    
## To be copied in the UI
# mod_taxaboxplot_ui("taxaboxplot_ui_1")
    
## To be copied in the server
# callModule(mod_taxaboxplot_server, "taxaboxplot_ui_1")
 
