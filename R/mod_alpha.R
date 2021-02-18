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
#' @import plotly
mod_alpha_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      infoBox("", 
              "Use phyloseq object without taxa merging step.", 
              icon = icon("info-circle"), fill=TRUE, width = 10),
      
      box(
        radioButtons(ns("metrics"), "Choose one index:", inline = TRUE,
                     choices =
                       list("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                            "InvSimpson"),
                     selected = c("Shannon")
        ),
        numericInput(ns("minAb"), "Minimum raw abundance:", 1, min = 1, max = NA),
        selectInput(
          ns("Fact1"),
          label = "Select factor to test: ",
          choices = ""
        ),
        checkboxInput(ns("checkbox1"), label = "Automatic order factor", value = TRUE), 
        
        actionButton(ns("go1"), "Run Alpha Diversity", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),

      box(
        dataTableOutput(ns("alphaout")),
        downloadButton(outputId = ns("alpha_download"), label = "Download Table", icon = icon("download"), class = "butt",
                       style="background-color: #3b9ef5"),
        width=12, status = "primary", solidHeader = TRUE, title = "Alpha indexes table", collapsible = TRUE, collapsed = TRUE
      ),
      box(
        plotlyOutput(ns("plot2")),
        width=12, status = "primary", solidHeader = TRUE, title = "Boxplot"
      ),
      box(
        downloadButton(outputId = ns("boxtab_download"), label = "Download Table", icon = icon("download")), 
        dataTableOutput(ns("boxstats")),
        box(verbatimTextOutput(ns("testalpha")), width=12, status = "primary"),
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
#' @importFrom DT datatable
#' @import plotly
#' @importFrom agricolae HSD.test
#' @importFrom gtools mixedsort

mod_alpha_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  observe({
    updateSelectInput(session, "Fact1",
                      choices = r$data16S()@sam_data@names)
  })
  
  
  alpha1 <- eventReactive(input$go1, {
    if(is.null(r$subdata())){return(NULL)}
    print(input$metrics)
    print(input$minAb)
    
    data <- prune_taxa(taxa_sums(r$subdata()) > input$minAb, r$subdata())
    if(r$RankGlom() == "ASV"){
      data <- prune_taxa(r$asvselect(), data)
    }
    print(data)
    
    alphatab <- estimate_richness(data, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                                                     "InvSimpson") )
    row.names(alphatab) = sample_names(data)
      
    LL=list()
    LL$alphatab = as.data.frame(alphatab)
    LL$data = data
    
    LL
  })
  
  output$alphaout <- renderDataTable({
    if(is.null(alpha1())){return(NULL)}
    LL = alpha1()
    LL$alphatab
  }, filter="top",options = list(pageLength = 5, scrollX = TRUE)) ##filter = "top",

  output$alpha_download <- downloadHandler(
    filename = "alpha_index.csv",
    content = function(file) {
      LL = alpha1()
      write.table(LL$alphatab, file, sep="\t", col.names=NA)}
  )
  

boxtab <- reactive(
  {
    print("plotAlpha")
    LL = alpha1()
    sdat = tibble::rownames_to_column(as.data.frame(as.matrix(r$subdata()@sam_data)))
    alphatab =  tibble::rownames_to_column(LL$alphatab)
    
    print(head(sdat))
    print(head(alphatab))
    boxtab <- dplyr::left_join(sdat, alphatab, by = "rowname")
    print(names(boxtab))
    
    print(input$checkbox1)

    if(input$checkbox1){
      print("ORDER factor")
      fun = glue::glue( "boxtab${input$Fact1} = factor( boxtab${input$Fact1}, levels = gtools::mixedsort(levels(boxtab${input$Fact1})) ) ")
      eval(parse(text=fun))
    }

    if( !any(names(boxtab)=="sample.id") ) { print("change rowname to sample.id"); dplyr::rename(boxtab, sample.id = rowname) }
    print(head(boxtab))
    
    boxtab$Depth <- sample_sums(r$subdata())
    
    boxtab  
  }
)

  
output$plot2 <- renderPlotly({
   plot_ly(boxtab(), x = as.formula(glue("~{input$Fact1}")), y = as.formula(glue("~{input$metrics}")),
           color = as.formula(glue("~{input$Fact1}")), type = 'box') %>% #, name = ~variable, color = ~variable) %>% #, color = ~variable
     layout(title=input$metrics, yaxis = list(title = glue('{input$metrics}')), xaxis = list(title = 'Samples'), barmode = 'stack') %>%
    config(toImageButtonOptions = list(format = "svg"))
 })

  
  reacalpha <- reactive({
   anova_data = boxtab()
   form1 = glue::glue("{input$metrics} ~ Depth + {input$Fact1}")
   anova_res1 <- aov( as.formula(form1), anova_data)
   outhsd <- HSD.test(anova_res1,input$Fact1)

   LL = list()
   LL$form1 = form1
   LL$aov1 = summary(anova_res1)
   LL$groups1 = outhsd$groups[levels(anova_data[,input$Fact1]),]
   LL$stats1 = outhsd$means[levels(anova_data[,input$Fact1]),]

   LL
 })

 output$testalpha <- renderPrint({
   req(reacalpha)
   LL = reacalpha()
   
   cat("ANOVA\n##########\n")
   print(LL$form1)
   print(LL$aov1)
   cat("\nHSD.test groups\n##########\n")
   print(LL$groups1)
   
 })
 
 output$boxstats <- renderDataTable({
   req(reacalpha)
   LL = reacalpha()
   LL$stats1
 }, filter="top",options = list(pageLength = 5, scrollX = TRUE))
 
 output$boxtab_download <- downloadHandler(
   filename = "alpha_boxplot_stats.csv",
   content = function(file) {
     LL = reacalpha()
     write.table(LL$stats1, file, sep="\t", col.names=NA)}
 )
  
}
    
## To be copied in the UI
# mod_alpha_ui("alpha_ui_1")
    
## To be copied in the server
# callModule(mod_alpha_server, "alpha_ui_1")
 
