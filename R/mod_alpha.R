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
      radioButtons(ns("metrics"), "Choose one index:", inline = TRUE,
                         choices =
                           list("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                                "InvSimpson", "Fisher"),
                         selected = c("Shannon")
      ),
      numericInput(ns("minAb"), "Minimum raw abundance:", 1, min = 1, max = NA),
      selectInput(
        ns("Fact1"),
        label = "Select factor to test: ",
        choices = ""
      ),
      
      verbatimTextOutput(ns("print1")),
      dataTableOutput(ns("alphaout")),
      # plotOutput(ns("plotalpha1")),
      plotlyOutput(ns("plot2"))
      # plotOutput(ns("plot3"))
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
    
mod_alpha_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  observe({
    updateSelectInput(session, "Fact1",
                      choices = r$data16S()@sam_data@names)
  })
  
  alpha1 <- reactive({
    if(is.null(r$subdata())){return(NULL)}
    print(input$metrics)
    print(input$minAb)
    
    data <- prune_taxa(taxa_sums(r$subdata()) > input$minAb, r$subdata())
    print(data)
    
    alphatab <- estimate_richness(data, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                                                     "InvSimpson", "Fisher") )
    
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

  output$print1 <- renderPrint({
    LL = alpha1()
    LL$data
  })

 plot1 <- reactive({
   LL = alpha1()
   p <- plot_richness(LL$data,x=input$Fact1, color=input$Fact1, measures=input$metrics)
   p$layers <- p$layers[-1]
   p1 <- p + ggtitle('Alpha diversity indexes') +  geom_boxplot(position = position_dodge(width = 0.5),alpha = 0.7, outlier.shape = NA) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=11), plot.title = element_text(hjust = 0.5)) + theme_bw()
   
   p1
 })

 output$plotalpha1 <-renderPlot({
   print("renderAlphaplot")
   if (is.null(plot1()))
     return(NULL)
   ggplotly(plot1())
 })
 
 output$plot2 <- renderPlotly({
   LL = alpha1()
   sdat = tibble::rownames_to_column(as.data.frame(as.matrix(r$subdata()@sam_data)))
   alphatab =  tibble::rownames_to_column(LL$alphatab)
   
   print(sdat)
   print(alphatab)
   
   boxtab <- dplyr::left_join(sdat, alphatab, by = "rowname") %>%
     dplyr::rename(sample.id = rowname)
   boxtab
   print(dim(boxtab))
   
   plot_ly(boxtab, x = as.formula(glue("~{input$Fact1}")), y = as.formula(glue("~{input$metrics}")),
           color = as.formula(glue("~{input$Fact1}")), type = 'box') %>% #, name = ~variable, color = ~variable) %>% #, color = ~variable
     layout(title=input$metrics, yaxis = list(title = glue('{input$metrics}')), xaxis = list(title = 'Samples'), barmode = 'stack')
 })

 output$plot3 <- renderPlot({
   LL = alpha1()
   sdat = tibble::rownames_to_column(as.data.frame(as.matrix(r$subdata()@sam_data)))
   alphatab =  tibble::rownames_to_column(LL$alphatab)
   boxtab <- dplyr::left_join(sdat, alphatab, by = "rowname") %>%
     dplyr::rename(sample.id = rowname)
   
   boxplot(as.formula(glue("{input$metrics} ~ {input$Fact1}")), color=input$Fact1, data=boxtab)
 })

  
}
    
## To be copied in the UI
# mod_alpha_ui("alpha_ui_1")
    
## To be copied in the server
# callModule(mod_alpha_server, "alpha_ui_1")
 
