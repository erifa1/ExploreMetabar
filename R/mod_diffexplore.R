# Module UI
  
#' @title   mod_diffexplore_ui and mod_diffexplore_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_diffexplore
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList 
#' @import ggplot2

mod_diffexplore_ui <- function(id){
  ns <- NS(id)
  tagList(
    fileInput(ns("diffData"),
              label = "Select an aggregate_diff_*.csv file : ",
              placeholder = "yourfile.csv"),
    h3("Sélectionner une comparaison à afficher:"),
    DT::dataTableOutput(ns("DTdiff")),
    # verbatimTextOutput(ns("print1")),
    plotOutput(ns("DiffPlot"), height="600px")
  
  )
}
    
# Module Server
    
#' @rdname mod_diffexplore
#' @export
#' @keywords internal

    
mod_diffexplore_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  diff1 <- reactive({
    print("InputTable")
    ne <- new.env() ## new env to store RData content and avoid border effects
    if (!is.null(input$diffData)){
      print(input$diffData$datapath)
      ne$A1 <- read.table(input$diffData$datapath, sep="\t", h=TRUE) 
    } else {return(NULL)}
    if (class(ne$A1) == "data.frame")
      return(ne$A1)
  })
  
  output$DTdiff <- DT::renderDataTable({
    print("RenderTable")
    # print(head(diff1()[,1:10]))
    if (is.null(diff1()))
      return(NULL)
    diff1()
  },filter = "top", options = list(pageLength = 5, scrollX = TRUE),
  server = FALSE)
  
  
  # output$print1 <- renderPrint({
  #   print(head(diff1()))
  # })
  
  p1 <- reactive({
    if (is.null(diff1()))
      return(NULL)
    
    print("plot")
    TABbar=diff1()[input$DTdiff_rows_all,] %>% arrange(desc(absDESeqLFC))
    
    if (length(unique(TABbar$Condition))>2)
      return(NULL)
    
    TABbar$tax = paste( substr(TABbar[,"Species"],1,40),"...","_", TABbar[,"seqid"], sep="")
    # print(head(TABbar))
    
    p<-ggplot(data=TABbar[1:50,], aes(x=reorder(tax, -abs(DESeqLFC)), y=DESeqLFC, fill=Condition ) ) +
      geom_bar(stat="identity", alpha = 0.7) + ggtitle("") + labs(x='Features') +
      coord_flip() + theme_bw() +
      scale_y_continuous(minor_breaks = seq(-1E4 , 1E4, 1), breaks = seq(-1E4, 1E4, 5))
    print(p)
  })
  
  output$DiffPlot <- renderPlot({
    if (is.null(p1()))
      return(NULL)
    p1()
  })
  
}
    
## To be copied in the UI
# mod_diffexplore_ui("diffexplore_ui_1")
    
## To be copied in the server
# callModule(mod_diffexplore_server, "diffexplore_ui_1")
 
