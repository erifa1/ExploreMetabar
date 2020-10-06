#' asvenn UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_asvenn_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    fluidPage(
      
      infoBox("", 
              "Select conditions to highlight shared taxa", 
              icon = icon("info-circle"), fill=TRUE, width = 10),
      box(
        selectInput(
          ns("Fact1"),
          label = "Select factor to test: ",
          choices = ""
        ),
        uiOutput(ns("lvls1")),
        numericInput(ns("minAb"), "Minimum raw abundance:", 1, min = 1, max = NA),
        actionButton(ns("go1"), "Run/Update ASVenn", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(
        plotOutput(ns("venn1"), height = "800px"),
        title = "Venn Diagram:", width = 12, status = "primary", solidHeader = TRUE, height = 900
      ),
      box(
        downloadButton(outputId = ns("otable_download"), label = "Download Table"),
        dataTableOutput(ns("tabvenn1")),
        title = "Venn table:", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE
      )
    )
  )
}
    
#' asvenn Server Function
#'
#' @noRd 
#' 
#' @importFrom futile.logger flog.threshold
#' @importFrom grid grid.draw
#' @importFrom grDevices rainbow
#' @importFrom VennDiagram calculate.overlap
#' @importFrom VennDiagram venn.diagram
#' @importFrom ranomaly ASVenn_fun
#' 

mod_asvenn_server <- function(input, output, session, r=r){
  ns <- session$ns
 
  observe({
    updateSelectInput(session, "Fact1",
                      choices = r$data16S()@sam_data@names)
  })
  
# output$lvls1 <- reactive({
#   level1 <- na.omit(levels(as.factor(sample_data(data)[,input$Fact1]@.Data[[1]])) )
#   level1
# })
  
output$lvls1 = renderUI({
  req(input$Fact1)
  level1 <- na.omit(levels(as.factor(sample_data(r$data16S())[,input$Fact1]@.Data[[1]])) )
  checkboxGroupInput(ns("lvls1"), label = "Select up to 5 levels :", 
                     choices = level1, inline = TRUE, selected = level1[1:3])
  
})
  
resVenn <- eventReactive(input$go1, {   #TF
  req(r$subglom(), r$asvselect())
  
  data <- prune_taxa(taxa_sums(r$subglom()) > input$minAb, r$subglom())
  data <- prune_taxa(r$asvselect(), data)
  print(data)
  
  resVenn = ASVenn_fun(data = data, output = "./ASVenn/", rank = "ASV", column1 = input$Fact1, lvls = input$lvls1, shared = TRUE)
  
  resVenn
  
  
  # #Nombre d'espèce par matrice
  # print('Parsing factor ...')
  # # level1 <- na.omit(levels(as.factor(sample_data(data)[,input$Fact1]@.Data[[1]])) )
  # TFdata <- list()
  # TFtax <- list()
  # databak <- data
  # for(i in 1:length(input$lvls1)){
  #   databak -> data
  #   LOC=as.character(input$lvls1[i])
  #   print(LOC)
  #   fun <- paste("data <- subset_samples(data, ",input$Fact1," %in% '",LOC,"')",sep="")
  #   eval(parse(text=fun))
  #   data <- prune_taxa(taxa_sums(data) > 0, data)
  # 
  #   ttable <- data@tax_table@.Data
  #   otable <- as.data.frame(otu_table(data))
  #   # print(nrow(ttable))
  # 
  #   if(!any(rownames(ttable) == rownames(otable))){print("Different order in otu table and tax table")} #;quit()
  # 
  #   TT = cbind(otable,ttable)
  #   TFdata[[i]] <- TT
  #   TFtax[[i]] <- cbind(row.names(TT), as.character(apply(TT[,colnames(ttable)], 1, paste, collapse=";") ) ) #c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  #   row.names(TFtax[[i]]) = TFtax[[i]][,1]
  #   # write.table(TT, paste(output,"/otu_table_sp_",LOC,".csv",sep=""), sep="\t", quote=FALSE, col.names=NA)
  # 
  # }
  # 
  # 
  # 
  # TF <- sapply(TFtax, row.names)
  # names(TF) = input$lvls1
  # 
  # print('Defining unique taxa ...')
  # alltax <- do.call(rbind, TFtax)
  # alltax <- alltax[!duplicated(alltax[,1]),]
  # row.names(alltax)=alltax[,1]
  
  ## END asvennfunction 
  
  # print('return TF')
  # TFlist = list()
  # TFlist$TF = TF
  # TFlist$alltax = alltax
  # print(str(TFlist))
  # TFlist
})


# PLOT  
# venn1 <-eventReactive(input$go1,{
#   print("venn react")
#   req(TF())
#   invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
#   if(length(input$lvls1) > 1 & length(input$lvls1) <= 5){
#     venn.plot <- venn.diagram(TF()$TF, filename = NULL, col = "black",
#                               alpha = 0.50, fill = rainbow(length(input$lvls1)),
#                               cex = 1.5, cat.col = 1, lty = "blank",
#                               cat.cex = 1.8, cat.fontface = "bold",
#                               margin = 0.1, main=glue::glue("{input$Fact1}"), main.cex=2.5,
#                               fontfamily ="Arial",main.fontfamily="Arial",cat.fontfamily="Arial",
#                               height = 3000, width = 5333) #cat.dist = 0.09,
#     grid.draw(venn.plot)
#   }else{
#     showNotification("Choose 1 to 5 levels...", type="error", duration = 5)
#     return(NULL)
#     }
# 
# })

output$venn1 <- renderPlot({
  invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
  if(length(input$lvls1) <= 5){
    grid.draw( resVenn()$venn_plot )
  }else{
    resVenn()$venn_plot
  }
})
  
# tabvenn1 <- eventReactive(input$go1, {
#   print("renderDT")
#   TF<- TF()$TF
#   print(names(TF))
#   alltax <- TF()$alltax
#   
#   if(length(input$lvls1) > 1 & length(input$lvls1) <= 5){
#     ov <- VennDiagram::calculate.overlap(TF)
#     print(sapply(ov, length))
#     
#     print('Calculating lists ...')
#     uniqTax = TABf = unique(do.call(c,TF))
#     for (j in 1:length(input$lvls1)){
#       TABtest = TF[[j]]
#       TABtest_filt=rep(0, length(uniqTax))
#       for (i in 1:length(uniqTax)) {
#         featureI = uniqTax[i]
#         res=grep( paste('^',featureI,'$', sep="" ) , TABtest)
#         if(length(res)>0){TABtest_filt[i]=length(res)
#         }
#       }
#       TABf=cbind.data.frame( TABtest_filt, TABf )
#       names(TABf)[1] = names(TF)[j]
#     }
#     
#     if(!is.null(refseq(r$data16S(), errorIfNULL=FALSE))){
#       refseq1 <- as.data.frame(refseq(r$data16S())); names(refseq1)="seq"
#     }else{print('No Sequences ...')}
#     
#     if(exists("refseq1")){
#       TABf <- cbind(TABf,alltax[as.character(TABf$TABf),2], refseq1[as.character(TABf$TABf),])
#       names(TABf) <- c(rev(names(TF)), "ASV", "taxonomy", "seq")
#     }else{
#       TABf <- cbind(TABf,alltax[as.character(TABf$TABf),2])
#       names(TABf) <- c(rev(names(TF)), "ASV", "taxonomy")
#     }
#     
#     TABf
#   }else{return(NULL)}
# }
# )
  
output$tabvenn1 <-DT::renderDataTable({
  resVenn()$TABf #tabvenn1()
}, filter="top", options = list(scrollX = TRUE))

output$otable_download <- downloadHandler(
  filename = "venn_table.csv",
  content = function(file) {
    write.table(tabvenn1(), file, sep="\t", row.names=FALSE)
  }
)
  
  
  
}
    
## To be copied in the UI
# mod_asvenn_ui("asvenn_ui_1")
    
## To be copied in the server
# callModule(mod_asvenn_server, "asvenn_ui_1")
 
