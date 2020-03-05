# Module UI
  
#' @title   mod_export_asvtaxtable_ui and mod_export_asvtaxtable_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_export_asvtaxtable
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList 
#' @importFrom DT dataTableOutput

mod_export_asvtaxtable_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      h1("ASV Taxonomy Table"),
      selectInput(
        ns("RankGlom"),
        label = "Select rank to glom : ",
        choices = "",
        selected = 1
      ),
      radioButtons(
        ns("NORM"),
        label = "Normalization : ",
        inline = TRUE,
        choices = list(
          `No Norm` = "Raw",
          `CLR` = "CLR norm.",
          `TSS` = "TSS norm.",
          `VST` = "VST norm."
        ), selected = "CLR norm."
      ),
      verbatimTextOutput(ns("print1")),
      h3("Glom/Subset Object:"),
      verbatimTextOutput(ns("print2")),
      verbatimTextOutput(ns("print3")),
      h1("ASV table"),
      fluidRow(column(dataTableOutput(ns("otable1")) , width = 12)),
      downloadButton(outputId = ns("otable_download"), label = "Download Table")
    
    )
  )
}
    
# Module Server
    
#' @rdname mod_export_asvtaxtable
#' @export
#' @keywords internal
#' @import dplyr 
#' @import tibble
#' @importFrom DT renderDataTable
#' @importFrom DESeq2 varianceStabilizingTransformation
    
mod_export_asvtaxtable_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  output$print1 <- renderPrint({
    r$data16S()
  })

  observe({
    updateSelectInput(session, "RankGlom",
                      choices = c( rank_names(r$data16S()), "ASV" ),
                      selected = "ASV")
  })

  glom <- reactive({
    if (is.null(r$data16S()))
      return(NULL)
    Fdata <- r$data16S()
    # print(head(otu_table(Fdata)))
    # Glom
    print("Glom")
    withProgress({
      
      if(input$RankGlom != "ASV"){
        FGdata <- tax_glom(Fdata, input$RankGlom)
        # FGotab <- otu_table(FGdata); 
        FGnames <- tax_table(FGdata)[,input$RankGlom]
        nnames <- paste(substr(FGnames, 1, 20), taxa_names(FGdata), sep="_")
        taxa_names(FGdata) <- nnames
      }else{FGdata <- Fdata}
      FGdata
    }, message = "Glom step on global table, please wait...")
    
  })
  
  output$print2 <- renderPrint({
    print(glom())
    r$rowselect()
  })

  subglom <- reactive({
    print("subset")
    Fdata <- prune_samples(sample_names(glom())[r$rowselect()], glom())
    Fdata <- prune_taxa(taxa_sums(Fdata) > 0, Fdata)
    Fdata
  })
  
  output$print3 <- renderPrint({
    subglom()
  })

  
  dat <- reactive({
    print("normalize")
    if (is.null(r$data16S()))
      return(NULL)
    FGdata <- subglom()
    print("Norm")
    if(input$NORM=="TSS norm."){
      normf = function(x){ x/sum(x) }
      FNGdata <- transform_sample_counts(FGdata, normf)
    }

    if(input$NORM=="Raw"){
      FNGdata <- FGdata
    }

    if(input$NORM=="CLR norm."){
      clr = function(x){log(x+1) - rowMeans(log(x+1))}
      otable <- otu_table(FGdata)
      otableCLR <- clr(otable)
      FNGdata <- FGdata; otu_table(FNGdata) <- otableCLR
    }

    #VST deseq2
    if(input$NORM=="VST norm."){
      withProgress({
        otable <- FGdata@otu_table@.Data+1
        otableVST <- varianceStabilizingTransformation(otable, fitType='local')
        FNGdata <- FGdata; FNGdata@otu_table@.Data <- otableVST
      },message = "VST normalization, please wait...")
    }

    FNGdata
  })

  merge1 <- reactive({
    print("Merge tables")
    FNGdata <- dat()
    if(input$RankGlom=="ASV"){rank1 = "Species"}else{rank1 = input$RankGlom}
    ttable <- FNGdata %>%
      tax_table() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      dplyr::select(1:rank1) %>%
      tibble::rownames_to_column() %>%
      as.matrix() %>% as.data.frame()
    otable <- FNGdata %>%
      otu_table() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      tibble::rownames_to_column()
    joinGlom <-
      dplyr::left_join(ttable, otable, by = "rowname") %>%
      dplyr::rename(asvname = rowname)

    print(str(as.data.frame(as.matrix(ttable))))
    as.data.frame(joinGlom)

  })

  output$otable1 <- DT::renderDataTable({
    merge1()
  }, filter="top", options = list(scrollX = TRUE))
  

  output$otable_download <- downloadHandler(
    filename = "asv_taxtable.csv",
    content = function(file) {write.table(merge1(), file, sep="\t", row.names=FALSE)}
  )

  #Saving variable for other modules.
  # Object only glom
  r$subglom <- reactive(
    subglom()
  )
  
  # Object glom + norm
  r$dat <- reactive(
    dat()
  )
  r$RankGlom <- reactive(
    input$RankGlom
  )
  
  
}
    
## To be copied in the UI
# mod_export_asvtaxtable_ui("export_asvtaxtable_ui_1")
    
## To be copied in the server
# callModule(mod_export_asvtaxtable_server, "export_asvtaxtable_ui_1")
 
