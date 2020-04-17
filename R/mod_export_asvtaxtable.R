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
      box(
        h3("Glom and normalized object will be used for next modules"),
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
            `TSS` = "TSS norm.",
            `CLR` = "CLR norm.",
            `VST` = "VST norm."
          ), selected = "TSS norm."
        ),
        numericInput(ns("minAb"), "Minimum raw abundance:", 1, min = 0, max = NA),
        title = "Settings:", width = 12, status = "primary", solidHeader = TRUE
      ),
      

      box(
        h3("Raw object:"),
        verbatimTextOutput(ns("print1")),
        h3("Glom/Subset object:"),
        verbatimTextOutput(ns("print2")),
        verbatimTextOutput(ns("print3")),
        title = "Details", collapsible = TRUE, collapsed = TRUE, width = 10, status = "primary", solidHeader = TRUE),
      box(
        downloadButton(outputId = ns("otable_download"), label = "Download Table"),
        downloadButton(outputId = ns("refseq_download"), label = "Download FASTA sequences"),
        downloadButton(outputId = ns("rdata_download"), label = "Download transformed Phyloseq object"),
        dataTableOutput(ns("otable1")),
        title = "ASV table", width = 12, status = "primary", solidHeader = TRUE
      )
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
#' @importFrom Biostrings writeXStringSet
    
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
    Fdata@phy_tree <- NULL # to improve speed remove TREE 
    # print(head(otu_table(Fdata)))
    # Glom
    print("Glom object")
    withProgress({
      
      if(input$RankGlom != "ASV"){
        FGdata <- tax_glom(Fdata, input$RankGlom)
        # FGotab <- otu_table(FGdata); 
        FGnames <- tax_table(FGdata)[,input$RankGlom]
        nnames <- paste(substr(FGnames, 1, 50), taxa_names(FGdata), sep="_")
        taxa_names(FGdata) <- nnames
      }else{FGdata <- Fdata}
      FGdata
    }, message = "Glom step on global table, please wait...")
    
  })
  
  output$print2 <- renderPrint({
    print(glom())
    # r$rowselect()
  })

  subglom <- reactive({
    req(input$minAb, glom())
    print("Subset object")
    Fdata <- prune_samples(sample_names(glom())[r$rowselect()], glom())
    Fdata <- prune_taxa(taxa_sums(Fdata) > input$minAb, Fdata)
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
    
    rawtaxasum1 <-  subglom() %>%
      taxa_sums() %>%
      as.data.frame %>%
      tibble::rownames_to_column()
    names(rawtaxasum1)[2] <- "RawAbundanceSum"
    
    joinGlom <-
      dplyr::left_join(ttable, rawtaxasum1, by = "rowname") %>%
      mutate(RawFreq = RawAbundanceSum / sum(RawAbundanceSum)) %>%
      dplyr::left_join(otable, by = "rowname") 
    
    if(input$RankGlom=="ASV" & !is.null(refseq(dat(), errorIfNULL=FALSE)) ){
      print("add sequence to dataframe")
      showNotification("Sequences added to dataframe.", type="message", duration = 5)
      refseq1 <- FNGdata %>%
        refseq %>%
        as.data.frame %>%
        tibble::rownames_to_column() %>%
        rename(sequences = x)
        
      joinGlom2 <- dplyr::left_join(joinGlom, refseq1, by = "rowname") %>%
        dplyr::rename(asvname = rowname)
      FTAB = as.data.frame(joinGlom2)
    }else{
      showNotification("No refseq in object.", type="error", duration = 5)
      dplyr::rename(joinGlom, asvname = rowname)
      # print(str(as.data.frame(as.matrix(ttable))))
      FTAB = as.data.frame(joinGlom)
    }
    
    FTAB

  })

  output$otable1 <- DT::renderDataTable({
    merge1()
  }, filter="top", options = list(scrollX = TRUE))
  

  output$otable_download <- downloadHandler(
    filename = "asv_taxtable.csv",
    content = function(file) {
        write.table(merge1(), file, sep="\t", row.names=FALSE)
    }
  )

  output$refseq_download <- downloadHandler(
    filename = "ref-seq.fasta",
    content = function(file) {
      req(dat())
      if(!is.null(refseq(dat(), errorIfNULL=FALSE))){
        writeXStringSet(refseq(dat()), file)
      }else(showNotification("FASTA Download failed. No refseq in object.", type="error", duration = 5))
    }
  )
  
  output$rdata_download <- downloadHandler(
    filename = "robject.rdata",
    content = function(file) {
      req(subglom())
        data = subglom()
        save(data, file = file)
    }
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
 
