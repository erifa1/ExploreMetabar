# Module UI
  
#' @title   mod_metadata_subset_ui and mod_metadata_subset_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_metadata_subset
#'
#' @keywords internal
#' @export 
#' @importFrom shiny NS tagList 
#' @import phyloseq
#' @importFrom DT dataTableOutput

mod_metadata_subset_ui <- function(id){
  ns <- NS(id)
  
  tagList(
    fluidPage(h1("Metadatas / Subset")),
    
    
    tags$div(
      title = "RData where 'data' is a phyloseq object.",
      fileInput(ns("fileRData"),
                label = "RData : ",
                placeholder = "data.RData")
    ),
    
    verbatimTextOutput(ns("print1")),
    h3("Metadata table"),
    fluidPage(h4("Use filters to subset phyloseq object")),
    dataTableOutput(ns("sdata3")),
    h3("Selected samples names"),
    verbatimTextOutput(ns("sids"))
  )
}
    
# Module Server
    
#' @rdname mod_metadata_subset
#' @export
#' @keywords internal
#' @importFrom DT renderDataTable
    
mod_metadata_subset_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  data16S <- reactive({
    ne <- new.env() ## new env to store RData content and avoid border effects
    if (!is.null(input$fileRData)){
      load(input$fileRData$datapath, envir = ne) 
    } else {load("~/Bureau/robjects_600.Rdata", envir = ne)}
    if (class(ne$data) == "phyloseq")
      return(ne$data)
  })
  
  sdat <- reactive({
    as.data.frame(as.matrix(sample_data(data16S())))
  })
  
  output$sdata3 <- renderDataTable({
    sdat()
  }, filter="top",options = list(pageLength = 5, scrollX = TRUE))
  
  output$print1 <- renderPrint({
    data16S()
  })
  
  sids <- output$sids <- reactive({
    stab <- as.data.frame(as.matrix(sample_data(data16S())))
    select  <- row.names(stab[input$sdata3_rows_all,])
    return(select)
  })
  
  subdata <- reactive({
    print("subset")
    Fdata <- prune_samples(sample_names(data16S())[input$sdata3_rows_all], data16S())
    Fdata <- prune_taxa(taxa_sums(Fdata) > 0, Fdata)
    Fdata
  })
  

#Saving variable for other modules.
  r$data16S <- reactive(
    data16S()
  )
  
  r$rowselect <- reactive({
    input$sdata3_rows_all
  })
  
  r$subdata <- reactive({
    subdata()
  })
}


# To improve
#add download button for metadata table
#add option to fix minimum abundance for ASV

    
## To be copied in the UI
# mod_metadata_subset_ui("metadata_subset_ui_1")
    
## To be copied in the server
# callModule(mod_metadata_subset_server, "metadata_subset_ui_1")
 
