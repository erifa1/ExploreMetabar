#' data_loading UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_data_loading_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        box(
          title = "Input phyloseq object", status = "warning", solidHeader = TRUE,
          tags$div(
            title = "RData where 'data' is a phyloseq object.",
            fileInput(ns("fileRData"),
                      label = "RData with phyloseq object : ",
                      placeholder = "data.RData")
          )
        ),
        box(
          title = 'Phyloseq preview', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
          verbatimTextOutput(ns("phy_prev"))
        )
      ),
      fluidRow(
        box(
          solidHeader = TRUE, status = "primary", title ="Metadata table", collapsible=TRUE, collapsed=FALSE, width=12,
          fluidPage(
            h3(icon("diagnoses"), "Use table filters to subset your dataset based on your metadata.")
          ),
          dataTableOutput(ns("metadata_table")),
          actionButton(ns('update_metadata'), "Update Sample", class='butt2')
        )
      ),
      fluidRow(
        box(
          title = "Select Your Taxonomy Rank", width = 12, solidHeader = TRUE, status = "warning", collapsible=FALSE, collapsed=FALSE,
          selectInput(
            ns("rank_glom"),
            label='Select rank to merge taxonomy table',
            choices='',
            selected=1
          ),
          numericInput(ns("minAb"), "Minimum taxa overall raw abundance:", 1, min = 0, max = NA),
          numericInput(ns("minPrev"), "Minimum taxa prevalence in samples:", 1, min = 0, max = NA),
          actionButton(ns('update_taxo'), "Update Taxonomy", class='butt2')
        )
      ),
      fluidRow(
        box(
          solidHeader = TRUE, status = "primary", title ="Taxonomy table", collapsible=TRUE, collapsed=FALSE, width=12,
          fluidPage(
            h3(icon("diagnoses"), "Use table filters to subset your dataset based on your taxonomy.")
          ),
          dataTableOutput(ns("taxonomy_table")),
          actionButton(ns('subset_taxo'), "Update Taxonomy", class='butt2')
        )
      ),
      fluidRow(
        box(
          title = 'Phyloseq final object', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
          verbatimTextOutput(ns("phy_after"))
        )
      )
    )
  )
}

#' data_loading Server Function
#'
#' @noRd
mod_data_loading_server <- function(input, output, session, r=r){
  ns <- session$ns
  phyloseq_data <- reactive({
    ne <- new.env()
    if (!is.null(input$fileRData)){
      load(input$fileRData$datapath, envir = ne)
    }
    else{
      load(system.file("data_test", "robjects_600.Rdata", package="ExploreMetabar"), envir = ne)
    }
    classes1 = sapply(ne, class)
    obj = classes1[classes1 == "phyloseq"]
    fun = glue::glue("return(ne${names(obj)})")
    eval(parse(text = fun))
  })

  output$phy_prev <- renderPrint({
    phyloseq_data()
  })

  sdat <- reactive({
    as.data.frame(as.matrix(phyloseq::sample_data(phyloseq_data())))
  })

  output$metadata_table <- DT::renderDataTable({
    sdat()
  }, filter="top",options = list(pageLength = 10, scrollX = TRUE))

  subset_samples <- reactive({
    cat(file=stderr(), 'subset samples...', "\n")
    cat(file=stderr(), 'initial number of samples before',phyloseq::nsamples(phyloseq_data()), "\n")
    physeq <- phyloseq::prune_samples(phyloseq::sample_names(phyloseq_data())[input$metadata_table_rows_all],phyloseq_data())
    physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq)>0, physeq)
    cat(file=stderr(), 'initial number of samples after',phyloseq::nsamples(physeq), "\n")
    physeq
  })

  observeEvent(input$update_metadata, {
    subset_samples()
  })

  observe({
    updateSelectInput(session, "rank_glom",
                      choices = c( rank_names(r$phyloseq_data()), "ASV" ),
                      selected = "ASV")
  })

  filter_taxonomy <- reactive({
    req(input$minAb, input$minPrev, input$rank_glom)
    cat(file=stderr(), 'filtering taxonomy...', "\n")
    tmp <- subset_samples()
    withProgress({
      if(input$rank_glom != 'ASV'){
        tmp <- tax_glom(tmp, input$rank_glom)
        FGnames <- tax_table(tmp)[,input$rank_glom]
        nnames <- paste(substr(FGnames, 1, 50), taxa_names(tmp), sep="_")
        taxa_names(tmp) <- nnames
      }
    }, message = 'Taxonomy agglomeration, please wait.')

    tmp <- prune_taxa(taxa_sums(tmp) > input$minAb, tmp)
    prevdf <- apply(X = otu_table(tmp), MARGIN = ifelse(taxa_are_rows(tmp), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
    taxToKeep1 <- names(prevdf)[(prevdf >= input$minPrev)]
    tmp <- prune_taxa(taxToKeep1, tmp)

    tmp
  })

  observeEvent(input$update_taxo, {
    filter_taxonomy()
  })


  render_taxonomy_table <- reactive({
    cat(file=stderr(), 'Rendering taxonomy table', "\n")
    phyloseq_obj <- filter_taxonomy()

    if(input$rank_glom=="ASV"){
      rank1 = "Species"
    }
    else{
      rank1 = input$rank_glom
    }

    ttable <- phyloseq_obj %>%
      tax_table() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      dplyr::select(1:rank1) %>%
      tibble::rownames_to_column() %>%
      as.matrix() %>% as.data.frame()

    otable <- phyloseq_obj %>%
      otu_table() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      tibble::rownames_to_column()

    rawtaxasum1 <-  phyloseq_obj %>%
      taxa_sums() %>%
      as.data.frame %>%
      tibble::rownames_to_column()
    names(rawtaxasum1)[2] <- "RawAbundanceSum"

    joinGlom <-
      dplyr::left_join(ttable, rawtaxasum1, by = "rowname") %>%
      mutate(RawFreq = RawAbundanceSum / sum(RawAbundanceSum)) %>%
      dplyr::left_join(otable, by = "rowname")

    if(input$rank_glom=="ASV" & !is.null(refseq(phyloseq_obj, errorIfNULL=FALSE)) ){
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

  output$taxonomy_table <- DT::renderDataTable({
    render_taxonomy_table()
  }, filter="top", options = list(pageLength = 10, scrollX = TRUE))



  subset_taxa <- reactive({
    selected <- render_taxonomy_table()[input$taxonomy_table_rows_all, 1]
    phy_obj <- prune_taxa(selected, filter_taxonomy())
    phy_obj
  })

  observeEvent(input$subset_taxo, {
    subset_taxa()
  })

  phyloseq_after <- reactive({
    subset_taxa()
  })

  output$phy_after <- renderPrint({
    phyloseq_after()
  })

  # Saving variable for other modules.
  # Raw object loaded from file.
  r$phyloseq_data <- reactive(
    phyloseq_data()
  )

  # sample + glom + taxa filtered
  r$data <- reactive({
    subset_taxa()
  })

  # Chosen rank to glom taxa
  r$rank_glom <- reactive(
    input$rank_glom
  )

}

## To be copied in the UI
# mod_data_loading_ui("data_loading_ui_1")

## To be copied in the server
# callModule(mod_data_loading_server, "data_loading_ui_1")
