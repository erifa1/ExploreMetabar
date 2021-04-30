#' data_loading UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import phyloseq
#' @import DT
#' @import Biostrings
#' @import phyloseq.extended
#'
mod_data_loading_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      infoBox("",
        HTML(paste("New interface to select your data.", br(), "You must validate each step by clicking each button, even if you did not make any modification.", br())),
        icon = icon("info-circle"), fill=TRUE, width = 10
      ),
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
          solidHeader = TRUE, status = "primary", title ="STEP 1: Select your samples", collapsible=TRUE, collapsed=FALSE, width=12,
          fluidPage(
            h3(icon("diagnoses"), "Use table filters to subset your dataset based on your metadata.")
          ),
          DT::dataTableOutput(ns("metadata_table")),
          # actionButton(ns('update_metadata'), "Update Sample", class = "btn-warning")
          shinyBS::bsButton(ns('update_metadata'), "Update Sample", type="action")
        )
      ),
      fluidRow(
        box(
          title = "STEP 2: Select Your Taxonomy Rank and filtering options", solidHeader = TRUE, status = "primary", collapsible=FALSE, collapsed=FALSE,
          selectInput(
            ns("rank_glom"),
            label='Select rank to merge taxonomy table',
            choices='',
            selected = 1,
          ),
          numericRangeInput(ns("minAb"), "Minimum taxa overall raw abundance:", c(1,1), width = NULL, separator = " to "),
          # numericInput(ns("minAb"), "Minimum taxa overall raw abundance:", 1, min = 0, max = NA),
          # numericInput(ns("minPrev"), "Minimum taxa prevalence in samples:", 1, min = 0, max = NA),
          numericRangeInput(ns("minPrev"), "Minimum taxa prevalence in samples:", c(1,1), width = NULL, separator = " to "),
          actionButton(ns('update_taxo'), "Update Taxonomy", class='butt2')
        )
      ),
      fluidRow(
        box(
          solidHeader = TRUE, status = "primary", title ="STEP 3: Select your taxa, preview abundances & representative sequences", collapsible=TRUE, collapsed=FALSE, width=12,
          fluidPage(
            h3(icon("diagnoses"), "Use table filters to subset your dataset based on your taxonomy.")
          ),
          DT::dataTableOutput(ns("taxonomy_table")),
          actionButton(ns('subset_taxo'), "Update Taxonomy", class='butt2')
        )
      ),
      fluidRow(
        box(
          title = 'STEP 4: Normalisation options', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
          radioButtons(
            ns("norm_method"),
            label = "Normalization : ",
            inline = TRUE,
            choices = list(
              "Raw" = 0 ,
              "TSS (total-sum normalization)" = 1,
              "CLR (center log-ration)" = 2,
              "VST (variance stabilizing transformation)" = 3
            ), selected = 1
          ),
          actionButton(ns('norm'), "Normalize", class='butt2')
        ),
        box(
          title = 'Phyloseq final object', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
          verbatimTextOutput(ns("phy_after"))
        ),
        box(
          title = 'Phyloseq normalized object', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
          verbatimTextOutput(ns("phy_norm"))
        )
      ),
      fluidRow(
        box(
          title = 'Download RAW tables', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          downloadButton(outputId = ns("raw_otable_download"), label = "Download raw ASV table"),
          downloadButton(outputId = ns("raw_refseq_download"), label = "Download raw FASTA sequences")
        ),
        box(
          title = 'Download filtered tables', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
          downloadButton(outputId = ns("filt_otable_download"), label = "Download filtered ASV table"),
          downloadButton(outputId = ns("filt_norm_otable_download"), label = "Download filtered and normalized ASV table"),
          downloadButton(outputId = ns("filt_rdata_download"), label = "Download filtered Phyloseq object"),
          downloadButton(outputId = ns("filt_rdata_norm_download"), label = "Download filtered and normalized Phyloseq object"),
          downloadButton(outputId = ns("filt_refseq_download"), label = "Download filtered FASTA sequences")
        )
      )
    )
  )
}


merge_table <- function(rank, table){
  # print("Merge tables")
  FNGdata <- table
  if(rank=="ASV"){
    rank1 = "Species"
  }
  else{
    rank1 = rank
  }
  # print("ttable")
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
  # print("otable")

  rawtaxasum1 <-  table %>%
    taxa_sums() %>%
    as.data.frame %>%
    tibble::rownames_to_column()
  names(rawtaxasum1)[2] <- "RawAbundanceSum"
  # print("taxsum")
  joinGlom <-
    dplyr::left_join(ttable, rawtaxasum1, by = "rowname") %>%
    mutate(RawFreq = RawAbundanceSum / sum(RawAbundanceSum)) %>%
    dplyr::left_join(otable, by = "rowname")

  if(rank=="ASV" & !is.null(refseq(table, errorIfNULL=FALSE)) ){
    # print("add sequence to dataframe")
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
  return(FTAB)
}





#' data_loading Server Function
#'
#' @noRd
mod_data_loading_server <- function(input, output, session, r=r){
  ns <- session$ns
  # r_values <- reactiveValues(phyobj_initial=phyloseq_data(), phyobj_sub_samples=NULL, phyobj_taxglom=NULL, phyobj_final=NULL)
  r_values <- reactiveValues(phyobj_initial=NULL, phyobj_sub_samples=NULL, phyobj_norm=NULL, phyobj_taxglom=NULL, phyobj_final=NULL, phyobj_tmp=NULL)

  phyloseq_data <- reactive({
    cat(file=stderr(), 'phyloseq_data fun', "\n")
    ne <- new.env()
    if (!is.null(input$fileRData)){
      load(input$fileRData$datapath, envir = ne)
    }
    else{
      load(system.file("data_test", "robjects_600.Rdata", package="ExploreMetabar"), envir = ne)
    }
    classes1 = sapply(ne, class)
    obj = classes1[classes1 == "phyloseq"]
    fun = glue::glue("r_values$phyobj_initial <- ne${names(obj)}")
    eval(parse(text = fun))
    r_values$phyobj_tmp <- r_values$phyobj_initial
    # fun = glue::glue("return(ne${names(obj)})")
    # eval(parse(text = fun))

    print(r_values$phyobj_initial)
  })

  output$phy_prev <- renderPrint({
    cat(file=stderr(), 'rendering phy_prev', "\n")
    phyloseq_data()
  })

  sdat <- reactive({
    as.data.frame(as.matrix(phyloseq::sample_data(r_values$phyobj_initial)))
  })
  
  rowCallback <- c(
    "function(row, data){",
    "  for(var i=0; i<data.length; i++){",
    "    if(data[i] === null){",
    "      $('td:eq('+i+')', row).html('NA')",
    "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
    "    }",
    "  }",
    "}"  
  )

  output$metadata_table <- DT::renderDataTable({
    sdat()
  }, filter="top",options = list(pageLength = 10, scrollX = TRUE, rowCallback = JS(rowCallback)), server=TRUE)

  subset_samples <- reactive({
    req(r_values$phyobj_initial)
    cat(file=stderr(), 'subset samples...', "\n")
    cat(file=stderr(), 'initial number of samples before',phyloseq::nsamples(r_values$phyobj_initial), "\n")
    physeq <- phyloseq::prune_samples(phyloseq::sample_names(r_values$phyobj_initial)[input$metadata_table_rows_all],r_values$phyobj_initial)
    physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq)>0, physeq)
    cat(file=stderr(), 'initial number of samples after',phyloseq::nsamples(physeq), "\n")
    r_values$phyobj_sub_samples <- physeq
    r_values$phyobj_tmp <- physeq

  })

  # shinyjs::onclick("update_metadata",{
  #   print("click")
  #   shinyBS::updateButton(session, "update_metadata", style="success")
  # })


  observeEvent(input$update_metadata, {
    cat(file=stderr(), 'button update_metadata', "\n")
    # shinyBS::updateButton(session, "update_metadata", style="success")
    subset_samples()

    # print(r_values$phyobj_sub_samples)
  },
  ignoreNULL = TRUE, ignoreInit = TRUE)


  observe({
    cat(file=stderr(), 'updating rank_glom selectInput...', "\n")
    updateSelectInput(session, "rank_glom",
                      choices = c( rank_names(phyloseq_data()), "ASV" ),
                      selected = "ASV")
  }) #updateSelectInput

  observe({
    cat(file=stderr(), 'updating minAb numericInput...', "\n")
    updateNumericRangeInput(session, 'minAb',"Minimum taxa overall raw abundance:", value=c(1,max(taxa_sums(r_values$phyobj_tmp))))
  }) #updateNumericRangeInput

  observe({
    cat(file=stderr(), 'updating minPrev numericInput...', "\n")
    updateNumericRangeInput(session, 'minPrev',"Minimum taxa prevalence in samples:", value=c(1,max(nsamples(r_values$phyobj_tmp))))
  }) #updateNumericRangeInput

  glom_taxo <- reactive({
    withProgress({

      req(input$minAb, input$minPrev, input$rank_glom, r_values$phyobj_sub_samples)
      cat(file=stderr(), 'filter_taxonomy...', "\n")
      tmp <- r_values$phyobj_sub_samples
      # print(rank_names(tmp))
      withProgress({
        if(input$rank_glom != 'ASV'){
          if(nsamples(tmp)>1000){
            showNotification("Phylogentic tree removed, too much samples...", type="message", duration = 5)
            tmp <- phyloseq.extended::fast_tax_glom(tmp, input$rank_glom)
          }else{
            tmp <- tax_glom(tmp, input$rank_glom)
          }
          FGnames <- tax_table(tmp)[,input$rank_glom]
          nnames <- paste(substr(FGnames, 1, 50), taxa_names(tmp), sep="_")
          taxa_names(tmp) <- nnames
        }
      }, message = 'Taxonomy agglomeration, please wait.')

      tmp <- prune_taxa(taxa_sums(tmp) >= input$minAb[1], tmp)
      tmp <- prune_taxa(taxa_sums(tmp) <= input$minAb[2], tmp)
      prevdf <- apply(X = otu_table(tmp), MARGIN = ifelse(taxa_are_rows(tmp), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
      taxToKeep1 <- names(prevdf)[(prevdf >= input$minPrev[1] & prevdf <= input$minPrev[2])]
      tmp <- prune_taxa(taxToKeep1, tmp)
      if(input$rank_glom != 'ASV'){
        tax_table(tmp) <- tax_table(tmp)[,1:match(input$rank_glom, rank_names(tmp))]
      }
      
      r_values$phyobj_taxglom <- tmp

      r_values$phyobj_tmp <- tmp

      cat(file=stderr(), 'filter_taxonomy done.', "\n")
    },message = "Update taxonomy, please wait...")
  })

  observeEvent(input$update_taxo, {
    glom_taxo()
  },ignoreInit = TRUE)

  render_taxonomy_table <- reactive({
    withProgress({

      req(r_values$phyobj_tmp, input$rank_glom)
      cat(file=stderr(), 'render_taxonomy_table fun', "\n")

      phyloseq_obj <- r_values$phyobj_tmp

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
        # print("add sequence to dataframe")
        showNotification("Sequences added to dataframe.", type="message", duration = 5)
        refseq1 <- phyloseq_obj %>%
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
        FTAB = as.data.frame(joinGlom)
      }
      cat(file=stderr(), 'render_taxonomy_table done.', "\n")
      return(FTAB)
    },message = "Processing, please wait...")

  })

  output$taxonomy_table <- DT::renderDataTable({
    if(ncol(render_taxonomy_table()) > 100){
      showNotification("Truncated abundances for preview...", type="message", duration = 5)
      render_taxonomy_table()[,c(1:20, ncol(render_taxonomy_table()))]
    }else{
      render_taxonomy_table()
    }
  }, filter="top", options = list(pageLength = 10, scrollX = TRUE), server=TRUE)


  subset_taxa <- reactive({
    withProgress({
      req(r_values$phyobj_taxglom)
      cat(file=stderr(), 'subset_taxa fun', "\n")
      selected <- render_taxonomy_table()[input$taxonomy_table_rows_all, 1]
      phy_obj <- prune_taxa(selected, r_values$phyobj_taxglom)
      r_values$phyobj_final <- phy_obj
      r_values$phyobj_tmp <- phy_obj
      cat(file=stderr(), 'subset_taxa fun done.', "\n")
      # phy_obj

    },message = "Subset taxonomy, please wait...")
  })


  observeEvent(input$subset_taxo, {
    cat(file=stderr(), 'button subset_taxo', "\n")
    subset_taxa()
  },ignoreNULL = TRUE, ignoreInit = TRUE)

  normalize <- reactive({
    req(r_values$phyobj_final, input$norm_method)
    FGdata <- r_values$phyobj_final

    if(input$norm_method == 0){
      FNGdata <- FGdata
    }

    if(input$norm_method == 1){
      normf = function(x){ x/sum(x) }
      FNGdata <- transform_sample_counts(FGdata, normf)
    }

    if(input$norm_method == 2){
      clr = function(x){log(x+1) - rowMeans(log(x+1))}
      otable <- otu_table(FGdata)
      otableCLR <- clr(otable)
      FNGdata <- FGdata; otu_table(FNGdata) <- otableCLR
    }

    #VST deseq2
    if(input$norm_method == 3){
      withProgress({
        otable <- FGdata@otu_table@.Data+1
        otableVST <- DESeq2::varianceStabilizingTransformation(otable, fitType='local')
        FNGdata <- FGdata; FNGdata@otu_table@.Data <- otableVST
      },message = "VST normalization, please wait...")
    }
    r_values$phyobj_norm <- FNGdata
    # print(r_values$phyobj_norm)
  })

  observeEvent(input$norm, {
    cat(file=stderr(), 'button normalize', "\n")
    normalize()
  },ignoreNULL = TRUE, ignoreInit = TRUE)


  output$phy_after <- renderPrint({
    cat(file=stderr(), 'rendering phyloseq_after...', "\n")
    print(r_values$phyobj_tmp)
    cat(file=stderr(), 'rendering phyloseq_after done.', "\n")
  })

  output$phy_norm <- renderPrint({
    req(r_values$phyobj_norm)
    cat(file=stderr(), 'rendering phyloseq_norm...', "\n")
    print(r_values$phyobj_norm)
    cat(file=stderr(), 'rendering phyloseq_norm done.', "\n")
  })


  output$raw_otable_download <- downloadHandler(
    filename = "raw_asv_taxtable.csv",
    content = function(file) {
      req(r_values$phyobj_initial)
      write.table(merge_table(input$rank_glom, r_values$phyobj_initial), file, sep="\t", row.names=FALSE)
    }
  )

  output$raw_refseq_download <- downloadHandler(
    filename = "raw_ref-seq.fasta",
    content = function(file) {
      req(r_values$phyobj_initial)
      if(!is.null(refseq(r_values$phyobj_initial, errorIfNULL=FALSE))){
        writeXStringSet(refseq(r_values$phyobj_initial), file)
      }else(showNotification("FASTA Download failed. No refseq in object.", type="error", duration = 5))
    }
  )

  output$filt_otable_download <- downloadHandler(
    filename = "filt_asv_table.csv",
    content = function(file) {
      req(r_values$phyobj_final)
      write.table(merge_table(input$rank_glom, r_values$phyobj_initial), file, sep="\t", row.names=FALSE)
    }
  )

  output$filt_norm_otable_download <- downloadHandler(
    filename = "filt_norm_asv_table.csv",
    content = function(file) {
      req(r_values$phyobj_norm)
      write.table(merge_table(input$rank_glom, r_values$phyobj_norm), file, sep="\t", row.names=FALSE)
    }
  )

  output$filt_rdata_download <- downloadHandler(
    filename = "filt_robject.rdata",
    content = function(file) {
      req(r_values$phyobj_final)
      data = r_values$phyobj_final
      save(data, file = file)
    }
  )

  output$filt_rdata_norm_download  <- downloadHandler(
    filename = "filt_norm_robject.rdata",
    content = function(file) {
      req(r_values$phyobj_norm)
      data = r_values$phyobj_norm
      save(data, file = file)
    }
  )

  output$filt_refseq_download <- downloadHandler(
    filename = "filt_ref-seq.fasta",
    content = function(file) {
      req(r_values$phyobj_final)
      if(!is.null(refseq(r_values$phyobj_final, errorIfNULL=FALSE))){
        Biostrings::writeXStringSet(refseq(r_values$phyobj_final), file)
      }else(showNotification("FASTA Download failed. No refseq in object.", type="error", duration = 5))
    }
  )

  # Saving variable for other modules.
  # Raw object loaded from file.
  r$phyloseq_data <- reactive({
    req(r_values$phyobj_initial)
    r_values$phyobj_initial
  })

  # final filtered object
  r$phyloseq_filtered <- reactive({
    r_values$phyobj_final
    # r_values$phyobj_initial
  })


  # final filtered object normalize
  r$phyloseq_filtered_norm <- reactive({
    req(r_values$phyobj_norm)
    r_values$phyobj_norm
  })

  r$norm_method <- reactive({
    input$norm_method
  })

  # Chosen rank to glom taxa
  r$rank_glom <- reactive({
    input$rank_glom
  })

  r$sdat <- reactive({
    req(r_values$phyobj_final)
    as.data.frame(as.matrix(phyloseq::sample_data(r_values$phyobj_final)))
    # as.data.frame(as.matrix(phyloseq::sample_data(r_values$phyobj_initial)))
  })

}

## To be copied in the UI
# mod_data_loading_ui("data_loading_ui_1")

## To be copied in the server
# callModule(mod_data_loading_server, "data_loading_ui_1")