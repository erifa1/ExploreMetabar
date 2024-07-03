#' data_loading UI Function
#'
#' @description Module for loading phyloseq object from rdata file. This module allows to filter and select samples and taxa prior to analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom phyloseq sample_data nsamples prune_samples prune_taxa taxa_sums
#' @importFrom DT dataTableOutput renderDataTable JS
#' @importFrom Biostrings writeXStringSet
#' @importFrom shinyBS bsButton updateButton
#' @importFrom glue glue
#' @importFrom futile.logger flog.info flog.debug
#' @import datamods
#'
mod_data_loading_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(infoBox("",
        HTML(paste("You must validate each step (filtering, normalization) by clicking each button, even if you did not make any modification.")),
        HTML(paste("Otherwise you just need to click 'Launch all' button, then you can use others modules.")),
        icon = icon("info-circle"), fill=TRUE, width = 6
      )),

      fluidRow(
        box(
          title = "Input phyloseq object", status = "warning", solidHeader = TRUE,
          tags$div(
            title = "RData where 'data' is a phyloseq object.",
            fileInput(ns("fileRData"),
                      label = "RData with phyloseq object : ",
                      placeholder = "data.RData")
          ),
          shinyBS::bsButton(inputId = ns('launch_all'), label = "Launch all", block = F, style = 'danger', type='action')
        ),
        box(
          title = 'Phyloseq preview', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
          verbatimTextOutput(ns("phy_prev"))
        )
      ),
      
      fluidRow(
        uiOutput(ns('ui_combine_columns'))
      ),

      fluidRow(box(title = "STEP 1: Metadata table",solidHeader = TRUE, status = "warning", width=12,

        tabBox(width=12,

          tabPanel("Metadata / Filters",
            tags$h3(icon("diagnoses"), "Use filters to subset your dataset based on your metadata :"),

              fluidRow(
                column(
                  width = 3,
                  datamods::filter_data_ui(ns("filtering"), max_height = "500px")
                ),
                column(
                  width = 9,
                  # progressBar(
                  #   id = ns("pbar"), value = 100,
                  #   total = 100, display_pct = TRUE
                  # ),
                  DT::dataTableOutput(outputId = ns("metadata_table"))
                )
              )

            ),
          tabPanel("Update variables",
            fluidRow(
                column(
                  width = 12,
                  datamods::update_variables_ui(ns("vars"))
                )
              )
            )
          ),
          actionButton(ns('update_metadata'), " Update Sample", icon("paper-plane"),
                  style="color: #fff; background-color: #D73925; border-color: #DB4836")
        )),

      fluidRow(
        box(
          title = "STEP 2: Taxonomy rank and filtering options", solidHeader = TRUE, status = "primary", collapsible=FALSE, collapsed=FALSE,
          selectInput(
            ns("rank_glom"),
            label='Select rank to merge taxonomy table',
            choices='',
            selected = 1,
          ),
          shinyBS::bsButton(inputId = ns('update_taxo0'), label = "Launch glom", block = F, style = 'danger', type='action'),
          numericRangeInput(ns("minAb"), "Minimum taxa overall raw abundance:", c(1,1), width = NULL, separator = " to "),
          numericRangeInput(ns("minPrev"), "Minimum taxa prevalence in samples:", c(1,1), width = NULL, separator = " to "),
          shinyBS::bsButton(inputId = ns('update_taxo'), label = "Update Filters", block = F, style = 'danger', type='action')
        )
      ),

      fluidRow(box(title = "STEP 3: Taxa filtering, preview abundances & representative sequences",solidHeader = TRUE, status = "warning", width=12,

        # tabBox(width=12,

          # tabPanel("Metadata / Filters",
            h3(icon("diagnoses"), "Use filters to subset your dataset based on taxonomy."),

              fluidRow(
                column(
                  width = 3,
                  datamods::filter_data_ui(ns("filtering_taxo"), max_height = "500px")
                ),
                column(
                  width = 9,
                  # progressBar(
                  #   id = ns("pbar"), value = 100,
                  #   total = 100, display_pct = TRUE
                  # ),
                  DT::dataTableOutput(outputId = ns("table_taxoFILT"))
                )
              ),
          # actionButton(inputId = ns('subset_taxo'), label = "Update Taxo"),
          actionButton(ns('subset_taxo'), " Update Taxo", icon("paper-plane"),
                style="color: #fff; background-color: #D73925; border-color: #DB4836")
        # )
      )),


      fluidRow(
        box(
          title = 'STEP 4: Abundance normalization', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
          radioButtons(
            ns("norm_method"),
            label = "Normalization : ",
            inline = TRUE,
            choices = list(
              "Raw" = 0 ,
              "TSS (total-sum normalization)" = 1,
              "CLR (centered log-ratio)" = 2,
              "VST (variance stabilizing transformation)" = 3
            ), selected = 1
          ),
          shinyBS::bsButton(inputId = ns('norm'), label = "Normalize", block = F, style = 'danger', type='action')
          # actionButton(ns('norm'), "Normalize", class='butt2')
        ),
        box(
          title = 'Final phyloseq object', status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
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
      ),
      fluidRow(
        box(title = "Color selection", width = 12, status = "warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            infoBox("",
                    HTML(paste("Accepted formats for color assignation.")),
                    HTML(paste("When the variable is a <b>factor</b>, modalities must be written in the first column of the Excel file and colors associated in hexadecimal code form in the second column.<br/> When the variable is <b>numeric</b>, a palette must be chosen with the form package::palette written in the Excel file (only packages ggthemes and viridis are accepted).")),
                    icon = icon("info-circle"), fill = FALSE, width = 12
            ),
            uiOutput(ns("var_color")),
            downloadButton(ns("var_color_download"), label = "Download xlsx file"),
            fileInput(ns("file_var_color"), label = "var-color file", placeholder = "var_colors.xlsx"),
            verbatimTextOutput(ns("fact_color_file_log")),
            uiOutput(ns("taxa_color")),
            downloadButton(ns("taxa_color_download"), label = "Download xlsx file"),
            fileInput(ns("file_taxa_color"), label = "taxa-color file", placeholder = "taxa_colors.xlsx"),
            verbatimTextOutput(ns("taxa_color_file_log"))
        )
      )
    )
  )
}


merge_table <- function(rank, table){
  FNGdata <- table
  rnames <- phyloseq::rank_names(FNGdata)
  if(rank=="ASV"){
    rank1 = rnames[length(rnames)]
  }
  else{
    rank1 = rank
  }
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


  rawtaxasum1 <-  table %>%
    taxa_sums() %>%
    as.data.frame %>%
    tibble::rownames_to_column()
  names(rawtaxasum1)[2] <- "RawAbundanceSum"

  joinGlom <-
    dplyr::left_join(ttable, rawtaxasum1, by = "rowname") %>%
    mutate(RawFreq = RawAbundanceSum / sum(RawAbundanceSum)) %>%
    dplyr::left_join(otable, by = "rowname")

  if(rank=="ASV" & !is.null(refseq(table, errorIfNULL=FALSE)) ){
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
    showNotification("No refseq in object.", type="error", duration = 3)
    dplyr::rename(joinGlom, asvname = rowname)
    FTAB = as.data.frame(joinGlom)
  }
  return(FTAB)
}





#' data_loading Server Function
#'
#' @noRd
mod_data_loading_server <- function(input, output, session, r=r){
  ns <- session$ns
  r_values <- reactiveValues(phyobj_initial=NULL, phyobj_sub_samples=NULL, phyobj_norm=NULL, phyobj_taxglom=NULL, phyobj_final=NULL, phyobj_tmp=NULL)

  ###Filtering metadata

  res_filter <- datamods::filter_data_server(
    id = "filtering",
    # data = data,
    data = reactive({
      req(sdat_initial())
      if(is.null(updated_data())){
        sdat_initial()
      }else{
      req(updated_data())
        updated_data()
      }
    }),
    name = reactive("feature_table"),
    vars = reactive(NULL),
    widget_num = "slider",
    widget_date = "slider",
    label_na = "Missing"
  )

  output$metadata_table <- DT::renderDT({
    res_filter$filtered()
  }, options = list(pageLength = 10, scrollX = TRUE))

  ## update tab
  updated_data <- datamods::update_variables_server(
    id = "vars",
    data = reactive({
        req(sdat_initial())
        sdat_initial()  #data()
    })
  )


  phyloseq_data <- reactive({
    ne <- new.env()
    if (!is.null(input$fileRData)){
      load(input$fileRData$datapath, envir = ne)
    }
    else{
      load(system.file("data_test", "robjects.Rdata", package="ExploreMetabar"), envir = ne)
      #load(system.file("data_test", "phy_test_numeric.rdata", package="ExploreMetabar"), envir = ne)
    }
    classes1 = sapply(ne, class)
    obj = classes1[classes1 == "phyloseq"]
    fun = glue::glue("r_values$phyobj_initial <- ne${names(obj)}")
    eval(parse(text = fun))
    if(is.null(refseq(r_values$phyobj_initial, errorIfNULL=FALSE)) ){
      showNotification("No refseq in object.", type="error", duration = 3)
    }
    r_values$phyobj_tmp <- r_values$phyobj_initial
    return(r_values$phyobj_initial)
  })

  output$phy_prev <- renderPrint({
    cat(file=stderr(), 'rendering phy_prev', "\n")
    cat('Running ExploreMetabar v2.0.1\n')
    phyloseq_data()
  })

  sdat_initial <- reactive({
    req(r_values$phyobj_initial)
    phyobj <- r_values$phyobj_initial
    sdat <- do.call(cbind.data.frame, phyobj@sam_data)
    sdat[sdat == ""] <- NA
    
    if( !"sample.id" %in% colnames(sdat) ){
      sdat <- sdat %>% dplyr::mutate(sample.id = sample_names(phyobj), .before = 1)
    }
    sdat <- sdat %>% select(where(~ n_distinct(.) > 1))
    
    for(i in seq(1, 6, 2)){
      if(!is.null(input[[paste("combin", i+1, sep = "_")]])){
        sdat <- combine_columns_sdat(metadata = sdat, list_columns = c(input[[paste("combin", i, sep = "_")]], input[[paste("combin", i+1, sep = "_")]]))
      }
    }
    return(sdat)
  })
  
  combine_columns_sdat <- function(metadata, list_columns){
    if(!all(sapply(metadata[, list_columns], is.numeric))){
      meta_col <- paste0(list_columns, collapse = '_')
      metadata <- tidyr::unite(metadata, !!meta_col, list_columns, na.rm = TRUE, remove = FALSE)
      metadata[, meta_col] <- as.factor(metadata[, meta_col])
    }
    return(metadata)
  }


  subset_samples <- reactive({
    req(r_values$phyobj_initial, res_filter$filtered)
    filt_sdata <- res_filter$filtered()
    physeq0 <- r_values$phyobj_initial

    flog.info('subset samples() starting...')
    flog.info(paste0('number of samples before ', phyloseq::nsamples(physeq0)))
    if(!any(sample_names(physeq0) %in% as.vector(filt_sdata$sample.id))){
      print("No match")
      shinyalert(title = "Oops", text="Sample names do not match with names in metadata. Check the phyloseq object.", type='error')
      return()
    }
    physeq <- phyloseq::prune_samples(as.vector(filt_sdata$sample.id),physeq0)
    physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq)>0, physeq)

    flog.info(paste0('number of samples after',phyloseq::nsamples(physeq)))
    rownames(filt_sdata) <- filt_sdata  %>% pull('sample.id')
    sample_data(physeq) <- filt_sdata
    r_values$phyobj_sub_samples <- r_values$phyobj_tmp <- physeq

  })

  #update button color when clicked
  observeEvent(input$update_metadata,{
    shinyBS::updateButton(session = session, ns('update_metadata'), block = F, style = 'success')
    shinyBS::updateButton(session = session, ns('update_taxo'), block = F, style = 'danger')
    shinyBS::updateButton(session = session, ns('subset_taxo'), block = F, style = 'danger')
    shinyBS::updateButton(session = session, ns('norm'), block = F, style = 'danger')
  })
  observeEvent(input$update_taxo,{
    shinyBS::updateButton(session = session, ns('update_taxo'), block = F, style = 'success')
    shinyBS::updateButton(session = session, ns('subset_taxo'), block = F, style = 'danger')
    shinyBS::updateButton(session = session, ns('norm'), block = F, style = 'danger')
  })
  observeEvent(input$subset_taxo,{
    shinyBS::updateButton(session = session, ns('subset_taxo'), block = F, style = 'success')
    shinyBS::updateButton(session = session, ns('norm'), block = F, style = 'danger')
  })
  observeEvent(input$norm,{
    shinyBS::updateButton(session = session, ns('norm'), block = F, style = 'success')
  })


  observeEvent(input$launch_all, {
    subset_samples()
    glom_taxo0()
    glom_taxo()
    subset_taxa()
    normalize()
  })

  observeEvent(input$update_metadata, {
    flog.info('button update_metadata')
    subset_samples()
  },
  ignoreNULL = TRUE, ignoreInit = TRUE)


  observe({
    flog.info('updating rank_glom selectInput...')
    updateSelectInput(session, "rank_glom",
                      choices = c( rank_names(phyloseq_data()), "ASV" ),
                      selected = "ASV")
  }) #updateSelectInput


  observe({
    flog.info('updating minAb numericInput...')
    updateNumericRangeInput(session, 'minAb',"Minimum taxa overall raw abundance:", value=c(1,max(taxa_sums(r_values$phyobj_tmp))))
  }) #updateNumericRangeInput


  observe({
    flog.info('updating minPrev numericInput...')
    updateNumericRangeInput(session, 'minPrev',"Minimum taxa prevalence in samples:", value=c(1,max(nsamples(r_values$phyobj_tmp))))
  }) #updateNumericRangeInput


  glom_taxo0 <- reactive({
    req(input$minAb, input$minPrev, input$rank_glom, r_values$phyobj_sub_samples)
    flog.info('filter_taxonomy...')
    tmp <- r_values$phyobj_sub_samples
    withProgress({
      if(input$rank_glom != 'ASV'){
        if(nsamples(tmp)>1000){
          showNotification("Phylogentic tree removed, too much samples...", type="message", duration = 5)
          tmp <- fast_tax_glom(tmp, input$rank_glom)
        }else{
          tmp <- tax_glom(tmp, input$rank_glom)

        }
        FGnames <- tax_table(tmp)[,input$rank_glom]
        nnames <- paste(substr(FGnames, 1, 50), taxa_names(tmp), sep="_")
        taxa_names(tmp) <- nnames
      }
    showNotification("Taxonomy agglomeration done...", type="message", duration = 1)
    }, message = 'Processing, please wait.')
    r_values$phyobj_taxglom0 <- r_values$phyobj_tmp <- tmp
    flog.info('done.')
  })


  glom_taxo <- reactive({
    req(input$minAb, input$minPrev, input$rank_glom, r_values$phyobj_sub_samples)
      tmp <- r_values$phyobj_taxglom0

      tmp <- prune_taxa(taxa_sums(tmp) >= input$minAb[1], tmp)

      tmp <- prune_taxa(taxa_sums(tmp) <= input$minAb[2], tmp)

      prevdf <- apply(X = otu_table(tmp), MARGIN = ifelse(taxa_are_rows(tmp), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
      taxToKeep1 <- names(prevdf)[(prevdf >= input$minPrev[1] & prevdf <= input$minPrev[2])]
      tmp <- prune_taxa(taxToKeep1, tmp)
      if(input$rank_glom != 'ASV'){
        tax_table(tmp) <- tax_table(tmp)[,1:match(input$rank_glom, rank_names(tmp))]
      }

      flog.info('glom object')

      r_values$phyobj_taxglom <- r_values$phyobj_tmp <- tmp

      flog.info('filter_taxonomy done.')
      showNotification("Filter taxonomy done...", type="message", duration = 1)
  })


  observeEvent(input$update_taxo0, {
    glom_taxo0()
  },ignoreInit = TRUE)


  observeEvent(input$update_taxo, {
    glom_taxo()
  },ignoreInit = TRUE)


  render_taxonomy_table <- reactive({
    withProgress({

      req(r_values$phyobj_tmp, input$rank_glom)
      flog.info('render_taxonomy_table fun')

      phyloseq_obj <- r_values$phyobj_tmp
      rnames <- phyloseq::rank_names(phyloseq_obj)
      if(input$rank_glom=="ASV"){
        rank1 = rnames[length(rnames)]
      }
      else{
        rank1 = input$rank_glom
      }
      ttable <- phyloseq_obj %>%
      tax_table() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      dplyr::select(1:rank1) %>%
      tibble::rownames_to_column() %>%
      as.matrix() %>% as.data.frame(stringsAsFactors = TRUE)

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
        showNotification("Sequences added to dataframe.", type="message", duration = 5)
        refseq1 <- phyloseq_obj %>%
        refseq %>%
        as.data.frame %>%
        tibble::rownames_to_column() %>%
        rename(sequences = x)

        joinGlom2 <- dplyr::left_join(joinGlom, refseq1, by = "rowname") %>%
        dplyr::rename(asvname = rowname)
        FTAB = as.data.frame(joinGlom2, stringsAsFactors = TRUE)
      }else{
        dplyr::rename(joinGlom, asvname = rowname)
        FTAB = as.data.frame(joinGlom, stringsAsFactors = TRUE)
      }
      flog.info('render_taxonomy_table done.')
      showNotification("Render taxonomy table ...", type="message", duration = 1)
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
      req(r_values$phyobj_taxglom, res_filter_taxo$filtered)
      flog.info('subset_taxa fun')
      filttax <- res_filter_taxo$filtered()
      selected <- filttax[,1]

      # selected <- render_taxonomy_table()[input$taxonomy_table_rows_all, 1]
      phy_obj <- prune_taxa(selected, r_values$phyobj_taxglom)
      r_values$phyobj_final <- phy_obj
      r_values$phyobj_tmp <- phy_obj
      flog.info('subset_taxa fun done.')
      # phy_obj

    }, message = "Subset taxonomy, please wait...")
  })


  observeEvent(input$subset_taxo, {
    flog.info('button subset_taxo')
    subset_taxa()
  },ignoreNULL = TRUE, ignoreInit = TRUE)

  ## Filter taxo

    res_filter_taxo <- datamods::filter_data_server(
      id = "filtering_taxo",
      data = reactive({
        req(render_taxonomy_table())
        render_taxonomy_table()
      }),
      name = reactive("tax_table"),
      vars = reactive({
        req(render_taxonomy_table())
        s_names <- phyloseq::sample_names(r_values$phyobj_tmp)
        col_names <- colnames(render_taxonomy_table())
        filt <- dplyr::setdiff(col_names, s_names)
      return(filt)
      }),
      widget_num = "slider",
      widget_date = "slider",
      label_na = "Missing"
    )

    output$table_taxoFILT <- DT::renderDT({
      res_filter_taxo$filtered()
    }, options = list(pageLength = 10, scrollX = TRUE))




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
    showNotification("Dataset ready !", type="message", duration = 5)
    r_values$phyobj_norm <- FNGdata
  })

  observeEvent(input$norm, {
    flog.info('button normalize')
    normalize()
  },ignoreNULL = TRUE, ignoreInit = TRUE)


  output$phy_after <- renderPrint({
    print(r_values$phyobj_tmp)
  })

  output$phy_norm <- renderPrint({
    req(r_values$phyobj_norm)
    print(r_values$phyobj_norm)
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
    # r_values$phyobj_initial #dev
  })


  # final filtered object normalize
  r$phyloseq_filtered_norm <- reactive({
    req(r_values$phyobj_norm)
    r_values$phyobj_norm
    # r_values$phyobj_initial #dev
  })

  r$norm_method <- reactive({
    input$norm_method
  })

  # Chosen rank to glom taxa
  r$rank_glom <- reactive({
    input$rank_glom
  })

  # Export metadata
  r$sdat <- reactive({
    req(r_values$phyobj_final)
    sdat <- sample_data(r_values$phyobj_final)
    sdat <- sdat[,which(unlist(lapply(sdat, function(x)!all(is.na(x))))),with=F]
    sdat <- as(sdat, "data.frame")
    if(! 'sample.id' %in% colnames(sdat)){
      sdat$sample.id <- rownames(sdat)
    }
    return(sdat)
  })

  r$var_list <- reactive({
    req(r_values$phyobj_final, r$sdat)
    sdat <- r$sdat()
    var_list <- colnames(sdat)
    if('sample.id' %in% var_list){
      var_list <- sort(var_list[! var_list %in% 'sample.id'])
    }
    return(var_list)
  })
  
  # combine variables of sample_data
  choice_combine_var <- reactive({
    sdat <- sample_data(r_values$phyobj_initial)
    sdat <- do.call(cbind.data.frame, sdat)
    sdat <- select(sdat, where(~ n_distinct(.) > 1))
    choices <- colnames(sdat)
    choices <- choices[choices != 'sample.id']
    return(choices)
  })
 
  output$ui_combine_columns <- renderUI({
    box(title = "Combine columns", width = 12, status = "warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
        column(width = 6,
               shinyWidgets::pickerInput(
                 ns("combin_1"),
                 label = "Variable to combine",
                 choices = choice_combine_var(),
                 selected = NULL,
                 multiple = TRUE,
                 options = pickerOptions(maxOptions = 1)
               )),
        lapply(2:8, FUN = function(i){
          column(width = 6, uiOutput(ns(paste0("ui_combin_", i))))
        })
    )
  })
  
  output$ui_combin_2 <- renderUI({
    choice <- choice_combine_var()
    shinyWidgets::pickerInput(
      ns("combin_2"),
      label = "Variable to combine",
      choices = choice[choice != input$combin_1],
      selected = NULL,
      multiple = TRUE,
      options = pickerOptions(maxOptions = 1)
    )
  })
  
  get_ui_combin <- function(i){
    choice <- choice_combine_var()
    if(i%%2 == 0){
      prec <- i - 2
      choice <- choice[choice != input[[paste0("combin_", i - 1)]]]
    }else{
      prec <- i - 1
    }
    output[[paste0("ui_combin_", i)]] <- renderUI({
      if(!is.null(input[[paste0("combin_", prec)]])){
        shinyWidgets::pickerInput(
          ns(paste0("combin_", i)),
          label = "Variable to combine",
          choices = choice,
          selected = NULL,
          multiple = TRUE,
          options = pickerOptions(maxOptions = 1)
        )
      }
    })
  }
  
  observe({
    lapply(3:6, FUN = function(i){
      if(is.null(input[[paste0("combin_", i)]])){
        get_ui_combin(i)
      }
    })
  })
  
  
  # for numeric variable colors
  palette_numeric_var <- reactive({
    pal <- paletteer::palettes_c_names
    pal <- pal[pal$package %in% c("viridis", "ggthemes"), ]
    return(pal)
  })
  
  
  # for variable colors
  output$var_color <- renderUI({
    req(r$var_list())
    shinyWidgets::pickerInput(ns("list_var_color"),
                              label = "Variables whose associated colors are exported in xlsx file",
                              choices = r$var_list(),
                              selected = r$var_list()[1],
                              multiple = TRUE,
                              options = pickerOptions(
                                actionsBox = TRUE,
                                liveSearch = TRUE,
                                showContent = FALSE
                              ),
                              choicesOpt = list(
                                content = unlist(lapply(
                                  X = r$var_list(),
                                  FUN = function(x){
                                    htmltools::doRenderTags(
                                      tags$div(
                                        splitLayout(cellWidths = 200,
                                                    tags$div(
                                                      style = htmltools::css(fontWeight = "bold"),
                                                      x
                                                    ),
                                                    tags$div(
                                                      style = htmltools::css(color = 'grey'),
                                                      class(r$sdat()[, x])
                                                    )
                                        )
                                      )
                                    )
                                  }
                                ))
                              )
    )
  })
  
  r$factor_list <- reactive({
    req(r$sdat(), r$var_list())
    num <- sapply(r$sdat()[, r$var_list()], is.numeric)
    fact <- r$var_list()[!num]
    return(fact)
  })
  
  modality <- reactive({
    req(r$phyloseq_filtered(), r$factor_list())
    mod <- get_cat_colors(type = "fact", list_fact = r$factor_list(), phy_object = r$phyloseq_filtered())
    return(mod)
  })
  
  palette_choices <- reactive({
    pal <- paletteer::palettes_d_names
    pal <- pal[pal$type == "qualitative", ]
    pal <- pal[pal$package %in% c("ggsci", "ggthemes", "jcolors", "pals", "Polychrome", "RColorBrewer"), ]
    pal <- pal[pal$package != "ggthemes" | ! pal$palette %in% c("fivethirtyeight", "Seattle_Grays", "Classic_Gray_5", "excel_Grayscale", "stata_mono", "stata_economist"), ]
    pal <- pal[pal$package != "Polychrome" | ! pal$palette %in% c("glasbey", "kelly"), ]
    pal <- pal[pal$package != "pals" | pal$palette != "kelly", ] # remove palettes with white or grey
    return(pal)
  })
  
  modality_colors <- reactive({
    req(r$var_list(), r$factor_list(), modality())
    colors <- get_modality_colors(type = "fact", list_fact = r$factor_list(), modality = modality(), palettes = palette_choices())
    num <- dplyr::setdiff(r$var_list(), r$factor_list())
    if(length(num) > 0){
      available_pal <- palette_numeric_var()
      pal <- lapply(1:length(num), FUN = function(i){
        paste(available_pal$package[i], available_pal$palette[i], sep = "::")
      })
      names(pal) <- num
      colors <- c(colors, pal)
    }
    return(colors)
  })
  
  modality_file <- reactive({
    req(r$var_list(), input$file_var_color)
    num <- sapply(r$sdat()[, r$var_list()], is.numeric)
    palettes <- palette_numeric_var()
    pal <- sapply(1:dim(palettes)[1], FUN = function(i){
      paste(palettes$package[i], palettes$palette[i], sep = "::")
    })
    colors <- get_modality_file(path_file = input$file_var_color$datapath, type = "fact", list_fact = r$var_list(), num_fact = num, num_palettes = pal, phy_obj = r$phyloseq_filtered())
    return(colors)
  })
  
  complete_modality_file <- reactive({
    req(modality_file(), missing_fact_colors())
    color <- get_complete_color_file(mod_file = modality_file(), missing_colors = missing_fact_colors())
    return(color)
  })
  
  missing_fact_colors <- reactive({
    req(modality_colors(), modality_file())
    fact_loaded <- names(modality_file())
    missing_fact <- lapply(1:length(fact_loaded), FUN = function(i){
      miss_fact <- dplyr::setdiff(names(modality_colors()[[fact_loaded[i]]]), names(modality_file()[[fact_loaded[i]]]))
      if(length(miss_fact) == 0){
        miss_fact <- c("No color is missing.")
      }
      return(miss_fact)
    })
    names(missing_fact) <- fact_loaded
    return(missing_fact)
  })
  
  observe({
    req(modality_file())
    length_modality_file <- length(modality_file())
    names_modality_file <- names(modality_file())
    output$fact_color_file_log <- renderPrint({
      if(length_modality_file > 0){
        cat(names(modality_file()), "have been loaded.\n", sep = " ")
        cat("Warning ! If some modalities are missing, colors associated to these modalities will be randomly chosen.\n")
        cat("Missing factor-color values for each factor loaded :\n")
        missing_fact_colors()
      }else{
        cat("No colors have been loaded.")
      }
    })
  })
  
  output$var_color_download <- downloadHandler(
    filename = "var_colors.xlsx",
    content = function(file){
      if(length(input$list_var_color) < 50){
        download_color <- r$factor_colors()[input$list_var_color]
        num <- sapply(r$sdat()[, input$list_var_color], is.numeric)
        for(i in 1:length(input$list_var_color)){
          xlsx::write.xlsx(data.frame(download_color[[i]]), file = file, sheetName = input$list_var_color[i], col.names = FALSE, row.names = !num[i] , append = TRUE)
        }
      }else{
        showNotification("Too much variables to export.", type = "warning", duration = 10)
      }
    }
  )
  
  r$factor_colors <- reactive({
    req(modality_colors())
    if(is.null(input$file_var_color)){
      factors <- modality_colors()
    }else{
      factors <- c(complete_modality_file(), modality_colors()[setdiff(r$var_list(), names(complete_modality_file()))])
    }
    return(factors)
  })
  
  
  # for taxa colors
  output$taxa_color <- renderUI({
    req(rank_list())
    shinyWidgets::pickerInput(ns("list_taxa_color"),
                              label = "Taxonomic ranks whose associated colors are exported in xlsx file",
                              choices = rank_list(),
                              selected = rank_list()[1],
                              multiple = TRUE,
                              options = pickerOptions(
                                actionsBox = TRUE,
                                liveSearch = TRUE,
                                showContent = FALSE
                              )
    )
  })
  
  rank_list <- reactive({
    phyloseq::rank_names(r$phyloseq_filtered())
  })
  
  taxa_modality <- reactive({
    req(rank_list(), r$var_list())
    mod <- get_cat_colors(type = "taxa", list_fact = rank_list(), phy_object = r$phyloseq_filtered())
    return(mod)
  })
  
  observe({
    req(input$file_taxa_color)
    taxa_modality_file()
    })
  
  taxa_modality_colors <- reactive({
    req(rank_list(), taxa_modality())
    colors <- get_modality_colors(type = "taxa", list_fact = rank_list(), modality = taxa_modality(), palettes = palette_choices())
    return(colors)
  })
  
  taxa_modality_file <- reactive({
    req(rank_list(), input$file_taxa_color)
    colors <- get_modality_file(path_file = input$file_taxa_color$datapath, type = "taxa", list_fact = rank_list(), phy_obj = r$phyloseq_filtered())
    return(colors)
  })
  
  missing_taxa_colors <- reactive({
    req(taxa_modality_colors(), taxa_modality_file())
    rank_loaded <- names(taxa_modality_file())
    missing_taxa <- lapply(1:length(rank_loaded), FUN = function(i){
      miss_taxa <- dplyr::setdiff(names(taxa_modality_colors()[[rank_loaded[i]]]), names(taxa_modality_file()[[rank_loaded[i]]]))
      return(miss_taxa)
    })
    names(missing_taxa) <- rank_loaded
    return(missing_taxa)
  })
  
  complete_taxa_file <- reactive({
    req(taxa_modality_file(), missing_taxa_colors())
    color <- get_complete_color_file(mod_file = taxa_modality_file(), missing_colors = missing_taxa_colors())
    return(color)
  })
  
  nb_missing_taxa_colors <- reactive({
    req(taxa_modality_colors(), taxa_modality_file())
    rank_loaded <- names(taxa_modality_file())
    missing_taxa <- sapply(1:length(rank_loaded), FUN = function(i){
      miss_taxa <- dplyr::setdiff(names(taxa_modality_colors()[[rank_loaded[i]]]), names(taxa_modality_file()[[rank_loaded[i]]]))
      return(length(miss_taxa))
    })
    names(missing_taxa) <- rank_loaded
    return(missing_taxa)
  })
  
  output$taxa_color_file_log <- renderPrint({
    if(length(taxa_modality_file()) > 0){
      cat(names(taxa_modality_file()), "have been loaded.\n", sep = " ")
      cat("Warning ! If some taxa are missing, colors associated to these taxa will be randomly chosen.\n")
      cat("Number of missing taxa-color values for each rank loaded :\n")
      nb_missing_taxa_colors()
    }else{
      cat("No colors have been loaded.")
    }
  })
  
  output$taxa_color_download <- downloadHandler(
    filename = "taxa_colors.xlsx",
    content = function(file){
      download_color <- taxa_modality_colors()[input$list_taxa_color]
      for(i in 1:length(input$list_taxa_color)){
        xlsx::write.xlsx(data.frame(download_color[[i]]), file = file, sheetName = input$list_taxa_color[i], col.names = FALSE, append = TRUE)
      }
    }
  )
  
  r$taxa_colors <- reactive({
    req(taxa_modality_colors())
    if(is.null(input$file_taxa_color)){
      colors <- taxa_modality_colors()
    }else{
      colors <- c(complete_taxa_file(), taxa_modality_colors()[dplyr::setdiff(rank_list(), names(complete_taxa_file()))])
    }
    return(colors)
  })
  
}

## To be copied in the UI
# mod_data_loading_ui("data_loading_ui_1")

## To be copied in the server
# callModule(mod_data_loading_server, "data_loading_ui_1")
