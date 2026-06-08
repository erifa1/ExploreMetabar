#' data_loading UI Function
#'
#' @description Module for loading phyloseq object from rdata file. This module allows to filter and select samples and taxa prior to analysis.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import DT
#' @importFrom Biostrings writeXStringSet
#' @importFrom glue glue
#' @importFrom futile.logger flog.info flog.debug
#' @import datamods
#' @import phyloseq
#' @import bslib
#' @importFrom bsicons bs_icon
#'
mod_data_loading_ui <- function(id){
  ns <- NS(id)
  tagList(
    layout_sidebar(
      fillable = TRUE,
      sidebar = sidebar(
        title = "Settings",
        width = "350px",
        open = "desktop",
        accordion(
          id = ns("config_accordion"),
          open = "Data input",
          multiple = TRUE,

          accordion_panel(
            "Data input",
            icon = bs_icon("upload"),
            tooltip(
              fileInput(ns("fileRData"),
                        label = "RData with phyloseq object :",
                        placeholder = "data.RData"),
              "RData where 'data' is a phyloseq object."
            )
          ),

          accordion_panel(
            "Combine columns",
            icon = bs_icon("columns-gap"),
            uiOutput(ns("ui_combine_columns"))
          ),

          accordion_panel(
            "Sample filtering",
            icon = bs_icon("funnel"),
            datamods::filter_data_ui(ns("filtering"), max_height = "400px")
          ),

          accordion_panel(
            "Taxonomy rank",
            icon = bs_icon("diagram-3"),
            selectInput(
              ns("rank_glom"),
              label = "Rank to merge taxonomy",
              choices = "",
              selected = 1
            ),
            layout_columns(
              col_widths = c(6, 6),
              autonumericInput(ns("minAb"), "Min. abundance:", value = 0, width = NULL, decimalPlaces = 6),
              autonumericInput(ns("minPrev"), "Min. prevalence:", value = 0, width = NULL, decimalPlaces = 6)
            )
          ),

          accordion_panel(
            "Taxa filtering",
            icon = bs_icon("filter"),
            datamods::filter_data_ui(ns("filtering_taxo"), max_height = "400px"),
            tooltip(
              actionButton(
                ns("apply_taxa_filter"),
                label = "Apply taxa filter",
                icon = bs_icon("check2-circle"),
                class = "btn-outline-secondary w-100 mt-2"
              ),
              "Optional. Process Data keeps all taxa; click to restrict the dataset to the taxa selected above."
            )
          ),

          accordion_panel(
            "Normalization",
            icon = bs_icon("sliders"),
            radioButtons(
              ns("norm_method"),
              label = "Method :",
              choices = list(
                "Raw" = 0,
                "TSS" = 1,
                "CLR" = 2,
                "VST" = 3
              ), selected = 1
            )
          )
        ),
        tags$hr(),
        input_task_button(
          ns("process_data"),
          label = "Process Data",
          icon = bs_icon("play-fill"),
          class = "btn-danger w-100",
          label_busy = "Processing..."
        )
      ),

      # ── Main content area ──
      card(
        fill = FALSE,
        card_header("Phyloseq preview", class = "bg-primary"),
        card_body(fillable = FALSE, verbatimTextOutput(ns("phy_prev")))
      ),

      navset_card_underline(
        id = ns("main_tabs"),

        nav_panel(
          title = "Metadata",
          icon = bs_icon("table"),
          navset_underline(
            nav_panel(
              "Table",
              DT::dataTableOutput(outputId = ns("metadata_table"))
            ),
            nav_panel(
              "Update variables",
              datamods::update_variables_ui(ns("vars"))
            )
          )
        ),

        nav_panel(
          title = "Taxonomy",
          icon = bs_icon("diagram-3"),
          DT::dataTableOutput(outputId = ns("table_taxoFILT"))
        ),

        nav_panel(
          title = "Downloads",
          icon = bs_icon("download"),
          layout_column_wrap(
            width = 1/2,
            heights_equal = "row",
            card(
              fill = FALSE,
              card_header("RAW tables"),
              card_body(
                fillable = FALSE,
                downloadButton(outputId = ns("raw_otable_download"), label = "Download raw ASV table"),
                downloadButton(outputId = ns("raw_refseq_download"), label = "Download raw FASTA sequences")
              )
            ),
            card(
              fill = FALSE,
              card_header("Filtered tables"),
              card_body(
                fillable = FALSE,
                downloadButton(outputId = ns("filt_otable_download"), label = "Download filtered ASV table"),
                downloadButton(outputId = ns("filt_norm_otable_download"), label = "Download filtered and normalized ASV table"),
                downloadButton(outputId = ns("filt_rdata_download"), label = "Download filtered Phyloseq object"),
                downloadButton(outputId = ns("filt_rdata_norm_download"), label = "Download filtered and normalized Phyloseq object"),
                downloadButton(outputId = ns("filt_refseq_download"), label = "Download filtered FASTA sequences")
              )
            )
          )
        ),

        nav_panel(
          title = "Colors",
          icon = bs_icon("palette"),
          mod_color_ui(ns("color_ui_1"))
        ),

        nav_panel(
          title = "Pipeline",
          icon = bs_icon("terminal"),
          layout_column_wrap(
            width = 1/3,
            heights_equal = "row",
            card(
              fill = FALSE,
              card_header("Filtered object"),
              card_body(fillable = FALSE, verbatimTextOutput(ns("phy_after")))
            ),
            card(
              fill = FALSE,
              card_header("Normalized object"),
              card_body(fillable = FALSE, verbatimTextOutput(ns("phy_norm")))
            ),
            card(
              fill = FALSE,
              card_header("Pipeline log"),
              card_body(fillable = FALSE, verbatimTextOutput(ns("pipeline_log")))
            )
          )
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
    phyloseq::tax_table() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::select(1:rank1) %>%
    tibble::rownames_to_column() %>%
    as.matrix() %>% as.data.frame()


  otable <- FNGdata %>%
    phyloseq::otu_table() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column()


  rawtaxasum1 <-  table %>%
    phyloseq::taxa_sums() %>%
    as.data.frame %>%
    tibble::rownames_to_column()
  names(rawtaxasum1)[2] <- "RawAbundanceSum"

  joinGlom <-
    dplyr::left_join(ttable, rawtaxasum1, by = "rowname") %>%
    dplyr::mutate(RawFreq = RawAbundanceSum / sum(RawAbundanceSum)) %>%
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


#' Correct sample-metadata column types
#'
#' Best-effort type repair for a sample-metadata data frame. Character/factor
#' columns whose values are really *decimal* numbers — including French
#' decimal-comma notation (`"3,5"`) and (non-breaking) space thousand
#' separators (`"1 000,5"`) — are converted to numeric; the remaining
#' non-numeric character columns are declared factors.
#'
#' Integer-only text columns (e.g. coded grouping variables `"1","2","3"`) are
#' deliberately **kept categorical** (turned into factors): conversion to
#' numeric only happens when at least one value carries a decimal separator, so
#' a coded group is never silently turned into a continuous variable. Columns
#' already numeric/logical, and `exclude_cols` (sample identifiers), are left
#' untouched.
#'
#' @param metadata A data frame (typically `phyloseq::sample_data` coerced to a
#'   data frame).
#' @param exclude_cols Column names never coerced (sample identifiers).
#'
#' @return `metadata` with corrected column types, carrying a `"coercions"`
#'   attribute: a character vector describing each change (length 0 if nothing
#'   changed).
#' @noRd
coerce_metadata_types <- function(metadata, exclude_cols = "sample.id") {
  changes <- character(0)
  cols <- setdiff(names(metadata), exclude_cols)

  for (col in cols) {
    x <- metadata[[col]]
    if (is.numeric(x) || is.logical(x)) next     # already typed, leave as-is

    xc <- trimws(as.character(x))
    xc[xc == ""] <- NA
    non_na <- !is.na(xc)
    if (!any(non_na)) next                        # nothing to infer from

    # French-locale numeric candidate: drop (thousand) spaces, comma -> dot.
    cand <- gsub("[[:space:]]", "", xc)
    cand <- gsub("[   ]", "", cand)   # nbsp / narrow / thin spaces
    cand <- gsub(",", ".", cand, fixed = TRUE)
    num <- suppressWarnings(as.numeric(cand))

    all_numeric <- all(!is.na(num[non_na]))
    has_decimal <- any(grepl("[.,]", xc[non_na]))

    if (all_numeric && has_decimal) {
      metadata[[col]] <- num
      changes <- c(changes, sprintf("'%s': text -> numeric (decimal)", col))
    } else if (is.character(x)) {
      metadata[[col]] <- as.factor(x)
      changes <- c(changes, sprintf("'%s': character -> factor", col))
    }
  }

  attr(metadata, "coercions") <- changes
  metadata
}




#' data_loading Server Function
#'
#' @noRd
#' @import metagMisc
#' @import phyloseq
#' @import dplyr
#' @import tibble
#' @import speedyseq
#' 
mod_data_loading_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns
  r_values <- reactiveValues(phyobj_initial = NULL)

  # `final` holds the user-facing filtered object: processed() optionally refined
  # by the manual taxa filter. It is a reactiveVal — never a stale one-shot
  # snapshot — so the taxa filter (calibrated on the post-subset table) can update
  # it after the browser round-trip without re-running the heavy pipeline.
  final <- reactiveVal(NULL)

  ###Filtering metadata

  res_filter <- datamods::filter_data_server(
    id = "filtering",
    # data = data,
    data = reactive({
      req(sdat_initial())
      if(is.null(updated_data())){
        d1 <- sdat_initial()
        if( any(sapply(d1, is.factor))){
          d1$sample.id <- as.factor(sample_names(phyloseq_data()))
          d1
        }else{
          # add dummy factor
          showNotification("Dummy factor added to display metadata.", type="warning", duration = 5)
          d1$sample.id <- as.factor(sample_names(phyloseq_data()))
          d1$dummy_fact <- as.factor(rep(LETTERS, each = 2, length.out = nrow(d1)))
          d1
        }
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
    r_values$phyobj_initial <- ne[[names(obj)[1]]]
    if(is.null(refseq(r_values$phyobj_initial, errorIfNULL=FALSE)) ){
      showNotification("No refseq in object.", type="error", duration = 3)
    }
    return(r_values$phyobj_initial)
  })

  output$phy_prev <- renderPrint({
    flog.info('rendering phy_prev')
    cat(paste0('Running ExploreMetabar v', as.character(utils::packageVersion("ExploreMetabar")), '\n'))
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

    # Repair column types: French decimal-comma numerics -> numeric, remaining
    # categorical character columns -> factor (sample.id is left untouched).
    sdat <- coerce_metadata_types(sdat, exclude_cols = "sample.id")
    coercions <- attr(sdat, "coercions")
    if (length(coercions)) {
      flog.info(paste0("Metadata type coercion: ", paste(coercions, collapse = "; ")))
      showNotification(
        ui = tagList(
          tags$b("Metadata column types corrected:"),
          tags$ul(lapply(coercions, tags$li))
        ),
        id = "metadata_type_coercion",
        type = "warning",
        duration = 8
      )
    }

    for(i in seq(1, 3, 2)){
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


  # Pipeline log
  pipeline_log <- reactiveVal("")
  log_msg <- function(msg) {
    pipeline_log(paste0(pipeline_log(), msg, "\n"))
  }

  output$pipeline_log <- renderPrint({
    cat(pipeline_log())
  })

  # ── Heavy pipeline, gated behind the "Process Data" button ──
  # subset samples -> drop empty taxa -> agglomerate -> abundance/prevalence filter.
  # Returns a single phyloseq object. The manual taxa filter is applied separately
  # (via `final`), so this stage never reads the lagging datamods taxa-filter
  # widgets — which is what previously dropped taxa using stale, full-dataset
  # slider thresholds.
  processed <- eventReactive(input$process_data, {
    req(r_values$phyobj_initial, res_filter$filtered, input$rank_glom)
    pipeline_log("")
    log_msg(paste0("[", Sys.time(), "] Pipeline started"))

    filt_sdata <- res_filter$filtered()
    physeq0 <- r_values$phyobj_initial

    # --- Step 1: subset samples, drop now-empty taxa (matches microViz::ps_filter) ---
    log_msg("Step 1/3: Subsetting samples...")
    flog.info('subset samples() starting...')
    flog.info(paste0('number of samples before ', phyloseq::nsamples(physeq0)))
    if(!any(sample_names(physeq0) %in% as.vector(filt_sdata$sample.id))){
      shinyalert(title = "Oops", text = "Sample names do not match with names in metadata. Check the phyloseq object.", type = 'error')
      return(NULL)
    }
    physeq <- phyloseq::prune_samples(as.vector(filt_sdata$sample.id), physeq0)
    physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq) > 0, physeq)
    flog.info(paste0('number of samples after ', phyloseq::nsamples(physeq)))
    rownames(filt_sdata) <- filt_sdata %>% pull('sample.id')
    sample_data(physeq) <- filt_sdata
    log_msg(paste0("  -> ", phyloseq::nsamples(physeq), " samples, ", phyloseq::ntaxa(physeq), " taxa retained"))

    # --- Step 2: agglomerate to the chosen rank (ASV = no-op) ---
    log_msg("Step 2/3: Taxonomy agglomeration & filtering...")
    withProgress({
      if(input$rank_glom != 'ASV'){
        physeq <- speedyseq::tax_glom(physeq, input$rank_glom)
        taxa_names(physeq) <- tax_table(physeq)[, input$rank_glom]
      }
    }, message = 'Agglomerating taxonomy…')

    # --- Step 3: abundance + prevalence filtering ---
    # phyloseq_filter_taxa_tot_fraction uses a strict '>' (frac is exclusive); at
    # frac = 0 it merely re-drops the empty taxa already removed above (no-op).
    physeq <- metagMisc::phyloseq_filter_taxa_tot_fraction(physeq, frac = input$minAb)
    physeq <- metagMisc::phyloseq_filter_prevalence(physeq, prev.trh = input$minPrev)
    if(input$rank_glom != 'ASV'){
      tax_table(physeq) <- tax_table(physeq)[, 1:match(input$rank_glom, rank_names(physeq))]
    }
    log_msg(paste0("  -> ", phyloseq::ntaxa(physeq), " taxa after glom + filters"))
    showNotification("Filter taxonomy done...", type = "message", duration = 1)
    physeq
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  # A fresh processed object resets the manual taxa filter to "keep all" (the
  # default, matching microViz). The user can refine it afterwards via the
  # "Apply taxa filter" button, once the table is calibrated on the subset.
  observeEvent(processed(), {
    final(processed())
    log_msg(paste0("[", Sys.time(), "] Pipeline complete!"))
    showNotification("Dataset ready!", type = "message", duration = 5)
    r$data_ready(TRUE)
  }, ignoreNULL = TRUE)


  observe({
    flog.info('updating rank_glom selectInput...')
    updateSelectInput(session, "rank_glom",
                      choices = c( rank_names(phyloseq_data()), "ASV" ),
                      selected = "ASV")
  })


  # minAb / minPrev sliders are calibrated on the loaded object (they are *inputs*
  # to processed(), set before clicking Process Data) — not on a pipeline
  # intermediate, which previously created a mutable-cursor dependency.
  observe({
    req(phyloseq_data())
    flog.info('updating minAb numericInput...')
    updateAutonumericInput(session,
                           'minAb',
                           paste0("Minimum taxa overall percent abundance (max: ",
                                  round(max(microbiome::abundances(phyloseq_data(), transform = 'compositional'))), "):"),
                           value = 0,
                           options = list(maximumValue = 1,
                                          minimumValue = 0))
  })


  observe({
    req(phyloseq_data())
    flog.info('updating minPrev numericInput...')
    updateAutonumericInput(session,
                           'minPrev',
                           paste0("Minimum taxa prevalence in percent of samples (min:",
                                  round(min(microbiome::prevalence(phyloseq_data())),4),
                                  " max:",
                                  round(max(microbiome::prevalence(phyloseq_data())),4),")"),
                           value = 0,
                           options = list(maximumValue = 1,
                                          minimumValue = 0))
  })


  # Taxonomy preview source: the loaded object with the *current sample filter*
  # applied, then agglomerated to the chosen rank. Available immediately (NOT
  # gated on Process Data) so the taxonomy table + its datamods filter render on
  # load and refresh live as samples are selected/deselected.
  #
  # res_filter$filtered() lags by one datamods round-trip — harmless for a live
  # preview, and this is never the committed pipeline: processed() reads the
  # sample filter itself at Process-Data time. Naming mirrors processed()
  # (glommed taxa renamed to the rank value) and processed() applies the same
  # sample subset (plus abundance/prevalence filters), so processed()'s taxa are
  # always a same-named subset of this preview — keeping the "Apply taxa filter"
  # intersect against taxa_names(processed()) valid.
  tax_glom_obj <- reactive({
    req(phyloseq_data(), input$rank_glom)
    ps <- phyloseq_data()

    # Track the live sample selection (drop now-empty taxa, like processed()).
    filt_sdata <- res_filter$filtered()
    if(!is.null(filt_sdata) && "sample.id" %in% colnames(filt_sdata)){
      keep <- as.vector(filt_sdata$sample.id)
      if(any(sample_names(ps) %in% keep)){
        ps <- phyloseq::prune_samples(keep, ps)
        ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
      }
    }

    if(input$rank_glom != 'ASV'){
      withProgress({
        ps <- speedyseq::tax_glom(ps, input$rank_glom)
        taxa_names(ps) <- tax_table(ps)[, input$rank_glom]
      }, message = 'Agglomerating taxonomy…')
    }
    ps
  })

  # Layer the abundance + prevalence thresholds on top of tax_glom_obj(), in the
  # same order as processed(), so the taxonomy preview is WYSIWYG of what Process
  # Data will keep. Kept as a separate (cheap) reactive so moving the minAb /
  # minPrev sliders re-runs only these two taxa filters, never the heavy
  # tax_glom. Sliders default to 0 (autonumericInput initial value), at which
  # both filters are no-ops — so the preview is unchanged until a slider moves.
  tax_preview_obj <- reactive({
    req(tax_glom_obj())
    ps <- tax_glom_obj()
    minAb <- if(is.null(input$minAb)) 0 else input$minAb
    minPrev <- if(is.null(input$minPrev)) 0 else input$minPrev
    ps <- metagMisc::phyloseq_filter_taxa_tot_fraction(ps, frac = minAb)
    ps <- metagMisc::phyloseq_filter_prevalence(ps, prev.trh = minPrev)
    ps
  })

  render_taxonomy_table <- reactive({
    withProgress({
      req(tax_preview_obj(), input$rank_glom)
      flog.info('render_taxonomy_table fun')

      phyloseq_obj <- tax_preview_obj()
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
    },message = "Building taxonomy table…")

  })

  ## Filter taxo — calibrated on the loaded object (tax_glom_obj) so the table
  ## and filter widgets render before Process Data; default ranges keep all taxa.

  res_filter_taxo <- datamods::filter_data_server(
    id = "filtering_taxo",
    data = reactive({
      req(render_taxonomy_table())
      render_taxonomy_table()
    }),
    name = reactive("tax_table"),
    vars = reactive({
      req(render_taxonomy_table())
      s_names <- phyloseq::sample_names(phyloseq_data())
      col_names <- colnames(render_taxonomy_table())
      # Exclude per-sample abundance columns and the DNA `sequences` column — a
      # sequence-as-filter-widget is meaningless and very heavy in datamods.
      filt <- dplyr::setdiff(col_names, c(s_names, "sequences"))
    return(filt)
    }),
    widget_num = "slider",
    widget_date = "slider",
    label_na = "Missing"
  )

  output$table_taxoFILT <- DT::renderDT({
    res_filter_taxo$filtered()
  }, options = list(pageLength = 10, scrollX = TRUE))

  # Manual taxa refinement: prune processed() to the taxa surviving the (now
  # correctly calibrated) datamods taxa filter. A dedicated button keeps the heavy
  # normalization off the slider-drag path and avoids reading a stale snapshot.
  observeEvent(input$apply_taxa_filter, {
    req(processed())
    withProgress({
      filttax <- res_filter_taxo$filtered()
      selected <- intersect(as.character(filttax[[1]]), phyloseq::taxa_names(processed()))
      final(phyloseq::prune_taxa(selected, processed()))
      flog.info(paste0('apply_taxa_filter -> ', length(selected), ' taxa'))
      log_msg(paste0("Manual taxa filter applied -> ", length(selected), " taxa retained"))
      showNotification(paste0(length(selected), " taxa retained after manual filter."), type = "message", duration = 3)
    }, message = "Subset taxonomy, please wait...")
  }, ignoreInit = TRUE)




  # Normalization is gated on final() (which changes only on Process Data or an
  # explicit manual taxa-filter apply), so heavy VST is not re-run when the user
  # merely toggles the Normalization radio — honouring the eventReactive
  # performance pattern. norm_method is read at trigger time, matching the
  # original behaviour where normalization ran inside the Process-Data handler.
  normalized <- eventReactive(final(), {
    req(final(), input$norm_method)
    FGdata <- final()

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
    FNGdata
  })


  output$phy_after <- renderPrint({
    req(final())
    print(final())
  })

  output$phy_norm <- renderPrint({
    req(normalized())
    print(normalized())
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
      req(final())
      write.table(merge_table(input$rank_glom, final()), file, sep="\t", row.names=FALSE)
    }
  )

  output$filt_norm_otable_download <- downloadHandler(
    filename = "filt_norm_asv_table.csv",
    content = function(file) {
      req(normalized())
      write.table(merge_table(input$rank_glom, normalized()), file, sep="\t", row.names=FALSE)
    }
  )

  output$filt_rdata_download <- downloadHandler(
    filename = "filt_robject.rdata",
    content = function(file) {
      req(final())
      data = final()
      save(data, file = file)
    }
  )

  output$filt_rdata_norm_download  <- downloadHandler(
    filename = "filt_norm_robject.rdata",
    content = function(file) {
      req(normalized())
      data = normalized()
      save(data, file = file)
    }
  )

  output$filt_refseq_download <- downloadHandler(
    filename = "filt_ref-seq.fasta",
    content = function(file) {
      req(final())
      if(!is.null(refseq(final(), errorIfNULL=FALSE))){
        Biostrings::writeXStringSet(refseq(final()), file)
      }else(showNotification("FASTA Download failed. No refseq in object.", type="error", duration = 5))
    }
  )

  # Saving variable for other modules.
  # Flag flipped TRUE once "Process Data" pipeline completes; gates navigation.
  r$data_ready <- reactiveVal(FALSE)

  # Raw object loaded from file.
  r$phyloseq_data <- reactive({
    req(r_values$phyobj_initial)
    r_values$phyobj_initial

  })

  # final filtered object
  r$phyloseq_filtered <- reactive({
    req(final())
    final()
  })


  # final filtered object normalize
  r$phyloseq_filtered_norm <- reactive({
    req(normalized())
    normalized()
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
    req(final())
    sdat <- sample_data(final())
    sdat <- sdat[,which(unlist(lapply(sdat, function(x)!all(is.na(x))))),with=F]
    sdat <- as(sdat, "data.frame")
    if(! 'sample.id' %in% colnames(sdat)){
      sdat$sample.id <- rownames(sdat)
    }
    return(sdat)
  })

  r$var_list <- reactive({
    req(final(), r$sdat)
    sdat <- r$sdat()
    var_list <- colnames(sdat)
    # if('sample.id' %in% var_list){
    #   var_list <- sort(var_list[! var_list %in% 'sample.id'])
    # }
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
    tagList(
      shinyWidgets::pickerInput(
        ns("combin_1"),
        label = "Variable to combine",
        choices = choice_combine_var(),
        selected = NULL,
        multiple = TRUE,
        options = pickerOptions(maxOptions = 1)
      ),
      lapply(2:4, FUN = function(i){
        uiOutput(ns(paste0("ui_combin_", i)))
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
    lapply(3:4, FUN = function(i){
      if(is.null(input[[paste0("combin_", i)]])){
        get_ui_combin(i)
      }
    })
  })
  
  
  r$factor_list <- reactive({
    req(r$sdat(), r$var_list())
    num <- sapply(r$sdat()[, r$var_list()], is.numeric)
    fact <- r$var_list()[!num]
    return(fact)
  })

  # Color reactives (r$factor_colors, r$numeric_palettes, r$taxa_colors) are
  # owned by mod_color, mounted below as a child of the "Colors" nav panel.
  mod_color_server("color_ui_1", r = r)

  })
}

## To be copied in the UI
# mod_data_loading_ui("data_loading_ui_1")

## To be copied in the server
# mod_data_loading_server("data_loading_ui_1")
