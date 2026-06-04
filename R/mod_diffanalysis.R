# Module UI

#' @title   mod_diffanalysis_ui and mod_diffanalysis_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_diffanalysis
#'
#' @keywords internal
#' @noRd
#' @importFrom shiny NS tagList
#' @importFrom bslib layout_sidebar sidebar accordion accordion_panel navset_card_underline nav_panel card card_header tooltip input_task_button
#' @importFrom bsicons bs_icon
mod_diffanalysis_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "Settings",
      open = "desktop",
      width = "350px",

      htmltools::p(
        "Configure the comparison, then launch DESeq2, MetaGenomeSeq and ",
        "MetaCoder in one click. Results and merged outputs appear on the right.",
        style = "font-size: 0.9em; color: grey;"
      ),

      tags$hr(),

      accordion(
        id = ns("config_accordion"),
        open = "Analysis Setup",
        multiple = TRUE,

        accordion_panel(
          "Analysis Setup",
          icon = bs_icon("sliders"),
          tooltip(
            uiOutput(ns("factor1")),
            "Factor used to build the contrast for all methods.",
            placement = "right"
          ),
          uiOutput(ns("cond1")),
          uiOutput(ns("cond2")),
          tooltip(
            numericInput(ns("pval"),
                         label = "Adjusted p-value threshold:",
                         min = 0, max = 1,
                         value = 0.05),
            "Threshold applied to adjusted p-values in every method and in the merge step.",
            placement = "right"
          )
        ),

        accordion_panel(
          "Plot Options",
          icon = bs_icon("bar-chart"),
          tooltip(
            sliderInput(ns("Nmeth"), "Number of diff. methods:",
                        min = 1, max = 3, value = 1),
            "Minimum number of methods (DESeq2 / MGseq / MetaCoder) flagging a feature as significant.",
            placement = "right"
          ),
          tooltip(
            numericInput(ns("minAb"), "Minimum mean relative abundance:",
                         value = 0.001, min = 0, max = 1, step = 1000),
            "Minimum mean relative abundance in at least one of the two conditions.",
            placement = "right"
          ),
          tooltip(
            sliderInput(ns("Nfeat"), "Number of features to plot:",
                        min = 0, max = 100, value = 50),
            "Top features (by |DESeq2 log2FC|) kept for the merged barplot.",
            placement = "right"
          )
        )
      ),

      tags$hr(),

      div(
        style = "margin: 1rem 0;",
        tooltip(
          input_task_button(ns("launch_diff"), "Run Differential Analyses",
                       icon = bs_icon("play-circle-fill"),
                       class = "btn-primary w-100 btn-lg",
                       label_busy = "Running…"),
          "Runs DESeq2, MetaGenomeSeq and MetaCoder with the current settings.",
          placement = "top"
        )
      )
    ),

    navset_card_underline(
      title = "Differential Analysis Results",
      full_screen = TRUE,

      nav_panel(
        "DESeq2",
        icon = bs_icon("table"),
        DT::dataTableOutput(ns("deseqTab"))
      ),

      nav_panel(
        "MetaGenomeSeq",
        icon = bs_icon("table"),
        DT::dataTableOutput(ns("MGseqTab"))
      ),

      nav_panel(
        "MetaCoder — table",
        icon = bs_icon("table"),
        DT::dataTableOutput(ns("mtcoderTab"))
      ),

      nav_panel(
        "Merged — table",
        icon = bs_icon("layers"),
        div(
          downloadButton(outputId = ns("merge_download"), label = "Download Table"),
          downloadButton(outputId = ns("fasta_download"), label = "Download FASTA")
        ),
        DT::dataTableOutput(ns("mergeTab"))
      ),

      nav_panel(
        "Merged — barplot",
        icon = bs_icon("bar-chart-line"),
        plotlyOutput(ns("barplot1"), height = "800px")
      ),

      nav_panel(
        "Merged — heat tree",
        icon = bs_icon("diagram-3"),
        div(
          downloadButton(outputId = ns("merged_heattree_download"), label = "Download plot")
        ),
        plotOutput(ns("plot_merged_heattree"), height = "1000px")
      )
    )
  )
}

# Module Server

#' @rdname mod_diffanalysis
#' @noRd
#' @keywords internal
#' @import DESeq2
#' @import metacoder
#' @import phyloseq
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @importFrom DT renderDataTable
#' @importFrom broom tidy
#' @import metagenomeSeq
#' @importFrom Biobase pData
#' @importFrom VennDiagram venn.diagram
#' @import ggplot2


mod_diffanalysis_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns

  factor_input   <- reactive({ input$diff_factor })
  get_meta_col   <- make_get_meta_col(factor_input, r)
  local_metadata <- make_local_metadata(factor_input, get_meta_col, r)
  isNumFactor    <- make_is_num_factor(get_meta_col, local_metadata)


  local_physeq <- reactive({
    req(local_metadata(), input$diff_factor, r$phyloseq_filtered())
    phy <- r$phyloseq_filtered()
    sample_data(phy) <- sample_data(local_metadata())
    if(isNumFactor()){
      s.list <- na.omit(local_metadata()[,c('sample.id', input$diff_factor)])[, 'sample.id']
      phy <- prune_samples(s.list, phy)
    } else {
      req(input$Cond1, input$Cond2)
      keep <- sample_data(phy)[[input$diff_factor]] %in% c(input$Cond1, input$Cond2)
      phy <- prune_samples(keep, phy)
    }
    phy <- prune_taxa(taxa_sums(phy) >= 1, phy)
    phy <- prune_samples(sample_sums(phy) >=1, phy)
    return(phy)
  })


  output$factor1 = renderUI({
    req(r$phyloseq_filtered(), r$var_list())
    shinyWidgets::pickerInput(ns("diff_factor"),
                              label = "Select factor to test: ",
                              choices = r$var_list(),
                              selected = r$var_list()[1],
                              multiple = FALSE,
                              options = pickerOptions(
                                actionsBox = TRUE,
                                liveSearch = TRUE,
                                showContent = FALSE
                              ),
                              choicesOpt = list(
                                content = unlist(lapply(
                                  X = r$var_list(),
                                  FUN = function(x) {
                                    htmltools::doRenderTags(
                                      tags$div(
                                        splitLayout(cellWidths = 200,
                                                    tags$div(
                                                      style = htmltools::css(fontWeight = "bold"),
                                                      x
                                                    ),
                                                    tags$div(
                                                      style = htmltools::css(color = 'grey'),
                                                      class(r$sdat()[,x])
                                                    ),
                                                    tags$div(
                                                      style = htmltools::css(color = 'grey'),
                                                      paste0(sum(is.na(r$sdat()[,x])), '/', nrow(r$sdat()), ' NAs')
                                                    )
                                        )
                                      )
                                    )
                                  }
                                ))
                              )
    )
  })


  output$cond1 = renderUI({
    req(input$diff_factor, r$phyloseq_filtered())
    if(! isNumFactor()){
      selectInput(ns("Cond1"),
                  label = "Condition 1 to compare: ",
                  choices = unique(local_metadata()[,input$diff_factor])
      )
    }
  })


  output$cond2 = renderUI({
    req(input$Cond1, input$diff_factor, r$phyloseq_filtered())
    if(! isNumFactor()){
      Conds <- unique(local_metadata()[,input$diff_factor])
      choices2 <- Conds[Conds != input$Cond1]
      selectInput(ns("Cond2"),
                  label = "Condition 2 to compare: ",
                  choices = choices2
      )
    }
  })


  deseqDA = eventReactive(input$launch_diff, {
    withProgress({
      req(input$diff_factor, r$phyloseq_filtered())
      if(! isNumFactor()) req(input$Cond1, input$Cond2)

      deseq <- phyloseq_to_deseq2(local_physeq(), as.formula(paste0("~", input$diff_factor)))
      gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      }
      geoMeans <- apply(counts(deseq), 1, gm_mean)
      deseq <- estimateSizeFactors(deseq, geoMeans = geoMeans)
      deseq <- try(DESeq(deseq, test="Wald", fitType="parametric"))

      if(class(deseq) == "try-error"){
        validate(
          need(FALSE, "Not enough replicates in these conditions")
        )
      }
      if( isNumFactor() ){
        res <- results(deseq, cooksCutoff = FALSE)
      } else {
        res <- results(deseq, cooksCutoff = FALSE, contrast = c(input$diff_factor, input$Cond1 , input$Cond2))
      }
      return(res)
    }, message = "Performing DESeq2...")
  })


  output$deseqTab <- DT::renderDataTable({
    req(deseqDA(), input$pval, snap$physeq)
    res <- deseqDA()
    phy_snap <- snap$physeq
    # Construct table
    ttable1 <- as.data.frame(phyloseq::tax_table(phy_snap)) %>%
              rownames_to_column()

    resDESeq <- as.data.frame(res) %>%
      rownames_to_column() %>%
      mutate(absLFC = abs(log2FoldChange)) %>%
      # filter(padj<=input$pval) %>%
      left_join(ttable1, by="rowname")


    if(!is.null(refseq(phy_snap, errorIfNULL=FALSE))){
      sseq1 <- as.data.frame(phyloseq::refseq(phy_snap)) %>%
        rownames_to_column()

      if(nrow(sseq1) != 0){
        sseq1 <- rename(sseq1, sequence = 2)
      }

      resDESeq <- resDESeq %>% left_join(sseq1, by="rowname")
    }
    round_df(as.data.frame(resDESeq), 4)
  }, filter="top", options = list(scrollX = TRUE))


  #### METAGENOMESEQ
  mgSeqDA = eventReactive(input$launch_diff, {
    withProgress({
      req(input$diff_factor, local_physeq())
      if(! isNumFactor()) req(input$Cond1, input$Cond2)
      flog.info('metaGseq...')
      MGdata <- phyloseq_to_metagenomeSeq(local_physeq())
      flog.info("phyloseq to metaGseq")

      #FitFeature model : zero-inflated log-normal model
      pd <- pData(MGdata)
      mod <- model.matrix(as.formula(paste("~", input$diff_factor)), data = pd)

      res1 = NULL
      tryCatch( {res1 = fitFeatureModel(MGdata, mod)} ,
                error=function(e){e;cat("ERROR :",conditionMessage(e), "\n")})
      if(is.null(res1)){
        validate(
          need(FALSE, "Not enough replicates in these conditions")
        )
      }
      TAB = MRcoefs(res1) #fdr adjustment
      return(TAB)
      }, message="Performing metagenomeSeq...")
  })


  output$MGseqTab <- DT::renderDataTable({
    round_df(as.data.frame(mgSeqDA()), 4)
  }, filter="top", options = list(scrollX = TRUE))


    ### MEtacoder
  mtcoderDA <- eventReactive(input$launch_diff, {
    withProgress({
      req(input$diff_factor, r$phyloseq_filtered())
      if(! isNumFactor()) req(input$Cond1, input$Cond2)
      mean_ratio <- function(abund_1, abund_2) {
        log_ratio <- log2(mean(abund_1) / mean(abund_2))
        if (is.nan(log_ratio)) {
          log_ratio <- 0
        }
        list(log2_mean_ratio = log_ratio,
             median_diff = median(abund_1) - median(abund_2),
             mean_diff = mean(abund_1) - mean(abund_2),
             wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
      }
      table <- NULL
      if( isNumFactor()){
        table <- NULL
      } else{
        psobj <- local_physeq()
        normf = function(x, tot=max(phyloseq::sample_sums(psobj))){ tot*x/sum(x) }
        psobj <- phyloseq::transform_sample_counts(local_physeq(), normf)

        incProgress(amount = 0.2, message = 'Zeroing low counts...')
        # metacoder::parse_phyloseq() (v0.3.9) references a bare `ranks_ref`
        # symbol that lives in metacoder's lazy-data env, which is *not* on its
        # namespace lookup chain. It only resolves when metacoder is attached
        # via library() (we only import it), so the bare lookup falls through to
        # .GlobalEnv. Seed it there so parse_phyloseq finds it. See git history:
        # data/ranks_ref.rda was an earlier, deploy-broken attempt at this.
        if (!exists("ranks_ref", envir = .GlobalEnv, inherits = FALSE)) {
          utils::data("ranks_ref", package = "metacoder", envir = .GlobalEnv)
        }
        flog.info('metacoder - parse_phyloseq')
        obj <- metacoder::parse_phyloseq(psobj, class_regex = "(.*)", class_key = "taxon_name")
        flog.info('metacoder - zero_low_counts')
        obj$data$otu_table <- metacoder::zero_low_counts(obj, "otu_table", min_count = 1000, use_total = TRUE, cols = obj$data$sample_data$sample_id)
        no_reads <- rowSums(obj$data$otu_table[, obj$data$sample_data$sample_id]) == 0
        flog.info('metacoder - filter_obs')
        obj <- metacoder::filter_obs(obj, "otu_table", ! no_reads, drop_taxa = TRUE)
        if(nrow(obj$data$otu_table)==0){return(NULL)}

        incProgress(amount = 0.2, message = 'Calculating taxon abundance...')
        flog.info('metacoder - calc_taxon_abund')
        obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "otu_table",  cols = obj$data$sample_data$sample_id)
        obj$data$tax_abund$total <- rowSums(obj$data$tax_abund[, -1]) # -1 = taxon_id column
        flog.info('metacoder - calc_n_samples')
        # cols = sample_id, otherwise calc_n_samples auto-detects all numeric
        # columns and counts the manually-added `total` column as a sample,
        # triggering metacoder's "groups without cols" NOTE and an off-by-one.
        obj$data$n_samples <- metacoder::calc_n_samples(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id)


        incProgress(amount = 0.4, message = 'Comparing groups...')
        flog.info('metacoder - compare_groups')
        obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = obj$data$sample_data[[input$diff_factor]], func = mean_ratio)
        flog.info('metacoder - wilcox_p_value')
        obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
        table <- merge(obj$data$diff_table, obj$data$tax_data,by='taxon_id')
      }
      }, message = "Performing metacoder...", min = 0, max = 1)
    return_obj <- list(table = table)
      return(return_obj)
    })


    # ---- Run-on-click orchestration -----------------------------------------
    # Snapshot the contrast inputs (factor + conditions + filtered phyloseq)
    # at the moment the user clicks the Run button. Downstream displays
    # (deseqTab, mergeTab, barplot) read from this snapshot so they stay
    # consistent with the data that was actually analysed, even if the user
    # later tweaks Cond1 / Cond2 / diff_factor without re-clicking Run.
    snap <- reactiveValues(
      factor = NULL,
      cond1  = NULL,
      cond2  = NULL,
      is_num = NULL,
      physeq = NULL
    )

    observeEvent(input$launch_diff, {
      req(input$diff_factor, r$phyloseq_filtered())
      if (! isNumFactor()) req(input$Cond1, input$Cond2)

      # Capture state of the run
      snap$factor <- input$diff_factor
      snap$cond1  <- input$Cond1
      snap$cond2  <- input$Cond2
      snap$is_num <- isNumFactor()
      snap$physeq <- local_physeq()

      # Force all three differential analyses to evaluate now,
      # regardless of which result tab is currently visible. Their
      # eventReactive caches are populated here, so subsequent tab
      # navigation just reads cached results without re-running.
      deseqDA()
      mgSeqDA()
      mtcoderDA()
    }, ignoreInit = TRUE)


    output$mtcoderTab <- DT::renderDataTable({
      round_df(as.data.frame(mtcoderDA()$table), 4)
    }, filter="top", options = list(scrollX = TRUE))


    mergeList <- reactive({
      withProgress({
        # Use the snapshotted contrast (factor / Cond1 / Cond2 / is_num)
        # captured at the last button click. This guarantees the merged
        # output is always consistent with the DA results currently cached
        # in deseqDA() / mgSeqDA() / mtcoderDA(), even if the user has
        # since changed the contrast inputs without re-clicking Run.
        req(input$pval, snap$factor)
        if(! snap$is_num) req(snap$cond1, snap$cond2)
        flog.info("mergeList")

        flog.info("mergeList - metacoder")
        mtList <- NULL
        if(! snap$is_num){
          mtTab <- mtcoderDA()$table
          if(!is.null(mtTab)){
            mtList <- mtTab[mtTab$wilcox_p_value <= input$pval, "otu_id"]
          }
        }

        if(all(is.na(mtList))){
          mtList <- NULL
        }

        flog.info("mergeList - deseq")
        deTab <- deseqDA()
        deTab <- deTab[!is.na(deTab$padj),]
        deList <- row.names(deTab[deTab$padj <= input$pval,])
        if(all(is.na(deList))){
          deList <- NULL
        }

        flog.info("mergeList - metagenomeSeq")
        mgTab <- mgSeqDA()
        mgList <- row.names(na.omit(mgTab[mgTab$adjPvalues <= input$pval,]))
        if(all(is.na(mgList))){
          mgList <- NULL
        }

        TF <- list(x1=deList, x2=mgList, x3=mtList)
        names(TF) <- c("DESeq", "metagenomeSeq", "metacoder")

        ListAllOtu = unique(unlist(TF))

        #Construction de la table
        flog.info('mergeList - Building table...')
        comp1 <- paste(snap$cond1, '_vs_' , snap$cond2, sep='')
        col_comp <- rep(comp1, length(ListAllOtu))
        TABf <- cbind.data.frame(ListAllOtu, col_comp)


        # Test si chaque ASV est diff dans les methdes.
        flog.info('mergeList - Check methods ...')
        for (j in 1:length(TF)){
          TABtest <- TF[[j]]

          TABtest_signif <- rep(0, length(ListAllOtu))
          names(TABtest_signif) <- ListAllOtu
          TABtest_signif[TABtest] <- 1

          TABf <- cbind.data.frame(TABf, TABtest_signif)
          names(TABf)[ncol(TABf)] <- names(TF)[j]
        }

        TABfbak <- TABf

        flog.info('mergeList - Adding LFC...')
        TABf <- cbind( TABf, sumMethods = apply(TABf[3:5], 1, sum, na.rm=TRUE),
                       DESeqLFC = deTab[as.character(TABf[,1]),"log2FoldChange"],
                       absDESeqLFC = abs(deTab[as.character(TABf[,1]),"log2FoldChange"]) )

        flog.info('mergeList - Calculating mean relative abundance...')
        data <- r$phyloseq_filtered()
        normf <- function(x){ x/sum(x) }
        data.norm <- transform_sample_counts(data, normf)
        otableNORM <- otu_table(data.norm)
        ssample <- as.matrix(sample_data(data.norm))
        ttax <- tax_table(data.norm)

        seqs = NULL
        if(!is.null(refseq(data.norm, errorIfNULL=FALSE)) ){
          seqs <- refseq(data.norm)
        }

        flog.info("mergeList - mean1")
        Gtab <- cbind(as.data.frame(ssample), t(otableNORM))
        MeanRelAbcond1 <- NULL
        for(i in TABf$ListAllOtu){
          tt <- mean(Gtab[Gtab[, snap$factor] == snap$cond1, i], na.rm=TRUE)
          MeanRelAbcond1 <- c(MeanRelAbcond1,tt)
        }
        flog.info("mergeList - mean2")
        MeanRelAbcond2 <- NULL
        for(i in TABf$ListAllOtu){
          tt <- mean(Gtab[Gtab[, snap$factor] == snap$cond2, i], na.rm=TRUE)
          MeanRelAbcond2 <- c(MeanRelAbcond2,tt)
        }
        TABfbak <- TABf <- cbind(TABf, MeanRelAbcond1, MeanRelAbcond2)

        #Adjust table
        flog.info('mergeList - Adjusting table...')
        TABf <- TABf[!is.na(TABf$DESeqLFC),]
        TABf$Condition <- rep(NA, nrow(TABf))
        TABf[TABf$DESeqLFC>0, "Condition"] <- as.character(snap$cond1)
        TABf[TABf$DESeqLFC<0, "Condition"] <- as.character(snap$cond2)
        TABf$Condition <- factor(TABf$Condition,  levels = c(as.character(snap$cond1), as.character(snap$cond2)) )

        flog.info("mergeList - Adding taxonomy and sequences...")

        TABf <- cbind.data.frame(TABf, ttax[as.character(TABf[,1]),]@.Data)
        if(!is.null(seqs)){
          TABf <- cbind.data.frame(TABf, sequences = seqs[as.character(TABf[,1])])
        }



        LL <- list()
        LL$TABf <- TABf
        LL$ttax <- ttax
        LL$signifseqs <- seqs[as.character(TABf[,1])]
        flog.info("mergeList - done")
      }, message = "Merging results...")
      return(LL)
    })


    output$mergeTab <- DT::renderDataTable({
      round_df(mergeList()$TABf, 4)
      }, filter="top", options = list(scrollX = TRUE), rownames = FALSE)


    output$merge_download <- downloadHandler(
      filename = "aggregate_table.csv",
      content = function(file) {
        req(mergeList())
        write.table(mergeList()$TABf, file, sep="\t", row.names=FALSE)
      }
    )

    output$fasta_download <- downloadHandler(
      filename = "signif_seqs.fasta",
      content = function(file) {
        req(mergeList(), r$phyloseq_filtered())
        if(!is.null(refseq(r$phyloseq_filtered(), errorIfNULL=FALSE))){
          writeXStringSet(mergeList()$signifseqs, file)
        }else(showNotification("FASTA Download failed. No refseq in object.", type="error", duration = 5))
      }
    )


    reacbarplot1 <- reactive({
      req(mergeList(), input$Nmeth, input$minAb, input$Nfeat)
      flog.info('reacbarplot1')
      TABf <- mergeList()$TABf
      ttax <- mergeList()$ttax
      # Barplot
      ## differentialy abundant features on 2 methods
      ## Top abs(LogFolchange)
      ## feature with relative abondance > 0.1% (0.001)

      TABbar <- TABf[TABf$sumMethods >= input$Nmeth, ]
      TABbar <- TABbar[TABbar$MeanRelAbcond1 >= input$minAb | TABbar$MeanRelAbcond2 >= input$minAb, ]
      TABbar <- tail(TABbar[order(abs(TABbar$DESeqLFC)),], input$Nfeat)

      if(nrow(TABbar)!=0){
        if(r$rank_glom() == "ASV"){
          TABbar$tax <- paste( substr(ttax[as.character(TABbar$ListAllOtu),"Species"],1,20),"...","_", TABbar$ListAllOtu, sep="")
        }else{
          TABbar$tax <- TABbar$ListAllOtu
        }
        flog.info('reacbarplot1 - ggplot2')

        validate(need(!snap$factor %in% names(r$numeric_palettes()),
                      "Bar plot requires a categorical factor."))
        p <- ggplot2::ggplot(data = TABbar, aes(x = reorder(tax, -abs(DESeqLFC)), y = DESeqLFC, fill = Condition ) ) +
                geom_bar(stat="identity", alpha = 0.7) + ggtitle(glue("{snap$cond1} vs. {snap$cond2}")) + labs(x='Features') +
                coord_flip() + theme_bw() +
                scale_y_continuous(minor_breaks = seq(-1E4 , 1E4, 1), breaks = seq(-1E4, 1E4, 5)) + scale_fill_manual(values = r$factor_colors()[[snap$factor]])

        p <- plotly::ggplotly(p)

      }else{
        showNotification("No ASV to plot... (all features show summethods <=2 ? too high mean relative abundance filters ?)", type="error", duration = 10)
      }
      return(p)
    })

    output$barplot1 <- renderPlotly({
        reacbarplot1()
    })


    # Heat-tree view of the same merged consensus table as the barplot.
    # Builds a metacoder taxmap straight from TABf's taxonomy columns via
    # parse_tax_data() (no parse_phyloseq, so no `ranks_ref` workaround needed)
    # and colours nodes by DESeqLFC, sized by feature count.
    reac_merged_heattree <- reactive({
      req(mergeList(), input$Nmeth, input$minAb, input$Nfeat)
      validate(need(!snap$factor %in% names(r$numeric_palettes()),
                    "Merged heat tree requires a categorical factor."))
      flog.info("reac_merged_heattree")

      TABf  <- mergeList()$TABf
      ttax  <- mergeList()$ttax

      # Identical filtering to reacbarplot1 so tree and barplot share features.
      TABbar <- TABf[TABf$sumMethods >= input$Nmeth, ]
      TABbar <- TABbar[TABbar$MeanRelAbcond1 >= input$minAb |
                       TABbar$MeanRelAbcond2 >= input$minAb, ]
      TABbar <- tail(TABbar[order(abs(TABbar$DESeqLFC)), ], input$Nfeat)
      validate(need(nrow(TABbar) > 0,
                    "No features to plot with the current Plot Options filters."))

      # Build a taxmap from the merged table's taxonomy columns. Drop ranks
      # that are entirely NA in the filtered subset (tax_glom case), then
      # placeholder remaining NA/"" so lineages stay well-formed and unrelated
      # features don't collapse mid-lineage.
      tax_cols <- colnames(ttax)
      keep <- vapply(tax_cols, function(cn) !all(is.na(TABbar[[cn]])), logical(1))
      tax_cols <- tax_cols[keep]
      tax_input <- TABbar[, c(tax_cols, "DESeqLFC")]
      tax_input[tax_cols] <- lapply(tax_input[tax_cols], function(v) {
        v <- as.character(v); v[is.na(v) | v == ""] <- "unknown"; v
      })

      obj <- metacoder::parse_tax_data(tax_input, class_cols = tax_cols,
                                       named_by_rank = TRUE)

      # Aggregate the per-feature DESeqLFC to each node (mean of descendant
      # features); obj$obs() is recursive, so internal nodes average all
      # features beneath them while leaves keep their own value.
      obj$data$taxon_lfc <- data.frame(
        taxon_id = obj$taxon_ids(),
        mean_lfc = vapply(obj$obs("tax_data"),
                          function(i) mean(obj$data$tax_data$DESeqLFC[i], na.rm = TRUE),
                          numeric(1))
      )

      # heat_tree() labels via ggfittext, which warns "Ignoring unknown
      # aesthetics: xmin/xmax/ymin/ymax" under current ggplot2 (metacoder
      # internals we can't change). Mute only that, like the metacoder tab.
      withCallingHandlers(
        heat_tree(obj,
                  node_label = taxon_names,
                  node_size  = n_obs,
                  node_color = mean_lfc,
                  node_color_range = c(r$factor_colors()[[snap$factor]][[snap$cond2]],
                                       "gray",
                                       r$factor_colors()[[snap$factor]][[snap$cond1]]),
                  node_size_axis_label  = "Feature count",
                  node_color_axis_label = "DESeq2 log2FoldChange",
                  layout = "davidson-harel",
                  initial_layout = "reingold-tilford"),
        warning = function(w) {
          if (grepl("Ignoring unknown aesthetics", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
      )
    })

    output$plot_merged_heattree <- renderPlot({
      reac_merged_heattree()
    })

    output$merged_heattree_download <- downloadHandler(
      filename = "merged_heat_tree.svg",
      content = function(file){
        req(reac_merged_heattree())
        grDevices::svg(filename = file, width = 12, height = 12)
        print(reac_merged_heattree())
        dev.off()
      }
    )

  })
}

## To be copied in the UI
# mod_diffanalysis_ui("diffanalysis_ui_1")

## To be copied in the server
# mod_diffanalysis_server("diffanalysis_ui_1", r = r)
