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
#' @export
#' @importFrom shiny NS tagList
mod_diffanalysis_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      box(title = "Settings:", width = 6, status = "warning", solidHeader = TRUE,
        uiOutput(ns("factor1")),
        fluidRow(column(3, uiOutput(ns("cond1")) ),
                 column(3, uiOutput(ns("cond2")) )
        ),
        uiOutput(ns("alpha1"))
      ),

      box(title = "Differential analysis:", width = 12, status = "primary", solidHeader = TRUE,
          tabsetPanel(
            tabPanel("DESeq2",
                     verbatimTextOutput(ns("print1")),
                     actionButton(ns("go1"), "Run DESeq2", icon = icon("play-circle"),
                                  style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
                     DT::dataTableOutput(ns("deseqTab"))

                     ),
            tabPanel("MetaGenomeSeq",
                     h1("Run MGseq with same settings:"),
                     actionButton(ns("go2"), "Run MGSeq", icon = icon("play-circle"),
                                  style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
                     DT::dataTableOutput(ns("MGseqTab"))
            ),
            tabPanel("MetaCoder",
                     h1("Run Metacoder with same settings:"),
                     actionButton(ns("go4"), "Run Metacoder", icon = icon("play-circle"),
                                  style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
                     DT::dataTableOutput(ns("mtcoderTab"))
            ),
            # tabPanel("Wilcox non parametric test",
            #          h1("Run Wilcox tests with same settings:"),
            #          actionButton(ns("go3"), "Run Wilcox"),
            #          DT::dataTableOutput(ns("WilcoxTab"))
            # ),
            tabPanel("Merge results",
                     h1("Merge results of differential analysis:"),
                     # verbatimTextOutput(ns("mergePrint")),
                     box(
                     downloadButton(outputId = ns("merge_download"), label = "Download Table"),
                     downloadButton(outputId = ns("fasta_download"), label = "Download FASTA"),
                     DT::dataTableOutput(ns("mergeTab")),
                     title = "Aggregate table:", width = 12, status = "primary", solidHeader = TRUE),

                     box(
                       sliderInput(ns("Nmeth"), "Number of diff. methods :",
                                   min = 1, max = 3, value = 1
                       ),
                     numericInput(ns("minAb"), "Minimum mean relative abundance:", value = 0.001, min = 0, max = 1, step = 1000),
                     sliderInput(ns("Nfeat"), "Number of features to plot:",
                                 min = 0, max = 100, value = 50
                     ),
                     plotlyOutput(ns("barplot1")), #, height = "1000px"
                     title = "Plotting features:", width = 12, status = "primary", solidHeader = TRUE)
            )

        )
      )
    )
  )
}

# Module Server

#' @rdname mod_diffanalysis
#' @export
#' @keywords internal
#' @importFrom DESeq2 estimateSizeFactors DESeq results counts
#' @importFrom metacoder parse_phyloseq zero_low_counts filter_obs calc_taxon_abund calc_n_samples compare_groups
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
#' @importFrom ggplot2 ggplot


mod_diffanalysis_server <- function(input, output, session, r = r){
  ns <- session$ns

  observeEvent(r$tabs$tabselected, {
    if(r$tabs$tabselected=='tab_diff' && !isTruthy(r$phyloseq_filtered())){
      shinyalert(title = "Oops", text="Phyloseq object not present. Return to input data and validate all steps.", type='error')
    }
  })

  output$factor1 = renderUI({
    req(r$phyloseq_filtered())
    selectInput(
      ns("Fact1"),
      label = "Select factor to test: ",
      choices = r$phyloseq_filtered()@sam_data@names
    )
  })

  output$cond1 = renderUI({
    req(input$Fact1, r$phyloseq_filtered())
    selectInput(ns("Cond1"),
      label = "Select Condition 1 to compare: ",
      choices = unique(r$phyloseq_filtered()@sam_data[,input$Fact1])
      )
  })

  output$cond2 = renderUI({
    req(input$Cond1, input$Fact1, r$phyloseq_filtered())
    Conds = unique(r$phyloseq_filtered()@sam_data[,input$Fact1])
    if(length(Conds) <= 2){
      print(Conds)
      choices2 = Conds
    }else{
      choices2 = Conds[Conds != input$Cond1]
    }
print(choices2)


    selectInput(ns("Cond2"),
                label = "Select Condition 2 to compare: ",
                choices = choices2
    )

  })

  output$alpha1 = renderUI({
    numericInput(ns("Alpha1"),
                 label = "Define pvalue threshold:",
                 min = 0, max = 1,
                 value = 0.05
                 )
  })


    output$print1 <- renderPrint({
      # req(input$Fact1)
      # print(input$Fact1)
      # print(unique(r$phyloseq_filtered()@sam_data[,input$Fact1]))
      cat("\n")
      print( glue("Compare {input$Cond1} and {input$Cond2} ") )

      # print(head(tax_table(r$phyloseq_filtered())))
    })


    #Taxonomy subset (asvselect)
    # data1 <- reactive({
    #   req(r$asvselect(), r$phyloseq_filtered())
    #   Fdata <- prune_taxa(r$asvselect(), r$phyloseq_filtered())
    #   Fdata
    # })


    deseqDA = eventReactive(input$go1, {

      withProgress({

      req(input$Fact1, input$Cond1, input$Cond2, r$phyloseq_filtered())

      print("deseqtophy")
      fun <- glue(" tmp <- subset_samples(r$phyloseq_filtered(), {input$Fact1} %in% c('{input$Cond1}','{input$Cond2}')) ")
      eval(parse(text=fun))
      # print("prunetaxa")
      tmp <- prune_taxa(taxa_sums(tmp) >= 1, tmp)
      tmp <- prune_samples(sample_sums(tmp) >=1, tmp)
      fun = glue ("deseq <- phyloseq_to_deseq2(tmp, ~ {input$Fact1})")
      eval(parse(text=fun))
      print(deseq)

      gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
      }

      print('geoMeans')
      geoMeans = apply(counts(deseq), 1, gm_mean)
      print(geoMeans)
      deseq = estimateSizeFactors(deseq, geoMeans = geoMeans)
      print("deseq")
      deseq = DESeq(deseq, test="Wald", fitType="parametric")

      print("res")
      res = results(deseq, cooksCutoff = FALSE, contrast = c(input$Fact1, input$Cond1 , input$Cond2))
      res

      }, message = "Performing DESeq2...")
    })


    output$deseqTab <- DT::renderDataTable({
      req(deseqDA(), input$Alpha1, r$phyloseq_filtered())

      res <- deseqDA()
      # Construct table
      ttable1 <- as.data.frame(r$phyloseq_filtered()@tax_table@.Data) %>%
        rownames_to_column()

      sseq1 <- as.data.frame(r$phyloseq_filtered()@refseq) %>%
        rownames_to_column()

      if(nrow(sseq1) != 0){sseq1 <- rename(sseq1, sequence = 2)}

      print(head(sseq1))

      resDESeq <- as.data.frame(res) %>%
        rownames_to_column() %>%
        mutate(absLFC = abs(log2FoldChange)) %>%
        # filter(padj<=input$Alpha1) %>%
        left_join(ttable1, by="rowname") %>%
        left_join(sseq1, by="rowname")

      print(head(resDESeq))

      as.data.frame(resDESeq)

    }, filter="top", options = list(scrollX = TRUE))


    #### METAGENOMESEQ
    mgSeqDA = eventReactive(input$go2, {

      withProgress({

        req(input$Fact1, input$Cond1, input$Cond2, r$phyloseq_filtered())

        fun <- glue(" tmp <- subset_samples(r$phyloseq_filtered(), {input$Fact1} %in% c('{input$Cond1}','{input$Cond2}')) ")
        eval(parse(text=fun))
        tmp <- prune_taxa(taxa_sums(tmp) >= 1, tmp)
        tmp <- prune_samples(sample_sums(tmp) >=1, tmp)
        # print(tmp)
        print(head(otu_table(tmp)))
        print("MGseqtophy")
        tax_table(tmp) <- NULL #problem with taxonomy table conversion
        MGdata <- phyloseq_to_metagenomeSeq(tmp)
        print(MGdata)

        print(head(rowSums(MGdata@assayData$counts)))

        featuresToKeep = which(rowSums(MGdata@assayData$counts) > 0)
        print(featuresToKeep)

        samplesToKeep = which(pData(MGdata)[,input$Fact1] == input$Cond1 | pData(MGdata)[,input$Fact1] == input$Cond2)
        print(samplesToKeep)
        obj_f = MGdata[featuresToKeep, samplesToKeep]

        #FitFeature model : zero-inflated log-normal model
        print('Fitzig Model')
        pd <- pData(obj_f)
        mod <- model.matrix(as.formula(paste("~", input$Fact1)), data = pd)

        res1 = NULL
        tryCatch( {res1 = fitFeatureModel(obj_f, mod)} ,
                  error=function(e){e;cat("ERROR :",conditionMessage(e), "\n")})

        TAB = MRcoefs(res1) #fdr adjustment
        # TAB = TAB[TAB$adjPvalues<=0.05,]

        print(head(TAB))

        TAB

      }, message="Performing metagenomeSeq...")

    })

    output$MGseqTab <- DT::renderDataTable({
      mgSeqDA()
      }, filter="top", options = list(scrollX = TRUE))


    ### MEtacoder
    mtcoderDA = eventReactive(input$go4, {

      withProgress({

      req(input$Fact1, input$Cond1, input$Cond2, r$phyloseq_filtered())

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

      fun <- glue(" psobj <- subset_samples(r$phyloseq_filtered(), {input$Fact1} %in% c('{input$Cond1}','{input$Cond2}')) ")
      eval(parse(text=fun))
      psobj <- prune_taxa(taxa_sums(psobj) >= 1, psobj)
      psobj <- prune_samples(sample_sums(psobj) >=1, psobj)


      normf = function(x, tot=max(sample_sums(psobj))){ tot*x/sum(x) }
      psobj <- transform_sample_counts(psobj, normf)

      print('Zeroing low counts...')
      obj <- parse_phyloseq(psobj, class_regex = "(.*)", class_key = "taxon_name")
      obj$data$otu_table <- zero_low_counts(obj, "otu_table", min_count = 1000, use_total = TRUE)
      no_reads <- rowSums(obj$data$otu_table[, obj$data$sample_data$sample_id]) == 0
      obj <- metacoder::filter_obs(obj, "otu_table", ! no_reads, drop_taxa = TRUE)
      if(nrow(obj$data$otu_table)==0){return(NULL)}

      print('Calculating taxon abundance...')
      obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table",  cols = obj$data$sample_data$sample_id)
      obj$data$tax_abund$total <- rowSums(obj$data$tax_abund[, -1]) # -1 = taxon_id column
      obj$data$n_samples <- calc_n_samples(obj,data="tax_abund")

      print('Comparing groups...')
      fun <- paste('obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = obj$data$sample_data$', input$Fact1, ',func = mean_ratio)', sep='')
      eval(parse(text=fun))

      table <- merge(obj$data$diff_table, obj$data$tax_data,by='taxon_id')
      table

      }, message = "Performing metacoder...")
    })

    output$mtcoderTab <- DT::renderDataTable({
      mtcoderDA()
    }, filter="top", options = list(scrollX = TRUE))




    ### Wilcox test M. Mariadassou

    wilcoxDA = eventReactive(input$go3, {

      withProgress({
        req(r$phyloseq_filtered(), input$Fact1, input$Cond1, input$Cond2)
        fun <- glue(" tmp <- subdata <- subset_samples(r$phyloseq_filtered(), {input$Fact1} %in% c('{input$Cond1}','{input$Cond2}')) ")
        eval(parse(text=fun))

        print("format")
        print(tmp)
        tax_table(tmp) <- NULL
        wilcoxon_data <- tmp %>%
          transform_sample_counts(function(x) {x / sum(x)}) %>%
          psmelt()


        print("perform test")
        log_fold_change <- function(data, variable = input$Fact1, cond1 = "Feces", cond2 = "Soil") {
          geom_mean <- function(x) {
            x <- x[x > 0]
            mean(log(x))
          }
          data %>%
            group_by(!!enquo(variable)) %>%
            summarize(geom_mean = geom_mean(Abundance)) %>%
            filter(!!enquo(variable) %in% c(cond1, cond2)) %>%
            pull(geom_mean) -> ratios
          if (is.na(ratios[1])) return(-10)
          if (is.na(ratios[2])) return(10)
          ratios[1] - ratios[2]
        }

        res <- wilcoxon_data %>%
          group_by(OTU) %>%
          nest() %>%
          mutate(
            wilcox_fit     = map(data, ~ wilcox.test(as.formula(glue("Abundance ~ {input$Fact1}")), data = .x)),
            log2FoldChange = purrr::map_dbl(data, log_fold_change),
            tidied         = map(wilcox_fit, tidy)
          ) %>%
          select(-data, -wilcox_fit) %>%
          unnest(tidied) %>%
          select(-alternative, -method) %>%
          rename(pvalue = p.value) %>%
          ungroup() %>%
          mutate(padj = p.adjust(pvalue, method = "fdr"))

        #Out
        LL=list()
        LL$res = res
        LL$subdata = subdata

        LL

      }, message="Performing Wilcox test")

    })

    output$WilcoxTab <- DT::renderDataTable({
      req(input$Alpha1)
      LL <-wilcoxDA()
      res = LL$res
      subdata = LL$subdata

      print("format res")
      wilcoxon_results <- res %>%
        rename(ASV = OTU) %>%
        filter(padj <= input$Alpha1) %>%
        inner_join(tax_table(subdata) %>% as("matrix") %>%
                     as_tibble() %>% mutate(ASV = taxa_names(subdata)),
                   by = "ASV") %>%
        mutate(ASV = forcats::fct_reorder(ASV, log2FoldChange))

      if(nrow(wilcoxon_results) != 0){
        as.data.frame(wilcoxon_results)
      }else(return(NULL))
      })


    output$mergePrint <- renderPrint({
      print("PRINT")
      # req(wilcoxDA(), mtcoderDA(), deseqDA(), mgSeqDA())
      print("MGSEQ")
      print(head(mgSeqDA()))

      # print("Others")
      # print(head(wilcoxDA()$res))

      print("metacoder")
      print(head(mtcoderDA()))
      print("deseq")
      print(head(deseqDA()))


    })




    mergeList <- reactive({

      withProgress({
      req(input$Alpha1, input$Cond1, input$Cond2, input$Fact1)
      print("merge")
      # req(wilcoxDA(), mtcoderDA(), deseqDA(), mgSeqDA(), input$Alpha1)

        # wTab = wilcoxDA()$res
        # wList = wTab[wTab$padj<=input$Alpha1, "OTU"]

        mtTab = mtcoderDA()
        mtList = mtTab[mtTab$wilcox_p_value <= input$Alpha1, "otu_id"]
        if(all(is.na(mtList))){
          mtList <- NULL
        }
        print("mtList")
        print(head(mtList))

        deTab = deseqDA()
        deTab = deTab[!is.na(deTab$padj),]
        deList = row.names(deTab[deTab$padj <= input$Alpha1,])
        if(all(is.na(deList))){
          deList <- NULL
        }

        mgTab = mgSeqDA()
        mgList = row.names(na.omit(mgTab[mgTab$adjPvalues <= input$Alpha1,]))
        if(all(is.na(mgList))){
          mgList <- NULL
        }

        TF = list(x1=deList, x2=mgList, x3=mtList)
        names(TF) <- c("DESeq", "metagenomeSeq", "metacoder")

        ListAllOtu = unique(unlist(TF))

        #Construction de la table
        print('Building table...')
        comp1 = paste(input$Cond1, '_vs_' , input$Cond2,sep='')
        col_comp = rep(comp1, length(ListAllOtu))
        TABf = cbind.data.frame(ListAllOtu, col_comp)


        # Test si chaque ASV est diff dans les méthdes.
        print('Check methods ...')
        for (j in 1:length(TF)){
          TABtest = TF[[j]]
          # TABtest=gsub("\\[|\\]", "", TF[[j]]) # cherche les ASVids

          TABtest_signif=rep(0, length(ListAllOtu))
          names(TABtest_signif) = ListAllOtu
          TABtest_signif[TABtest] = 1

          TABf=cbind.data.frame(TABf, TABtest_signif)
          names(TABf)[ncol(TABf)] = names(TF)[j]
        }

        TABfbak = TABf

        print('Add LFC...')
        TABf <- cbind( TABf, sumMethods = apply(TABf[3:5], 1, sum, na.rm=TRUE),
                       DESeqLFC = deTab[as.character(TABf[,1]),"log2FoldChange"],
                       absDESeqLFC = abs(deTab[as.character(TABf[,1]),"log2FoldChange"]) )

        # input$Fact1="SampleType"
        # input$Cond1="Feces"
        # input$Cond2="Soil"

        print('Calculating mean relative abundance...')
        data = r$phyloseq_filtered()
        normf = function(x){ x/sum(x) }
        data.norm <- transform_sample_counts(data, normf)
        otableNORM <- otu_table(data.norm)
        ssample <- as.matrix(sample_data(data.norm))
        ttax <- tax_table(data.norm)
        head(ttax)

        seqs = NULL
        if(!is.null(refseq(data.norm, errorIfNULL=FALSE)) ){
          seqs <- refseq(data.norm)
        }

        print("mean1")
        Gtab <- cbind(as.data.frame(ssample), t(otableNORM))
        MeanRelAbcond1 = NULL
        for(i in TABf$ListAllOtu){
          tt=mean(Gtab[Gtab[,input$Fact1]==input$Cond1,i], na.rm=TRUE)
          MeanRelAbcond1=c(MeanRelAbcond1,tt)
        }
        print("mean2")
        MeanRelAbcond2=NULL
        for(i in TABf$ListAllOtu){
          tt=mean(Gtab[Gtab[,input$Fact1]==input$Cond2,i], na.rm=TRUE)
          MeanRelAbcond2=c(MeanRelAbcond2,tt)
        }
        TABfbak <- TABf <- cbind(TABf, MeanRelAbcond1, MeanRelAbcond2)

        #Adjust table
        print('Adjusting table...')
        TABf <- TABf[!is.na(TABf$DESeqLFC),]
        TABf$Condition = rep(NA, nrow(TABf))
        TABf[TABf$DESeqLFC>0, "Condition"] = as.character(input$Cond1)
        TABf[TABf$DESeqLFC<0, "Condition"] = as.character(input$Cond2)
        TABf$Condition = factor(TABf$Condition,
                                levels=c(as.character(input$Cond1),as.character(input$Cond2)) )

        print("Adding taxonomy and sequences...")
        # Debug
        # LL = list()
        # LL$TABf = TABf; LL$ttax =ttax[as.character(TABf[,1]),] ;  LL$seq = seqs[as.character(TABf[,1])]
        # save(LL, file = "~/Téléchargements/debug_explore.rdata")

        TABf <- cbind.data.frame(TABf, ttax[as.character(TABf[,1]),]@.Data)
        if(!is.null(seqs)){TABf <- cbind.data.frame(TABf, sequences = seqs[as.character(TABf[,1])])}
        print("done")


        LL = list()
        LL$TABf = TABf
        LL$ttax = ttax
        LL$signifseqs = seqs[as.character(TABf[,1])]
        LL
      }, message = "Merging results...")

    })









    output$mergeTab <- DT::renderDataTable({
      mergeList()$TABf
      }, filter="top", options = list(scrollX = TRUE))

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
      print("plot")
      req(mergeList(), input$Nmeth, input$minAb, input$Nfeat)
      TABf = mergeList()$TABf
      ttax = mergeList()$ttax
      # Barplot
      ## differentialy abundant features on 2 methods
      ## Top abs(LogFolchange)
      ## feature with relative abondance > 0.1% (0.001)
      print("Filters...")
      print(input$minAb)
      print(input$Nfeat)
      TABbar = TABf[TABf$sumMethods >= input$Nmeth, ] #TABf$DESeq ==1 | TABf$metagenomeSeq ==1 &
      TABbar = TABbar[TABbar$MeanRelAbcond1 >= input$minAb | TABbar$MeanRelAbcond2 >= input$minAb, ]    # min mean Abundance to choose
      TABbar = tail(TABbar[order(abs(TABbar$DESeqLFC)),], input$Nfeat)      # number of features to plot

      print(TABbar)


      if(nrow(TABbar)!=0){
        print("Plot")
        print(r$rank_glom())
        if(r$rank_glom() == "ASV"){
          TABbar$tax = paste( substr(ttax[as.character(TABbar$ListAllOtu),"Species"],1,20),"...","_", TABbar$ListAllOtu, sep="")
        }else{
          TABbar$tax = TABbar$ListAllOtu
        }

          p<-ggplot2::ggplot(data=TABbar, aes(x=reorder(tax, -abs(DESeqLFC)), y=DESeqLFC, fill=Condition ) ) +
          geom_bar(stat="identity", alpha = 0.7) + ggtitle(glue("{input$Cond1} vs. {input$Cond2}")) + labs(x='Features') +
          coord_flip() + theme_bw() +
          scale_y_continuous(minor_breaks = seq(-1E4 , 1E4, 1), breaks = seq(-1E4, 1E4, 5))
        plotly::ggplotly(p)

      }else{
        print("No ASV to plot")
        showNotification("No ASV to plot... (all features show summethods <=2 ? too high mean relative abundance filters ?)", type="error", duration = 10)
        }
    })

    output$barplot1 <- renderPlotly({
        reacbarplot1()
    })


    # output$mergePlot <- renderPlot({
    #   TF = mergeList()
    #   venn.plot <- venn.diagram(TF, filename = NULL, col = "black",
    #                             fill = rainbow(length(TF)), alpha = 0.50,
    #                             cex = 2, cat.col = 1, lty = "blank",
    #                             cat.cex = 2.5, cat.fontface = "bold",
    #                             margin = 0.07, main="test", main.cex=2.5,
    #                             fontfamily ="Arial",main.fontfamily="Arial",cat.fontfamily="Arial");
    #
    #   grid.draw(venn.plot)
    #
    # })

}

## To be copied in the UI
# mod_diffanalysis_ui("diffanalysis_ui_1")

## To be copied in the server
# callModule(mod_diffanalysis_server, "diffanalysis_ui_1")
