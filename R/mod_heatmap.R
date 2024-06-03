#' heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import phyloseq
#' @import ggplot2
#' @import shinycustomloader
#' @importFrom DT renderDataTable
#' @import pheatmap
#' @importFrom indicspecies multipatt
#' @import DESeq2
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly ggplotly

mod_heatmap_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        box(title = "Settings", width = 6, status = "warning", solidHeader = TRUE,
            selectInput(
              ns("Rank"),
              label = "Rank to agglomerate",
              choices = ""
            ),
            selectInput(
              ns("norm"),
              label = "Normalization",
              choices = list(
                "Raw" = 0 ,
                "TSS (total-sum normalization)" = 1,
                "CLR (center log-ration)" = 2,
                "VST (variance stabilizing transformation)" = 3,
                "hellinger" = 4,
                "log10" = 5),
              selected = 0
            ),
            selectInput(
              ns("sample_label"),
              label = "Sample label",
              choices = ""
            ),
            checkboxInput(ns("select_features"), label = "Features selection", value = FALSE),
            checkboxInput(ns("clust_taxa"), label = "Taxa clustering", value = FALSE),
            checkboxInput(ns("clust_samp"), label = "Sample clustering", value = FALSE),
            selectInput(
              ns("fact_annot"),
              label = "Factor annotation",
              choices = "",
              multiple = TRUE
            ),
            selectInput(
              ns("taxa_annot"),
              label = "Taxa annotation",
              choices = "",
              multiple = TRUE
            ),
            actionButton(ns("launch_heatmap"), "Launch heatmap", icon = icon("play-circle"),
                         style="color: #fff; background-color: #3b9ef5; border-color: #1a4469")
        ),
        box(title = "Display settings", width = 6, status = "warning", solidHeader = TRUE,
            sliderInput(ns('plot_height'), label = 'Plot Height', min = 300, max = 2000, step = 100, value = 300),
            sliderInput(ns('plot_width'), label = 'Plot Width', min = 300, max = 2000, step = 100, value= 600),
            checkboxInput(ns("print_taxa"), label = "Taxa labels", value = TRUE),
            checkboxInput(ns("print_sample"), label = "Sample labels", value = TRUE),
            checkboxInput(ns("print_nb"), label = "Frequencies", value = FALSE),
            selectInput(
              ns("color_map"),
              label = "Heatmap color palette",
              choices = ""
            )
        ),
        uiOutput(ns("ui_sample_clustering")),
        uiOutput(ns("ui_select_features")),
        uiOutput(ns("ui_factor_annot_file")),
        uiOutput(ns("ui_taxa_annot_file")),
        uiOutput(ns("ui_fact_colors")),
        uiOutput(ns("ui_taxa_colors"))
      ),
      fluidRow(
        uiOutput(ns("ui_box_heatmap"))
      ),
      fluidRow(
        box(title = "Selected features", width = 12, status = "primary", solidHeader = TRUE,
            downloadButton(ns("table_download"), label = "Download table"),
            shinycustomloader::withLoader(
              DT::dataTableOutput(ns("feat")),
              type = "html", loader = "loader2"
            )
        )
      ),
    ),
  )
}


#' heatmap Server Function
#'
#' @noRd
mod_heatmap_server <- function(input, output, session, r){
  ns <- session$ns
  
  observe({
    req(r$phyloseq_filtered(), r$phyloseq_filtered_norm())
    updateSelectInput(session, "sample_label",
                      choices = names(sample_data(r$phyloseq_filtered())))
    updateSelectInput(session, "fact_annot",
                      choices = r$var_list())
    ranks <- phyloseq::rank_names(r$phyloseq_filtered())
    updateSelectInput(session, "Rank",
                      choices = ranks,
                      selected = ranks[length(ranks)])
    color_choices <- RColorBrewer::brewer.pal.info
    updateSelectInput(session, "color_map",
                      choices = rownames(color_choices[color_choices$category != "qual" & color_choices$colorblind == TRUE , ]),
                      selected = "RdYlBu")
  })
  
  observe({
    req(r$phyloseq_filtered(), input$Rank)
    ranks <- phyloseq::rank_names(r$phyloseq_filtered())
    updateSelectInput(session, "taxa_annot",
                      choices = ranks[1:which(ranks == input$Rank)])
  })
  
  output$ui_box_heatmap <- renderUI({
    box(title = "Heatmap", width = 12, status = "primary", solidHeader = TRUE, height = plot_height() + 100,
        downloadButton(ns("heatmap_download"), label = "Download plot"),
        shinycustomloader::withLoader(
          plotOutput(ns('heatmap_t')),
          type = "html", loader = "loader1"
        )
    )
  })
  
  output$ui_sample_clustering <- renderUI({
    req(input$clust_samp)
    box(title = "Sample clustering", width = 6, status = "warning", solidHeader = TRUE,
        selectInput(
          ns("dist"),
          label = "Distance method",
          choices = c("euclidean", "bray", "jaccard", "dpcoa", "unifrac", "wunifrac"),
          selected = "euclidean"
        ),
        selectInput(
          ns("clust_method"),
          label = "Clustering method",
          choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
          selected = "ward.D2"
        )
    )
  })
  
  output$ui_select_features <- renderUI({
    req(input$select_features)
    box(title = "Features selection", width = 6, status = "warning", solidHeader = TRUE,
        radioButtons(ns("selection_method"),
                     label = "Selection method",
                     inline = TRUE,
                     choices = c("indicspecies", "DESeq2"),
                     selected = "indicspecies"),
        selectInput(ns("test_fact"),
                    label = "Factor to test",
                    choices = names(sample_data(r$phyloseq_filtered())),
                    selected = "Lot"
        ),
        fluidRow(
        column(width = 6, uiOutput(ns("ui_test_cond1"))),
        column(width = 6, uiOutput(ns("ui_test_cond2")))
        ),
        numericInput(ns("pval"),
                     label = "p-value threshold",
                     min = 0, max = 1,
                     value = 0.05
        )
    )
  })
  
  output$ui_test_cond1 <- renderUI({
    if(input$selection_method == "DESeq2"){
      selectInput(ns("cond1"),
                  label = "Condition 1 to compare",
                  choices = unique(sample_data(r$phyloseq_filtered())[, input$test_fact]))
    }
  })
  
  output$ui_test_cond2 <- renderUI({
    if(input$selection_method == "DESeq2"){
      choi <- data.frame(unique(sample_data(r$phyloseq_filtered())[, input$test_fact]))
      choi <- choi[choi != input$cond1]
      selectInput(ns("cond2"),
                  label = "Condition 2 to compare",
                  choices = choi)
    }
  })
  
  plot_height <- reactive({
    return(input$plot_height)
  })
  
  plot_width <- reactive({
    return(input$plot_width)
  })
  
  # agglomerate taxa at chosen taxonomic rank
  agglom_data <- reactive({
    req(r$phyloseq_filtered(), input$Rank)
    data_glom <- phyloseq::tax_glom(r$phyloseq_filtered(), input$Rank)
    taxa_names(data_glom) <- tax_table(data_glom)[, input$Rank]
    tax_table(data_glom) <- tax_table(data_glom)[, 1:match(input$Rank, rank_names(data_glom))]
    return(data_glom)
  })
  
  # agglomerate taxa at chosen taxonomic rank then normalize data with the chosen method
  agglom_normalized_data <- reactive({
    req(r$phyloseq_filtered(), input$Rank, input$norm)
    FGdata <- agglom_data()
    
    if(input$norm == 0){
      FNGdata <- FGdata
    }
    
    if(input$norm == 1){
      normf = function(x){ x/sum(x) }
      FNGdata <- transform_sample_counts(FGdata, normf)
    }
    
    if(input$norm == 2){
      clr = function(x){log(x+1) - rowMeans(log(x+1))}
      otable <- otu_table(FGdata)
      otableCLR <- clr(otable)
      FNGdata <- FGdata
      otu_table(FNGdata) <- otableCLR
    }

    if(input$norm == 3){
      withProgress({
        otable <- as(otu_table(FGdata) + 1, "matrix")
        otableVST <- DESeq2::varianceStabilizingTransformation(otable, fitType = 'local')
        FNGdata <- FGdata
        otu_table(FNGdata) <- otu_table(otableVST, taxa_are_rows = TRUE)
      }, message = "VST normalization, please wait...")
    }
    
    if(input$norm == 4){
      spe <- veganifyOTU(FGdata)
      spe.hell <- vegan::decostand(spe, method = 'hell')
      FNGdata <- FGdata
      otu_table(FNGdata) <- otu_table(t(spe.hell), taxa_are_rows = TRUE)
    }
    
    if(input$norm == 5){
      log = function(x){log10(x+1)}
      otable <- otu_table(FGdata)
      otableLOG <- log(otable)
      FNGdata <- FGdata
      otu_table(FNGdata) <- otableLOG
    }
    return(FNGdata)
  })
  
  # use indicspecies package to select features
  indic_species_results <- reactive({
    req(agglom_data(), input$select_features, input$test_fact)
    return(indicspecies::multipatt(t(otu_table(agglom_data())), t(sample_data(agglom_data())[, input$test_fact]), control = permute::how(nperm = 999), duleg = TRUE))
  })
  
  # use DESeq2 package to select features
  deseq2_results <- reactive({
    req(agglom_data(), input$selection_method, input$test_fact, input$cond1, input$cond2)
    fun <- glue("deseq <- phyloseq_to_deseq2(agglom_data(), ~ {input$test_fact})")
    # fun <- glue("deseq <- phyloseq_to_deseq2(r$phyloseq_filtered(), ~ {input$test_fact})")
    eval(parse(text = fun))
    gm_mean <- function(x, na.rm = TRUE){
      exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
    }
    geoMeans <- apply(counts(deseq), 1, gm_mean)
    deseq <- DESeq2::estimateSizeFactors(deseq, geoMeans = geoMeans)
    deseq <- DESeq2::DESeq(deseq, test = "Wald", fitType = "parametric")
    res <-  DESeq2::results(deseq, cooksCutoff = FALSE, contrast = c(input$test_fact, input$cond1 , input$cond2))
    return(res)
  })
  
  selected_data <- reactive({
    req(agglom_normalized_data())
    if(input$select_features == FALSE){
      return(agglom_normalized_data())
    }else if(input$selection_method == "indicspecies"){
      slct <- indic_species_results()$sign
      slct <- slct[slct$p.value <= input$pval, ]
      return(phyloseq::prune_taxa(rownames(slct), agglom_normalized_data()))
    }else{
      return(phyloseq::prune_taxa(rownames(deseq2_results()[deseq2_results()$pvalue <= input$pval, ]), agglom_normalized_data()))
    }
  })

  
  # taxa annotation colors
  taxa_colors <- reactive({
    req(input$taxa_annot, selected_data())
    colors <- get_cat_colors(type = "taxa", list_fact = input$taxa_annot, phy_object = selected_data())
    return(colors)
  })
  
  annot_taxa_colors <- reactive({
    req(input$taxa_annot, r$taxa_colors(), taxa_colors())
    if(length(input$taxa_annot) == 0){
      annot <- NA
    }else{
      annot_colors <- r$taxa_colors()[input$taxa_annot]
      annot <- lapply(1:length(input$taxa_annot), FUN = function(i){
        annot_colors[[i]][lapply(taxa_colors()[[i]], FUN = as.character)[[1]]]
      })
      names(annot) <- input$taxa_annot
    }
    return(annot)
  })
  
  
  # factor annotation colors
  modality_colors <- reactive({
    req(input$fact_annot, selected_data())
    get_cat_colors(type = "fact", list_fact = input$fact_annot, phy_object = selected_data())
  })

  annot_fact_colors <- reactive({
    req(input$fact_annot, r$factor_colors(), modality_colors())
    if(length(input$fact_annot) == 0){
      annot <- NA
    }else{
      annot_colors <- r$factor_colors()[input$fact_annot]
      annot <- lapply(1:length(input$fact_annot), FUN = function(i){
        annot_colors[[i]][lapply(modality_colors()[[i]], FUN = as.character)[[1]]]
      })
      names(annot) <- input$fact_annot
    }
    return(annot)
  })
  
  
  annot_colors <- reactive({
    if(length(input$fact_annot) >= 1 && length(input$taxa_annot) >= 1){
      return(c(annot_fact_colors(), annot_taxa_colors()))
    }else if(length(input$fact_annot) >= 1){
      return(annot_fact_colors())
    }else if(length(input$taxa_annot) >= 1){
      return(annot_taxa_colors())
    }else{
      return(NA)
    }
  })
  
  heatmap <- eventReactive(input$launch_heatmap, {
    req(selected_data(), input$Rank)
    t_table <- tax_table(selected_data())
    samp_table <- sample_data(selected_data())
    
    if(length(input$fact_annot) == 0){
      annot_col <- NA
    }else{
      annot_col <- data.frame(samp_table[, input$fact_annot])
    }
    
    if(length(input$taxa_annot) == 0){
      annot_row <- NA
    }else{
      annot_row <- data.frame(t_table[, input$taxa_annot])
    }
    
    if(input$clust_samp == FALSE){
      clust_sample <- FALSE
    }else{
      dist_samp <- phyloseq::distance(selected_data(), method = input$dist, type = "sample")
      clust_sample <- hclust(dist_samp, method = input$clust_method)
    }
    
    heatmap <- pheatmap::pheatmap(otu_table(selected_data()),
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = input$color_map)))(100),
                       cluster_cols = clust_sample,
                       cluster_rows = input$clust_taxa,
                       show_rownames = input$print_taxa,
                       show_colnames = input$print_sample,
                       labels_row = as.expression(lapply(
                         stringr::str_replace(gsub("[a-z]__", "", t_table[, input$Rank]), "_", " "),
                         function(x) bquote(italic(.(x))))),
                       labels_col = sapply(samp_table[, input$sample_label], FUN = as.character),
                       angle_col = 90,
                       annotation_col = annot_col,
                       annotation_row = annot_row,
                       annotation_names_col = TRUE,
                       annotation_names_row = TRUE,
                       annotation_colors = annot_colors(),
                       display_numbers = input$print_nb,
                       number_format = "%.2f"
    )
    return(heatmap)
  })
  
  selection_features_results <- reactive({
    req(input$select_features, input$selection_method)
    if(input$selection_method == "indicspecies"){
      result <- indic_species_results()$sign
    }else if(input$selection_method == "DESeq2"){
      result <- data.frame(deseq2_results())
    }
    t_table <- as.data.frame(phyloseq::tax_table(agglom_data())) %>%
      rownames_to_column()
    res <- result %>%
      rownames_to_column() %>%
      left_join(t_table, by = "rowname")
    return(res)
  })
  
  output$feat <- DT::renderDataTable({
    req(selection_features_results())
    selection_features_results()
  }, filter = "top", options = list(scrollX = TRUE))
  
  observe({
    output$heatmap_t <- renderPlot({
      withProgress(message = 'Computing heatmap...',{
        heatmap()
      })
    }, height = plot_height(), width = plot_width())
  })
  
  output$heatmap_download <- downloadHandler(
    filename = "heatmap.svg",
    content = function(file){
      req(heatmap())
      svg(filename = file)
      print(heatmap())
      dev.off()
    }
  )
  
  output$table_download <- downloadHandler(
    filename = "heatmap_table.csv",
    content = function(file){
      req(selection_features_results())
      write.table(selection_features_results(), file, sep = ",", row.names = FALSE)
    }
  )
}

## To be copied in the UI
# mod_heatmap_ui("heatmap_ui_1")

## To be copied in the server
# callModule(mod_heatmap_server, "heatmap_ui_1")
