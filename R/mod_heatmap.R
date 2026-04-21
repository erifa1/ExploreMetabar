#' heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom DT dataTableOutput renderDataTable
#' 
#' 
#' @import phyloseq
#' @import ggplot2
#' @import shinycustomloader
#' @importFrom DT renderDataTable
#' @import pheatmap
#' @importFrom indicspecies multipatt
#' @import DESeq2
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly ggplotly
#' @importFrom bslib layout_sidebar sidebar accordion accordion_panel navset_card_underline nav_panel input_task_button tooltip
#' @importFrom bsicons bs_icon

mod_heatmap_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "Heatmap configuration",
      open = "desktop",
      width = "350px",

      accordion(
        id = ns("config_accordion"),
        open = "Heatmap Setup",
        multiple = TRUE,

        # ── Panel 1: Core heatmap setup ──
        accordion_panel(
          "Heatmap Setup",
          icon = bs_icon("gear"),
          tooltip(
            selectInput(ns("rank"), label = "Rank to agglomerate", choices = ""),
            "Taxonomic rank used to agglomerate ASV counts before computing the heatmap.",
            placement = "right"
          ),
          tooltip(
            selectInput(ns("norm"), label = "Normalization",
              choices = list(
                "Raw" = 0,
                "TSS (total-sum normalization)" = 1,
                "CLR (centered log-ratio)" = 2,
                "VST (variance stabilizing transformation)" = 3,
                "hellinger" = 4,
                "log10" = 5
              ),
              selected = 0
            ),
            "Method used to normalize/transform abundance values displayed in the heatmap.",
            placement = "right"
          ),
          tooltip(
            selectInput(ns("sample_label"), label = "Sample label", choices = ""),
            "Metadata column used as label for samples (columns of the heatmap).",
            placement = "right"
          )
        ),

        # ── Panel 2: Features selection ──
        accordion_panel(
          "Features Selection",
          icon = bs_icon("funnel"),
          tooltip(
            checkboxInput(ns("select_features"), label = "Enable features selection", value = FALSE),
            "Restrict the heatmap to a subset of taxa selected by indicator species analysis or by abundance.",
            placement = "right"
          ),
          uiOutput(ns("ui_select_features"))
        ),

        # ── Panel 3: Clustering (samples + taxa) ──
        accordion_panel(
          "Clustering",
          icon = bs_icon("diagram-3"),
          tooltip(
            checkboxInput(ns("clust_taxa"), label = "Taxa clustering", value = FALSE),
            "Cluster taxa (rows of the heatmap) based on the displayed values.",
            placement = "right"
          ),
          tooltip(
            checkboxInput(ns("clust_samp"), label = "Sample clustering", value = FALSE),
            "Cluster samples (columns of the heatmap) using a chosen distance and linkage method.",
            placement = "right"
          ),
          uiOutput(ns("ui_sample_clustering"))
        ),

        # ── Panel 4: Annotations & Display ──
        accordion_panel(
          "Annotations & Display",
          icon = bs_icon("brush"),
          tooltip(
            selectInput(ns("fact_annot"), label = "Factor annotation", choices = "", multiple = TRUE),
            "Metadata variables to annotate samples (top of the heatmap).",
            placement = "right"
          ),
          tooltip(
            selectInput(ns("taxa_annot"), label = "Taxa annotation", choices = "", multiple = TRUE),
            "Taxonomic ranks to annotate taxa (left side of the heatmap).",
            placement = "right"
          ),
          tooltip(
            selectInput(ns("color_map"), label = "Heatmap color palette", choices = ""),
            "Color palette used for the heatmap cells (RColorBrewer).",
            placement = "right"
          ),
          tags$hr(),
          tooltip(
            sliderInput(ns('plot_height'), label = 'Plot Height', min = 300, max = 2000, step = 100, value = 300),
            "Adjust the rendered plot height (pixels).",
            placement = "right"
          ),
          tooltip(
            sliderInput(ns('plot_width'), label = 'Plot Width', min = 300, max = 2000, step = 100, value = 600),
            "Adjust the rendered plot width (pixels).",
            placement = "right"
          ),
          tooltip(
            checkboxInput(ns("print_taxa"), label = "Show taxa labels", value = TRUE),
            "Display taxa names on the rows of the heatmap.",
            placement = "right"
          ),
          tooltip(
            checkboxInput(ns("print_sample"), label = "Show sample labels", value = TRUE),
            "Display sample names on the columns of the heatmap.",
            placement = "right"
          ),
          tooltip(
            checkboxInput(ns("print_nb"), label = "Show frequencies", value = FALSE),
            "Print numerical values inside each cell of the heatmap.",
            placement = "right"
          )
        )
      ),

      tags$hr(),

      # Always-visible launch button placed AFTER the accordion
      tooltip(
        input_task_button(
          ns("launch_heatmap"),
          label = "Launch heatmap",
          icon = bs_icon("play-fill"),
          class = "btn-primary w-100",
          label_busy = "Computing..."
        ),
        "Compute the heatmap with the current settings.",
        placement = "top"
      )
    ),

    # ── Main content area ──
    navset_card_underline(
      id = ns("main_tabs"),
      full_screen = TRUE,
      title = "Heatmap results",

      nav_panel(
        title = "Heatmap",
        icon = bs_icon("grid-3x3"),
        downloadButton(ns("heatmap_download"), label = "Download plot"),
        plotOutput(ns('heatmap_t'))
      ),

      nav_panel(
        title = "Selected Features",
        icon = bs_icon("table"),
        uiOutput(ns("ui_selected_features"))
      )
    )
  )
}


#' heatmap Server Function
#'
#' @noRd
#' @import phyloseq
#' @importFrom pheatmap pheatmap
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @import vegan
#' @importFrom indicspecies multipatt
#' @import shinyalert
#' @importFrom grDevices svg
mod_heatmap_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns
  
  observe({
    req(r$phyloseq_filtered(), r$phyloseq_filtered_norm())
    num <- sapply(r$sdat()[, names(r$sdat())], is.numeric)
    samp_lab <- names(r$sdat())[!num]
    updateSelectInput(session, "sample_label",
                      choices = samp_lab)
    updateSelectInput(session, "fact_annot",
                      choices = r$var_list())
    ranks <- phyloseq::rank_names(r$phyloseq_filtered())
    if(r$rank_glom() == "ASV"){
      ranks <- c(ranks, "ASV")
    }
    updateSelectInput(session, "rank",
                      choices = ranks)
    color_choices <- RColorBrewer::brewer.pal.info
    updateSelectInput(session, "color_map",
                      choices = rownames(color_choices[color_choices$category != "qual" & color_choices$colorblind == TRUE , ]),
                      selected = "RdYlBu")
  })
  
  observe({
    req(r$phyloseq_filtered(), r$phyloseq_filtered_norm(), input$rank)
    ranks <- phyloseq::rank_names(r$phyloseq_filtered())
    if(input$rank == "ASV"){
      updateSelectInput(session, "taxa_annot",
                        choices = ranks)
    }else{
      updateSelectInput(session, "taxa_annot",
                        choices = ranks[1:which(ranks == input$rank)])
    }
  })
  
  output$ui_sample_clustering <- renderUI({
    req(input$clust_samp)
    if(is.null(phyloseq::phy_tree(r$phyloseq_filtered(), errorIfNULL = FALSE))){
      choice = list("euclidean", "bray", "jaccard")
    }else{
      choice = list("euclidean", "bray", "jaccard", "unifrac", "wunifrac", "dpcoa")
    }
    tagList(
      selectInput(
        ns("dist_method"),
        label = "Distance method",
        choices = choice,
        selected = "bray"
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
    tagList(
      radioButtons(ns("selection_method"),
                   label = "Selection method",
                   inline = TRUE,
                   choices = c("indicspecies", "abundance"),
                   selected = "indicspecies"),
      uiOutput(ns("ui_test_fact")),
      uiOutput(ns("ui_pval_indicspecies")),
      uiOutput(ns("ui_type_features"))
    )
  })
  
  output$ui_test_fact <- renderUI({
    req(input$selection_method)
    if(input$selection_method == "indicspecies"){
        selectInput(ns("test_fact"),
                    label = "Factor to test",
                    choices = r$factor_list()
        )
    }
  })
  
  output$ui_pval_indicspecies <- renderUI({
    req(input$selection_method)
    if(input$selection_method == "indicspecies"){
      numericInput(ns("pval"),
                   label = "p-value threshold",
                   min = 0, max = 1,
                   value = 0.05
      )
    }
  })
  
  output$ui_type_features <- renderUI({
    req(input$selection_method)
    if(input$selection_method == "abundance"){
      fluidRow(
        column(width = 8,
               radioButtons(ns("type_features"),
                   label = paste0("Display features"),
                   inline = FALSE,
                   choices = list(
                     "the most abundant" = "top",
                     "the least abundant" = "bottom"),
                   selected = "top")
               ),
        column(width = 4,
               numericInput(ns("nb_feat"),
                            label = "Number of features",
                            min = 0, max = dim(otu_table(agglom_data()))[1],
                            value = dim(otu_table(agglom_data()))[1]
               ))
      )
      
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
    req(r$phyloseq_filtered(), input$rank)
    if(input$rank != "ASV"){
      data_glom <- speedyseq::tax_glom(r$phyloseq_filtered(), input$rank)
      taxa_names(data_glom) <- tax_table(data_glom)[, input$rank]
      tax_table(data_glom) <- tax_table(data_glom)[, 1:match(input$rank, rank_names(data_glom))]
    }else{
      data_glom <- r$phyloseq_filtered()
    }
    return(data_glom)
  })
  
  # normalize data with the chosen method
  agglom_normalized_data <- reactive({
    req(r$phyloseq_filtered(), input$norm)
    FGdata <- selected_data()
    
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
  
  observeEvent(input$test_fact, {
    req(r$phyloseq_filtered, r$phyloseq_filtered_norm(), input$test_fact, input$select_features)
    if(anyNA(sample_data(r$phyloseq_filtered())[, input$test_fact])){
      shinyalert::shinyalert(title = "Oops", text = paste0("The factor ", input$test_fact, " has missing values. Features selection with indicspecies will not work. Please remove missing values."), type = "error")
      req(FALSE)
    }
  })
  
  # use indicspecies package to select features
  indic_species_results <- reactive({
    req(agglom_data(), input$select_features, input$test_fact)
    result <- indicspecies::multipatt(t(otu_table(agglom_data())), t(sample_data(agglom_data())[, input$test_fact]), control = permute::how(nperm = 999), duleg = TRUE)
    return(result)
  })
  
  selected_data <- reactive({
    req(agglom_data())
    if(input$select_features == FALSE){
      selected_data <- agglom_data()
    }else if(input$selection_method == "indicspecies"){
      slct <- indic_species_results()$sign
      slct <- slct[slct$p.value <= input$pval, ]
      if(dim(slct)[1] == 0){
        shinyalert::shinyalert(title = "Oops", text = paste0("No taxa have been selected by indicspecies. Taxonomic rank chosen to agglomerate taxa may be too high.\n"), type = "error")
        req(FALSE)
      }
      selected_data <- phyloseq::prune_taxa(rownames(slct), agglom_data())
    }else{
      sum_taxa <- phyloseq::taxa_sums(agglom_data())
      decreasing_order <- ifelse(input$type_features == "top", TRUE, FALSE)
      sort_taxa <- sum_taxa[order(sum_taxa, decreasing = decreasing_order)]
      sort_taxa <- head(sort_taxa, n = input$nb_feat)
      selected_data <- phyloseq::prune_taxa(names(sort_taxa), agglom_data())
    }
    return(selected_data)
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
    colors <- get_cat_colors(type = "fact", list_fact = input$fact_annot, phy_object = selected_data())
    return(colors)
  })

  annot_fact_colors <- reactive({
    req(input$fact_annot, r$factor_colors(), modality_colors(), agglom_normalized_data())
    if(length(input$fact_annot) == 0){
      annot <- NA
    }else{
      annot <- list()
      var_nb_na <- list()
      for(i in 1:length(input$fact_annot)){
        colors <- r$factor_colors()[input$fact_annot[i]]
        modality <- lapply(modality_colors()[[i]], FUN = as.character)[[1]]
        if(is.numeric(r$sdat()[, input$fact_annot[i]])){
          result_colors <- paletteer::paletteer_c(colors[[1]], n = 30)
          if(NA %in% modality){
            var_values <- sample_data(agglom_normalized_data())[, input$fact_annot[i]]
            var_nb_na[[input$fact_annot[i]]] <- sum(is.na(var_values))
          }
        }else{
          result_colors <- colors[[1]][modality[!is.na(modality)]]
          if(NA %in% modality){
            result_colors["NA"] <- "#000000"
          }
        }
        annot[[input$fact_annot[i]]] <- result_colors
      }
      if(length(var_nb_na) > 0){
        nb_na <- paste(sapply(1:length(var_nb_na), FUN = function(i){
          paste0(names(var_nb_na)[i], " (", var_nb_na[[i]][1], " NA)")
        }), collapse = "\n")
        shinyalert::shinyalert(title = "NA values in heatmap annotation", text = paste0("The following numeric variables have NA values which will be displayed in blank on heatmap annotation.\n", nb_na), type = "warning")
      }
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
    req(selected_data(), input$rank, agglom_normalized_data())
    withProgress(message = 'Computing heatmap...',{
      t_table <- tax_table(agglom_normalized_data())
      samp_table <- sample_data(agglom_normalized_data())
      
      if(input$rank != "ASV"){
        row_labels <- as.expression(lapply(
          stringr::str_replace(gsub("[a-z]__", "", t_table[, input$rank]), "_", " "),
          function(x) bquote(italic(.(x)))))
      }else{
        row_labels <- NULL
      }
      
      if(length(input$fact_annot) == 0){
        annot_col <- NA
      }else{
        num <- sapply(r$sdat()[, input$fact_annot], is.numeric)
        annot_col <- data.frame(samp_table[, input$fact_annot])
        for(i in length(input$fact_annot)){
          if(!is.numeric(r$sdat()[, input$fact_annot[i]])){
            ann <- sapply(annot_col[, input$fact_annot[i]], FUN = as.character)
            ann <- data.frame(replace(ann, list = which(is.na(ann)), "NA"))
            annot_col[input$fact_annot[i]] <- ann
          }
        }
      }
      
      if(length(input$taxa_annot) == 0){
        annot_row <- NA
      }else{
        annot_row <- data.frame(t_table[, input$taxa_annot])
      }
      
      if(input$clust_samp == FALSE){
        clust_sample <- FALSE
      }else{
        no_empty_samples <- colnames(otu_table(selected_data())[, colSums(otu_table(selected_data())) > 0])
        select_data <- phyloseq::prune_samples(no_empty_samples, selected_data())
        dist_samp <- phyloseq::distance(select_data, method = input$dist_method, type = "sample")
        clust_sample <- stats::hclust(dist_samp, method = input$clust_method)
      }
      
      heatmap <- pheatmap::pheatmap(otu_table(agglom_normalized_data()),
                                    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = input$color_map)))(100),
                                    cluster_cols = clust_sample,
                                    cluster_rows = input$clust_taxa,
                                    show_rownames = input$print_taxa,
                                    show_colnames = input$print_sample,
                                    labels_row = row_labels,
                                    labels_col = sapply(samp_table[, input$sample_label], FUN = as.character),
                                    angle_col = 90,
                                    annotation_col = annot_col,
                                    annotation_row = annot_row,
                                    annotation_names_col = TRUE,
                                    annotation_names_row = TRUE,
                                    annotation_colors = annot_colors(),
                                    border_color = "grey60",
                                    display_numbers = input$print_nb,
                                    number_format = "%.2f"
      )
      return(heatmap)
    })
  })
  
  selection_features_results <- eventReactive(input$launch_heatmap, {
    req(input$select_features, input$selection_method)
    if(input$selection_method == "indicspecies"){
      result <- indic_species_results()$sign
      t_table <- as.data.frame(phyloseq::tax_table(agglom_data())) %>%
        rownames_to_column()
      res <- result %>%
        rownames_to_column() %>%
        left_join(t_table, by = "rowname")
      return(res)
    }
  })

  output$ui_selected_features <- renderUI({
    req(input$selection_method)
    if(input$selection_method == "indicspecies"){
      tagList(
        downloadButton(ns("table_download"), label = "Download table"),
        DT::dataTableOutput(ns("feat"))
      )
    }
  })

  output$feat <- DT::renderDataTable({
    req(selection_features_results())
    selection_features_results()
  }, filter = "top", options = list(scrollX = TRUE))

  observe({
    output$heatmap_t <- renderPlot({
      req(heatmap(), plot_height(), plot_width())
      withProgress(message = 'Computing heatmap...',{
        heatmap()
      })
    }, height = plot_height(), width = plot_width())
  })

  output$heatmap_download <- downloadHandler(
    filename = "heatmap.svg",
    content = function(file){
      req(heatmap())
      grDevices::svg(filename = file, width = plot_width()/96, height = plot_height()/96)
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
  })
}


## To be copied in the UI
# mod_heatmap_ui("heatmap_ui_1")

## To be copied in the server
# mod_heatmap_server("heatmap_ui_1", r = r)
