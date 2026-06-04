#' asvenn UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session,r Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom bslib layout_sidebar sidebar accordion accordion_panel navset_card_underline nav_panel input_task_button
#' @importFrom bsicons bs_icon
mod_asvenn_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "Settings",
      open = "desktop",
      width = "350px",
      htmltools::p(
        "Select conditions to highlight shared taxa. Venn diagrams show ASV overlap across groups.",
        style = "font-size: 0.9em; color: grey;"
      ),
      tags$hr(),
      accordion(
        id = ns("config_accordion"),
        open = "Venn settings",
        multiple = TRUE,
        accordion_panel(
          "Venn settings",
          icon = bs_icon("diagram-2"),
          selectInput(ns("Fact1"), label = "Select factor:", choices = ""),
          uiOutput(ns("lvls1")),
          numericInput(
            ns("minAb"),
            "Minimum raw abundance to detect a taxon in a group:",
            value = 1, min = 1, max = NA
          )
        ),
        accordion_panel(
          "Shared taxa chart",
          icon = bs_icon("bezier2"),
          htmltools::p(
            "Available after running the Venn diagram.",
            style = "font-size: 0.9em; color: grey;"
          ),
          uiOutput(ns("ui_alluvial_ranks"))
        )
      ),
      tags$hr(),
      input_task_button(
        ns("go1"),
        "Run / Update Venn",
        icon = bs_icon("play-fill"),
        class = "btn-primary w-100 btn-lg",
        label_busy = "Computing…"
      )
    ),
    navset_card_underline(
      title = "Results",
      full_screen = TRUE,
      nav_panel(
        "Venn Classic",
        icon = bs_icon("diagram-2"),
        plotOutput(ns("venn2"), height = "600px")
      ),
      nav_panel(
        "Venn VennR",
        icon = bs_icon("bezier2"),
        imageOutput(ns("venn1"), width = "100%", height = "600px")
      ),
      nav_panel(
        "Table",
        icon = bs_icon("table"),
        DT::dataTableOutput(ns("tabvenn1")),
        downloadButton(outputId = ns("otable_download"), label = "Download Table")
      ),
      nav_panel(
        "Boxplot",
        icon = bs_icon("bar-chart-line"),
        htmltools::p("Click on a row in the Table tab to generate the boxplot for that taxon."),
        plotly::plotlyOutput(ns("boxplot_chart"), width = "100%", height = "500px")
      ),
      nav_panel(
        "Shared Taxa",
        icon = bs_icon("bezier2"),
        plotOutput(ns("alluvial_plot"), height = "600px")
      )
    )
  )
}




#' asvenn Server Function
#'
#' @noRd
#'
#' @importFrom futile.logger flog.threshold
#' @importFrom grDevices rainbow
#' @importFrom venn venn
#' @importFrom nVennR plotVenn
#' @import ggpolypath
#' @import ggalluvial


mod_asvenn_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns


  observe({
    req(r$phyloseq_filtered(), r$sdat())
    metadata <- r$sdat()
    num_col_names <- metadata %>% dplyr::select_if(is.numeric) %>% colnames
    tmp <- dplyr::setdiff(colnames(metadata), num_col_names)
    updateSelectInput(session, "Fact1",
                      choices = tmp)
  })


  output$lvls1 = renderUI({
    req(input$Fact1, r$sdat())
    metadata <- r$sdat()
    level1 <- as.character(na.omit(unique(metadata[,input$Fact1])))
    shinyWidgets::pickerInput(ns("lvls1"), label = "Select up to 5 levels :",
          choices = level1, selected = level1[1:3], multiple = TRUE
      )

  })

  output$ui_alluvial_ranks <- renderUI({
    req(r$phyloseq_filtered())
    ranks <- phyloseq::rank_names(r$phyloseq_filtered())
    selectInput(
      ns("alluvial_rank"),
      label = "Represent taxa at rank:",
      choices = ranks,
      selected = ranks[ceiling(length(ranks) / 2)]
    )
  })

  # The rank selectInput lives in a collapsed accordion panel (display:none),
  # which Shiny would otherwise suspend — leaving input$alluvial_rank NULL and
  # the Shared Taxa plot blank until the panel is manually expanded.
  outputOptions(output, "ui_alluvial_ranks", suspendWhenHidden = FALSE)

  resVenn <- eventReactive(input$go1, {
    req(r$phyloseq_filtered(), input$lvls1)
    flog.info('compute Venn diagram...')
    if(length(input$lvls1) < 2 || length(input$lvls1) > 5){
      shinyalert("Oops!", "You need to choose between 2 and 5 factors...", type = "error")
      req(FALSE)
    }
    else{
      withProgress(message = "Computing Venn overlaps…", {
        res <- list()
        TFdata <- list()
        TFtax <- tibble(taxa = character(), taxo = character())

        for(lvl in input$lvls1){
          flog.info(lvl)
          keep <- sample_data(r$phyloseq_filtered())[[input$Fact1]] %in% lvl
          data.tmp <- prune_samples(keep, r$phyloseq_filtered())
          sp_data <- prune_taxa(taxa_sums(data.tmp) > 0, data.tmp)

          abund_to_zero = function(x){
            x[x < input$minAb] <- 0
            return(x)
          }
          sp_data <- transform_sample_counts(sp_data, fun = abund_to_zero)
          sp_data <- prune_taxa(taxa_sums(sp_data) > 0, sp_data)

          TT = cbind(otu_table(sp_data),tax_table(sp_data))

          TFdata[[lvl]] <- TT
          TFtax <- dplyr::full_join(TFtax, as_tibble(cbind(taxa = row.names(TT), taxo =  as.character(apply(TT[,colnames(tax_table(sp_data))], 1, paste, collapse=";") ) )), by = c("taxa", "taxo"))
          row.names(TFtax[[lvl]]) = TFtax[[lvl]][,1]
          incProgress(1/length(input$lvls1), detail = lvl)
        }

        TF <- sapply(TFdata, row.names, simplify = FALSE)
        names(TF) = input$lvls1
        res$TF <- TF
        all_taxa <- unique(unlist(TF))
        mtab <- sapply(TF, function(x) as.integer(all_taxa %in% x))
        rownames(mtab) <- all_taxa
        v.table <- as_tibble(mtab, rownames = "taxa")
        v.table <- full_join(v.table, TFtax, by = 'taxa')
        res$v.table <- v.table
      })
      return(res)
    }
  })

  getVenn1 <- reactive({
    req(resVenn)
    outfile <- tempfile(fileext='.svg')
    nVennR::plotVenn(resVenn()$TF, showPlot = T, labelRegions = T, systemShow=F, outFile = outfile, setColors = getPal())
    list(
      src = normalizePath(outfile),
      contentType = "image/svg+xml",
      # nVennR writes a square SVG; object-fit keeps its 1:1 aspect ratio
      # inside the wide/short container instead of stretching it horizontally.
      style = "display:block; width:100%; height:100%; object-fit:contain; margin:0 auto;"
    )
  })
  
  
  output$venn1 <- renderImage(
    getVenn1()
  , deleteFile=TRUE)
  
  output$venn2 <- renderPlot({
    invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
    venn::venn(resVenn()$TF, zcolor = getPal(), ilcs = 1.5, sncs = 2,
                      ggplot = TRUE, ilabels = 'counts')
  })


  getPal <- reactive({
    req(input$Fact1, resVenn())
    validate(need(!input$Fact1 %in% names(r$numeric_palettes()),
                  "Venn diagram requires a categorical variable."))
    pal <- r$factor_colors()[[input$Fact1]]
    pal <- pal[names(pal) %in% names(resVenn()$TF)]
    if (is.null(pal)){
      pal <- RColorBrewer::brewer.pal(n= length(names(resVenn()$TF)), name = 'Set1')
    }
    return(pal)
  })
  
  output$tabvenn1 <-DT::renderDataTable({
    resVenn()$v.table
  }, filter="top", selection = "single", options = list(scrollX = TRUE))

  output$otable_download <- downloadHandler(
    filename = "venn_table.csv",
    content = function(file) {
      req(resVenn()$v.table)
      write.table(resVenn()$v.table, file, sep="\t", row.names=FALSE)
    }
  )

  get_boxplot_data <- reactive({
    req(input$tabvenn1_row_last_clicked, input$lvls1)
    keep <- sample_data(r$phyloseq_filtered())[[input$Fact1]] %in% input$lvls1
    data.tmp <- prune_samples(keep, r$phyloseq_filtered())

    obj <- prune_taxa(pull(resVenn()$v.table[input$tabvenn1_row_last_clicked,1]), data.tmp)
    ot <- as.data.frame(otu_table(obj))
    ot <- as.data.frame(t(ot))
    mt <- as.data.frame(as.matrix(sample_data(obj)))
    ot[input$Fact1] <- as.vector(mt[rownames(ot), input$Fact1])

    return(ot)
  })

  get_boxplot <- reactive({
    dt <- get_boxplot_data()
    fig <- plotly::plot_ly(x =~dt[,2], y=~dt[,1], type = "box",color = ~dt[,2], colors = getPal()) %>%
      plotly::layout(xaxis = list(title = colnames(dt)[2]),
      yaxis = list(title = colnames(dt)[1]),
      title = colnames(dt)[1] )

    return(fig)
  })


  output$boxplot_chart <- plotly::renderPlotly({
    get_boxplot()
  })

  get_alluvial_data <- reactive({
    req(resVenn(), input$lvls1, input$alluvial_rank)

    vtab       <- resVenn()$v.table
    level_cols <- intersect(input$lvls1, colnames(vtab))

    # ASVs present in ALL selected levels → shared; rest → unique to their group(s)
    shared_asvs <- vtab$taxa[
      rowSums(vtab[, level_cols, drop = FALSE]) == length(level_cols)
    ]

    tax_at_rank <- as.data.frame(tax_table(r$phyloseq_filtered())) %>%
      dplyr::select(all_of(input$alluvial_rank)) %>%
      tibble::rownames_to_column("taxa") %>%
      dplyr::rename(taxon = !!input$alluvial_rank)

    # Work with the full union of taxa across all selected levels
    phy <- prune_taxa(vtab$taxa, r$phyloseq_filtered())

    rows <- dplyr::bind_rows(lapply(input$lvls1, function(lvl) {
      keep    <- phyloseq::sample_data(phy)[[input$Fact1]] %in% lvl
      phy_lvl <- prune_samples(keep, phy)
      phy_lvl <- prune_samples(sample_sums(phy_lvl) > 0, phy_lvl)
      if (phyloseq::nsamples(phy_lvl) == 0) return(NULL)

      data.frame(
        taxa      = taxa_names(phy_lvl),
        abundance = taxa_sums(phy_lvl),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(
          share_type = ifelse(taxa %in% shared_asvs, "shared", "unique")
        ) %>%
        dplyr::left_join(tax_at_rank, by = "taxa") %>%
        dplyr::group_by(taxon, share_type) %>%
        dplyr::summarise(abundance = sum(abundance), .groups = "drop") %>%
        dplyr::mutate(
          group       = lvl,
          # shared alluvia use taxon as ID → ribbon flows across groups
          # unique alluvia use a per-group ID → bar segment only, no ribbon
          alluvium_id = ifelse(
            share_type == "shared", taxon,
            paste(taxon, lvl, "uniq", sep = "_")
          )
        ) %>%
        dplyr::filter(abundance > 0)
    }))

    # Normalise to relative abundance (%) within each group
    rows %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(pct = abundance / sum(abundance) * 100) %>%
      dplyr::ungroup()
  })

  output$alluvial_plot <- renderPlot({
    req(get_alluvial_data())
    df       <- get_alluvial_data()
    df$group <- factor(df$group, levels = input$lvls1)
    pal      <- r$taxa_colors()[[input$alluvial_rank]]

    p <- ggplot2::ggplot(df,
      ggplot2::aes(x = group, y = pct,
                   stratum = taxon, alluvium = alluvium_id,
                   fill = taxon, alpha = share_type, label = taxon)) +
      ggalluvial::geom_alluvium(width = 1/3) +
      ggalluvial::geom_stratum(width = 1/3, color = "white") +
      ggplot2::geom_text(stat = "stratum", size = 3, check_overlap = TRUE) +
      ggplot2::scale_alpha_manual(
        values = c(shared = 0.9, unique = 0.3),
        labels = c(shared = "Shared (all groups)", unique = "Group-specific"),
        name   = NULL
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        x     = input$Fact1,
        y     = "Relative abundance (%)",
        title = paste("Taxa at rank:", input$alluvial_rank)
      ) +
      ggplot2::theme(legend.position = "right")

    if (!is.null(pal)) {
      p <- p + ggplot2::scale_fill_manual(
        values = pal, na.value = "grey70", guide = "none"
      )
    }
    p
  })
  })
}

## To be copied in the UI
# mod_asvenn_ui("asvenn_ui_1")

## To be copied in the server
# mod_asvenn_server("asvenn_ui_1", r = r)
