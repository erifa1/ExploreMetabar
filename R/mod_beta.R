#' mod_beta UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom DT dataTableOutput
#' @importFrom plotly plot_ly
#' @importFrom plotly add_trace

mod_beta_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        infoBox("",
                "Use phyloseq object without taxa merging step.",
                icon = icon("info-circle"), fill=TRUE, width = 10
        )
      ),
      fluidRow(
        box(
          fluidPage(
            fluidRow(
              radioButtons(
                ns("beta_norm_bool"),
                label = "Use normalized data (prefer TSS normalization)",
                inline = TRUE,
                choices = list(
                  "Raw" = 0 ,
                  "Normalized" = 1
                ), selected = 1
              )
            ),
            fluidRow(
              radioButtons(ns('ordi_type'), 'Choose your ordination type:',
                           inline = T,
                           choices = c('Unconstrained', 'Constrained', 'Distance-based'),
                           selected = 'Unconstrained')
            ),
            fluidRow(
              uiOutput(ns('ui_metrics'))
            ),
            fluidRow(
              radioButtons(ns("ordination"), "Choose one ordination:", inline = TRUE,
                           choices = '',
                           selected = ''
              )
            ),
            fluidRow(
              radioButtons(ns("plot_type"), "Choose plot type:", inline = TRUE,
                           choices =
                             list("samples", "taxa", "biplot"),
                           selected = c("samples")
              )
            ),
            fluidRow(
              uiOutput(ns('ui_beta_fact1')),
              uiOutput(ns('ui_beta_fact2'))
            ),
            fluidRow(
              uiOutput(ns('rank_select'))
            ),
            fluidRow(
              materialSwitch(
                ns('envfit_switch'),
                label = 'Try envfit',
                value = FALSE,
                status = 'primary'
              ),
              uiOutput(ns('ui_taxa'))
            ),
            fluidRow(
              actionButton(ns("launch_beta"), "Run Beta Plot", icon = icon("play-circle"),
                           style="color: #fff; background-color: #3b9ef5; border-color: #1a4469")
            )
          ),  title = "Settings:", width = 6, status = "warning", solidHeader = TRUE
        ),
        uiOutput(ns('envfit_box')) 
      ),
      fluidRow(
        box(
          shinycustomloader::withLoader(
            plotly::plotlyOutput(ns("plot1"), height = "730px"),
            type = "html", loader = "loader4"
          ),
          title = "Ordination plot:", width = 12, height = "800px", status = "primary", solidHeader = TRUE
        ), style = "height:800px;"
      ),
      fluidRow(
        box(
          title = "Permanova with adonis:", width = 12, status = "primary", solidHeader = TRUE,
          uiOutput(ns("factor2")),
          uiOutput(ns("interac_factor")),
          actionButton(ns("go1"), "Update Test", style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
          h3('ADONIS formula:'),
          verbatimTextOutput(ns("adonis_formula")),
          h2('Permanova Adonis Test Result: '),
          DT::dataTableOutput(ns('adonistest')),
          uiOutput(ns('pairwise_res')),
          uiOutput(ns('disper_res'))
        )
      )
    )
  )
}


phyloseq_to_ampvis2 <- function(physeq) {
  #check object for class
  if(!any(class(physeq) %in% "phyloseq"))
    stop("physeq object must be of class \"phyloseq\"", call. = FALSE)
  
  #ampvis2 requires taxonomy and abundance table, phyloseq checks for the latter
  if(is.null(physeq@tax_table))
    stop("No taxonomy found in the phyloseq object and is required for ampvis2", call. = FALSE)
  
  #OTUs must be in rows, not columns
  if(phyloseq::taxa_are_rows(physeq))
    abund <- as.data.frame(phyloseq::otu_table(physeq)@.Data)
  else
    abund <- as.data.frame(t(phyloseq::otu_table(physeq)@.Data))
  
  #tax_table is assumed to have OTUs in rows too
  tax <- phyloseq::tax_table(physeq)@.Data
  
  #merge by rownames (OTUs)
  otutable <- merge(
    abund,
    tax,
    by = 0,
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  colnames(otutable)[1] <- "OTU"
  
  #extract sample_data (metadata)
  if(!is.null(physeq@sam_data)) {
    metadata <- data.frame(
      phyloseq::sample_data(physeq),
      row.names = phyloseq::sample_names(physeq), 
      stringsAsFactors = FALSE, 
      check.names = FALSE
    )
    
    #check if any columns match exactly with rownames
    #if none matched assume row names are sample identifiers
    samplesCol <- unlist(lapply(metadata, function(x) {
      identical(x, rownames(metadata))}))
    
    if(any(samplesCol)) {
      #error if a column matched and it's not the first
      if(!samplesCol[[1]])
        stop("Sample ID's must be in the first column in the sample metadata, please reorder", call. = FALSE)
    } else {
      #assume rownames are sample identifiers, merge at the end with name "SampleID"
      if(any(colnames(metadata) %in% "SampleID"))
        stop("A column in the sample metadata is already named \"SampleID\" but does not seem to contain sample ID's", call. = FALSE)
      metadata$SampleID <- rownames(metadata)
      
      #reorder columns so SampleID is the first
      metadata <- metadata[, c(which(colnames(metadata) %in% "SampleID"), 1:(ncol(metadata)-1L)), drop = FALSE]
    }
  } else
    metadata <- NULL
  
  #extract phylogenetic tree, assumed to be of class "phylo"
  if(!is.null(physeq@phy_tree)) {
    tree <- phyloseq::phy_tree(physeq)
  } else
    tree <- NULL
  
  #extract OTU DNA sequences, assumed to be of class "XStringSet"
  if(!is.null(physeq@refseq)) {
    #convert XStringSet to DNAbin using a temporary file (easiest)
    fastaTempFile <- tempfile(pattern = "ampvis2_", fileext = ".fa")
    Biostrings::writeXStringSet(physeq@refseq, filepath = fastaTempFile)
  } else
    fastaTempFile <- NULL
  
  #load as normally with amp_load
  ampvis2::amp_load(
    otutable = otutable,
    metadata = metadata,
    tree = tree,
    fasta = fastaTempFile
  )
}

veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#' mod_beta Server Function
#'
#' @importFrom vegan vegdist
#' @importFrom vegan adonis2
#' @importFrom plotly ggplotly
#' @importFrom DT renderDataTable
#'
#' @noRd
mod_beta_server <- function(input, output, session, r = r){
  ns <- session$ns

  # observe({
  #   req(r$phyloseq_filtered())
  #   updateSelectInput(session, "beta_fact1",
  #                     choices = r$phyloseq_filtered()@sam_data@names)
  # })
  
  isNumFactor <- reactive({
    req(input$beta_fact1, r$sdat())
    metadata <- as(r$sdat(), "data.frame")
    if(is.numeric(metadata[, input$beta_fact1])){
      return(TRUE)
    } else{
      return(FALSE)
    }
  })
  
  
  output$ui_metrics <- renderUI({
    if(input$ordi_type == 'Distance-based'){
      radioButtons(ns("metrics"), "Choose one distance metric:", inline = TRUE,
                   choices =c('bray', 'jaccard', 'none'),
                   selected = c("bray")
      )
    }
  })
  
  
  output$rank_select <- renderUI({
    if(input$ordination == 'NMDS'){
      tags$div(
        hr(style = "border-top: 1px solid #000000;"),
        h4('Taxa ordination options: '),
        column(4,
           selectInput(
             ns("rank_color"),
             label = "Select rank to color taxa points: ",
             choices = rank_names(physeq()),
             selected = rank_names(physeq())[length(rank_names(physeq()))]
           )
        )
      )
    }
  })
  
  
  output$pairwise_res <- renderUI({
    if(! isNumFactor() && get_meta_col() != 'sample.id'){
      box(
        title = "Pairwise Adonis Test", width = 12, status = "primary", solidHeader = TRUE,
        DT::dataTableOutput(ns("adonispairwisetest"))
      )
    }
  })
  
  
  output$disper_res <- renderUI({
    if(! isNumFactor() && get_meta_col() != 'sample.id'){
      box(
        title = "Dispersion results:", width = 12, status = "primary", solidHeader = TRUE,
        h3('Boxplots distance to centroid for each group:'),
        checkboxInput(ns("order1"), label = "Automatic order factor", value = TRUE),
        plotlyOutput(ns("dispersionPlot")),
        h3('Anova on dispersion:'),
        DT::dataTableOutput(ns("dispersionTable")),
        h3('TukeyHSD test on dispersion'),
        DT::dataTableOutput(ns("dispersionTukey"))
      )
    }
  })
  
  
  output$ui_beta_fact1 <- renderUI({
    req(r$sdat())
    
    metadata <- as(r$sdat(), "data.frame")
    tags$div(
      hr(style = "border-top: 1px solid #000000;"),
      h4('Sample ordination options: '),
      column(4,
             selectInput(
               ns("beta_fact1"),
               label = "Select main factor to test + color plot: ",
               choices = colnames(metadata)
             )
      )
    )
    
  })
  
  
  output$ui_beta_fact2 <- renderUI({
    req(r$sdat())
    metadata <- as(r$sdat(), "data.frame")
    # if(! isNumFactor() && input$plot_type == 'samples'){
        num_col_names <- metadata %>% dplyr::select_if(is.numeric) %>% colnames
        tmp <- dplyr::setdiff(colnames(metadata), num_col_names)
        tags$div(
          column(4,
             selectInput(
               ns("beta_fact2"),
               label = "Select second factor to combine: ",
               choices = c('none',dplyr::setdiff(tmp, input$beta_fact1))
             )
          ),
        )
    # }
  })
  
  
  output$envfit_box <- renderUI({
    if(input$ordination %in% c('NMDS', 'PCOA') && input$envfit_switch){
      box(
        multiInput(
          ns('envfit_param'),
          label = "Select variables: ",
          choices = NULL,
          choiceValues = colnames(r$sdat()),
          choiceNames = colnames(r$sdat())
        ),
        numericInput(ns('envfit_pval'), 'p-value', value = 1, min = 0, max = 1, step = 0.01),
        title = 'VEGAN envfit', status = 'primary'
      )
    }
  })
  
  output$ui_taxa <- renderUI({
    if(input$ordination == 'NMDS'){
      fluidRow(
        materialSwitch(
          ns('taxa_switch'),
          label = 'Add taxa to plot',
          value = FALSE,
          status = 'primary'
        ),
        materialSwitch(
          ns('plotly_switch'),
          label = 'Plotly on sample or taxa',
          value = FALSE,
          status = 'primary'
        )
      )
      
      
    }
  })
  
  
  local_metadata <- reactive({
    req(input$beta_fact1)
    metadata <- as(r$sdat(), "data.frame")
    if(!isNumFactor()){
      if(input$beta_fact2 != 'none' && ! is.numeric(metadata[, input$beta_fact1])){
        metadata <- tidyr::unite(metadata, !!get_meta_col(), input$beta_fact1, input$beta_fact2, na.rm=TRUE)
        metadata[, get_meta_col()] <- as.factor(metadata[, get_meta_col()])
        metadata <- select(metadata, "sample.id", get_meta_col())
      }
    }
    else{
      metadata <- select(metadata, "sample.id", input$beta_fact1)
    }
    return(metadata)
  })
  
  
  get_meta_col <- reactive({
    if(!isNumFactor()){
      if(input$beta_fact2 != 'none'){
        meta.col <- paste0(input$beta_fact1, '_', input$beta_fact2)
      }
      else{
        meta.col <- input$beta_fact1
      }
    } else{
      meta.col <- input$beta_fact1
    }
    return(meta.col)
  })
  
  
  observeEvent(input$ordi_type, {
    if(input$ordi_type == 'Unconstrained'){
      ch <- c('PCA', 'CA', 'DCA')
    } else if (input$ordi_type == 'Constrained'){
      ch <- c('RDA', 'CCA')
    } else{
      ch <- c('PCOA', 'NMDS')
    }
    updateRadioButtons(session,
                       "ordination",
                       choices = ch ,
                       inline = T)
  })
  

  observe({
    req(r$phyloseq_filtered())
    if(is.null(phy_tree(r$phyloseq_filtered(), errorIfNULL=FALSE))){
      flog.info("no phytree beta metrics update")
      ch1 = list("bray", "jaccard")
    }else{
      ch1 = list("bray", "jaccard", "unifrac", "wunifrac")
    }
    updateRadioButtons(session, "metrics",
                      choices = ch1, inline = TRUE)
  })


  output$factor2 = renderUI({
    req(get_meta_col(), r$sdat())
    facts = phyloseq::sample_variables(r$sdat())
    Fchoices = facts[facts != get_meta_col()]

    checkboxGroupInput(
      ns("covariate_fact"),
      label = "Select covariable(s) to test: ",
      choices = Fchoices,
      inline = TRUE
    )
  })

  
  output$interac_factor <- renderUI({
    req(get_meta_col(), r$sdat())
    facts = phyloseq::sample_variables(r$sdat())
    Fchoices = facts[facts != get_meta_col()]

    checkboxGroupInput(
      ns("interFactor"),
      label = "Select interaction factor(s) to test: ",
      choices = Fchoices,
      inline = TRUE
    )
  })

  
  physeq <- reactive({
    req(r$phyloseq_filtered, r$phyloseq_filtered_norm, input$beta_norm_bool)
    if(input$beta_norm_bool==0){
      data <- phyloseq::rarefy_even_depth(r$phyloseq_filtered(), rngseed = 20210225, verbose = FALSE)
    }
    if(input$beta_norm_bool==1){
      data <- r$phyloseq_filtered_norm()
    }
    sample_data(data) <- sample_data(local_metadata())
    return(data)
  })

  
  physeq_dist <- reactive({
    req(input$metrics)
    res <- phyloseq::distance(physeq(), method = input$metrics)
    return(res)
  })

  
  ord <- reactive({
    req(input$ordination)
    if(input$ordination == 'NMDS'){
      validate(
        need(input$metrics %in% c('bray', 'jaccard', 'none'), 'For NMDS ordination only bray and jaccard distances allowed.')
        )
      res <- vegan::metaMDS(veganifyOTU(physeq()), distance = input$metrics, wascores=TRUE, trace=FALSE, autotransform = FALSE)
    } else if(input$ordination == 'PCOA'){
      validate(
        need(input$metrics %in% c('bray', 'jaccard', 'none'), 'For PCoA ordination only bray and jaccard distances allowed.')
      )
      res <- vegan::wcmdscale(physeq_dist(), k=2)
    } else{
      res <- phyloseq::ordinate(physeq= physeq(), distance = physeq_dist(), method= input$ordination)
    }
    return(res)
  })

  
  get_sites_nmds_coord <- reactive({
    nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='sites') %>% 
                    as_tibble(rownames="sample.id")
    nmds_coord <- nmds_coord %>% 
                    inner_join(., local_metadata() %>% 
                                 select(sample.id, !!get_meta_col()), by="sample.id")
    return(nmds_coord)
  })
  
  
  get_species_nmds_coord <- reactive({
    # browser()
    nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='specie') %>% as_tibble(rownames=r$rank_glom())
    
    if(r$rank_glom() == 'ASV'){
      nmds_coord <- nmds_coord %>% 
                      inner_join(., tax_table(physeq()) %>% 
                                   as.data.frame() %>% 
                                   as_tibble(rownames="ASV") %>% 
                                   select(ASV, input$rank_color), by='ASV')
    } else{
      nmds_coord <- nmds_coord %>% 
                      inner_join(., tax_table(physeq()) %>% 
                                   as.data.frame() %>% 
                                   as_tibble(rownames=r$rank_glom()) %>% 
                                   select(r$rank_glom(), input$rank_color), by=r$rank_glom())
    }
    return(nmds_coord)
  })
  
  
  get_dist_plot <- reactive({
    plot_str <- "p <- ampvis2::amp_ordinate(phyloseq_to_ampvis2(physeq()), type = input$ordination, distmeasure = input$metrics, transform = 'none', filter_species = 0, sample_color_by = get_meta_col(), sample_colorframe = T"
    # browser()
    if(input$envfit_switch){
      validate(
        need(length(input$envfit_param) > 0, "Select one metadata value for envfit." )
      )
      env_num <- names(r$sdat()[,input$envfit_param])[which(sapply(r$sdat()[,input$envfit_param], is.numeric))]
      env_fact <- names(r$sdat()[,input$envfit_param])[which(sapply(r$sdat()[,input$envfit_param], is.character))]
      
      if(length(env_num) > 0 || length(env_fact) > 0){
        env_str <- ", envfit_show = T"
        
      }
      if(length(env_num) > 0){
        env_num <- gsub('\\.','_',env_num)
        env_str <- paste0(env_str, ", envfit_numeric = env_num")
      }
      if(length(env_fact) > 0){
        validate(
          need(length(env_fact) == 1, 'Choose only one factor.')
        )
        env_fact <- gsub('\\.','_',env_fact)
        env_str <- paste0(env_str, ", envfit_factor = env_fact")
      }
      
      plot_str <- paste0(plot_str, env_str, glue::glue(", envfit_signif_level = {input$envfit_pval}"))
    }
    if('taxa_switch' %in% names(input)){
      if(input$taxa_switch && input$ordination == 'NMDS'){
        taxa_str <- glue::glue(", species_plot = T, species_label_taxonomy = '{input$rank_color}'")
        plot_str <- paste0(plot_str, taxa_str)
      }
      if('plotly_switch' %in% names(input)){
        if(!input$plotly_switch){
          plotly_str <- ", sample_plotly = 'all', species_plotly = F"
        } else{
          plotly_str <- ", species_plotly = T"
        }
        plot_str <- paste0(plot_str, plotly_str)
      }
    } else{
      plotly_str <- ", sample_plotly = 'all'"
      plot_str <- paste0(plot_str, plotly_str)
    }
    
    plot_str <- paste0(plot_str, ')')
    # browser()
    eval(parse(text=plot_str))
    return(p)
  })

  
  base_plot <- reactive({
    # p <- phyloseq::plot_ordination(physeq = physeq(), type = input$plot_type, ordination = ord(), axes = c(1, 2))
    # if(input$ordination %in% c('PCOA', 'NMDS')){
    #   p <- get_dist_plot()
    # }
    fig <- plotly::plot_ly()
    if(input$plot_type == 'samples'){
      nmds_coord <- get_sites_nmds_coord()

      p <- ggplot2::ggplot() +
            geom_point(data = nmds_coord, mapping = aes(x=NMDS1, y=NMDS2, color=.data[[get_meta_col()]], text = paste('sample.id:',phyloseq::sample_names(physeq())))) +
            stat_ellipse(data = nmds_coord, mapping = aes(x=NMDS1, y=NMDS2, group = !!sym(get_meta_col()), color = !!sym(get_meta_col())))
    } else if (input$plot_type == 'taxa'){
      nmds_coord <- get_species_nmds_coord()
      if(r$rank_glom() == 'ASV'){
        taxa <- tax_table(physeq()) %>% as.data.frame() %>% as_tibble(rownames="ASV") %>% select(ASV) %>% pull
      } else{
        taxa <- tax_table(physeq()) %>% as.data.frame() %>% as_tibble(rownames=r$rank_glom()) %>% select(r$rank_glom()) %>% pull
      }
      p <- ggplot2::ggplot() +
            geom_point(data = nmds_coord, aes(x=NMDS1, y=NMDS2, color=.data[[input$rank_color]], taxa = taxa))
    } else if(input$plot_type == 'biplot'){
      nmds_coord_species <- get_species_nmds_coord()
      nmds_coord_sites <- get_sites_nmds_coord()
      if(r$rank_glom() == 'ASV'){
        taxa <- tax_table(physeq()) %>% as.data.frame() %>% as_tibble(rownames="ASV") %>% select(ASV) %>% pull
      } else{
        taxa <- tax_table(physeq()) %>% as.data.frame() %>% as_tibble(rownames=r$rank_glom()) %>% select(r$rank_glom()) %>% pull
      }
      p <- ggplot() +
        geom_point(data = nmds_coord_species, aes(x=NMDS1, y=NMDS2, color=.data[[input$rank_color]], taxa=taxa), size=4) +
        geom_point(data = nmds_coord_sites, aes(x=NMDS1, y=NMDS2, fill=.data[[get_meta_col()]], sample.id = phyloseq::sample_names(physeq())), shape=23, size=6)
    }

    if(input$envfit_switch){
      # browser()
      env <- as(r$sdat(), 'data.frame')
      env <- env[, input$envfit_param]
      en <- vegan::envfit(ord(), env)
      if(length(env %>% select_if(is.numeric) %>% colnames())>0){
        en_coord_cont <- as.data.frame(vegan::scores(en, "vectors")) * vegan::ordiArrowMul(en)
        p <- p + geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                              data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
                geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30",
                    fontface = "bold", label = row.names(en_coord_cont))
      }
      if(length(env %>% select_if(is.character) %>% colnames())>0){
        en_coord_cat <- as.data.frame(vegan::scores(en, "factors")) * vegan::ordiArrowMul(en)
        p <- p + geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2),
                            shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
          geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2),
                    label = row.names(en_coord_cat), colour = "navy", fontface = "bold")
      }

    }
    # p$layers[[1]] <- NULL
    # 
    # xrange <- c()
    # xrange[1] <- layer_scales(p)$x$range$range[1] - abs(layer_scales(p)$x$range$range[1])*3
    # xrange[2] <- layer_scales(p)$x$range$range[2] + abs(layer_scales(p)$x$range$range[2])*3
    # 
    # yrange <- c()
    # yrange[1] <- layer_scales(p)$y$range$range[1] - abs(layer_scales(p)$y$range$range[1])*3
    # yrange[2] <- layer_scales(p)$y$range$range[2] + abs(layer_scales(p)$y$range$range[2])*3
    # return(list('plot'=p, 'xrange'=xrange, 'yrange'=yrange))
    return(list('plot'=p))
  })


  output$plot1 <- plotly::renderPlotly({
    beta_plot()
  })


  beta_plot <- eventReactive(input$launch_beta, {
    withProgress({
      p <- base_plot()$plot
      # p <- p + xlim(base_plot()$xrange) + ylim(base_plot()$yrange)
      # p <- p + geom_point() + theme_bw()
      # browser()
    #   if(input$plot_type == 'samples'){
    #     p <- ggplotly(p, tooltip=c("x", "y", "sample.id"))
    #   } else if(input$plot_type == 'taxa'){
    #     p <- ggplotly(p, tooltip=c("x", "y", "taxa"))
    #   }
    #   else if(input$plot_type == 'biplot'){
    #     p <- ggplotly(p, tooltip=c("x", "y", "taxa", "sample.id"))
    #   }else{
    #     
    #   }
    #   p <- p %>% config(toImageButtonOptions = list(format = "svg"))
    #   
    }, message = "Plot Beta...")
    return(p)
  })

  
  get_formula <- reactive({
    req(input$metrics, get_meta_col())
    form <- glue::glue('dist ~ Depth + ')
    if(!is.null(input$covariate_fact)){
      cov1 = paste(input$covariate_fact, collapse = " + ")
      form <- paste(form, glue::glue('{cov1} + {get_meta_col()}'), sep='')
    }
    else if(!is.null(input$interFactor)){
      cov1 = paste(input$interFactor, collapse = "*")
      form <- paste(form, glue::glue('{get_meta_col()}*{cov1}'), sep='')
    }
    else{
      form <- paste(form, glue::glue('{get_meta_col()}'), sep='')
    }
    return(form)
  })
  
  
  get_dispersion_res <- reactive({
    req(physeq_dist(), get_meta_col())
    res <- vegan::betadisper(physeq_dist(), local_metadata()[,get_meta_col()])
    return(res)
  })
  
  
  get_dispersion_anova <- reactive({
    req(get_dispersion_res())
    res <- anova(get_dispersion_res())
    return(res)
  })
  
  
  get_dispersion_tukey <- reactive({
    req(get_dispersion_res())
    res <- TukeyHSD(get_dispersion_res())
    return(res)
  })
  
  
  get_adonis_res <- reactive({
    req(physeq_dist(), get_formula())
    dist <- physeq_dist()
    mdata <- local_metadata()
    mdata$Depth <- sample_sums(physeq())
    # Filter NA value in metadata
    mdata <- mdata %>% filter(!is.na(get_meta_col()))
    res <- vegan::adonis2(as.formula(get_formula()), data = mdata, permutations = 1000)
    return(data.frame(res))
  })
  
  
  get_pairwise_res <- reactive({
    req(physeq_dist(), get_meta_col(), local_metadata())
    res <- pairwise.adonis(physeq_dist(), local_metadata()[,get_meta_col()], p.adjust.m = "fdr")
    return(res)
  })
  

  output$adonis_formula <- renderText({
    print(get_formula())
  })


  output$adonistest <- DT::renderDataTable({
    get_adonis_res()
  })

  output$adonispairwisetest <- DT::renderDataTable({
    get_pairwise_res()
  })


  dfdisper <- reactive({
    cat(file=stderr(),'dfdisper ...',"\n")
    
    df1 = cbind.data.frame(distances = get_dispersion_res()$distances, group = get_dispersion_res()$group)

    if(input$order1){
      print("ORDER factor")
      df1$group = factor( df1$group, levels = gtools::mixedsort(levels(df1$group)) ) 
    }
    cat(file=stderr(),'Done ...',"\n")
    
    return(df1)
  })


  output$dispersionPlot <- renderPlotly({
   df1 <- dfdisper()
   plot_ly(df1, x = ~group, y = ~distances,
           color = ~group, type = 'box') %>%
     layout(title="", yaxis = list(title = "Distance to centroid"), xaxis = list(title = 'Group'), barmode = 'stack') %>%
    config(toImageButtonOptions = list(format = "svg"))
 })

  output$dispersionTable <- DT::renderDataTable({
    get_dispersion_anova()
  })

  output$dispersionTukey <- DT::renderDataTable({
    get_dispersion_tukey()$group
  })
}

## To be copied in the UI
# mod_mod_beta_ui("mod_beta_ui_1")

## To be copied in the server
# callModule(mod_mod_beta_server, "mod_beta_ui_1")
