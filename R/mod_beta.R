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
#' @importFrom DT dataTableOutput renderDataTable JS
#' @importFrom plotly plot_ly
#' @importFrom plotly add_trace
#' @importFrom bslib navset_card_underline nav_panel accordion accordion_panel input_task_button
#' @import PCAmixdata
#' @import shinycustomloader
#' @import shinyWidgets
#' 

mod_beta_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "Settings",
      open = "desktop",
      width = "350px",
      
      # Section 1: Workflow Status (smaller font)
      htmltools::p(
        'ASVs normalized by Hellinger transformation. Environmental variables centered/scaled.',
        style = "font-size: 0.9em; color: grey;"
      ),
      
      tags$hr(),
      
      # Single clean accordion like data_loading.R
      accordion(
        id = ns("config_accordion"),
        open = "Ordination Setup",
        multiple = TRUE,
        
        # Core controls → first panel
        accordion_panel(
          "Ordination Setup",
          icon = bs_icon("tablet"),
          uiOutput(ns('ui_metrics')),
          tooltip(
            radioButtons(ns("ordination"), "Method:", inline = TRUE,
                         choices = c('PCOA', 'NMDS', 'dbRDA'),
                         selected = 'PCOA'),
            "Choose ordination technique based on your data type and research question.",
            placement = "right"
          ),
          uiOutput(ns('ui_beta_factor')),
          
          # Primary action button
          div(
            style = "margin: 1.5rem 0;",
            tooltip(
              input_task_button(ns("launch_beta"), "Run Ordination Analysis",
                          icon = bs_icon("play-fill"),
                          class = "btn-primary w-100 btn-lg",
                          label_busy = "Computing…"),
              "Computes ordination + tests. Progress shown below.",
              placement = "top"
            )
          )
        ),
        
        # Consolidated model parameters
        accordion_panel(
          "Model Parameters", 
          icon = bs_icon("sliders"),
          tooltip(bs_icon("question-circle"), "Advanced model configuration", placement = "right"),
          uiOutput(ns('ui_constrain')),
          uiOutput(ns('ui_envfit_box')),
          uiOutput(ns('ui_envfit_box_res'))
        ),
        
        # Plot options
        accordion_panel(
          "Plot Options", 
          icon = bs_icon("brush"), 
          shinyWidgets::prettyCheckboxGroup(ns("plot_type"), "Elements:", inline = TRUE,
                                            choices = list("Samples" = "samples", "Taxa" = "taxa", "Env" = "env"),
                                            selected = c("samples")),
          uiOutput(ns('ui_taxa_rank')),
          uiOutput(ns('ui_axe_x')),
          uiOutput(ns('ui_axe_y'))
        )
      )
    ),
    
    # Main content area with PERMANOVA/Dispersion integrated
    navset_card_underline(
      title = "Ordination Results",
      full_screen = TRUE,
      nav_panel("Ordination Plot", 
        full_screen = TRUE,
        plotOutput(ns('beta_ggplot'), height = "80vh")
      ),
      nav_panel("Screeplot", plotOutput(ns('screeplot'), height = "400px")),
      
      # PERMANOVA (moved from sidebar ui_permanova_box)
      nav_panel(
        title = "PERMANOVA", 
        icon = bs_icon("calculator"),
        h4("About PERMANOVA"),
        htmltools::p(
          "PERMANOVA (Permutational Multivariate Analysis of Variance) tests whether ",
          "groups of samples differ in their multivariate composition. It partitions the ",
          "variation of a dissimilarity matrix among explanatory factors and assesses ",
          "significance via permutations. The formula below uses the factors selected ",
          "for colors and shapes in the sidebar as explanatory variables, with sequencing ",
          "depth included as a covariate."
        ),
        hr(),
        h3("Adonis formula:"),
        verbatimTextOutput(ns("adonis_formula")),
        hr(),
        h3("Adonis results:"),
        DT::dataTableOutput(ns('adonistest')),
        hr(),
        h3("Pairwise Adonis results:"),
        DT::dataTableOutput(ns('adonispairwisetest'))
      ),
      
      # Dispersion (flat layout, full_screen only on boxplot)
      nav_panel(
        title = "Dispersion", 
        icon = bs_icon("funnel"),
        h4("Dispersion Analysis"),
        htmltools::p(
          "Multivariate homogeneity of group dispersions (variances) computed with ",
          "vegan::betadisper(). Boxplots show the distance of each sample to its group ",
          "centroid; ANOVA and TukeyHSD test for differences in dispersion among groups."
        ),
        hr(),
        card(
          full_screen = TRUE,
          card_header("Boxplots — Distance to centroid"),
          plotlyOutput(ns("dispersionPlot"))
        ),
        card(
          card_header("ANOVA on dispersion"),
          DT::dataTableOutput(ns("dispersionTable"))
        ),
        card(
          card_header("TukeyHSD on dispersion"),
          DT::dataTableOutput(ns("dispersionTukey"))
        )
      )
    )
  )
}


#' Convert a phyloseq OTU table to a vegan-compatible matrix
#'
#' Transposes the OTU table if taxa are rows, then coerces to matrix.
#' Used by mod_beta and mod_heatmap.
#'
#' @param physeq A phyloseq object
#' @return A numeric matrix with samples as rows and taxa as columns
#' @noRd
veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}


#' Build pickerInput choicesOpt content with column class info
#'
#' Generates the HTML content for pickerInput choices showing column name,
#' class, and optionally NA counts. Avoids repeating the same splitLayout
#' pattern across multiple pickerInputs.
#'
#' @param df A data.frame to inspect column classes from
#' @param cols Character vector of column names to display
#' @param show_na Logical; if TRUE, show NA count per column
#' @return A character vector of rendered HTML tags
#' @noRd
make_picker_choices <- function(df, cols, show_na = FALSE) {
  unlist(lapply(cols, function(x) {
    na_tag <- if (show_na) {
      tags$div(
        style = htmltools::css(color = "grey"),
        paste0(sum(is.na(df[, x])), "/", nrow(df), " NAs")
      )
    }
    htmltools::doRenderTags(
      tags$div(
        splitLayout(
          cellWidths = 200,
          tags$div(style = htmltools::css(fontWeight = "bold"), x),
          tags$div(style = htmltools::css(color = "grey"), class(df[, x])),
          na_tag
        )
      )
    )
  }))
}


#' mod_beta Server Function
#'
#' @import vegan
#' @importFrom plotly ggplotly
#' @importFrom DT renderDataTable
#' @importFrom permute how
#' @importFrom futile.logger flog.info flog.debug
#' @import htmltools
#' @import formula.tools
#' @import tidyr
#' @import glue
#' @import PCAmixdata
#' @import grid
#' @import gtools
#' @import goeveg
#' @import ggnewscale
#' @importFrom ggrepel geom_text_repel
#'
#' @noRd
mod_beta_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns

  # ---- Shared reactives setup ----
  factor_input   <- reactive({ input$beta_factor })
  get_meta_col   <- make_get_meta_col(factor_input, r)
  local_metadata <- make_local_metadata(factor_input, get_meta_col, r, keep_all_cols = TRUE)
  isNumFactor    <- make_is_num_factor(get_meta_col, local_metadata)
  local_physeq   <- make_local_physeq(local_metadata, r)

  # ---- Dynamic UI: distance metric ----
  output$ui_metrics <- renderUI({
    tooltip(
      radioButtons(ns("metrics"), "Distance metric:", inline = TRUE,
                   choices =c('Bray-Curtis' = 'bray', 'Jaccard' = 'jaccard', 
                             'Unifrac' = 'unifrac', 'Weighted Unifrac' = 'wunifrac'),
                   selected = c("bray")),
      "Distance metric for dissimilarity matrix. Unifrac requires phylogeny.",
      placement = "right"
    )
  })

  # ---- Dynamic UI: constrained model parameters (dbRDA only) ----
  output$ui_constrain <- renderUI({
    req(input$ordination)
    if(input$ordination == 'dbRDA'){
      tagList(
        tags$strong(bs_icon("calculator"), " dbRDA Constrained Model"),
        tags$hr(),
        htmltools::p(
          class = "text-info",
          "ℹ️ For dbRDA, the constrained model defines the explanatory variables ",
          "of the ordination itself (model parameters shown as biplot vectors). ",
          "This is distinct from envfit, which post-hoc maps environmental ",
          "variables onto an unconstrained PCoA/NMDS plot."
        ),
        htmltools::p(
          class = "text-warning",
          '⚠️ Samples with missing values in selected env variables will be omitted.'
        ),
        tooltip(
          radioButtons(inputId = ns('param_mode'),
                       label = 'Parameter selection:',
                       choices = c('Picker' = 'picker', 'Manual formula' = 'manual'),
                       selected = 'picker'),
          "Picker: Select variables for ~ formula. Manual: Write formula directly.",
          placement = "right"
        ),
        uiOutput(ns('constr_select')),
        verbatimTextOutput(ns('formula'))
      )
    }
  })

  output$constr_select <- renderUI({
    req(input$param_mode)
    if(input$param_mode == 'picker'){
      # Preselect the sidebar factor(s) used for color (beta_factor) when available
      preselected <- NULL
      if(isTruthy(input$beta_factor) && input$beta_factor %in% colnames(local_metadata())){
        preselected <- input$beta_factor
      }
      shinyWidgets::pickerInput(inputId = ns('constr_picker'),
                                label = 'Select terms to build the formula:',
                                choices = colnames(local_metadata()),
                                selected = preselected,
                                multiple = TRUE,
                                options = pickerOptions(
                                  actionsBox = TRUE,
                                  liveSearch = TRUE,
                                  showContent = FALSE
                                ),
                                choicesOpt = list(
                                  content = make_picker_choices(
                                    local_metadata(), colnames(local_metadata()), show_na = TRUE
                                  )
                                )
      )
    } else if(input$param_mode == 'manual'){
      textInput(ns('constr_manual_formula'), label = 'Enter your own formula:', placeholder = 'spe ~')
    }
  })

  output$formula <- renderPrint({
    get_constr_formula()
  })

  # ---- Constrained formula builder (dbRDA only) ----
  get_constr_formula <- reactive({
    mode <- input$param_mode %||% 'picker'
    if(mode == 'picker'){
      if(is.null(input$constr_picker) || length(input$constr_picker) == 0){
        f <- 'spe ~ 1'
      } else {
        f <- paste0('spe ~ ', paste(input$constr_picker, collapse = ' + '))
      }
    } else if(mode == 'manual'){
      f <- input$constr_manual_formula %||% 'spe ~ 1'
    } else {
      f <- 'spe ~ 1'
    }
    flog.debug('get_constr_formula(): formula=%s', f)
    return(f)
  })

  # ---- Dynamic UI: taxa rank selector + top-N contributors ----
  output$ui_taxa_rank <- renderUI({
    req(input$ordination, local_physeq())
    if(input$ordination %in% c('NMDS', 'PCOA', 'dbRDA', 'CCA', 'RDA')){
      tagList(
        selectInput(
         ns("rank_color"),
         label = "Select rank to color taxa arrows: ",
         choices = rank_names(local_physeq()),
         selected = rank_names(local_physeq())[length(rank_names(local_physeq()))]
        ),
        tooltip(
          numericInput(
            ns("n_top_taxa"),
            label = "Top N contributing taxa:",
            value = 10, min = 1, step = 1
          ),
          "Keep only the N taxa with the largest score vectors (strongest contribution) on the displayed axes.",
          placement = "right"
        )
      )
    }
  })


  # ---- Screeplot ----
  get_screeplot <- reactive({
    req(ord())
    ordination <- isolate(input$ordination)
    flog.debug('get_screeplot(): ordination=%s', ordination)
    if(ordination == 'NMDS'){
      validate(
        need(input$metrics %in% c('bray', 'jaccard'), 'Only bray distance supported for NMDS screeplot.')
      )
      return(goeveg::dimcheckMDS(get_species_table(), distance = input$metrics, k=5))
    } else {
      p <- screeplot(ord())
      return(p)
    }
  })


  output$screeplot <- renderPlot({
    get_screeplot()
  })

  # ---- Dynamic UI: factor & shape pickers ----
  output$ui_beta_factor <- renderUI({
    req(r$var_list())
    tagList(
    shinyWidgets::pickerInput(
       ns("beta_factor"),
       label = "Select factor(s) to color samples and ellipses:",
       choices = r$var_list(),
       selected = r$var_list()[1],
       multiple = FALSE,
       options = pickerOptions(
         actionsBox = TRUE,
         liveSearch = TRUE,
         showContent = FALSE
       ),
       choicesOpt = list(
         content = make_picker_choices(r$sdat(), r$var_list())
       )
     ),
    shinyWidgets::pickerInput(
      ns("beta_shape"),
      label = "Select factor(s) to shape samples:",
      choices = r$var_list(),
      selected = NULL,
      multiple = TRUE,
      options = pickerOptions(
        actionsBox = TRUE,
        liveSearch = TRUE,
        showContent = FALSE,
        maxOptions = 1
      ),
      choicesOpt = list(
        content = make_picker_choices(r$sdat(), r$var_list())
      )
    )
    )
  })

  # ---- Dynamic UI: envfit (PCoA / NMDS only) ----
  output$ui_envfit_box <- renderUI({
    req(input$ordination, local_metadata())
    if(input$ordination %in% c('NMDS', 'PCOA')){
      all_cols <- sort(colnames(local_metadata()))
      # Preselect the factors chosen in the sidebar for color and shape
      preselected <- intersect(c(input$beta_factor, input$beta_shape), all_cols)
      if(length(preselected) == 0) preselected <- NULL

      tagList(
        tags$strong(bs_icon("compass"), " VEGAN envfit (PCoA / NMDS)"),
        tags$hr(),
        htmltools::p(
          class = "text-info",
          "ℹ️ Envfit post-hoc fits environmental vectors/factors onto an ",
          "unconstrained ordination (PCoA or NMDS). This is distinct from a dbRDA ",
          "constrained model, where the explanatory variables are part of the ",
          "ordination itself."
        ),
        shinyWidgets::pickerInput(
          ns("envfit_param"),
          label = "Variables to map onto the ordination:",
          choices = all_cols,
          selected = preselected,
          multiple = TRUE,
          options = pickerOptions(
            actionsBox = TRUE,
            liveSearch = TRUE,
            showContent = FALSE
          ),
          choicesOpt = list(
            content = make_picker_choices(
              local_metadata(), all_cols
            )
          )
        ),
        numericInput(
          ns('envfit_pval'),
          'p-value threshold to be plotted', value = 1,
          min = 0,
          max = 1,
          step = 0.01
        )
      )
    }
  })


  output$ui_envfit_box_res <- renderUI({
    req(input$ordination)
    if(input$ordination %in% c('NMDS', 'PCOA')){
      tagList(
        tags$strong("envfit results"),
        verbatimTextOutput(ns('envfit_res'))
      )
    }
  })

  # ---- Distance matrix computation ----
  physeq_dist <- eventReactive(input$launch_beta, {
    req(input$metrics, local_physeq())
    flog.info('physeq_dist() starting...')
    flog.debug('physeq_dist(): metrics=%s, nsamples=%d',
               input$metrics, phyloseq::nsamples(local_physeq()))
      validate(
        need(!input$metrics %in% c("unifrac", "wunifrac") | !is.null(local_physeq()@phy_tree),
          message = 'Unifrac and wunifrac available if phylogenetic tree available.')
      )

    res <- phyloseq::distance(local_physeq(), method = input$metrics)
    flog.info('physeq_dist() end.')
    return(res)
  })

  # ---- Environmental data scaling ----
  get_env_scaled <- reactive({
    req(r$sdat())
    flog.info('get_env_scaled() starting...')
    env <- local_metadata()
    data.split <- PCAmixdata::splitmix(env)
    flog.debug('get_env_scaled(): ncol_quant=%d, ncol_qual=%d',
               length(data.split$col.quant), length(data.split$col.qual))
    env[,data.split$col.quant] <- scale(env[,data.split$col.quant], center = T, scale = T)
    flog.info('get_env_scaled() end.')
    return(env)
  })

  # ---- Species table (Hellinger-transformed) ----
  get_species_table <- reactive({
    req(local_physeq())
    flog.debug('get_species_table(): nsamples=%d, ntaxa=%d',
               phyloseq::nsamples(local_physeq()), phyloseq::ntaxa(local_physeq()))
    spe <- veganifyOTU(local_physeq())
    spe <- vegan::decostand(spe, method = 'hell')
  })

  # ---- Ordination computation ----
  ord <- eventReactive(input$launch_beta, {
    req(input$ordination, local_physeq())

    flog.info('ord() starting...')
    flog.debug('ord(): ordination=%s, metrics=%s', input$ordination, input$metrics)
    withProgress(message = 'Computing ordination...', min=0, max=10, value = 0,{
      setProgress(value = 4, detail = input$ordination)
      if(input$ordination == 'NMDS'){

        if(input$metrics %in% c('bray', 'jaccard')){
          res <- vegan::metaMDS(get_species_table(), k = 5, distance = input$metrics, wascores=TRUE, trace=FALSE, autotransform = FALSE)
        } else if(input$metrics %in% c('unifrac', 'wunifrac')){
          res <- vegan::metaMDS(comm = physeq_dist(), wascores=TRUE, trace=FALSE, autotransform = FALSE, k = 5)
        }
      } else if(input$ordination == 'PCOA'){

        spe <- get_species_table()
        if(input$metrics %in% c('bray', 'jaccard')){
          res <- vegan::capscale(spe ~ 1, spe, distance = input$metrics)
        } else if(input$metrics %in% c('unifrac', 'wunifrac')){
          dist <- physeq_dist()
          res <- vegan::capscale(dist ~ 1, comm = spe)
        }
      } else if(input$ordination == 'RDA'){
        env <- get_env_scaled()
        spe <- get_species_table()
        flog.debug('ord() RDA: formula=%s', get_constr_formula())
        res <- vegan::rda(as.formula(get_constr_formula()), data = env, na.action = 'na.omit')
      } else if(input$ordination == 'CCA'){
        env <- get_env_scaled()
        spe <- get_species_table()
        flog.debug('ord() CCA: formula=%s', get_constr_formula())
        res <- vegan::cca(as.formula(get_constr_formula()), data = env, na.action = 'na.omit')
      } else if(input$ordination == 'dbRDA'){
        env <- get_env_scaled()
        flog.debug('ord() dbRDA: formula=%s, metrics=%s', get_constr_formula(), input$metrics)
        if(input$metrics %in% c('bray', 'jaccard')){
          spe <- get_species_table()
          res <- vegan::capscale(as.formula(get_constr_formula()), data = env, na.action = 'na.omit', distance = input$metrics)
        } else if(input$metrics %in% c('unifrac', 'wunifrac')){
          spe <- physeq_dist()
          res <- vegan::capscale(as.formula(get_constr_formula()), data = env, na.action = 'na.omit', distance = input$metrics, comm = veganifyOTU(local_physeq()))
        }
      } else{
        res <- phyloseq::ordinate(physeq= local_physeq(), distance = physeq_dist(), method= input$ordination)
      }
      setProgress(value = 10, detail = 'done')
    })

    flog.info('ord() end.')
    return(res)
  })

  # ---- Site coordinates extraction ----
  get_sites_coord <- reactive({
    req(ord(), local_metadata(), get_meta_col())
    flog.info('get_sites_coord() starting...')
    ordination <- isolate(input$ordination)
    if(ordination %in% c('PCA', 'RDA', 'PCOA', 'dbRDA')){
      nmds_coord <- vegan::scores(ord(), display='sites', correlation = TRUE) %>%
        as_tibble(rownames="sample.id")
    } else if(ordination %in% c('CCA')){
      nmds_coord <- vegan::scores(ord(), display='sites', hill = TRUE) %>%
        as_tibble(rownames="sample.id")
    } else {
      nmds_coord <- vegan::scores(ord(), display='sites') %>%
        as_tibble(rownames="sample.id")
    }

    shape_cols <- if(is.null(input$beta_shape)) character(0) else input$beta_shape
    nmds_coord <- nmds_coord %>%
                    inner_join(., local_metadata() %>%
                                 dplyr::select(all_of(c('sample.id', get_meta_col(), shape_cols))), by="sample.id")
    flog.debug('get_sites_coord(): %d sites, %d columns', nrow(nmds_coord), ncol(nmds_coord))
    flog.info('get_sites_coord() end.')
    return(nmds_coord)
  })

  # ---- Species coordinates extraction ----
  get_species_coord <- reactive({
    req(local_physeq(), ord(), r$rank_glom(), input$rank_color)
    flog.info('get_species_coord() starting...')
    ordination <- isolate(input$ordination)
    flog.debug('get_species_coord(): ordination=%s, rank_glom=%s, rank_color=%s',
               ordination, r$rank_glom(), input$rank_color)
    if(ordination %in% c('PCA', 'RDA', 'PCOA', 'dbRDA')){
      if(ordination %in% c('dbRDA') && input$metrics %in% c('unifrac', 'wunifrac')){
        nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='specie', scalling = 2) %>% as_tibble(rownames=r$rank_glom())
      } else {
        nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='specie', correlation = TRUE) %>% as_tibble(rownames=r$rank_glom())
      }
    } else if(ordination %in% c('CCA')){
      nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='specie', hill = TRUE) %>% as_tibble(rownames=r$rank_glom())
    } else if(ordination %in% c('NMDS') && input$metrics %in% c('unifrac', 'wunifrac')){
      spe <- veganifyOTU(local_physeq())
      spe <- vegan::decostand(spe, method = 'hell')
      res <- ord()
      vegan::sppscores(res) <- spe
      nmds_coord <- vegan::scores(res, choices=c(1,2), display='specie') %>% as_tibble(rownames=r$rank_glom())
    } else {
      nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='specie') %>% as_tibble(rownames=r$rank_glom())
    }

    if(r$rank_glom() == 'ASV'){
      nmds_coord <- nmds_coord %>%
                      inner_join(., tax_table(local_physeq()) %>%
                                   as.data.frame() %>%
                                   as_tibble(rownames="ASV") %>%
                                   select(ASV, input$rank_color), by='ASV')
    } else{
      nmds_coord <- nmds_coord %>%
                      inner_join(., tax_table(local_physeq()) %>%
                                   as.data.frame() %>%
                                   as_tibble(rownames=r$rank_glom()) %>%
                                   select(r$rank_glom(), input$rank_color), by=r$rank_glom())
    }
    flog.debug('get_species_coord(): %d species extracted', nrow(nmds_coord))
    flog.info('get_species_coord() end.')
    return(nmds_coord)
  })

  # ---- Envfit computation ----
  output$envfit_res <- renderPrint({
    get_env_fit()
  })

  get_env_fit <- eventReactive( input$launch_beta, {
    flog.info('get_env_fit() starting...')
    req(ord())
    env <- get_env_scaled()
    flog.debug('get_env_fit(): envfit_param=[%s]',
               paste(input$envfit_param, collapse = ', '))
    env <- env[, input$envfit_param, drop=F]
    en <- vegan::envfit(ord(), env, na.rm = T)
    flog.info('get_env_fit() end.')
    return(en)
  })

  # ---- Axis names ----
  output$ui_axe_x <- renderUI({
    req(ord())
    shinyWidgets::pickerInput(inputId = ns('axe_x'), choices = get_axis_names(), selected = get_axis_names()[1], width = 100)
  })

  output$ui_axe_y <- renderUI({
    req(ord())
    shinyWidgets::pickerInput(inputId = ns('axe_y'), choices = get_axis_names(), selected = get_axis_names()[2], width = 100)
  })

get_axis_names <- reactive({
  flog.info('get_axis_names() starting...')
  ordination <- isolate(input$ordination)

  # Safe fallback before ord() available
  if (is.null(ord())) {
    fallback <- switch(ordination %||% "PCOA",
      "NMDS" = c("MDS1", "MDS2"),
      "PCOA" = c("PC1", "PC2"),
      "RDA" = c("RDA1", "RDA2"),
      "CCA" = c("CCA1", "CCA2"),
      "dbRDA" = c("dbRDA1", "dbRDA2"),
      c("PC1", "PC2")
    )
    flog.debug('get_axis_names(): ord NULL, fallback=%s', paste(fallback, collapse=', '))
    return(fallback)
  }

  axes <- colnames(vegan::scores(ord(), display = 'sites'))

  # Fallback defaults by ordination type
  if (length(axes) == 0) {
    fallback <- switch(ordination,
      "NMDS" = c("MDS1", "MDS2"),
      "PCOA" = c("PC1", "PC2"),
      "RDA" = c("RDA1", "RDA2"),
      "CCA" = c("CCA1", "CCA2"),
      "dbRDA" = c("dbRDA1", "dbRDA2"),
      c("PC1", "PC2")
    )
    axes <- fallback
    flog.debug('get_axis_names(): using fallback axes=%s', paste(axes, collapse = ', '))
  } else {
    flog.debug('get_axis_names(): axes=[%s]', paste(axes, collapse = ', '))
  }
  flog.info('get_axis_names() end.')
  return(axes)
})

  # ---- Ordination plot ----
  base_plot <- reactive({
    req(ord(), input$beta_factor)
    flog.info('base_plot() starting...')
    
    # Robust axes fallback chain
    axe_x <- input$axe_x %||% get_axis_names()[1] %||% "PC1"
    axe_y <- input$axe_y %||% get_axis_names()[2] %||% "PC2"
    
    # Validate axes ready
    validate(need(!is.null(axe_x) && axe_x != "" && !is.null(axe_y) && axe_y != "", 
                  "Ordination axes not ready. Click '🚀 Run Ordination Analysis' first."))
    
    # Default plot_type safe
    plot_types <- input$plot_type %||% c("samples")
    
    ordination <- isolate(input$ordination)
    flog.debug('base_plot(): plot_types=%s, axe_x=%s, axe_y=%s',
               paste(plot_types, collapse = ', '), axe_x, axe_y)

    withProgress(message = 'Plotting...', min=0, max=10, value = 0,{
      p <- ggplot2::ggplot()
      if('samples' %in% input$plot_type){
        flog.info('base_plot() plotting samples...')
        setProgress(value = 4, detail = 'plotting samples')
        sites_coord <- get_sites_coord()
        shape_mapping <- if(is.null(input$beta_shape)) {
          aes(x=!!sym(axe_x), y=!!sym(axe_y), 
              color=.data[[get_meta_col()]], 
              text = paste('sample.id:', sample.id))
        } else {
          aes(x=!!sym(axe_x), y=!!sym(axe_y), 
              color=.data[[get_meta_col()]], 
              text = paste('sample.id:', sample.id),
              shape = .data[[input$beta_shape]])
        }
        p <- p +
          geom_point(
            data = sites_coord,
            mapping = shape_mapping,
            size = 5
          )
        if (isNumFactor()) {
          p <- p + paletteer::scale_color_paletteer_c(
            palette = r$numeric_palettes()[[input$beta_factor]]
          ) + ggnewscale::new_scale_color()
        } else {
          p <- p +
            scale_color_manual(
              values = r$factor_colors()[[input$beta_factor]],
              drop = FALSE
            ) +
            ggnewscale::new_scale_color() +
            stat_ellipse(
              data = sites_coord,
              mapping = aes(x=!!sym(axe_x), y=!!sym(axe_y),
              group = !!sym(get_meta_col()),
              color = !!sym(get_meta_col()))
            ) +
            scale_color_manual(
              values = r$factor_colors()[[input$beta_factor]],
              drop = FALSE
            ) +
            ggnewscale::new_scale_color()
        }
      }

      if(ordination == "PCOA"){
        eig1 <- eigenvals(ord())
        percent1 <- eig1/sum(eig1)*100
        p <- p + xlab(glue::glue("{input$axe_x} ({round(percent1[input$axe_x], 2)} %)")) +
          ylab(glue::glue("{input$axe_y} ({round(percent1[input$axe_y], 2)} %)"))
      }

      if ('taxa' %in% input$plot_type){
        flog.info('base_plot() plotting taxa...')
        species_coord <- get_species_coord()
        rank_id <- r$rank_glom()
        n_top   <- input$n_top_taxa %||% 10
        # Label arrows with the taxon name: the agglomeration rank already holds
        # a name, except when glommed at ASV level — then fall back to the
        # chosen color rank so labels read as taxa, not ASV ids.
        label_col <- if(rank_id == 'ASV') input$rank_color else rank_id

        # Keep the N most contributing taxa: rank by the length of their score
        # vector in the displayed plane (sqrt(x^2 + y^2)).
        species_coord <- species_coord %>%
          dplyr::mutate(.norm = sqrt(.data[[axe_x]]^2 + .data[[axe_y]]^2)) %>%
          dplyr::arrange(dplyr::desc(.norm)) %>%
          utils::head(n_top)

        validate(need(nrow(species_coord) > 0, "No taxa scores available for this ordination."))

        # Biplot scaling: stretch taxa vectors to roughly fill the sample cloud
        # so the arrows are readable alongside the sample points.
        sites_coord <- get_sites_coord()
        denom_x <- max(abs(species_coord[[axe_x]]))
        denom_y <- max(abs(species_coord[[axe_y]]))
        multiplier <- if(denom_x > 0 && denom_y > 0){
          min(max(abs(sites_coord[[axe_x]])) / denom_x,
              max(abs(sites_coord[[axe_y]])) / denom_y) * 0.8
        } else 1
        species_coord <- species_coord %>%
          dplyr::mutate(.xend = .data[[axe_x]] * multiplier,
                        .yend = .data[[axe_y]] * multiplier)

        p <- p +
          geom_segment(
            data = species_coord,
            aes(x = 0, y = 0, xend = .xend, yend = .yend, color = .data[[input$rank_color]]),
            arrow = grid::arrow(length = grid::unit(0.02, "npc")),
            linewidth = 0.6, alpha = 0.8
          ) +
          scale_color_manual(values = r$taxa_colors()[[input$rank_color]], drop = FALSE) +
          ggrepel::geom_text_repel(
            data = species_coord,
            aes(x = .xend, y = .yend, label = .data[[label_col]], color = .data[[input$rank_color]]),
            size = 5, fontface = "italic", max.overlaps = 20, show.legend = FALSE
          )
      }

      if ('env' %in% input$plot_type){
        flog.info('base_plot() plotting env...')
        if(ordination %in% c('RDA', 'CCA', 'dbRDA')){
          scr <- vegan::scores(ord())
          if(!is.null(scr$biplot)){
            en_coord_cont <- as.data.frame(scr$biplot)
            en_coord_cont <- en_coord_cont[rownames(en_coord_cont) %in% colnames(local_metadata()),]
            if(nrow(en_coord_cont) > 0){
              p <- p + geom_segment(aes(x = 0, y = 0, xend = !!sym(input$axe_x), yend = !!sym(input$axe_y)),
                                    data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = grid::arrow()) +
              geom_text(data = en_coord_cont, aes(x = !!sym(axe_x), y = !!sym(axe_y)), colour = "grey30",
                        fontface = "bold", label = row.names(en_coord_cont))
            }
          }
          if(!is.null(scr$centroids)){
            en_coord_cat <- as.data.frame(scr$centroids)
            p <- p + geom_point(data = en_coord_cat, aes(x = !!sym(input$axe_x), y = !!sym(input$axe_y)),
                                shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
              geom_text(data = en_coord_cat, aes(x = !!sym(input$axe_x), y = !!sym(input$axe_y)),
                        label = row.names(en_coord_cat), colour = "navy", fontface = "bold")
          }
        } else if(ordination %in% c('PCOA', 'NMDS')){
          if(is.null(input$envfit_param)){
            shinyalert::shinyalert(title = "Oops", text="You need to use ENVFIT module to use env type.", type='error')
            return(NULL)
          }

          en <- get_env_fit()
          if(!is.null(en$vectors)){
            en_coord_cont <- as.data.frame(vegan::scores(en, "vectors"))
            p <- p + geom_segment(aes(x = 0, y = 0, xend = !!sym(input$axe_x), yend = !!sym(input$axe_y)),
                                  data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = grid::arrow()) +
              geom_text(data = en_coord_cont, aes(x = !!sym(axe_x), y = !!sym(axe_y)), colour = "grey30",
                        fontface = "bold", label = row.names(en_coord_cont))
          }
          if(!is.null(en$factors)){
            en_coord_cat <- as.data.frame(vegan::scores(en, "factors"))
            p <- p + geom_point(data = en_coord_cat, aes(x = !!sym(input$axe_x), y = !!sym(input$axe_y)),
                                shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
              geom_text(data = en_coord_cat, aes(x = !!sym(input$axe_x), y = !!sym(input$axe_y)),
                        label = row.names(en_coord_cat), colour = "navy", fontface = "bold")
          }
        }
      }
    })

    flog.info('base_plot() end.')
    return(p)
  })

  # ---- Plot rendering ----
  observeEvent(list(input$axe_x, input$axe_y, input$plot_type, ord()), {
    output$beta_ggplot <- renderPlot({
      base_plot()
    })
  }, ignoreNULL = FALSE)

  # ---- PERMANOVA formula (built from sidebar color + shape selections) ----
  get_formula <- reactive({
    req(input$metrics, get_meta_col())
    terms <- c('Depth', get_meta_col())
    # Add shape factor(s) only if set and distinct from the color factor
    if(!is.null(input$beta_shape) && length(input$beta_shape) > 0){
      extra <- setdiff(input$beta_shape, get_meta_col())
      terms <- c(terms, extra)
    }
    form <- paste('dist ~', paste(terms, collapse = ' + '))
    flog.debug('get_formula(): %s', form)
    return(form)
  })

  # ---- Dispersion tests ----
  get_dispersion_res <- eventReactive(input$launch_beta,{
    req(physeq_dist(), get_meta_col())
    validate(need(!isNumFactor(),
                  "Dispersion test requires a categorical factor."))
    withProgress(message = "Computing dispersion (betadisper)…", {
      flog.info('get_dispersion_res() starting...')
      dist <- physeq_dist()
      # Group vector ordered to the distance-matrix samples by name (betadisper
      # pairs the grouping to the dissimilarities positionally).
      groups <- local_metadata()[attr(dist, "Labels"), get_meta_col()]
      res <- vegan::betadisper(dist, groups)
      flog.info('get_dispersion_res() end.')
    })
    return(res)
  })

  get_dispersion_anova <- reactive({
    req(get_dispersion_res())
    flog.info('get_dispersion_anova() starting...')
    res <- anova(get_dispersion_res())
    flog.info('get_dispersion_anova() end.')
    return(res)
  })

  get_dispersion_tukey <- reactive({
    req(get_dispersion_res())
    flog.info('get_dispersion_tukey() starting...')
    res <- TukeyHSD(get_dispersion_res())
    flog.info('get_dispersion_tukey() end.')
    return(res)
  })

  # ---- PERMANOVA (adonis) ----
  # adonis2 pairs the LHS distance matrix to the RHS `data` rows *by position*,
  # so align both to the distance-matrix sample order (by name) and drop any
  # sample with a missing value in a model term -- subsetting the distance
  # matrix to the same samples so the two never disagree in size or order.
  get_adonis_res <- eventReactive(input$launch_beta, {
    req(physeq_dist(), get_formula(), get_meta_col(), ord())
    withProgress(message = "Computing PERMANOVA (adonis2)…", {
      flog.info('get_adonis_res() starting...')
      flog.debug('get_adonis_res(): formula=%s', get_formula())
      dist  <- physeq_dist()
      samp  <- attr(dist, "Labels")
      mdata <- local_metadata()[samp, , drop = FALSE]
      # Sequencing-depth covariate, matched to each sample by name (not position).
      mdata$Depth <- sample_sums(r$phyloseq_filtered())[samp]
      term_vars <- intersect(c("Depth", get_meta_col(), input$beta_shape), colnames(mdata))
      keep <- samp[stats::complete.cases(mdata[, term_vars, drop = FALSE])]
      flog.debug('get_adonis_res(): %d of %d samples kept after dropping NA terms',
                 length(keep), length(samp))
      validate(need(length(keep) >= 3,
                    "Not enough samples with complete data for the selected PERMANOVA terms."))
      dist  <- stats::as.dist(as.matrix(dist)[keep, keep])
      mdata <- mdata[keep, , drop = FALSE]
      # by = "terms": assess each term sequentially so every term's R2 is
      # reported on its own row (default by = NULL collapses to one model line).
      res <- vegan::adonis2(as.formula(get_formula()), data = mdata,
                            permutations = 1000, by = "terms")
      flog.info('get_adonis_res() end.')
    })
    return(data.frame(res))
  })

  get_pairwise_res <- eventReactive(input$launch_beta, {
    req(physeq_dist(), get_meta_col(), local_metadata())
    validate(need(!isNumFactor(),
                  "Pairwise PERMANOVA requires a categorical factor."))
    withProgress(message = "Computing pairwise PERMANOVA…", {
      dist <- physeq_dist()
      # Group vector ordered to the distance-matrix samples by name, so the
      # `factors %in% ...` row selection inside pairwise.adonis stays aligned.
      groups <- local_metadata()[attr(dist, "Labels"), get_meta_col()]
      flog.info('get_pairwise_res() starting...')
      flog.debug('get_pairwise_res(): factor=%s, nsamples=%d',
                 get_meta_col(), length(groups))
      res <- pairwise.adonis(dist, groups, p.adjust.m = "fdr")
      flog.info('get_pairwise_res() end.')
    })
    return(res)
  })

  # ---- PERMANOVA outputs ----
  output$adonis_formula <- renderText({
    get_formula()
  })

  output$adonistest <- DT::renderDataTable({
    round_df(get_adonis_res(), 4)
  })

  output$adonispairwisetest <- DT::renderDataTable({
    round_df(get_pairwise_res(), 4)
  })

  # ---- Dispersion outputs ----
  dfdisper <- eventReactive(input$launch_beta, {
    req(get_dispersion_res())
    flog.info('dfdisper() starting...')

    df1 = cbind.data.frame(distances = get_dispersion_res()$distances, group = get_dispersion_res()$group)

    # Always order factor levels with natural sort (previous optional behavior made default)
    flog.debug('dfdisper(): ordering factor levels with mixedsort')
    df1$group = factor(df1$group, levels = gtools::mixedsort(levels(df1$group)))

    flog.info('dfdisper() end.')
    return(df1)
  })

  # ---- Eager test population ----
  # Force PERMANOVA / pairwise / dispersion to compute on Run, inside the
  # input_task_button "Computing…" window, instead of silently the first time the
  # user opens the (hidden) PERMANOVA / Dispersion tabs. Each eventReactive caches,
  # so subsequent tab navigation just reads the result. Wrapped in try() because
  # get_pairwise_res()/get_dispersion_res() validate() on factor type and sample
  # count; those conditions still surface gracefully when the output renders.
  observeEvent(input$launch_beta, {
    req(physeq_dist())
    try(get_adonis_res(),     silent = TRUE)   # req(ord()) pulls the ordination too
    try(get_pairwise_res(),   silent = TRUE)
    try(get_dispersion_res(), silent = TRUE)
    try(dfdisper(),           silent = TRUE)
  }, ignoreInit = TRUE)

  output$dispersionPlot <- renderPlotly({
   df1 <- dfdisper()
   plot_ly(df1, x = ~group, y = ~distances,
           color = ~group, type = 'box', colors = r$factor_colors()[[input$beta_factor]]) %>%
     layout(title="", yaxis = list(title = "Distance to centroid"), xaxis = list(title = 'Group'), barmode = 'stack') %>%
    config(toImageButtonOptions = list(format = "svg"))
 })

  output$dispersionTable <- DT::renderDataTable({
    round_df(as.data.frame(get_dispersion_anova()), 4)
  })

  output$dispersionTukey <- DT::renderDataTable({
    round_df(as.data.frame(get_dispersion_tukey()$group), 4)
  })
  })
}

## To be copied in the UI
# mod_beta_ui("beta_ui_1")

## To be copied in the server
# mod_beta_server("beta_ui_1", r = r)
