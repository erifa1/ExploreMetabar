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
#' @importFrom bslib navset_card_underline nav_panel accordion accordion_panel
#' @import PCAmixdata
#' @import shinycustomloader
#' @import shinyWidgets
#' 

mod_beta_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "🎯 Guided Ordination Workflow",
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
          tooltip(
            radioButtons(ns('ordi_type'), 'Type:',
                         inline = TRUE,
                         choices = c('Constrained', 'Distance-based'),
                         selected = 'Distance-based'),
            "Constrained uses environmental variables. Distance-based uses dissimilarity matrices.",
            placement = "right"
          ),
          uiOutput(ns('ui_metrics')),
          tooltip(
            radioButtons(ns("ordination"), "Method:", inline = TRUE,
                         choices = '', selected = ''),
            "Choose ordination technique based on your data type and research question.",
            placement = "right"
          ),
          uiOutput(ns('ui_beta_factor')),
          
          # Primary action button
          div(
            style = "margin: 1.5rem 0;",
            tooltip(
              actionButton(ns("launch_beta"), "🚀 Run Ordination Analysis", 
                          icon = bs_icon("play-fill"), 
                          class = "btn-primary w-100 btn-lg"),
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
        
        # ANOVA only (RDA/CCA specific)
        accordion_panel(
          "Model ANOVA", 
          icon = bs_icon("bezier"),
          tooltip(bs_icon("question-circle"), "ANOVA results for constrained models (RDA/CCA)", placement = "right"),
          uiOutput(ns('ui_anova_box'))
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
      nav_panel("Screeplot", plotOutput(ns('screeplot'))),
      
      # PERMANOVA (moved from sidebar ui_permanova_box)
      nav_panel(
        title = "PERMANOVA", 
        icon = bs_icon("calculator"),
        h3("Adonis formula:"),
        verbatimTextOutput(ns("adonis_formula")),
        hr(),
        h3("Adonis results:"),
        DT::dataTableOutput(ns('adonistest')),
        hr(),
        h3("Pairwise Adonis results:"),
        DT::dataTableOutput(ns('adonispairwisetest'))
      ),
      
      # Dispersion (moved from sidebar ui_dispersion)
      nav_panel(
        title = "Dispersion", 
        icon = bs_icon("funnel"),
        card(
          full_screen = TRUE,
          card_header("Dispersion Analysis"),
          checkboxInput(ns("order1"), label = "Automatic order factor", value = TRUE),
          navset_card_underline(
            nav_panel("Boxplots", plotlyOutput(ns("dispersionPlot"))),
            nav_panel("Anova", DT::dataTableOutput(ns("dispersionTable"))),
            nav_panel("TukeyHSD", DT::dataTableOutput(ns("dispersionTukey")))
          )
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
#'
#' @noRd
mod_beta_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns

  # Status reactives for workflow value_boxes (minimal, no new deps)
  output$status_type <- renderText({
    if(isTruthy(input$ordi_type)) input$ordi_type else "Select type"
  })
  
  output$status_method <- renderText({
    if(isTruthy(input$ordination)) input$ordination else "Select method"
  })
  
  output$status_ready <- renderText({
    ready <- isTruthy(input$beta_factor) && 
             isTruthy(input$ordination) && 
             input$ordination != ""
    if(ready) "Ready ✓" else "Complete setup"
  })
  
  # Simple accordion hints (no shinyjs needed)
  observeEvent(input$ordination, {
    # Tooltip already guides users to relevant sections
    NULL
  })

  # ---- Shared reactives setup ----
  factor_input   <- reactive({ input$beta_factor })
  get_meta_col   <- make_get_meta_col(factor_input, r)
  local_metadata <- make_local_metadata(factor_input, get_meta_col, r, keep_all_cols = TRUE)
  isNumFactor    <- make_is_num_factor(get_meta_col, local_metadata)
  local_physeq   <- make_local_physeq(local_metadata, r)

  # ---- Dynamic UI: distance metric ----
  output$ui_metrics <- renderUI({
    req(input$ordi_type)
    tooltip(
      radioButtons(ns("metrics"), "Distance metric:", inline = TRUE,
                   choices =c('Bray-Curtis' = 'bray', 'Jaccard' = 'jaccard', 
                             'Unifrac' = 'unifrac', 'Weighted Unifrac' = 'wunifrac'),
                   selected = c("bray")),
      "Distance metric for dissimilarity matrix. Unifrac requires phylogeny.",
      placement = "right"
    )
  })

  # ---- Dynamic UI: constrained model parameters ----
  output$ui_constrain <- renderUI({
    req(input$ordination)
    if(input$ordination %in% c('RDA','CCA', 'dbRDA')){
      tooltip(
        card(card_header(bs_icon("calculator"), "Constrained Model"),
          htmltools::p(
            class = "text-warning", 
            '⚠️ Samples with missing environmental values will be omitted.'
          ),
          tooltip(
            radioButtons(inputId = ns('param_mode'),
                         label = 'Parameter selection:',
                         choices = c('Picker' = 'picker', 'Auto (ordiR2step)' = 'ordiR2step', 'Manual formula' = 'manual')),
            "Picker: Select variables manually. ordiR2step: Automatic stepwise selection. Manual: Write formula.",
            placement = "right"
          ),
          uiOutput(ns('constr_select')),
          uiOutput(ns('ordiR2_out')),
          verbatimTextOutput(ns('formula')),
          uiOutput(ns('ui_ordiR2_btn'))
        ),
        "Build constrained ordination model (~ environmental variables)",
        placement = "right"
      )
    }
  })

  output$ui_ordiR2_btn <- renderUI({
    if(input$param_mode == 'ordiR2step'){
      shinyWidgets::actionBttn(inputId = ns('ordiR2_btn'), label = 'launch', size = 'sm')
    }
  })

  output$ordiR2_out <- renderUI({
    req(input$param_mode)
    if(input$param_mode == 'ordiR2step'){
      verbatimTextOutput(ns('constr_ordiR2step_out'))
    }
  })

  output$constr_select <- renderUI({
    req(input$param_mode)
    if(input$param_mode %in% c('picker', 'ordiR2step')){
      shinyWidgets::pickerInput(inputId = ns('constr_picker'),
                                label = 'Select terms to create a formula',
                                choices = colnames(local_metadata()),
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

  # ---- ordiR2step computation ----
  get_ordiR2step <- eventReactive(input$ordiR2_btn, {
    req(local_metadata(), local_physeq())
    flog.info('get_ordiR2step() starting...')
    flog.debug('get_ordiR2step(): ordination=%s, selected_terms=[%s]',
               input$ordination, paste(input$constr_picker, collapse = ', '))

    spe <- veganifyOTU(local_physeq())
    spe <- vegan::decostand(spe, method = 'hell')

    env <- local_metadata()[, c(input$constr_picker, 'sample.id')]
    nsample_before <- nrow(env)
    env <- na.omit(env)
    nsample_after <- nrow(env)

    flog.debug('get_ordiR2step(): samples before NA removal=%d, after=%d',
               nsample_before, nsample_after)

    spe <- spe[env$sample.id,]
    if('sample.id' %in% colnames(env)){
      env[,'sample.id'] <- NULL
    }
    validate(
      need(nsample_after > 0, message = 'Too much NAs in your env variables. No sample left.')
    )

    if(input$ordination == 'RDA'){
      mod0 <- vegan::rda(spe ~ 1, data = env, na.action = 'na.omit')
      mod1 <- vegan::rda(spe ~ ., data = env, na.action = 'na.omit')
    } else if(input$ordination == 'CCA'){
      mod0 <- vegan::cca(spe ~ 1, data = env, na.action = 'na.omit')
      mod1 <- vegan::cca(spe ~ ., data = env, na.action = 'na.omit')
    } else if(input$ordination == 'dbRDA'){
      validate(
        need(input$metrics %in% c('bray', 'jaccard'), message = 'ordiR2step works only with bray and jaccard distances.')
      )
      if(input$metrics %in% c('bray', 'jaccard')){
        mod0 <- vegan::capscale(spe ~ 1, data = env, na.action = "na.omit", distance = input$metrics)
        mod1 <- vegan::capscale(spe ~ ., data = env, na.action = "na.omit", distance = input$metrics)
      }
    }
    res <- list()

    res$sel <- vegan::ordiR2step(mod0, scope = formula(mod1), R2scope = FALSE, trace = FALSE)
    res$lostsamples <- nsample_before - nsample_after
    flog.info('get_ordiR2step() end.')
    return(res)
  })

  output$constr_ordiR2step_out <- renderPrint({
    res <- get_ordiR2step()
    if(res$lostsamples != 0){
      print(paste0('warn: ', res$lostsamples, ' samples omitted due to NAs in envirnomental variables.'))
    }
    print(res$sel$anova)
  })

  output$formula <- renderPrint({
    get_constr_formula()
  })

  # ---- Constrained formula builder ----
  get_constr_formula <- reactive({
    req(input$param_mode)
    if(input$param_mode == 'picker'){
      if(length(input$constr_picker) == 0){
        f <- 'spe ~ 1'
      } else {
        f <- paste0('spe ~ ', paste(input$constr_picker, collapse = ' + '))
      }
    } else if(input$param_mode == 'ordiR2step'){
      f <- as.character(formula(get_ordiR2step()$sel))
    } else if(input$param_mode == 'manual'){
      f <- input$constr_manual_formula
    }
    flog.debug('get_constr_formula(): formula=%s', f)
    return(f)
  })

  # ---- Dynamic UI: taxa rank selector ----
  output$ui_taxa_rank <- renderUI({
    req(input$ordination, local_physeq())
    if(input$ordination %in% c('NMDS', 'PCOA', 'dbRDA', 'CCA', 'RDA')){
      selectInput(
       ns("rank_color"),
       label = "Select rank to color taxa points: ",
       choices = rank_names(local_physeq()),
       selected = rank_names(local_physeq())[length(rank_names(local_physeq()))]
      )
    }
  })


  output$pairwise_res <- renderUI({
    req(get_meta_col())
    if(! isNumFactor() && get_meta_col() != 'sample.id'){
      NULL  # Will be handled in main PERMANOVA navset
    }
  })

  # ---- Screeplot ----
  get_screeplot <- reactive({
    req(ord())
    flog.debug('get_screeplot(): ordination=%s', input$ordination)
    if(input$ordination == 'NMDS'){
      validate(
        need(input$metrics %in% c('bray', 'jaccard'), 'Only bray distance supported for NMDS screeplot.')
      )
      goeveg::dimcheckMDS(get_species_table(), distance = input$metrics, k=5)
    } else{
      p <- screeplot(ord())
    }
    return(p)
  })


  output$screeplot <- renderPlot({
    get_screeplot()
  })

  # ---- Dynamic UI: dispersion results ----
  output$ui_dispersion <- renderUI({
    req(get_meta_col())
    if(! isNumFactor() && get_meta_col() != 'sample.id'){
      tagList(
        card(
          card_header("Dispersion Analysis"),
          checkboxInput(ns("order1"), label = "Automatic order factor", value = TRUE)
        ),
        navset_card_underline(
          title = "Dispersion Results",
          full_screen = TRUE,
          nav_panel("Boxplots Distance", plotlyOutput(ns("dispersionPlot"))),
          nav_panel("Anova on Dispersion", DT::dataTableOutput(ns("dispersionTable"))),
          nav_panel("TukeyHSD Test", DT::dataTableOutput(ns("dispersionTukey")))
        )
      )
    }
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

  # ---- Dynamic UI: envfit ----
  output$ui_envfit_box <- renderUI({
    req(input$ordination, local_metadata())
    if(input$ordination %in% c('NMDS', 'PCOA')){
      card(card_header("VEGAN envfit"),
        htmltools::p('The envfit function fits environmental vectors or factors onto an ordination.'),
        shinyWidgets::pickerInput(
          ns("envfit_param"),
          label = "Select one or more factor to use in envfit module:",
          choices = sort(colnames(local_metadata())),
          multiple = TRUE,
          options = pickerOptions(
            actionsBox = TRUE,
            liveSearch = TRUE,
            showContent = FALSE
          ),
          choicesOpt = list(
            content = make_picker_choices(
              local_metadata(), sort(colnames(local_metadata()))
            )
          )
        ),
        numericInput(
          ns('envfit_pval'),
          'p-value threshold to be plotted', value = 1,
          min = 0,
          max = 1,
          step = 0.01
        ),
      )
    }
  })


  output$ui_envfit_box_res <- renderUI({
    req(input$ordination)
    if(input$ordination %in% c('NMDS', 'PCOA')){
      card(card_header("envfit results"),
        verbatimTextOutput(ns('envfit_res'))
      )
    }
  })

  # ---- Ordination type observer ----
  observeEvent(input$ordi_type, {
    if(input$ordi_type == 'Unconstrained'){
      ch <- c('PCA', 'CA', 'DCA')
    } else if (input$ordi_type == 'Constrained'){
      ch <- c('RDA', 'CCA')
    } else{
      ch <- c('PCOA', 'NMDS', 'dbRDA')
    }
    updateRadioButtons(session,
                       "ordination",
                       choices = ch ,
                       inline = T)
  })

  # ---- Dynamic UI: ANOVA box ----
  output$ui_anova_box <- renderUI({
    if(input$ordination %in% c('RDA', 'CCA')){
      card(
        card_header("Anova on RDA/CCA results"),
        h3("Anova results on model"),
        verbatimTextOutput(ns('anova_res')),
        h3("Anova results on axis"),
        verbatimTextOutput(ns('anova_axis_res')),
        h3("Anova results on terms"),
        verbatimTextOutput(ns('anova_term_res')),
        h3('Anova results on contrasts'),
        verbatimTextOutput(ns('anova_contrasts_res'))
      )
    }
  })

  # ---- ANOVA computations ----
  get_anova_contrast <- eventReactive(input$launch_beta, {
    flog.info('get_anova_contrast() starting...')
    res <- anova(ord(), permutations = permute::how(nperm = 999), by = 'onedf')
    flog.info('get_anova_contrast() end.')
    return(res)
  })

  get_anova_term <- eventReactive(input$launch_beta, {
    flog.info('get_anova_term() starting...')
    res <- anova(ord(), permutations = permute::how(nperm = 999), by = 'term')
    flog.info('get_anova_term() end.')
    return(res)
  })

  get_anova_model <- eventReactive(input$launch_beta, {
    flog.info('get_anova_model() starting...')
    res <- anova(ord(), permutations = permute::how(nperm = 999))
    flog.info('get_anova_model() end.')
    return(res)
  })

  output$anova_res <- renderPrint({
    get_anova_model()
  })

  get_anova_axis <- eventReactive(input$launch_beta, {
    flog.info('get_anova_axis() starting...')
    res <- anova(ord(), permutations = permute::how(nperm = 999), by = "axis")
    flog.info('get_anova_axis() end.')
    return(res)
  })

  output$anova_term_res <- renderPrint({
    get_anova_term()
  })

  output$anova_contrasts_res <- renderPrint({
    get_anova_contrast()
  })

  output$anova_axis_res <- renderPrint({
    get_anova_axis()
  })

  # ---- Dynamic UI: PERMANOVA box ----
  output$ui_permanova_box <- renderUI({
    if(input$ordination %in% c('NMDS', 'PCOA', 'dbRDA')){
      tagList(
        card(
          card_header("PERMANOVA Analysis"),
          htmltools::p(paste0('Permanova is done on the dissimilarity matrix computed with the selected index.', ' (here ', input$metrics, ' is used)')),
          uiOutput(ns("ui_adonis_factor")),
          actionButton(ns("update_test_btn"), "Update Test", class = "btn-primary"),
          h3('ADONIS formula:'),
          verbatimTextOutput(ns("adonis_formula"))
        ),
        navset_card_underline(
          title = "PERMANOVA Results",
          full_screen = TRUE,
          nav_panel("Adonis Test Result",
            DT::dataTableOutput(ns('adonistest'))
          ),
          nav_panel("Pairwise Adonis Test",
            DT::dataTableOutput(ns('adonispairwisetest'))
          )
        )
      )
    }
  })

  output$ui_adonis_factor = renderUI({
    req(get_meta_col(), r$sdat())
    facts = r$var_list()
    Fchoices = facts[facts != get_meta_col()]

    shinyWidgets::pickerInput(inputId = ns('adonis_factor'),
                              label = 'Select factor(s) to add as covariable: ',
                              choices = Fchoices,
                              multiple = TRUE
    )
  })

  # ---- Distance matrix computation ----
  physeq_dist <- eventReactive(input$launch_beta, {
    req(input$metrics, local_physeq())
    flog.info('physeq_dist() starting...')
    flog.debug('physeq_dist(): metrics=%s, nsamples=%d',
               input$metrics, phyloseq::nsamples(local_physeq()))
      validate(
        need(!input$metrics %in% c("unifrac", "wunifrac") & !is.null(local_physeq()@phy_tree),
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
    flog.debug('ord(): ordination=%s, metrics=%s', input$ordination,
               if (input$ordi_type == 'Distance-based') input$metrics else 'N/A')
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
    if(input$ordination %in% c('PCA', 'RDA', 'PCOA', 'dbRDA')){
      nmds_coord <- vegan::scores(ord(), display='sites', correlation = TRUE) %>%
        as_tibble(rownames="sample.id")
    } else if(input$ordination %in% c('CCA')){
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
    flog.debug('get_species_coord(): ordination=%s, rank_glom=%s, rank_color=%s',
               input$ordination, r$rank_glom(), input$rank_color)
    if(input$ordination %in% c('PCA', 'RDA', 'PCOA', 'dbRDA')){
      if(input$ordination %in% c('dbRDA') && input$metrics %in% c('unifrac', 'wunifrac')){
        nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='specie', scalling = 2) %>% as_tibble(rownames=r$rank_glom())
      } else {
        nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='specie', correlation = TRUE) %>% as_tibble(rownames=r$rank_glom())
      }
    } else if(input$ordination %in% c('CCA')){
      nmds_coord <- vegan::scores(ord(), choices=c(1,2), display='specie', hill = TRUE) %>% as_tibble(rownames=r$rank_glom())
    } else if(input$ordination %in% c('NMDS') && input$metrics %in% c('unifrac', 'wunifrac')){
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
  req(ord())
  flog.info('get_axis_names() starting...')
  axes <- colnames(vegan::scores(ord(), display = 'sites'))
  
  # Fallback defaults by ordination type
  if (length(axes) == 0) {
    fallback <- switch(input$ordination,
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
    
    # Default axes if UI not rendered yet
    axe_x <- input$axe_x %||% get_axis_names()[1]
    axe_y <- input$axe_y %||% get_axis_names()[2]
    
    # Default plot_type
    plot_types <- if (is.null(input$plot_type)) "samples" else input$plot_type
    
    flog.debug('base_plot(): plot_type=[%s], axe_x=%s, axe_y=%s',
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
          ) +
          scale_color_manual(
            values = r$factor_colors()[[input$beta_factor]], 
            drop = FALSE
          ) +
          ggnewscale::new_scale_color()
        if(!isNumFactor()){
          p <- p + 
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

      if(input$ordination == "PCOA"){
        eig1 <- eigenvals(ord())
        percent1 <- eig1/sum(eig1)*100
        p <- p + xlab(glue::glue("{input$axe_x} ({round(percent1[input$axe_x], 2)} %)")) +
          ylab(glue::glue("{input$axe_y} ({round(percent1[input$axe_y], 2)} %)"))
      }

      if ('taxa' %in% input$plot_type){
        flog.info('base_plot() plotting taxa...')
        species_coord <- get_species_coord()
        if(r$rank_glom() == 'ASV'){
          taxa <- tax_table(local_physeq()) %>% as.data.frame() %>% as_tibble(rownames="ASV") %>% select(ASV) %>% pull
        } else{
          taxa <- tax_table(local_physeq()) %>% as.data.frame() %>% as_tibble(rownames=r$rank_glom()) %>% select(r$rank_glom()) %>% pull
        }
        p <- p +
          geom_point(data = species_coord, aes(x=!!sym(axe_x), y=!!sym(axe_y), color=.data[[input$rank_color]], taxa = taxa)) + scale_color_manual(values=r$taxa_colors()[[input$rank_color]], drop = FALSE)
      }

      if ('env' %in% input$plot_type){
        flog.info('base_plot() plotting env...')
        if(input$ordination %in% c('RDA', 'CCA', 'dbRDA')){
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
        } else if(input$ordination %in% c('PCOA', 'NMDS')){
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
  observe({
    base_plot()
    output$beta_ggplot <- renderPlot({
      base_plot()
    })
  })

  # ---- PERMANOVA formula ----
  get_formula <- reactive({
    req(input$metrics, get_meta_col())
    form <- glue::glue('dist ~ Depth + ')
    if(!is.null(input$adonis_factor)){
      cov1 = paste(input$adonis_factor, collapse = " + ")
      form <- paste(form, glue::glue('{cov1} + {get_meta_col()}'), sep='')
    } else{
      form <- paste(form, glue::glue('{get_meta_col()}'), sep='')
    }
    flog.debug('get_formula(): %s', form)
    return(form)
  })

  # ---- Dispersion tests ----
  get_dispersion_res <- eventReactive(input$launch_beta,{
    req(physeq_dist(), get_meta_col())
    flog.info('get_dispersion_res() starting...')
    res <- vegan::betadisper(physeq_dist(), local_metadata()[,get_meta_col()])
    flog.info('get_dispersion_res() end.')
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
  get_adonis_res <- eventReactive(input$launch_beta | input$update_test_btn, {
    req(physeq_dist(), get_formula(), ord())
    flog.info('get_adonis_res() starting...')
    flog.debug('get_adonis_res(): formula=%s', get_formula())
    dist <- physeq_dist()
    mdata <- local_metadata()
    mdata$Depth <- sample_sums(r$phyloseq_filtered())
    mdata <- mdata %>% filter(!is.na(get_meta_col()))
    flog.debug('get_adonis_res(): nsamples=%d after NA filtering', nrow(mdata))
    res <- vegan::adonis2(as.formula(get_formula()), data = mdata, permutations = 1000)
    flog.info('get_adonis_res() end.')
    return(data.frame(res))
  })

  get_pairwise_res <- eventReactive(input$launch_beta | input$update_test_btn, {
    req(physeq_dist(), get_meta_col(), local_metadata())
    flog.info('get_pairwise_res() starting...')
    flog.debug('get_pairwise_res(): factor=%s, nsamples=%d',
               get_meta_col(), nrow(local_metadata()))
    res <- pairwise.adonis(physeq_dist(), local_metadata()[,get_meta_col()], p.adjust.m = "fdr")
    flog.info('get_pairwise_res() end.')
    return(res)
  })

  # ---- PERMANOVA outputs ----
  output$adonis_formula <- renderText({
    get_formula()
  })

  output$adonistest <- DT::renderDataTable({
    round_df <- function(x, digits = 4) {
      if (is.data.frame(x)) {
        x[] <- lapply(x, function(col) if (is.numeric(col)) round(col, digits) else col)
      }
      x
    }
    round_df(get_adonis_res(), 4)
  })

  output$adonispairwisetest <- DT::renderDataTable({
    round_df <- function(x, digits = 4) {
      if (is.data.frame(x)) {
        x[] <- lapply(x, function(col) if (is.numeric(col)) round(col, digits) else col)
      }
      x
    }
    round_df(get_pairwise_res(), 4)
  })

  # ---- Dispersion outputs ----
  dfdisper <- eventReactive(input$launch_beta | input$update_test_btn,{
    req(get_dispersion_res())
    flog.info('dfdisper() starting...')

    df1 = cbind.data.frame(distances = get_dispersion_res()$distances, group = get_dispersion_res()$group)

    if(input$order1){
      flog.debug('dfdisper(): ordering factor levels with mixedsort')
      df1$group = factor( df1$group, levels = gtools::mixedsort(levels(df1$group)) )
    }
    flog.info('dfdisper() end.')
    return(df1)
  })

  output$dispersionPlot <- renderPlotly({
   df1 <- dfdisper()
   plot_ly(df1, x = ~group, y = ~distances,
           color = ~group, type = 'box', colors = r$factor_colors()[[input$beta_factor]]) %>%
     layout(title="", yaxis = list(title = "Distance to centroid"), xaxis = list(title = 'Group'), barmode = 'stack') %>%
    config(toImageButtonOptions = list(format = "svg"))
 })

  output$dispersionTable <- DT::renderDataTable({
    round_df <- function(x, digits = 4) {
      if (is.data.frame(x)) {
        x[] <- lapply(x, function(col) if (is.numeric(col)) round(col, digits) else col)
      }
      x
    }
    round_df(as.data.frame(get_dispersion_anova()), 4)
  })

  output$dispersionTukey <- DT::renderDataTable({
    round_df <- function(x, digits = 4) {
      if (is.data.frame(x)) {
        x[] <- lapply(x, function(col) if (is.numeric(col)) round(col, digits) else col)
      }
      x
    }
    round_df(as.data.frame(get_dispersion_tukey()$group), 4)
  })
  })
}

## To be copied in the UI
# mod_beta_ui("beta_ui_1")

## To be copied in the server
# mod_beta_server("beta_ui_1", r = r)
