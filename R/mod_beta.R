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
#' @import PCAmixdata
#' @import shinycustomloader
#' @import shinyWidgets
#' 

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
            htmltools::p('In this module, asv table is normalized by hellinger (Legendre & Gallagher 2001) method and environmental variables are centered and scaled.'),
            fluidRow(
              radioButtons(ns('ordi_type'), 'Choose your ordination type:',
                           inline = T,
                           # choices = c('Unconstrained', 'Constrained', 'Distance-based'),
                           choices = c('Constrained', 'Distance-based'),
                           selected = 'Distance-based')
            ),
            fluidRow(
              uiOutput(ns('ui_metrics'))
            ),
            fluidRow(
              radioButtons(ns("ordination"), "Choose one ordination:", inline = TRUE,
                           choices = '',
                           selected = ''
              ),
              uiOutput(ns('ui_beta_factor'))
            ),
            fluidRow(
              actionButton(ns("launch_beta"), "Run Beta Plot", icon = icon("play-circle"),
                           style="color: #fff; background-color: #3b9ef5; border-color: #1a4469")
            )
          ),  title = "Settings:", width = 6, status = "warning", solidHeader = TRUE
        ),
        uiOutput(ns('ui_constrain')),
        uiOutput(ns('ui_ordistep')),
        uiOutput(ns('ui_envfit_box')),
        uiOutput(ns('ui_envfit_box_res'))
      ),
      fluidRow(
        box(
          shinyWidgets::prettyCheckboxGroup(ns("plot_type"), "Choose plot type:", inline = TRUE,
                                            choices =
                                              list("samples", "taxa", "env"),
                                            selected = c("samples")
          ),
          uiOutput(ns('ui_taxa_rank')),
          #shinyWidgets::materialSwitch(inputId = ns('ggplot_switch'), label = 'ggplot2 or plotly'),
          uiOutput(ns('ui_axe_x')),
          uiOutput(ns('ui_axe_y')),
          title = "Plot options", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE
        )
      ),
      fluidRow(
        box(
          shinycustomloader::withLoader(
            uiOutput(ns('ui_plot')),
            # plotly::plotlyOutput(ns("plot1"), height = "730px"),
            type = "html", loader = "loader4"
          ),
          title = "Ordination plot:", width = 12, height = "800px", status = "primary", solidHeader = TRUE
        ), style = "height:800px;"
      ),
      fluidRow(
        box(
          plotOutput(ns('screeplot')),
          title = 'Screeplot', width = 12, collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE)
      ),
      fluidRow(
        uiOutput(ns('ui_permanova_box')),
        uiOutput(ns('ui_anova_box'))
      )
    )
  )
}



veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#' mod_beta Server Function
#'
#' @import vegan
#' @importFrom plotly ggplotly
#' @importFrom DT renderDataTable
#' @importFrom permute how
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
mod_beta_server <- function(input, output, session, r = r){
  ns <- session$ns

  ### local metadata and phyloseq object functions
  
  isNumFactor <- reactive({
    req(get_meta_col(), local_metadata())
    metadata <- local_metadata()
    if(is.numeric(metadata[, get_meta_col()])){
      return(TRUE)
    } else{
      return(FALSE)
    }
  })
  
  get_meta_col <- reactive({
    req(input$beta_factor, r$sdat())
    metadata <- r$sdat()
    if(length(input$beta_factor) == 1){
      meta.col <- input$beta_factor
    } else if(length(input$beta_factor) > 1) {
      validate(
        need(!any(sapply(metadata[, input$beta_factor], is.numeric)), message = "You can't select multiple numeric factors")
      )
      meta.col <- paste0(input$beta_factor, collapse='_')
    }
    return(meta.col)
  })
  
  
  local_metadata <- reactive({
    req(input$beta_factor, r$sdat())
    metadata <- r$sdat()
    if(! all(sapply(metadata[, input$beta_factor], is.numeric))){
      metadata <- tidyr::unite(metadata, !!get_meta_col(), input$beta_factor, na.rm=TRUE)
      metadata[, get_meta_col()] <- as.factor(metadata[, get_meta_col()])
    }
    return(metadata)
  })
  
  ### dynamic UI rendering functions
  
  # distance based radioButton
  output$ui_metrics <- renderUI({
    req(input$ordi_type)
    if(input$ordi_type == 'Distance-based'){
      radioButtons(ns("metrics"), "Choose one distance metric:", inline = TRUE,
                   choices =c('bray', 'jaccard', 'unifrac', 'wunifrac'),
                   selected = c("bray")
      )
    }
  })
  
  # constrain based radioButton
  output$ui_constrain <- renderUI({
    req(input$ordination)
    if(input$ordination %in% c('RDA','CCA', 'dbRDA')){
      box(title = 'Model parameters', width=6, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
          htmltools::p('Warning: when integrating environmental variable, samples with missing values are omitted.'),
          radioButtons(inputId = ns('param_mode'),
                       label = 'methode to select parameters',
                       choices = c('picker', 'ordiR2step', 'manual')),
          uiOutput(ns('constr_select')),
          uiOutput(ns('ordiR2_out')),
          verbatimTextOutput(
            ns('formula')
          ),
          uiOutput(ns('ui_ordiR2_btn'))
      )
    }
  })
  
  output$ui_ordiR2_btn <- renderUI({
    if(input$param_mode == 'ordiR2step'){
      shinyWidgets::actionBttn(inputId = ns('ordiR2_btn'), label = 'launch', size = 'sm')
    }
  })
  
  
  # ordiR2step verbatim output
  output$ordiR2_out <- renderUI({
    req(input$param_mode)
    if(input$param_mode == 'ordiR2step'){
      verbatimTextOutput(ns('constr_ordiR2step_out'))
    }
  })

  
  ## constrained based box for model parameters
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
                                  content = unlist(lapply(
                                    X = colnames(local_metadata()),
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
                                                        class(local_metadata()[,x])
                                                      ),
                                                      tags$div(
                                                        style = htmltools::css(color = 'grey'),
                                                        paste0(sum(is.na(local_metadata()[,x])), '/', nrow(local_metadata()), ' NAs')
                                                      )
                                          )
                                        )
                                      )
                                    }
                                  ))
                                )
      )
    } else if(input$param_mode == 'manual'){
      textInput(ns('constr_manual_formula'), label = 'Enter your own formula:', placeholder = 'spe ~')
    }
  })
  
  
  ## ordiR2step fonction
  get_ordiR2step <- eventReactive(input$ordiR2_btn, {
    req(local_metadata(), local_physeq())
    # flog.info('get_ordiR2step() starting...')
    
    spe <- veganifyOTU(local_physeq())
    spe <- vegan::decostand(spe, method = 'hell')
    
    env <- local_metadata()[, c(input$constr_picker, 'sample.id')]
    nsample_before <- nrow(env)
    env <- na.omit(env)
    nsample_after <- nrow(env)
    
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
        fun <- glue::glue('mod0 <- vegan::capscale(spe ~ 1, data = env, na.action = "na.omit", distance = "{input$metrics}")')
        eval(parse(text=fun))
        fun <- glue::glue('mod1 <- vegan::capscale(spe ~ ., data = env, na.action = "na.omit", distance = "{input$metrics}")')
        eval(parse(text=fun))
      } 
      # else if(input$metrics %in% c('unifrac', 'wunifrac')){
      #   phy <- phyloseq::prune_samples(rownames(spe), local_physeq())
      #   phy <- phyloseq::prune_taxa(phyloseq::taxa_sums(phy)>0, phy)
      #   dist <- phyloseq::distance(physeq = phy, method = input$metrics)
      #   browser()
      #   fun <- glue::glue('mod0 <- vegan::capscale(dist ~ 1, data = env, na.action = "na.omit", comm = spe)')
      #   eval(parse(text=fun))
      #   fun <- glue::glue('mod1 <- vegan::capscale(dist ~ ., data = env, na.action = "na.omit", comm = spe)')
      #   eval(parse(text=fun))
      # }
    }
    res <- list()
    
    res$sel <- vegan::ordiR2step(mod0, scope = formula(mod1), R2scope = FALSE, trace = FALSE)
    res$lostsamples <- nsample_before - nsample_after
    # flog.info('get_ordiR2step() end.')
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
    return(f)
  })
  
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
      box(
        title = "Pairwise Adonis Test", width = 12, status = "primary", solidHeader = TRUE,
        DT::dataTableOutput(ns("adonispairwisetest"))
      )
    }
  })
  
  get_screeplot <- reactive({
    req(ord())
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
  
  output$disper_res <- renderUI({
    req(get_meta_col())
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
  
  
  output$ui_beta_factor <- renderUI({
    req(r$var_list())
    shinyWidgets::pickerInput(
       ns("beta_factor"),
       label = "Select factor(s) to color samples and ellipses (multiple selection allowed for qualitative variables):",
       choices = r$var_list(),
       selected = r$var_list()[2],
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
                  )
                )
              )
            )
           }
         ))
       )
     )
  })

  
  output$ui_envfit_box <- renderUI({
    req(input$ordination, local_metadata())
    if(input$ordination %in% c('NMDS', 'PCOA')){
      box(
        htmltools::p('The envfit function fits environmental vectors or factors onto an ordination.'),
        shinyWidgets::pickerInput(
          ns("envfit_param"),
          label = "Select one or more factor to use in envfit module:",
          choices = sort(colnames(local_metadata())),
          # selected = colnames(local_metadata())[2],
          multiple = TRUE,
          options = pickerOptions(
            actionsBox = TRUE,
            liveSearch = TRUE,
            showContent = FALSE
          ),
          choicesOpt = list(
            content = unlist(lapply(
              X = sort(colnames(local_metadata())),
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
                                  class(local_metadata()[,x])
                                )
                    )
                  )
                )
              }
            ))
          )
        ),
        numericInput(
          ns('envfit_pval'), 
          'p-value threshold to be plotted', value = 1,
          min = 0,
          max = 1,
          step = 0.01
        ),
        title = 'VEGAN envfit', width=6, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE
      )
    }
  })
  
  
  output$ui_envfit_box_res <- renderUI({
    req(input$ordination)
    if(input$ordination %in% c('NMDS', 'PCOA')){
      box(width=6, status = "primary", solidHeader = TRUE, title = 'envfit results', collapsible = TRUE, collapsed = TRUE,
        verbatimTextOutput(
          ns('envfit_res')
        ),
      )
    }
  })
  

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

  
  
  output$ui_anova_box <- renderUI({
    if(input$ordination %in% c('RDA', 'CCA')){
      box(
        title = "Anova on RDA/CCA results", width = 12, status = "primary", solidHeader = TRUE,
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
  
  
  
  
  get_anova_contrast <- eventReactive(input$launch_beta, {
    res <- anova(ord(), permutations = permute::how(nperm = 999), by = 'onedf')
    return(res)
  })
  
  
  get_anova_term <- eventReactive(input$launch_beta, {
    res <- anova(ord(), permutations = permute::how(nperm = 999), by = 'term')
    return(res)
  })
  
  
  get_anova_model <- eventReactive(input$launch_beta, {
    res <- anova(ord(), permutations = permute::how(nperm = 999))
    return(res)
  })
  
  
  output$anova_res <- renderPrint({
    get_anova_model()
  })
  
  
  get_anova_axis <- eventReactive(input$launch_beta, {
    res <- anova(ord(), permutations = permute::how(nperm = 999), by = "axis")
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
  
  
  output$ui_permanova_box <- renderUI({
    if(input$ordination %in% c('NMDS', 'PCOA', 'dbRDA')){
      box(
        title = "Permanova with adonis:", width = 12, status = "primary", solidHeader = TRUE,
        htmltools::p(paste0('Permanova is done on the dissimilarity matrix computed with the selected index.', ' (here ', input$metrics, ' is used)')),
        uiOutput(ns("ui_adonis_factor")),
        actionButton(ns("update_test_btn"), "Update Test", style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        h3('ADONIS formula:'),
        verbatimTextOutput(ns("adonis_formula")),
        h2('Permanova Adonis Test Result: '),
        DT::dataTableOutput(ns('adonistest')),
        uiOutput(ns('pairwise_res')),
        uiOutput(ns('disper_res'))
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


  local_physeq <- reactive({
    req(r$phyloseq_filtered(), local_metadata())
    data <- r$phyloseq_filtered()
    sample_data(data) <- sample_data(local_metadata())
    return(data)
  })

  
  physeq_dist <- eventReactive(input$launch_beta, {
    req(input$metrics, local_physeq())
    flog.info('physeq_dist() starting...')
      validate(
        need(!input$metrics %in% c("unifrac", "wunifrac") & !is.null(local_physeq()@phy_tree), 
          message = 'Unifrac and wunifrac available if phylogenetic tree available.')
      )

    res <- phyloseq::distance(local_physeq(), method = input$metrics)
    flog.info('physeq_dist() end.')
    return(res)
  })
  
  
  get_env_scaled <- reactive({
    req(r$sdat())
    # flog.info('get_env_scaled() starting...')
    env <- local_metadata()
    data.split <- PCAmixdata::splitmix(env)
    env[,data.split$col.quant] <- scale(env[,data.split$col.quant], center = T, scale = T)
    # flog.info('get_env_scaled() end.')
    return(env)
  })
  
  get_species_table <- reactive({
    req(local_physeq())
    spe <- veganifyOTU(local_physeq())
    spe <- vegan::decostand(spe, method = 'hell')
  })

  
  ord <- eventReactive(input$launch_beta, {
    req(input$ordination, local_physeq())
    
    flog.info('ord() starting...')
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
          res <- vegan::capscale(spe ~ 1, spe, distance = input$metrics, )
        } else if(input$metrics %in% c('unifrac', 'wunifrac')){
          dist <- physeq_dist()
          res <- vegan::capscale(dist ~ 1, comm = spe)
        }
      } else if(input$ordination == 'RDA'){
        env <- get_env_scaled()
        spe <- get_species_table()
        res <- vegan::rda(as.formula(get_constr_formula()), data = env, na.action = 'na.omit')
      } else if(input$ordination == 'CCA'){
        env <- get_env_scaled()
        spe <- get_species_table()
        res <- vegan::cca(as.formula(get_constr_formula()), data = env, na.action = 'na.omit')
      } else if(input$ordination == 'dbRDA'){
        env <- get_env_scaled()
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
    
    nmds_coord <- nmds_coord %>% 
                    inner_join(., local_metadata() %>% 
                                 select(sample.id, !!get_meta_col()), by="sample.id")
    flog.info('get_sites_coord() end.')
    return(nmds_coord)
  })
  
  
  get_species_coord <- reactive({
    req(local_physeq(), ord(), r$rank_glom(), input$rank_color)
    flog.info('get_species_coord() starting...')
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
    flog.info('get_species_coord() end.')
    return(nmds_coord)
  })
  
  
  output$envfit_res <- renderPrint({
    get_env_fit()
  })
  
  get_env_fit <- eventReactive( input$launch_beta, {
    # flog.info(msg = 'get_env_fit() starting...')
    req(ord())
    env <- get_env_scaled()
    env <- env[, input$envfit_param, drop=F]
    en <- vegan::envfit(ord(), env, na.rm = T)
    # flog.info(msg = 'get_env_fit() end.')
    return(en)
  })
  
  
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
    flog.info('get_axis_names() end.')
    return(axes)
  })
  
  # 
  # observe({
  #   req(input$metrics)
  #   if(input$metrics %in% c('unifrac', 'wunifrac') && input$ordination %in% c('NMDS')){
  #     ch <- c('samples', 'env')
  #   } else{
  #     ch <- c('samples', 'taxa', 'env')
  #   }
  #   shinyWidgets::updatePrettyCheckboxGroup(inputId = 'plot_type', choices = ch, selected = 'samples', inline = TRUE)
  # })
  # 

  base_plot <- reactive({
    req(ord(), get_axis_names(), input$axe_x, input$axe_y)
    flog.info('base_plot() starting...')
    withProgress(message = 'Plotting...', min=0, max=10, value = 0,{
      p <- ggplot2::ggplot()
      if('samples' %in% input$plot_type){
        flog.info('base_plot() plotting samples...')
        setProgress(value = 4, detail = 'plotting samples')
        sites_coord <- get_sites_coord()
        p <- p +
          geom_point(data = sites_coord, mapping = aes(x=!!sym(input$axe_x), y=!!sym(input$axe_y), fill=.data[[get_meta_col()]], text = paste('sample.id:',sites_coord$sample.id)), shape=23, size = 5) + scale_fill_manual(values = r$factor_colors()[[input$beta_factor]]) + ggnewscale::new_scale_fill()
        if(!isNumFactor()){
          p <- p + stat_ellipse(data = sites_coord, mapping = aes(x=!!sym(input$axe_x), y=!!sym(input$axe_y), group = !!sym(get_meta_col()), color = !!sym(get_meta_col()))) + scale_color_manual(values = r$factor_colors()[[input$beta_factor]]) + ggnewscale::new_scale_color()
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
          geom_point(data = species_coord, aes(x=!!sym(input$axe_x), y=!!sym(input$axe_y), color=.data[[input$rank_color]], taxa = taxa)) + scale_color_manual(values=r$taxa_colors()[[input$rank_color]])
      } 
      
      if ('env' %in% input$plot_type){
        if(input$ordination %in% c('RDA', 'CCA', 'dbRDA')){
          scr <- vegan::scores(ord())
          if(!is.null(scr$biplot)){
            en_coord_cont <- as.data.frame(scr$biplot)
            en_coord_cont <- en_coord_cont[rownames(en_coord_cont) %in% colnames(local_metadata()),]
            if(nrow(en_coord_cont) > 0){
              p <- p + geom_segment(aes(x = 0, y = 0, xend = !!sym(input$axe_x), yend = !!sym(input$axe_y)),
                                    data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = grid::arrow()) +
                geom_text(data = en_coord_cont, aes(x = !!sym(input$axe_x), y = !!sym(input$axe_y)), colour = "grey30",
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
              geom_text(data = en_coord_cont, aes(x = !!sym(input$axe_x), y = !!sym(input$axe_y)), colour = "grey30",
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
    flog.info('base_plot() end.')
    return(p)
  })

  output$ui_plot <- renderUI({
    # if(input$ggplot_switch){
    #   flog.info('plotly ui...')
    #   plotly::plotlyOutput(ns("beta_plotly"), height = "730px")
    # } else{
      plotOutput(ns('beta_ggplot'), height = "730px")
    # }
  })

  
  observe({
    req(base_plot())
    # if(input$ggplot_switch){
    #   flog.info('plotly render...')
    #   output$beta_plotly <- plotly::renderPlotly({
    #     base_plot()
    #   })
    # } else {
      output$beta_ggplot <- renderPlot({
        base_plot()
      })
    # }
  })
  

# 
#   beta_plot <- reactive({
#     withProgress({
#       p <- base_plot()$plot
#       # p <- get_plot()
#       # p <- p + xlim(base_plot()$xrange) + ylim(base_plot()$yrange)
#       # p <- p + geom_point() + theme_bw()
#       # browser()
#     #   if(input$plot_type == 'samples'){
#     #     p <- ggplotly(p, tooltip=c("x", "y", "sample.id"))
#     #   } else if(input$plot_type == 'taxa'){
#     #     p <- ggplotly(p, tooltip=c("x", "y", "taxa"))
#     #   }
#     #   else if(input$plot_type == 'biplot'){
#     #     p <- ggplotly(p, tooltip=c("x", "y", "taxa", "sample.id"))
#     #   }else{
#     #     
#     #   }
#     #   p <- p %>% config(toImageButtonOptions = list(format = "svg"))
#     #   
#     }, message = "Plot Beta...")
#     return(p)
#   })

  
  get_formula <- reactive({
    req(input$metrics, get_meta_col())
    form <- glue::glue('dist ~ Depth + ')
    if(!is.null(input$adonis_factor)){
      cov1 = paste(input$adonis_factor, collapse = " + ")
      form <- paste(form, glue::glue('{cov1} + {get_meta_col()}'), sep='')
    } else{
      form <- paste(form, glue::glue('{get_meta_col()}'), sep='')
    }
    return(form)
  })
  
  
  get_dispersion_res <- eventReactive(input$launch_beta,{
    req(physeq_dist(), get_meta_col())
    flog.info(msg = 'get_dispersion_res() starting...')
    res <- vegan::betadisper(physeq_dist(), local_metadata()[,get_meta_col()])
    flog.info(msg = 'get_dispersion_res() end.')
    return(res)
  })
  
  
  get_dispersion_anova <- reactive({
    req(get_dispersion_res())
    flog.info(msg = 'get_dispersion_anova() starting...')
    res <- anova(get_dispersion_res())
    flog.info(msg = 'get_dispersion_anova() end.')
    return(res)
  })
  
  
  get_dispersion_tukey <- reactive({
    req(get_dispersion_res())
    flog.info(msg = 'get_dispersion_tukey() starting...')
    res <- TukeyHSD(get_dispersion_res())
    flog.info(msg = 'get_dispersion_tukey() end.')
    return(res)
  })
  
  
  get_adonis_res <- eventReactive(input$launch_beta | input$update_test_btn, {
    req(physeq_dist(), get_formula(), ord())
    flog.info(msg = 'get_adonis_res() starting...')
    dist <- physeq_dist()
    mdata <- local_metadata()
    mdata$Depth <- sample_sums(r$phyloseq_filtered())
    # Filter NA value in metadata
    mdata <- mdata %>% filter(!is.na(get_meta_col()))
    res <- vegan::adonis2(as.formula(get_formula()), data = mdata, permutations = 1000)
    flog.info(msg = 'get_adonis_res() end.')
    return(data.frame(res))
  })
  
  
  get_pairwise_res <- eventReactive(input$launch_beta | input$update_test_btn, {
    req(physeq_dist(), get_meta_col(), local_metadata())
    flog.info(msg = 'get_pairwise_res() starting...')
    res <- pairwise.adonis(physeq_dist(), local_metadata()[,get_meta_col()], p.adjust.m = "fdr")
    flog.info(msg = 'get_pairwise_res() end.')
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


  dfdisper <- eventReactive(input$launch_beta | input$update_test_btn,{
    req(get_dispersion_res())
    flog.info('dfdisper() starting...')
    
    df1 = cbind.data.frame(distances = get_dispersion_res()$distances, group = get_dispersion_res()$group)

    if(input$order1){
      flog.info("ORDER factor")
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
