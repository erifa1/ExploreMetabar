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
                           selected = 'Distance-based')
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
                             list("samples", "taxa", "biplot", "triplot"),
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
              uiOutput(ns('ui_envfit_switch'))
              
            ),
            fluidRow(
              actionButton(ns("launch_beta"), "Run Beta Plot", icon = icon("play-circle"),
                           style="color: #fff; background-color: #3b9ef5; border-color: #1a4469")
            )
          ),  title = "Settings:", width = 6, status = "warning", solidHeader = TRUE
        ),
        uiOutput(ns('ui_constrain')),
        uiOutput(ns('ui_ordistep')),
        uiOutput(ns('envfit_box')),
        uiOutput(ns('envfit_box_res'))
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



veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#' mod_beta Server Function
#'
#' @import vegan
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
  
  
  
  output$ui_envfit_switch <- renderUI({
    if(input$ordination %in% c('NMDS', 'PCOA')){
      materialSwitch(
        ns('envfit_switch'),
        label = 'Try envfit',
        value = FALSE,
        status = 'primary'
      )
    }
  })
  
  
  output$ui_constrain <- renderUI({
    if(input$ordination %in% c('RDA','CCA')){
      box(title = 'Model parameters',
          radioButtons(inputId = ns('param_mode'),
                       label = 'méthode to select parameters',
                       choices = c('picker', 'ordiR2step', 'manual')),
          uiOutput(ns('constr_select')),
          verbatimTextOutput(
            ns('formula')
          )
      )
    }
  })
  
  output$constr_select <- renderUI({
    if(input$param_mode == 'picker'){
      shinyWidgets::pickerInput(inputId = ns('constr_picker'),
                                label = 'Select terms to create a formula',
                                choices = colnames(r$sdat()),
                                multiple = TRUE
      )
    } else if(input$param_mode == 'ordiR2step'){
      verbatimTextOutput(ns('constr_ordiR2step_out'))
    } else if(input$param_mode == 'manual'){
      textInput(ns('constr_manual_formula'), label = 'Enter your own formula:', placeholder = 'spe ~')
    }
  })
  
  get_ordiR2step <- reactive({
    env <- as(r$sdat(), 'data.frame')
    
    spe <- veganifyOTU(physeq())
    mod0 <- vegan::rda(spe ~ 1, data = env)
    mod1 <- vegan::rda(spe ~ ., data = env)
    sel <- vegan::ordiR2step(mod0, scope = formula(mod1), R2scope = FALSE)
    return(sel)
  })
  
  output$constr_ordiR2step_out <- renderPrint({
    sel <- get_ordiR2step()
    print(sel)
    print(sel$anova)
  })
  
  output$formula <- renderPrint({
    get_constr_formula()
  })
  
  
  get_constr_formula <- reactive({
    if(input$param_mode == 'picker'){
      f <- paste0('spe ~ ', paste(input$constr_picker, collapse = ' + '))
    } else if(input$param_mode == 'ordiR2step'){
      f <- formula(get_ordiR2step())
    } else if(input$param_mode == 'manual'){
      f <- input$constr_manual_formula
    }
    return(f)
  })
  
  output$rank_select <- renderUI({
    if(input$ordination %in% c('NMDS', 'PCOA')){
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
  
  
  output$envfit_box_res <- renderUI({
    if(input$ordination %in% c('NMDS', 'PCOA') && input$envfit_switch){
      box(
        verbatimTextOutput(
          ns('envfit_res')
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
  

  # observe({
  #   req(r$phyloseq_filtered())
  #   if(is.null(phy_tree(r$phyloseq_filtered(), errorIfNULL=FALSE))){
  #     flog.info("no phytree beta metrics update")
  #     ch1 = list("bray", "jaccard")
  #   }else{
  #     ch1 = list("bray", "jaccard", "unifrac", "wunifrac")
  #   }
  #   updateRadioButtons(session, "metrics",
  #                     choices = ch1, inline = TRUE)
  # })


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
      spe <- veganifyOTU(physeq())
      res <- vegan::capscale(spe ~ 1, spe, distance = input$metrics)
    } else if(input$ordination == 'RDA'){
      env <- as(r$sdat(), 'data.frame')
      spe <- veganifyOTU(physeq())
      res <- vegan::rda(as.formula(get_constr_formula()), data = env)
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
  
  
  output$envfit_res <- renderPrint({
    get_env_fit()
  })
  
  get_env_fit <- eventReactive( input$launch_beta, {
    env <- as(r$sdat(), 'data.frame')
    env <- env[, input$envfit_param, drop=F]
    en <- vegan::envfit(ord(), env)
    return(en)
  })
  
  get_axis_names <- reactive({
    if(input$ordination == 'NMDS'){
      axes <- c('NMDS1', 'NMDS2')
    } else if(input$ordination == 'PCOA'){
      axes <- c('MDS1', 'MDS2')
    } else if(input$ordination == 'RDA'){
      axes <- colnames(vegan::scores(ord(), display = 'sites'))
    }
    
  })
  
  base_plot <- reactive({

    p <- ggplot2::ggplot()
    axes <- get_axis_names()
    if(input$plot_type %in% c('samples', 'biplot', 'triplot')){
      nmds_coord <- get_sites_nmds_coord()
      # browser()
      
      p <- p +
            geom_point(data = nmds_coord, mapping = aes(x=!!sym(axes[1]), y=!!sym(axes[2]), color=.data[[get_meta_col()]], text = paste('sample.id:',phyloseq::sample_names(physeq()))), shape=23) +
            stat_ellipse(data = nmds_coord, mapping = aes(x=!!sym(axes[1]), y=!!sym(axes[2]), group = !!sym(get_meta_col()), color = !!sym(get_meta_col())))
    } 
    
    if (input$plot_type %in% c('taxa', 'biplot', 'triplot')){
      nmds_coord <- get_species_nmds_coord()
      # browser()
      if(r$rank_glom() == 'ASV'){
        taxa <- tax_table(physeq()) %>% as.data.frame() %>% as_tibble(rownames="ASV") %>% select(ASV) %>% pull
      } else{
        taxa <- tax_table(physeq()) %>% as.data.frame() %>% as_tibble(rownames=r$rank_glom()) %>% select(r$rank_glom()) %>% pull
      }
      p <- p +
            geom_point(data = nmds_coord, aes(x=!!sym(axes[1]), y=!!sym(axes[2]), color=.data[[input$rank_color]], taxa = taxa))
    } 
    
    if (input$plot_type == 'triplot'){
      if(input$ordination == 'RDA'){
        scr <- vegan::scores(ord())
        if(!is.null(scr$biplot)){
          en_coord_cont <- as.data.frame(scr$biplot)
          p <- p + geom_segment(aes(x = 0, y = 0, xend = !!sym(axes[1]), yend = !!sym(axes[2])),
                                data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = grid::arrow()) +
            geom_text(data = en_coord_cont, aes(x = !!sym(axes[1]), y = !!sym(axes[2])), colour = "grey30",
                      fontface = "bold", label = row.names(en_coord_cont))
        }
        if(!is.null(scr$centroids)){
          en_coord_cat <- as.data.frame(scr$centroids)
          p <- p + geom_point(data = en_coord_cat, aes(x = !!sym(axes[1]), y = !!sym(axes[2])),
                              shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
            geom_text(data = en_coord_cat, aes(x = !!sym(axes[1]), y = !!sym(axes[2])),
                      label = row.names(en_coord_cat), colour = "navy", fontface = "bold")
        }
      }
      
      if(input$envfit_switch){
        
        en <- get_env_fit()
        
        if(!is.null(en$vectors)){
          en_coord_cont <- as.data.frame(vegan::scores(en, "vectors")) * vegan::ordiArrowMul(en)
          p <- p + geom_segment(aes(x = 0, y = 0, xend = !!sym(axes[1]), yend = !!sym(axes[2])),
                                data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30", arrow = grid::arrow()) +
            geom_text(data = en_coord_cont, aes(x = !!sym(axes[1]), y = !!sym(axes[2])), colour = "grey30",
                      fontface = "bold", label = row.names(en_coord_cont))
        }
        if(!is.null(en$factors)){
          en_coord_cat <- as.data.frame(vegan::scores(en, "factors"))
          p <- p + geom_point(data = en_coord_cat, aes(x = !!sym(axes[1]), y = !!sym(axes[2])),
                              shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
            geom_text(data = en_coord_cat, aes(x = !!sym(axes[1]), y = !!sym(axes[2])),
                      label = row.names(en_coord_cat), colour = "navy", fontface = "bold")
        }
        
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
  
  
  get_dispersion_res <- eventReactive(input$launch_beta,{
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
  
  
  get_adonis_res <- eventReactive(input$launch_beta, {
    req(physeq_dist(), get_formula())
    dist <- physeq_dist()
    mdata <- local_metadata()
    mdata$Depth <- sample_sums(physeq())
    # Filter NA value in metadata
    mdata <- mdata %>% filter(!is.na(get_meta_col()))
    res <- vegan::adonis2(as.formula(get_formula()), data = mdata, permutations = 1000)
    return(data.frame(res))
  })
  
  
  get_pairwise_res <- eventReactive(input$launch_beta, {
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


  dfdisper <- eventReactive(input$launch_beta,{
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
