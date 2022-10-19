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

      infoBox("",
              "Use phyloseq object without taxa merging step.",
              icon = icon("info-circle"), fill=TRUE, width = 10
      ),

      box(
        radioButtons(
          ns("beta_norm_bool"),
          label = "Use normalized data (prefer TSS normalization)",
          inline = TRUE,
          choices = list(
            "Raw" = 0 ,
            "Normalized" = 1
          ), selected = 1
        ),

        radioButtons(ns("metrics"), "Choose one index:", inline = TRUE,
                     choices ='',
                     selected = c("bray")
        ),

        radioButtons(ns("ordination"), "Choose one ordination:", inline = TRUE,
                     choices =
                       list("MDS", "NMDS", "CCA", "RDA"),
                     selected = c("NMDS")
        ),
        selectInput(
          ns("beta_fact1"),
          label = "Select main factor to test + color plot: ",
          choices = ''
        ),
        uiOutput(ns('ui_beta_fact2')),
        actionButton(ns("launch_beta"), "Run Beta Plot", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(
        shinycustomloader::withLoader(
          plotly::plotlyOutput(ns("plot1")),
          type = "html", loader = "loader4"
        ),
        title = "Ordination plot:", width = 12, status = "primary", solidHeader = TRUE
      ),
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

  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "beta_fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })
  
  isNumFactor <- reactive({
    req(input$beta_fact1, r$sdat())
    metadata <- as(r$sdat(), "data.frame")
    if(is.numeric(metadata[, input$beta_fact1])){
      return(TRUE)
    } else{
      return(FALSE)
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
  
  output$ui_beta_fact2 <- renderUI({
    req(r$phyloseq_filtered(), input$beta_fact1)
    if(! isNumFactor()){
      metadata <- as(r$sdat(), "data.frame")
      num_col_names <- metadata %>% dplyr::select_if(is.numeric) %>% colnames
      tmp <- dplyr::setdiff(colnames(metadata), num_col_names)
      selectInput(
        ns("beta_fact2"),
        label = "Select second factor to combine: ",
        choices = c('none',dplyr::setdiff(tmp, input$beta_fact1))
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
    phyloseq::distance(physeq(), method = input$metrics)
  })

  ord <- reactive({
    req(input$ordination)
    phyloseq::ordinate(physeq= physeq(), distance = physeq_dist(), method= input$ordination)
  })


  base_plot <- reactive({
    p <- phyloseq::plot_ordination(physeq = physeq(), ordination = ord(), axes = c(1, 2))
    p$layers[[1]] <- NULL

    xrange <- c()
    xrange[1] <- layer_scales(p)$x$range$range[1] - abs(layer_scales(p)$x$range$range[1])*3
    xrange[2] <- layer_scales(p)$x$range$range[2] + abs(layer_scales(p)$x$range$range[2])*3

    yrange <- c()
    yrange[1] <- layer_scales(p)$y$range$range[1] - abs(layer_scales(p)$y$range$range[1])*3
    yrange[2] <- layer_scales(p)$y$range$range[2] + abs(layer_scales(p)$y$range$range[2])*3
    return(list('plot'=p, 'xrange'=xrange, 'yrange'=yrange))
  })


  output$plot1 <- plotly::renderPlotly({
    beta_plot()
  })


  beta_plot <- eventReactive(input$launch_beta, {
    withProgress({
      sample.id = sample_names(physeq())
      p <- base_plot()$plot
      p <- p + aes(color = !!sym(get_meta_col()), sample.id = sample.id)
      p <- p + stat_ellipse(aes(group = !!sym(get_meta_col())))
      p <- p + xlim(base_plot()$xrange) + ylim(base_plot()$yrange)
      p <- p + geom_point() + theme_bw()
      ggplotly(p, tooltip=c("x", "y", "sample.id")) %>% config(toImageButtonOptions = list(format = "svg"))
    }, message = "Plot Beta...")
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
  
  
  get_paitwise_res <- reactive({
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
    get_paitwise_res()
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
