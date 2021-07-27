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
#' @importFrom pairwiseAdonis pairwise.adonis
#' @importFrom DT dataTableOutput

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
        h2('Pairwise Adonis Test Results: '),
        DT::dataTableOutput(ns("adonispairwisetest"))
      ),
      box(
        title = "Dispersion results:", width = 12, status = "primary", solidHeader = TRUE,
        h3('Boxplots distance to centroid for each group:'),
        plotlyOutput(ns("dispersionPlot")),
        h3('Anova on dispersion:'),
        DT::dataTableOutput(ns("dispersionTable")),
        h3('TukeyHSD test on dispersion'),
        DT::dataTableOutput(ns("dispersionTukey"))
      )
    )
  )
}

#' mod_beta Server Function
#'
#' @importFrom vegan vegdist
#' @importFrom vegan adonis
#' @importFrom plotly ggplotly
#' @importFrom DT renderDataTable
#'
#' @noRd
mod_beta_server <- function(input, output, session, r = r){
  ns <- session$ns

  observeEvent(r$tabs$tabselected, {
    if(r$tabs$tabselected=='tab_beta' && !isTruthy(r$phyloseq_filtered())){
      shinyalert(title = "Oops", text="Phyloseq object not present. Return to input data and validate all steps.", type='error')
    }
  })


  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "beta_fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })

  observe({
    req(r$phyloseq_filtered())
    if(is.null(phy_tree(r$phyloseq_filtered(), errorIfNULL=FALSE))){
      print("no phytree beta metrics update")
      ch1 = list("bray", "jaccard")
    }else{
      ch1 = list("bray", "jaccard", "unifrac", "wunifrac")
    }
    updateRadioButtons(session, "metrics",
                      choices = ch1, inline = TRUE)
  })


  output$factor2 = renderUI({
    req(input$beta_fact1, r$phyloseq_filtered())
    facts = r$phyloseq_filtered()@sam_data@names
    Fchoices = facts[facts != input$beta_fact1]

    checkboxGroupInput(
      ns("Fact2"),
      label = "Select covariable(s) to test: ",
      choices = Fchoices,
      inline = TRUE
    )
  })

  output$interac_factor <- renderUI({
    req(input$beta_fact1, r$phyloseq_filtered())
    facts = r$phyloseq_filtered()@sam_data@names
    Fchoices = facts[facts != input$beta_fact1]

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
    data
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
      p <- base_plot()$plot
      p <- p + aes(color = !!sym(input$beta_fact1))
      p <- p + stat_ellipse(aes(group = !!sym(input$beta_fact1)))
      p <- p + xlim(base_plot()$xrange) + ylim(base_plot()$yrange)
      p <- p + geom_point() + theme_bw()
      ggplotly(p) %>% config(toImageButtonOptions = list(format = "svg"))
    }, message = "Plot Beta...")
  })

  get_formula <- reactive({
    req(input$metrics, input$beta_fact1)
    form <- glue::glue('{input$metrics}.dist ~ Depth + ')
    if(!is.null(input$Fact2)){
      cov1 = paste(input$Fact2, collapse = " + ")
      form <- paste(form, glue::glue('{cov1} + {input$beta_fact1}'), sep='')
    }
    else if(!is.null(input$interFactor)){
      cov1 = paste(input$interFactor, collapse = "*")
      form <- paste(form, glue::glue('{input$beta_fact1}*{cov1}'), sep='')
    }
    else{
      form <- paste(form, glue::glue('{input$beta_fact1}'), sep='')
    }
    return(form)
  })


  betatest <- eventReactive(input$go1, {
    req(input$metrics, input$beta_fact1, r$phyloseq_filtered, r$phyloseq_filtered_norm, input$beta_norm_bool)

    if(input$beta_norm_bool==0){
      data <- r$phyloseq_filtered()
    }
    if(input$beta_norm_bool==1){
      data <- r$phyloseq_filtered_norm()
    }
    otable = otu_table(data)
    mdata = data.frame(sample_data(data))
    mdata$Depth <- sample_sums(data)

    if(any(input$metrics == c("bray", "jaccard")) ){
      fun = glue::glue("{input$metrics}.dist <<- vegdist(t(otable), distance={input$metrics})")
    }else{
      fun = glue::glue("{input$metrics}.dist <<- phyloseq::distance(data, '{input$metrics}')")
    }
    eval(parse(text=fun))

    fun = glue::glue("res.disper <- vegan::betadisper({input$metrics}.dist, mdata${input$beta_fact1})")
    eval(parse(text=fun))

    disper.anova <- anova(res.disper)
    disper.tukey <- TukeyHSD(res.disper)

    res.adonis = vegan::adonis2(as.formula(get_formula()), data = mdata, permutations = 1000)

    fun <- glue::glue('res.pairwise = pairwiseAdonis::pairwise.adonis({input$metrics}.dist, mdata[,input$beta_fact1], p.adjust.m = "fdr")')
    # fun <- glue::glue('res.pairwise = TukeyHSD(res.adonis, \"{input$beta_fact1}\")') <- marche pas TukeyHSD ne prend pas en charge les résultats d'adonis.
    eval(parse(text=fun))

    return(list(form = get_formula(), res.adonis = data.frame(res.adonis), res.pairwise = res.pairwise, res.disper = res.disper, disper.anova = disper.anova, disper.tukey = disper.tukey ))
  })

  output$adonis_formula <- renderText({
    print(betatest()$form)
  })


  output$adonistest <- DT::renderDataTable({
    betatest()$res.adonis
  })

  output$adonispairwisetest <- DT::renderDataTable({
    betatest()$res.pairwise
  })

  output$dispersionPlot <- renderPlotly({
    df1 = cbind.data.frame(distances = betatest()$res.disper$distances, group = betatest()$res.disper$group)
    print(head(df1))
   plot_ly(df1, x = ~group, y = ~distances,
           color = ~group, type = 'box') %>%
     layout(title="", yaxis = list(title = "Distance to centroid"), xaxis = list(title = 'Group'), barmode = 'stack') %>%
    config(toImageButtonOptions = list(format = "svg"))
 })

  output$dispersionTable <- DT::renderDataTable({
    betatest()$disper.anova
  })

  output$dispersionTukey <- DT::renderDataTable({
    betatest()$disper.tukey$group
  })
}

## To be copied in the UI
# mod_mod_beta_ui("mod_beta_ui_1")

## To be copied in the server
# callModule(mod_mod_beta_server, "mod_beta_ui_1")
