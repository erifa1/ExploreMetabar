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
        actionButton(ns("go1"), "Update Test", style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        h2('Permanova Adonis Test Result: '),
        DT::dataTableOutput(ns('adonistest')),
        #verbatimTextOutput(ns("adonistest")),
        h2('Pairwise Adonis Test Results: '),
        DT::dataTableOutput(ns("adonispairwisetest")),
        h2('Dispersion results:'),
        h3('Boxplots distance to centroid for each group:'),
        plotOutput(ns("dispersionPlot")),
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
#'
#' @noRd
mod_beta_server <- function(input, output, session, r = r){
  ns <- session$ns

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
    xrange[1] <- layer_scales(p)$x$range$range[1] - abs(layer_scales(p)$x$range$range[1])*2
    xrange[2] <- layer_scales(p)$x$range$range[2] + abs(layer_scales(p)$x$range$range[2])*2
    
    yrange <- c()
    yrange[1] <- layer_scales(p)$y$range$range[1] - abs(layer_scales(p)$y$range$range[1])*2
    yrange[2] <- layer_scales(p)$y$range$range[2] + abs(layer_scales(p)$y$range$range[2])*2
    return(list('plot'=p, 'xrange'=xrange, 'yrange'=yrange))
  })

  # betaplot1 <- eventReactive(input$launch_beta,{
  #   cat(file=stderr(), "betaplot1 fun...", "\n")
  #   req(input$ordination, input$metrics, input$beta_norm_bool, input$beta_fact1, r$phyloseq_filtered, r$phyloseq_filtered_norm)
  #   if(input$beta_norm_bool==0){
  #     data <- r$phyloseq_filtered()
  #   }
  #   if(input$beta_norm_bool==1){
  #     data <- r$phyloseq_filtered_norm()
  #   }
  #   cat(file=stderr(), "betaplot1 fun ordinate", "\n")
  #   ord1 = ordinate(data, input$ordination, input$metrics)
  #   cat(file=stderr(), "betaplot1 fun plot_samples", "\n")
  #   p1 <- plot_samples(data, ord1 , color = input$beta_fact1 ) + theme_bw() +
  #     ggtitle(paste( input$ordination, input$metrics, sep = "+" )) + stat_ellipse()
  #   cat(file=stderr(), "betaplot1 fun done.", "\n")
  #   ggplotly(p1) %>% config(toImageButtonOptions = list(format = "svg"))
  # })


  # output$plot1 <- plotly::renderPlotly({
  #   withProgress({
  #   betaplot1()
  #   }, message = "Plot Beta...")
  # })
  
  output$plot1 <- plotly::renderPlotly({
    withProgress({
      p <- base_plot()$plot
      p <- p + aes(color = !!sym(input$beta_fact1))
      p <- p + stat_ellipse(aes(group = !!sym(input$beta_fact1)))
      p <- p + xlim(base_plot()$xrange) + ylim(base_plot()$yrange)
      p <- p + geom_point() + theme_bw()
      ggplotly(p) %>% config(toImageButtonOptions = list(format = "svg"))
    }, message = "Plot Beta...")
  })

  # output$download_svg <- downloadHandler(
  #   filename = 'beta_div.svg',
  #   content = function(file){
  #     req(betaplot1())
  #     p <- betaplot1()
  #     orca(p)
  #   }
  # )


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

    if(is.null(input$Fact2)){
      form1 = glue::glue('{input$metrics}.dist ~ Depth + {input$beta_fact1}')
    }else{
      cov1 = paste(input$Fact2, collapse = " + ")
      form1 = glue::glue('{input$metrics}.dist ~ Depth + {cov1} + {input$beta_fact1}')
    }

    res.adonis = adonis(as.formula(form1), data = mdata, permutations = 1000)

    fun <- glue::glue('res.pairwise = pairwiseAdonis::pairwise.adonis({input$metrics}.dist, mdata[,input$beta_fact1], p.adjust.m = "fdr")')
    # fun <- glue::glue('res.pairwise = TukeyHSD(res.adonis, \"{input$beta_fact1}\")') <- marche pas TukeyHSD ne prend pas en charge les résultats d'adonis.
    eval(parse(text=fun))

    return(list(res.adonis = res.adonis$aov.tab, res.pairwise = res.pairwise, res.disper = res.disper, disper.anova = disper.anova, disper.tukey = disper.tukey ))
  })


  output$adonistest <- renderDataTable({
    betatest()$res.adonis
  })

  output$adonispairwisetest <- renderDataTable({
    betatest()$res.pairwise
  })

  output$dispersionPlot <- renderPlot({boxplot(betatest()$res.disper, col=unique(betatest()$res.disper$group),las=2)})

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
