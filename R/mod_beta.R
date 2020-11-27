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
              icon = icon("info-circle"), fill=TRUE, width = 10),

      box(
        radioButtons(ns("norm1"), "Choose TSS normalization:", inline = TRUE,
                     choices =
                       list("TSS (recommended)", "raw"),
                     selected = c("TSS (recommended)")
        ),

        radioButtons(ns("metrics"), "Choose one index:", inline = TRUE,
                     choices =
                       list("bray", "jaccard", "unifrac", "wunifrac"),
                     selected = c("bray")
        ),

        radioButtons(ns("ordination"), "Choose one ordination:", inline = TRUE,
                     choices =
                       list("MDS", "NMDS", "CCA", "RDA"),
                     selected = c("NMDS")
        ),
        selectInput(
          ns("Fact1"),
          label = "Select main factor to test + color plot: ",
          choices = ""
        ),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(plotlyOutput(ns("plot1")),
          title = "Ordination plot:", width = 12, status = "primary", solidHeader = TRUE),
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
    updateSelectInput(session, "Fact1",
                      choices = r$data16S()@sam_data@names)
  })


  output$factor2 = renderUI({
    req(input$Fact1)
    facts = r$subglom()@sam_data@names
    Fchoices = facts[facts != input$Fact1]

    checkboxGroupInput(
      ns("Fact2"),
      label = "Select covariable(s) to test: ",
      choices = Fchoices,
      inline = TRUE
    )
  })


  Fdata <- reactive( {
    print("Beta")
      Fdata <- prune_samples(sample_names(r$data16S())[r$rowselect()], r$data16S())
      Fdata <- prune_taxa(taxa_sums(Fdata) > 0, Fdata)
      if(r$RankGlom() == "ASV"){
        Fdata <- prune_taxa(r$asvselect(), Fdata)
      }
      print(Fdata)

      if(input$norm1 != "raw"){
        print("PROPORTIONS")
        normf = function(x){ x/sum(x) }
        Fdata <- transform_sample_counts(Fdata, normf)
      }

      Fdata
  })

  depth1 <- reactive( {
    Fdata <- prune_samples(sample_names(r$data16S())[r$rowselect()], r$data16S())
    Fdata <- prune_taxa(taxa_sums(Fdata) > 0, Fdata)
    if(r$RankGlom() == "ASV"){
      Fdata <- prune_taxa(r$asvselect(), Fdata)
    }
    depth1 = sample_sums(Fdata)
  })


  ord1 <- reactive({
    req(input$ordination, input$metrics, Fdata())
    data = Fdata()
    ord1 = ordinate(data, input$ordination, input$metrics)
    ord1
  } )

  betaplot1 <- reactive({
    req(ord1(), input$metrics, input$Fact1, Fdata())
    data = Fdata()

    p1 <- plot_samples(data, ord1() , color = input$Fact1 ) + theme_bw() +
      ggtitle(paste( input$ordination, input$metrics, sep = "+" )) + stat_ellipse()
    ggplotly(p1)
  })


  output$plot1 <- renderPlotly({
    withProgress({
    betaplot1()
    }, message = "Plot Beta...")
  })


  betatest <- eventReactive(input$go1, {
    req(input$metrics, input$Fact1, Fdata())

    data = Fdata()
    otable = otu_table(data)
    mdata = data.frame(sample_data(data))
    mdata$Depth <- depth1()

    if(any(input$metrics == c("bray", "jaccard")) ){
      fun = glue::glue("{input$metrics}.dist <<- vegdist(t(otable), distance={input$metrics})")
    }else{
      fun = glue::glue("{input$metrics}.dist <<- phyloseq::distance(data, '{input$metrics}')")
    }
    eval(parse(text=fun))

    fun = glue::glue("res.disper <- vegan::betadisper({input$metrics}.dist, mdata${input$Fact1})")
    eval(parse(text=fun))

    disper.anova <- anova(res.disper)
    disper.tukey <- TukeyHSD(res.disper)

    if(is.null(input$Fact2)){
      form1 = glue::glue('{input$metrics}.dist ~ Depth + {input$Fact1}')
    }else{
      cov1 = paste(input$Fact2, collapse = " + ")
      form1 = glue::glue('{input$metrics}.dist ~ Depth + {cov1} + {input$Fact1}')
    }

    res.adonis = adonis(as.formula(form1), data = mdata, permutations = 1000)

    fun <- glue::glue('res.pairwise = pairwiseAdonis::pairwise.adonis({input$metrics}.dist, mdata[,input$Fact1], p.adjust.m = "fdr")')
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

  output$dispersionTable <- renderDataTable({
    betatest()$disper.anova
  })

  output$dispersionTukey <- renderDataTable({
    betatest()$disper.tukey$group
  })
}

## To be copied in the UI
# mod_mod_beta_ui("mod_beta_ui_1")

## To be copied in the server
# callModule(mod_mod_beta_server, "mod_beta_ui_1")
 
