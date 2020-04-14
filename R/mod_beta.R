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
mod_beta_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      h1("Beta diversity analysis"),
      box(
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
        title = "Settings:", width = 12, status = "primary", solidHeader = TRUE
      ),
      box(plotlyOutput(ns("plot1")),
          title = "Ordination plot:", width = 12, status = "primary", solidHeader = TRUE),
      box(
      uiOutput(ns("factor2")),
      actionButton(ns("go1"), "Update Test"),
      verbatimTextOutput(ns("testprint")), 
      title = "Permanova with adonis:", width = 12, status = "primary", solidHeader = TRUE)
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
      Fdata
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
    mdata$Depth <- sample_sums(data)
    
    if(any(input$metrics == c("bray", "jaccard")) ){
      fun = glue::glue("{input$metrics}.dist <<- vegdist(t(otable), distance={input$metrics})")
      eval(parse(text=fun))
    }else{
      fun = glue::glue("{input$metrics}.dist <<- phyloseq::distance(data, '{input$metrics}')")
      eval(parse(text=fun))
    }
    
    if(is.null(input$Fact2)){
      form1 = glue::glue('{input$metrics}.dist ~ Depth + {input$Fact1}')
    }else{
      cov1 = paste(input$Fact2, collapse = " + ")
      form1 = glue::glue('{input$metrics}.dist ~ Depth + {cov1} + {input$Fact1}')
    }
    
    
    print(form1)
    res1 = adonis(as.formula(form1), data = mdata, permutations = 1000)
    print(res1)
  })
  
  
  output$testprint <- renderPrint({
    # print(input$Fact2)
    # print(str(betatest()))
    betatest()
  })
  
 
}
    
## To be copied in the UI
# mod_mod_beta_ui("mod_beta_ui_1")
    
## To be copied in the server
# callModule(mod_mod_beta_server, "mod_beta_ui_1")
 
