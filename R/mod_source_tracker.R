#' source_tracker UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session,r Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' 
mod_source_tracker_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      infoBox("",
              HTML(
                paste(
                  h4("New SourceTracker module."),
                  h5("The code ran behing this module is available at ", tags$a(href="https://github.com/danknights/sourcetracker", "github")),
                  h5("Please cite: "), 
                  h5(tags$a(href="https://www.sciencedirect.com/science/article/abs/pii/S0043135416300884","doi:10.1016/j.watres.2016.02.029")) )
              ),
              icon = icon("info-circle"), fill=TRUE, width = 12
      ),
      box(
        fluidPage(
          fluidRow(
            selectInput(
              ns("src_fact1"),
              label = "Select factor to test: ",
              choices = ""
            )
          ),
          fluidRow(
            uiOutput(ns("sources_box")),
            uiOutput(ns("sink_radio")),
          ),
          fluidRow(
            actionButton(ns('src_go'), label='Launch', style = "material-circle", color = "primary", icon = icon('play'))
          )
        ),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(
        h4('WARNING: increasing those values will slow down analysis !'),
        sliderInput(
          inputId = ns("src_burnin"),
          label = "Burnin value. number of 'burn-in' passes for Gibbs sampling.", 
          min=2, max=100, value=2, step=2
        ),
        sliderInput(
          inputId = ns("src_nrestarts"),
          label = "nrestarts value. number of times to restart the Gibbs sampling process.", 
          min=1,max=10,value=2, step=1
        ),
        title = "Advanced Settings:", width = 12, status = "warning", solidHeader = TRUE, collapsed = TRUE, collapsible = TRUE
      ),
      box(
        plotOutput(ns('boxplot')),
        title = "Plots", width = 12, status = "primary", solidHeader = TRUE
      ),
      box(
        verbatimTextOutput(ns("text_out")),
        title = "Text output:", status = "primary", width = 12, solidHeader = TRUE, collapsible = TRUE)
    )
  )
}
    
#' source_tracker Server Function
#'
#' @noRd 
#' 
#' 
mod_source_tracker_server <- function(input, output, session, r = r){
  ns <- session$ns
  
  observeEvent(r$tabs$tabselected, {
    if(r$tabs$tabselected=='source_tracker' && !isTruthy(r$phyloseq_filtered())){
      shinyalert(title = "Oops", text="Phyloseq object not present. Return to input data and validate all steps.", type='error')
    }
  })
  
  r_values <- reactiveValues(factor_list=NULL)
  
  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "src_fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })
  
  output$sources_box <- renderUI({
    req(r$sdat(), input$src_fact1)
    levels <- na.omit(levels(r$sdat()[,input$src_fact1]))
    checkboxGroupButtons(ns('sources_box'), "Choose your sources:", choices=levels, justified = TRUE, checkIcon = list(yes = icon("ok", lib = "glyphicon")))
  })
  
  output$sink_radio <- renderUI({
    req(input$sources_box)
    levels <- na.omit(levels(r$sdat()[,input$src_fact1]))
    ch <- dplyr::setdiff(levels, input$sources_box )
    radioGroupButtons(ns('sink_radio'), "Choose your sink", choices=ch, justified = TRUE, checkIcon = list(yes = icon("ok", lib = "glyphicon")))
  })
  
  
  launch_sourceTracker <- reactive({
    phy_obj <- r$phyloseq_filtered()
    ff <- glue::glue("tt <- r$sdat()${input$src_fact1} %in% c(input$sources_box, input$sink_radio )")
    eval(parse(text=ff))

    fun <- glue::glue("lmax <- length( rownames( r$sdat()[r$sdat()${input$src_fact1} == '{input$sink_radio}',] ))")
    eval(parse(text=fun))
    
    withProgress(message = 'Computing SourceTracker...', min=0, max=lmax+2, value = 0,{
      cat(file=stderr(),'prune_sample...')

      setProgress(value = 0.2, detail = 'Pruning samples...')

      fun <- glue::glue("phy_obj <- phyloseq::prune_samples(tt, phy_obj)")
      eval(parse(text=fun))
      cat(file=stderr(),'done.', "\n")
      
      phy_obj <- prune_taxa(taxa_sums(phy_obj) > 0, phy_obj)
      
      metadata <- as.data.frame(as.matrix(sample_data(phy_obj)[,input$src_fact1]))
      o_table <- as.data.frame(as.matrix(otu_table(phy_obj)))
      o_table <- t(o_table)
      
      lst <- c()
      for (src in input$sources_box){
        lst[[src]] <- 'source'
      }
      
      metadata$sourceSink <- metadata[[input$src_fact1]]
      metadata$sourceSink <- dplyr::recode(metadata$sourceSink, !!!lst)
      
      setProgress(value = 0.4, detail = 'Formatting metadata...')

      fun <- glue::glue("metadata$sourceSink <- dplyr::recode(metadata$sourceSink, '{input$sink_radio}' = 'sink')")
      eval(parse(text = fun))
      
      
      if(length(rownames(metadata)) != length(rownames(o_table))){
        stop('Number of samples in metadata differs from otu table.')
      }
      
      metadata <- apply(metadata,2,function(x) gsub("[[:punct:]]",'_',x))
      metadata <- apply(metadata,2,function(x) gsub("[[:space:]]",'_',x))
      metadata <- as.data.frame(metadata)
      
      cat(file=stderr(),'train & test',"\n")
      train <- which(metadata$sourceSink=='source')
      test <- which(metadata$sourceSink=='sink')
      
      cat(file=stderr(),'envs',"\n")
      envs <- metadata[[input$src_fact1]]
      

      alpha1 <- alpha2 <- 0.001
      cat(file=stderr(),'sourceTracker function...')

      setProgress(value = 0.6, detail = 'Source Tracker function')

      st <- sourcetracker(o_table[train,], envs[train], rarefaction_depth=1000)
      cat(file=stderr(),'done.',"\n")
      
      cat(file=stderr(),'predict function...')
      setProgress(value = 1, detail = 'Predict function (slow)')

      res <- predict(st, o_table[test,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth=1000, burnin=input$src_burnin, nrestarts=input$src_nrestarts, verbosity = 0)
      cat(file=stderr(),'done.',"\n")
    })
    res$metadata=metadata
    return(res)
  })
  
  get_box_plot <- reactive({
    tmp <- reshape2::melt(srcTracker_reactive()$proportions, id.vars=input$src_fact1)
    ggplot2::ggplot(tmp, aes(x=variable, y=value, fill=variable)) + 
      geom_boxplot() +
      ggtitle("SourceTracker Proportion per sources found in sink") + 
      xlab("Sources") + 
      ylab("Proportions") + 
      labs(fill = "")
  })
  
  srcTracker_reactive <- eventReactive(input$src_go, {
    cat(file=stderr(), 'launch sourceTracker...')
    res <- launch_sourceTracker()
    prop <- as.data.frame(res$proportions)
    fun <- glue::glue("prop[, '{input$src_fact1}' ] <- as.vector(res$metadata[rownames(res$proportions), '{input$src_fact1}'])")

    eval(parse(text=fun))
    
    prop_sd <- as.data.frame(res$proportions_sd)
    fun <- glue::glue("prop_sd[, '{input$src_fact1}' ] <- as.vector(res$metadata[rownames(res$proportions_sd), '{input$src_fact1}'])")

    eval(parse(text=fun))
    
    return(list('proportions'=prop, 'proportions_sd'=prop_sd))
  })
  
  output$text_out <- renderPrint({
    cat("Source Tracker proportions \n\n")
    print(srcTracker_reactive()$proportions)
    cat("\n\n")
    cat("Source Tracker proportions standard deviation \n\n")
    print(srcTracker_reactive()$proportions_sd)
  })
  
  output$boxplot <- renderPlot(
    get_box_plot()
  )
}
    
## To be copied in the UI
# mod_source_tracker_ui("source_tracker_ui_1")
    
## To be copied in the server
# callModule(mod_source_tracker_server, "source_tracker_ui_1")
