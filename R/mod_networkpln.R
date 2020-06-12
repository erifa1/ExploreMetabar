#' networkpln UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_networkpln_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    fluidPage(
      uiOutput(ns("factor1")),
      fluidRow(column(3, uiOutput(ns("cond1")) ),
               column(3, uiOutput(ns("cond2")) )
      ),
      # numericInput(ns("minAb"), "Minimum overall raw abundance for each ASV:", value = 1000, min = 0, max = 1e100, step = 1),
      actionButton(ns("go1"), "Run PLN", icon = icon("play-circle")),
      verbatimTextOutput(ns("print1")),
      
      plotOutput(ns("plotnet"))
    )
  )
}
    
#' networkpln Server Function
#'
#' @noRd 
#' 
#' @importFrom igraph degree
#' @importFrom igraph delete.vertices
#' @importFrom igraph layout_in_circle
#' @import PLNmodels
#' 
mod_networkpln_server <- function(input, output, session = session, r = r){
  ns <- session$ns
 
  
  output$factor1 = renderUI({
    selectInput(
      ns("Fact1"),
      label = "Select factor to test: ",
      choices = r$subglom()@sam_data@names
    )
  })
  
  output$cond1 = renderUI({
    req(input$Fact1)
    selectInput(ns("Cond1"), 
                label = "Select Condition 1 to compare: ", 
                choices = unique(r$subglom()@sam_data[,input$Fact1])
    )
  })
  
  output$cond2 = renderUI({
    req(input$Cond1)
    Conds = unique(r$subglom()@sam_data[,input$Fact1])
    if(length(Conds) <= 2){
      print(Conds)
      choices2 = Conds
    }else{
      choices2 = Conds[Conds != input$Cond1]
    }
    print(choices2)
    
    
    selectInput(ns("Cond2"), 
                label = "Select Condition 2 to compare: ", 
                choices = choices2
    )
    
  })
  
  
  my_data <- eventReactive(input$go1, {
    # req(r$subtax())
    data0 <- r$subglom()
  
    prepdata = function(x, rank = r$RankGlom()){
      counts <- as(phyloseq::otu_table(x), "matrix")
      rownames(counts) = x@tax_table@.Data[,rank]
      ## extract covariates (or prepare your own)
      covariates <- phyloseq::sample_data(x)
      ## prepare data
      my_data <- prepare_data(counts = counts, covariates = covariates)
      # str(my_data)
      return(my_data)
    }


    # Cond1
    print("network 1")
    eval(parse(
      text= glue::glue("data1 <- subset_samples(data0, {input$Fact1} == '{input$Cond1}')")
      ))

    # select = taxa_sums(data1)>500   ## à mettre plus haut.
    # data1 <- prune_taxa(taxa_sums(data1)>input$minAb, data1)
    print("prep data")
    my_data <- prepdata(data1, rank = r$RankGlom())
    print(my_data)
    print(ntaxa(my_data))
  
    
    eval(parse(
      text= glue::glue("data2 <- subset_samples(data0, {input$Fact1} == '{input$Cond2}')")
    ))
    # data2 <- prune_taxa(taxa_sums(data2)>0, data2)
    
    my_data2 <- prepdata(data2, rank = r$RankGlom())
    
    LL=list()
    LL$my_data = my_data
    LL$my_data2 = my_data2 
    
    LL
    
  })
    
  plnproc <- reactive({
    network_models <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = my_data()$my_data)
    model_StARS <- getBestModel(network_models, "StARS") # if StARS is requested, stabiltiy selection is performed if needed

    print(model_StARS)
    
    # Cond2
    print("network 2")

    network_models <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = my_data()$my_data2)
    model_StARS2 <- getBestModel(network_models, "StARS") # if StARS is requested, stabiltiy selection is performed if needed

    #Plotting
    my_graph1 = plot(model_StARS, plot = FALSE) #? plot.PLNnetworkfit
    my_graph2 = plot(model_StARS2, plot = FALSE)

    NODES0= table(names(c(which(degree(my_graph1) == 0), which(degree(my_graph2) == 0))))
    NODES0 = names(NODES0[NODES0==2])

    LL = list()
    LL$my_graph1 = my_graph1
    LL$my_graph2 = my_graph2
    LL$NODES0 = NODES0

    print(LL)

    LL
  })
  
  
  # output$print1 <- renderPrint({
  #   
  #   print(my_data())
  #   
  #   print("toto")
  #   
  # })
  
  output$plotnet = renderPlot({
    # plot(plnproc())
    LL <- plnproc()
    
    par(mfrow=c(1,2))
    plot(delete.vertices(LL$my_graph1, v = LL$NODES0), layout = layout_in_circle)
    title(glue::glue("{input$Fact1}  {input$Cond1}"))
    
    plot(delete.vertices(LL$my_graph2, v = LL$NODES0), layout = layout_in_circle)
    title(glue::glue("{input$Fact1}  {input$Cond2}"))
    
  })
  
  
  
  # output$plotnet = plotOutput({
  #   LL <- plnproc()
  # 
  #   par(mfrow=c(1,2))
  #   plot(delete.vertices(LL$my_graph1, v = LL$NODES0), layout = layout_in_circle)
  #   title(glue::glue("{input$Fact1}  {input$Cond1}"))
  # 
  #   plot(delete.vertices(LL$my_graph2, v = LL$NODES0), layout = layout_in_circle)
  #   title(glue::glue("{ENV}  {input$Cond2}"))
  # })
  
  
}
    
## To be copied in the UI
# mod_networkpln_ui("networkpln_ui_1")
    
## To be copied in the server
# callModule(mod_networkpln_server, "networkpln_ui_1")
 
