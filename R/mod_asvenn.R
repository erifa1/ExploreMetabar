#' asvenn UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session,r Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_asvenn_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      infoBox("",
              "Select conditions to highlight shared taxa",
              icon = icon("info-circle"), fill=TRUE, width = 10),
      box(
        selectInput(
          ns("Fact1"),
          label = "Select factor to test: ",
          choices = ""
        ),
        uiOutput(ns("lvls1")),
        numericInput(ns("minAb"), "Minimum raw abundance:", 1, min = 1, max = NA),
        actionButton(ns("go1"), "Run/Update ASVenn", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
      ),
      box(
        plotOutput(ns("venn1"), height = "800px"),
        title = "Venn Diagram:", width = 12, status = "primary", solidHeader = TRUE, height = 900
      ),
      box(
        downloadButton(outputId = ns("otable_download"), label = "Download Table"),
        DT::dataTableOutput(ns("tabvenn1")),
        title = "Venn table:", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE
      ),
      box(
        uiOutput(ns("krona_select")),
        uiOutput(ns("krona_exclud")),
        uiOutput(ns("krona_glom")),
        actionButton(ns("launch_krona"), "Generate Krona", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
        # htmlwidgets::shinyWidgetOutput(ns("krona_plot"), "Krona"),
        htmlOutput(ns("krona_plot")),
        title = "Krona plot", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE
      )
    )
  )
}


plot_krona<-function(physeq,output,variable, trim=F){
  # Check if KronaTools are installed.
  if( system(command = "which ktImportText",
             intern = FALSE,
             ignore.stdout = TRUE)) {
    stop("KronaTools are not installed. Please see https://github.com/marbl/Krona/wiki/KronaTools.")
  }
  if( is.null(tax_table(physeq)) ){
    stop("No taxonomy table available.")
  }
  if( ! variable %in% colnames(sample_data(physeq))){
    stop(paste(variable, "is not a variable in the sample data."))
  }
  if (trim == FALSE) {
    spec.char<- grepl(" |\\(|\\)", as(sample_data(physeq),"data.frame")[,variable] )
    if(sum(spec.char > 0 )){
      message("The following lines contains spaces or brackets.")
      print(paste(which(spec.char)))
      stop("Use trim=TRUE to convert them automatically or convert manually before re-run")
    }
  }
  # Melt the OTU table and merge associated metadata
  df<-psmelt(physeq)
  # Fetch only Abundance, Description and taxonomic rank names columns

  phyla <- intersect(rank_names(physeq), colnames(df))
  df<-df[ ,c("Abundance", variable, phyla) ]
  # Make sure there are no spaces left
  df[,2]<-gsub(" |\\(|\\)","",df[,2])
  # Convert the field of interest as factor.
  df[,2]<-as.factor(df[,2])
  # Create a directory for krona files
  dir.create(output)
  
  # For each level of the Description variable
  # Abundance and taxonomic assignations for each OTU are fetched
  # and written to a file that would be processed by Krona.
  for( lvl in levels(df[,2])){
    write.table(
      df[which(df[, 2] == lvl & df[,1] != 0), -2],
      file = paste0(output,"/",lvl, "taxonomy.txt"),
      sep = "\t",row.names = F,col.names = F,na = "",quote = F)
  }
  # Arguments for Krona command
  # taxonomic file and their associated labels.
  krona_args<-paste(output,"/",levels(df[,2]),
                    "taxonomy.txt,",
                    levels(df[,2]),
                    sep = "", collapse = " ")
  # Add html suffix to output
  output<-paste(output,".html",sep = "")
  # Execute Krona command
  system(paste("ktImportText",
               krona_args,
               "-o", output,
               sep = " "))
  # Run the browser to visualise the output.
  # browseURL(output)
}




#' asvenn Server Function
#'
#' @noRd
#'
#' @importFrom futile.logger flog.threshold
#' @importFrom grid grid.draw
#' @importFrom grDevices rainbow
#' @importFrom VennDiagram calculate.overlap
#' @importFrom VennDiagram venn.diagram
#' @importFrom ranomaly ASVenn_fun
#'

mod_asvenn_server <- function(input, output, session, r=r){
  ns <- session$ns
  
  observeEvent(r$tabs$tabselected, {
    if(r$tabs$tabselected=='tab_asvenn' && !isTruthy(r$phyloseq_filtered())){
      shinyalert(title = "Oops", text="Phyloseq object not present. Return to input data and validate all steps.", type='error')
    }
  })
  
  
  observe({
    req(r$phyloseq_filtered())
    updateSelectInput(session, "Fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
  })


  output$lvls1 = renderUI({
    req(input$Fact1, r$phyloseq_filtered())
    level1 <- na.omit(levels(r$sdat()[,input$Fact1]))
    checkboxGroupInput(ns("lvls1"), label = "Select up to 5 levels :",
                       choices = level1, inline = TRUE, selected = level1[1:3])

  })
  
  output$krona_select <- renderUI({
    req(input$lvls1)
    checkboxGroupInput(ns('krona_select'), "Select your shared factors:", choices=input$lvls1, inline = TRUE)
  })
  
  output$krona_glom <- renderUI({
    req(input$krona_select)
    radioButtons(ns('krona_glom'), "Select if you want to agglomerate samples by factor or not:", choices = list('TRUE'=1, 'FALSE'=0), inline=TRUE, selected=1)
  })
  
  output$krona_exclud <- renderUI({
    req(input$krona_select)
    ch <- dplyr::setdiff(input$lvls1, input$krona_select)
    checkboxGroupInput(ns('krona_exclud'), "Select your excluded factors:", choices=ch, inline = TRUE)
  })

  resVenn <- eventReactive(input$go1, {   #TF
    req(r$phyloseq_filtered())
    resVenn = ASVenn_fun(data = r$phyloseq_filtered(), output = NULL, rank = "ASV", column1 = input$Fact1, lvls = input$lvls1, shared = TRUE)

    resVenn
  })

  output$venn1 <- renderPlot({
    invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
    if(length(input$lvls1) >= 2 & length(input$lvls1) <= 5){
        resVenn()$venn_plot
    }else{showNotification("Choose 2 to 5 levels...", type="error", duration = 5)
          return(NULL)
          }
  })

  output$tabvenn1 <-DT::renderDataTable({
    resVenn()$TABf #tabvenn1()
  }, filter="top", options = list(scrollX = TRUE))

  output$otable_download <- downloadHandler(
    filename = "venn_table.csv",
    content = function(file) {
      req(tabvenn1())
      write.table(tabvenn1(), file, sep="\t", row.names=FALSE)
    }
  )
  
  get_krona_plot <- reactive({
    req(input$krona_select)
    cat(file=stderr(),"Drawing krona...", "\n")
    df <- resVenn()$TABf
    
    cat(file=stderr(),'Selected list: ', input$krona_select, "\n")
    cat(file=stderr(),'Excluded list: ', input$krona_exclud, "\n")
    
    if(length(input$krona_exclud) > 0){
      df_ex <- dplyr::select(df, input$krona_exclud )
      df_ex$sum <- rowSums(df_ex)
      ex_asv <- rownames(dplyr::filter(df_ex, df_ex$sum >= 1))
      df <- df[!(rownames(df) %in% ex_asv),]
      cat(file=stderr(),'Excluding ', length(ex_asv), "\n")
    }
    df <- dplyr::select(df, input$krona_select)
    
    df$sum <- rowSums(df)
    
    dff <- dplyr::filter(df, df$sum == length(input$krona_select))
    phy_obj <- phyloseq::phyloseq(otu_table(r$phyloseq_filtered()), tax_table(r$phyloseq_filtered()), sample_data(r$phyloseq_filtered())) 
    print(phy_obj)
    cat(file=stderr(),'prune_taxa...')
    phy_obj <- prune_taxa(rownames(dff), r$phyloseq_filtered())
    cat(file=stderr(),'done.', "\n")
 
    cat(file=stderr(),'bool vector...')
    ff <- glue::glue("tt <- r$sdat()${input$Fact1} %in% input$krona_select")
    eval(parse(text=ff))
    cat(file=stderr(),'done.', "\n")
    
    cat(file=stderr(),'prune_sample...')
    fun <- glue::glue("phy_obj <- phyloseq::prune_samples(tt, phy_obj)")
    eval(parse(text=fun))
    cat(file=stderr(),'done.', "\n")
    
    phy_obj <- prune_samples(sample_sums(phy_obj) > 0, phy_obj)
    print(phy_obj)
    phy_obj@sam_data$sample.id <- rownames(sample_data(phy_obj))
    cat(file=stderr(),"plot_krona...")
    fun <- glue::glue("sample_data(phy_obj)${input$Fact1}")
    if(input$krona_glom==1){
      plot_krona(phy_obj, '/tmp/krona', variable = input$Fact1, trim=T)
    }
    else{
      plot_krona(phy_obj, '/tmp/krona', variable = 'sample.id', trim=T)
    }
    
    cat(file=stderr(),'done.', "\n")
    return('tmp/krona.html')

    cat(file=stderr(),"DONE.", "\n")
  })
  
  krona_reactive <- eventReactive(input$launch_krona, {
    get_krona_plot()
  })
  
  addResourcePath('tmp', '/tmp')
  output$krona_plot <- renderUI({
    print(krona_reactive())
    tags$iframe(
      seamless="seamless",
      src=krona_reactive(),
      width=800,
      height=800
    )
  })
}

## To be copied in the UI
# mod_asvenn_ui("asvenn_ui_1")

## To be copied in the server
# callModule(mod_asvenn_server, "asvenn_ui_1")
