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
      useShinyalert(),
      fluidRow(
        infoBox("",
                "Select conditions to highlight shared taxa",
                icon = icon("info-circle"), fill=TRUE, width = 8
                )
      ),
      fluidRow(
        box(
          selectInput(
            ns("Fact1"),
            label = "Select factor to test: ",
            choices = ""
          ),
          uiOutput(ns("lvls1")),
          numericInput(ns("minAb"), "Minimum raw abundance to detect a taxa in group of samples:", 1, min = 1, max = NA),
          actionButton(ns("go1"), "Run/Update ASVenn", icon = icon("play-circle"),
                       style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
          title = "Settings:", width = 12, status = "warning", solidHeader = TRUE
        )
      ),
      fluidRow(
        box(
        # plotOutput(ns("venn1"), height = "800px"),
        
          imageOutput(ns("venn1"), width = "100%", height = "100%"),
          title = "Venn Diagram VennR:", width = 12, status = "primary", solidHeader = TRUE,
          collapsible = TRUE, collapsed = FALSE
        )
      ),
      fluidRow(
        box(
          plotOutput(ns("venn2"), height = "800px"),
          title = "Venn Diagram classic:", width = 12, status = "primary", solidHeader = TRUE, 
          collapsible = TRUE, collapsed = TRUE
        )
      ),
      fluidRow(
        box(
          downloadButton(outputId = ns("otable_download"), label = "Download Table"),
          DT::dataTableOutput(ns("tabvenn1")),
          title = "Venn table:", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE
        )
      ),
      fluidRow(
        box(
          plotly::plotlyOutput(ns('radar_chart'), width = '100%', height = '100%'),
          title = "Radar Chart:", width = 12, status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE
        )
      ),
      fluidRow(
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
  )
}


plot_krona <- function(physeq,output,variable, trim=F){
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
#' @importFrom grDevices rainbow recordPlot replayPlot
#' @importFrom venn venn
#' @importFrom qdapTools mtabulate
#' @importFrom nVennR plotVenn
#' @import ggpolypath


mod_asvenn_server <- function(input, output, session, r=r){
  ns <- session$ns  
  
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

  
  addResourcePath('tmp', '/tmp')
  resVenn <- eventReactive(input$go1, {
    req(r$phyloseq_filtered(), input$lvls1)
    flog.info('compute Venn diagram...')
    if(length(input$lvls1) < 2 || length(input$lvls1) > 5){
      shinyalert("Oops!", "You need to choose between 2 and 5 factors...", type = "error")
    }
    else{
      res <-
      TFdata <- list()
      TFtax <- tibble(taxa = character(), taxo = character())
      for(lvl in input$lvls1){
        flog.info(lvl)
        fun <- paste("data.tmp <- subset_samples(r$phyloseq_filtered(), ",input$Fact1," %in% '",lvl,"')",sep="")
        eval(parse(text=fun))
        sp_data <- prune_taxa(taxa_sums(data.tmp) > 0, data.tmp)
  
        abund_to_zero = function(x){
          x[x < input$minAb] <- 0
          return(x)
        }
        sp_data <- transform_sample_counts(sp_data, fun = abund_to_zero)
        sp_data <- prune_taxa(taxa_sums(sp_data) > 0, sp_data)
        
        TT = cbind(otu_table(sp_data),tax_table(sp_data))
        
        TFdata[[lvl]] <- TT
        TFtax <- dplyr::full_join(TFtax, as_tibble(cbind(taxa = row.names(TT), taxo =  as.character(apply(TT[,colnames(tax_table(sp_data))], 1, paste, collapse=";") ) )), by = c("taxa", "taxo"))
        row.names(TFtax[[lvl]]) = TFtax[[lvl]][,1]
      }
      
      TF <- sapply(TFdata, row.names, simplify = FALSE)
      names(TF) = input$lvls1

      outfile <- tempfile(fileext='.svg')
      venn.res <- nVennR::plotVenn(TF, showPlot = T, labelRegions = T, systemShow=F, outFile = outfile)
      
      res$svg.obj <- list(src = normalizePath(outfile), width = "100%", height = "100%")
      v.table <- as_tibble(t(qdapTools::mtabulate(TF)), rownames = "taxa")
      v.table <- full_join(v.table, TFtax, by = 'taxa')
      res$v.table <- v.table
      res$TF <- TF

      return(res)
    }
  })

  
  output$venn1 <- renderImage(
    resVenn()$svg.obj
  , deleteFile=TRUE)

  output$venn2 <- renderPlot({
    invisible(flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
        # grid.draw
        # grDevices::replayPlot(resVenn()$venn.plot2)
        venn::venn(resVenn()$TF, zcol = rainbow(7), ilcs = 1.5, sncs = 2, 
                          ggplot = FALSE) 
  })


  output$tabvenn1 <-DT::renderDataTable({
    resVenn()$v.table
  }, filter="top", selection = "single", options = list(scrollX = TRUE))

  output$otable_download <- downloadHandler(
    filename = "venn_table.csv",
    content = function(file) {
      req(resVenn()$v.table)
      write.table(resVenn()$v.table, file, sep="\t", row.names=FALSE)
    }
  )
  
  get_radar_data <- reactive({
    req(input$tabvenn1_row_last_clicked, input$lvls1)
    
    fun <- paste("data.tmp <- subset_samples(r$phyloseq_filtered(), ",input$Fact1," %in% c('",paste(input$lvls1, collapse='\',\''),"'))",sep="")
    eval(parse(text=fun))
    
    obj <- prune_taxa(pull(resVenn()$v.table[input$tabvenn1_row_last_clicked,1]), data.tmp)
    ot <- as.data.frame(otu_table(obj))
    ot <- as.data.frame(t(ot))
    mt <- as.data.frame(as.matrix(sample_data(obj)))
    ot[input$Fact1] <- as.vector(mt[rownames(ot), input$Fact1])
    
    return(ot)
  })
  
  get_radar <- reactive({
    dt <- get_radar_data()
    fig <- plotly::plot_ly(x =~dt[,2], y=~dt[,1], type = "box" )
    return(fig)
  })
  
  
  output$radar_chart <- plotly::renderPlotly({
    get_radar()
  })

  get_krona_plot <- reactive({
    req(input$krona_select)
    cat(file=stderr(),"Drawing krona...", "\n")
    df <- resVenn()$v.table

    cat(file=stderr(),'Selected list: ', input$krona_select, "\n")
    cat(file=stderr(),'Excluded list: ', input$krona_exclud, "\n")

    if(length(input$krona_exclud) > 0){
      df_ex <- dplyr::select(df, c('taxa', input$krona_exclud) )
      df_ex$sum <- rowSums(select_if(df_ex, is.numeric ))
      ex_asv <- rownames(dplyr::filter(df_ex, df_ex$sum >= 1))
      df <- df[!(df$taxa %in% ex_asv),]
      cat(file=stderr(),'Excluding ', length(ex_asv), "\n")
    }
    df <- dplyr::select(df, c('taxa', input$krona_select))

    df$sum <- rowSums(select_if(df, is.numeric))
  
    dff <- dplyr::filter(df, df$sum == length(input$krona_select))
    phy_obj <- phyloseq::phyloseq(otu_table(r$phyloseq_filtered()), tax_table(r$phyloseq_filtered()), sample_data(r$phyloseq_filtered()))

    cat(file=stderr(),'prune_taxa...')
    phy_obj <- prune_taxa(dff$taxa, phy_obj)
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


  output$krona_plot <- renderUI({
    print(krona_reactive())
    tags$iframe(
      seamless="seamless",
      src=krona_reactive(),
      width="100%",
      height=800
    )
  })
}

## To be copied in the UI
# mod_asvenn_ui("asvenn_ui_1")

## To be copied in the server
# callModule(mod_asvenn_server, "asvenn_ui_1")
