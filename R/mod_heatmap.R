#' heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import phyloseq
#' @import ggplot2
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly ggplotly

mod_heatmap_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(

      tabsetPanel(
          tabPanel("Heatmap on raw abundance with ecological distances",
          fluidRow(
            infoBox("",
            HTML(
              paste(
                h4("New Heatmap module."),
                h5("Ecology organized heatmap. Row an columns are organized using ordination methods (NMDS, PCA) and computed on ecological distances (bray, unifrac, jaccard)."),
                h5("Please cite: "),
                h5(tags$a(href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-45", "doi:10.1186/1471-2105-11-45")) )
              ),
              icon = icon("info-circle"), fill=TRUE, width = 12
            )
          ),

              fluidRow(
                box(title = "Settings", width = 12, status = "warning", solidHeader = TRUE,
                  selectInput(
                    ns("src_fact1"),
                    label = "Sample label",
                    choices = ""
                  ),
                  selectInput(
                    ns("ord.method"),
                    label = "Ordination method",
                    choices = c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"),
                    selected = 'NMDS'
                  ),
                  selectInput(
                    ns("dist"),
                    label = "Ecological distance method ",
                    choices = c("bray", "jaccard", "dpcoa", "unifrac", "wunifrac"),
                    selected = 'NMDS'
                  ),
                  actionButton(ns("launch_heatmap"), "Launch heatmap", icon = icon("play-circle"),
                               style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
                )
              ),
              dropdown(
                tags$h3("Plot size"),
                sliderInput(ns('plot_height'), label = 'Plot Height', min = 300, max = 2000, step = 100, value = 300),
                sliderInput(ns('plot_width'), label = 'Plot Width', min = 300, max = 2000, step = 100, value= 600),
                style = "unite",
                icon = icon('cogs'),
                width = "300px"
              ),
              shinycustomloader::withLoader(
                plotOutput(ns('heatmap_plot')),
              type = "html", loader = "loader1"
              )

            ),

          tabPanel("Heatmap on normalized abundances",
                box(title = "Settings", width = 12, status = "warning", solidHeader = TRUE,
                  h5("The heatmap is generated with preformated (filtered, agglomerated to rank...) data and normalized with method chosen in the 'Input data' module (VST is recommended)."),
                  selectInput(
                    ns("fact2"),
                    label = "Factor to separate samples",
                    choices = ""
                  ),
                  checkboxInput(ns("clust1"), label = "Default clustering (euclidean, complete)", value = TRUE),
                  actionButton(ns("launch_heatmap2"), "Launch heatmap", icon = icon("play-circle"),
                               style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"),
                ),

                box(title = "Heatmap", width = 12, status = "warning", solidHeader = TRUE,
                  downloadButton(outputId = ns("heatmap_download"), label = "Download heatmap"),
                  downloadButton(outputId = ns("heatmap_table_download"), label = "Download heatmap tsv table"),
                  plotlyOutput(ns("heatmap2"))
                )

              )
          )
      )
  )
}


#' heatmap Server Function
#'
#' @noRd
mod_heatmap_server <- function(input, output, session, r){
  ns <- session$ns

  observe({
    req(r$phyloseq_filtered(), r$phyloseq_filtered_norm())
    updateSelectInput(session, "src_fact1",
                      choices = r$phyloseq_filtered()@sam_data@names)
    updateSelectInput(session, "fact2",
                  choices = r$phyloseq_filtered_norm()@sam_data@names)
  })

  get_heatmap <- eventReactive(input$launch_heatmap,{
    req(input$src_fact1, r$phyloseq_filtered())
    phyloseq::plot_heatmap(r$phyloseq_filtered(), method = input$ord.method, distance = input$dist, sample.label = input$src_fact1)
  })

  plot_height <- reactive({
    return(input$plot_height)
  })

  plot_width <- reactive({
    return(input$plot_width)
  })

  observe({
    output$heatmap_plot <- renderPlot({
      withProgress(message = 'Computing heatmap...',{
        get_heatmap()
      })
    }, height = plot_height(), width = plot_width())
  })


#Heatmap2
  table_heatmap1 <- reactive({
    req(r$phyloseq_filtered_norm())
    cat(file=stderr(), 'heatmap2', "\n")
    dataglom <- r$phyloseq_filtered_norm()
    otable <- otu_table(dataglom)

    if(input$clust1){
      h1 <- hclust(dist(otable))
      otable <- otable[h1$labels,]
    }

    sdata <- as.data.frame(as.matrix(sample_data(dataglom)))

    data.com <- reshape2::melt(otable)
    data.com$xlabel <- as.factor(sdata[as.character(data.com$Var2),match(input$fact2,names(sdata))])
    names(data.com) <- c("Tax", "Sample", "Abundance", "xlabel")

    #order in alphabetical order if no clustering.
    if(!input$clust1){
      data.com$Tax = factor(data.com$Tax, levels = sort(unique(as.character(data.com$Tax))))
    }

    cat(file=stderr(), 'done', "\n")
    print(head(data.com))
    data.com
  })

  heatmap_plot2 <- eventReactive(input$launch_heatmap2,{
      req(table_heatmap1())
      cat(file=stderr(), 'heatmap2_plot', "\n")
      data.com <- table_heatmap1()
      p.heat <- ggplot(data.com, aes(x = Sample, y = Tax)) + geom_tile(aes(fill = Abundance))
      p.heat <- p.heat + scale_fill_distiller("Normalized abundance", palette = "RdYlBu") + theme_bw()

      # Make bacterial names italics
      p.heat <- p.heat + theme(axis.text.y = element_text(colour = 'black',
      size = 10,
      face = 'italic'))
      # Make seperate samples based on main varaible
      p.heat <- p.heat + facet_grid(~xlabel, scales = "free")

      p.heat <- p.heat + ylab("Taxa")

      #Clean the x-axis
      p.heat <- p.heat + theme(axis.title.x=element_blank(),
      axis.text.x=element_text(angle = 90),
      axis.ticks.x=element_blank())

      # Clean the facet label box
      p.heat <- p.heat + theme(legend.key = element_blank(),
      strip.background = element_rect(colour="black", fill="white"))
      cat(file=stderr(), 'done', "\n")
      # print(str(p.heat))
      # save(p.heat, file = "~/Bureau/debug_explore.rdata")
    p.heat

    })

  output$heatmap2 <- renderPlotly({
    req(heatmap_plot2())
    cat(file=stderr(), 'heatmap2_render', "\n")
    showNotification("Rendering plot, please wait.", type="message", duration = 3)
    p2 <- heatmap_plot2()
    p2 %>% config(toImageButtonOptions = list(format = "svg"))
  })

  output$heatmap_download <- downloadHandler(
    filename = "heatmap_normalized.html",
    content = function(file) {
      req(heatmap_plot2())
      p2 <- heatmap_plot2()
      pltly.heat <- ggplotly(p2)
      saveWidget(pltly.heat, file= file)
    }
  )

  output$heatmap_table_download <- downloadHandler(
    filename = "heatmap_table.csv",
    content = function(file) {
      req(table_heatmap1())
      data.com <- table_heatmap1()
      write.table(data.com, file, sep="\t", row.names=FALSE)
    }
  )


}

## To be copied in the UI
# mod_heatmap_ui("heatmap_ui_1")

## To be copied in the server
# callModule(mod_heatmap_server, "heatmap_ui_1")
