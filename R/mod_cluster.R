#' cluster UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_cluster_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      box(
        radioButtons(ns("dist.meth"), "Choose distance method:", inline = TRUE,
                     choices ='',
                     selected = c("bray")
        ),
        radioButtons(ns("hclust.meth"), "Choose clustering method:", inline = TRUE,
                     choices =c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                     selected = c("ward.D2")
        ),
        radioButtons(ns("k.meth"), "Choose optimal number of cluster method:", inline = TRUE,
                     choices =c("silhouette", "pearson"),
                     selected = c("silhouette")
        ),
        actionButton(ns("launch_clust"), "Run Clustering", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"
        )
      ),
      box(
        plotOutput(ns('dendro.plot'))
      )
    )
  )
}
    
#' cluster Server Functions
#'
#' @noRd 
mod_cluster_server <- function(input, output, session, r = r){

    ns <- session$ns
    
    
    observe({
      req(r$phyloseq_filtered())
      if(is.null(phy_tree(r$phyloseq_filtered(), errorIfNULL=FALSE))){
        ch1 = list("bray", "jaccard")
      }else{
        ch1 = list("bray", "jaccard", "unifrac", "wunifrac")
      }
      updateRadioButtons(session, "dist.meth",
                         choices = ch1, inline = TRUE)
    })

    
    compute.dist <- reactive({
      dd <- phyloseq::distance(r$phyloseq_filtered_norm(), method = input$dist.meth, type="sample")
      return(dd)
    })
    
    compute.clust <- reactive({
      hc <- hclust(compute.dist(), method = input$hclust.meth)
      return(hc)
    })
    
    compute.k <- reactive({
      if(input$k.meth == "silhouette"){
        Si <- numeric(nrow(t(otu_table(r$phyloseq_filtered_norm()))))
        for (k in 2:(nrow(t(otu_table(r$phyloseq_filtered_norm()))) -1)){
          sil <- cluster::silhouette(cutree(compute.clust(), k = k), compute.dist())
          Si[k] <- summary(sil)$avg.width
        }
        k.best <- whiwh.max(Si)
      }
      return(k.best)
    })
    
    
    plot.dendro <- eventReactive(input$launch_clust,{
      dd <- as.dendrogram(compute.clust())
      colours <- colourvalues::colour_values(as.matrix(sample_data(r$phyloseq_filtered_norm()))[labels(dd),"SampleType"], palette="viridis")
      dd <- dendextend::color_branches(dd, col=colours)
      plot(dd, leaflab="none") 
    })
    
    output$dendro.plot <- renderPlot({
      plot.dendro()
    })
}
    
## To be copied in the UI
# mod_cluster_ui("cluster_ui_1")
    
## To be copied in the server
# mod_cluster_server("cluster_ui_1")
