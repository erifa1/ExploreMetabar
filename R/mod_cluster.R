#' cluster UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' @import dendextend
mod_cluster_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
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
        shinyWidgets::materialSwitch(inputId = ns("leaflabels"), label = "Leaf label"),
        shinyWidgets::materialSwitch(inputId = ns("branchcolor"), label = "Branch color cluster or metadata"),
        selectInput(
          ns("clust_fact1"),
          label = "Select metadata column to replace lables",
          choices = ''
        ),
        actionButton(ns("launch_clust"), "Run Clustering", icon = icon("play-circle"),
                     style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"
        )
      )
      ),
      fluidRow(
        box(
          fluidPage(
            plotOutput(ns('dendro.plot')),
            verbatimTextOutput(ns('nb_clstr'))
            
          ),width=12, title='Dendrogram'
        )
      ),
      fluidRow(
        box(
          plotOutput(ns('sample.by.clstr')),
          width=12, title = 'Number of samples by cluster'
        )
      ),
      fluidRow(
        box(
          uiOutput(ns('clstr_input')),
          uiOutput(ns('sub_input')),
          width = 12,
          height = 12
        )
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
      updateSelectInput(session, "clust_fact1",
                        choices = r$phyloseq_filtered()@sam_data@names)
    })
    
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
        k.best <- which.max(Si)
      }
      if(input$k.meth == "pearson"){
        kt <- data.frame(k=2:nrow(t(otu_table(r$phyloseq_filtered_norm()))), r=0)
        for (i in 2:(nrow(t(otu_table(r$phyloseq_filtered_norm())))-1)){
          gr <- cutree(compute.clust(), i)
          veg <- as.data.frame(as.factor(gr))
          distgr <- cluster::daisy(veg, "gower")
          mt <- cor(compute.dist(), distgr, method = "pearson")
          kt[i,2] <- mt
        }
        k.best <- which.max(kt$r)
      }
      return(k.best)
    })
    
    compute.cutree <- reactive({
      ct <- cutree(compute.clust(), k=compute.k())
      ct <- as.data.frame(ct)
      colnames(ct) <- c('clstr')
      ct$fact <- r$sdat()[rownames(ct),input$clust_fact1]
      return(ct)
    })
    
    plot.dendro <- eventReactive(input$launch_clust,{
      dd <- as.dendrogram(compute.clust())
      colours <- colourvalues::colour_values(r$sdat()[labels(dd),input$clust_fact1], palette="viridis")
      dd <- dendextend::color_labels(dd, col=colours)
      if(input$branchcolor){
        dd <- dendextend::color_branches(dd, col=colours)
      } else{
        dd <- dendextend::color_branches(dd, compute.k())
      }
      # dd <- dendextend::labels_cex(dd, 1.2)
      # browser()
      par(mar = c(2,2,2,15))
      if(input$leaflabels){
        labels(dd) <- r$sdat()[labels(dd),input$clust_fact1]
        plot(dd, horiz = TRUE)
      } else{
        plot(dd, leaflab="none", horiz = TRUE)
      }
        
      
    })
    
    
    output$dendro.plot <- renderPlot({
      plot.dendro()
    })
    
    
    
    plot.sample.by.clstr <- eventReactive(input$launch_clust, {
      ct <- compute.cutree()
      ct %>%
        group_by(clstr, fact) %>%
        summarise(count = n()) %>%
        ggplot(aes(x=(fact), y=count, fill=fact)) +
        geom_bar(stat="identity") +
        facet_grid(. ~ clstr) +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      
    })
    
    output$sample.by.clstr <- renderPlot({
      plot.sample.by.clstr()
    })
    
    
    plot.subtree <- reactive({
      k <- compute.k()
      dend <- compute.clust()
      dend_list <- get_subdendrograms(as.dendrogram(dend), k)
      dd <- dend_list[[as.numeric(input$clust_nb)]]
      sub.phy <- phyloseq::prune_samples(labels(dd), r$phyloseq_filtered_norm())
      # sub.phy <- phyloseq::prune_taxa(phyloseq::taxa_sums(sub.phy)>0, sub.phy)
      order.dendrogram(dd) <- as.integer(rank(order.dendrogram(dd)))
      gplots::heatmap.2(as.matrix(as.data.frame(phyloseq::otu_table(sub.phy))), Colv = dd, trace = "none", col = viridis::viridis(100))
      # if(input$leaflabels){
      #   labels(dd) <- r$sdat()[labels(dd),input$clust_fact1]
      #   plot(dd, horiz = TRUE)
      # } else{
      #   plot(dd, leaflab="none", horiz = TRUE)
      # }
    })
    
    output$plot.subtree <- renderPlot({
      plot.subtree()
    })
    
    
    observeEvent(input$launch_clust, {
      output$nb_clstr <- renderPrint(({
        k <- compute.k()
        cat('Number of clusters: ', k)
      }))
      
      
      output$clstr_input <- renderUI({
        selectInput(
          ns("clust_nb"),
          label = "Select cluster number",
          choices = 1:compute.k()
        )
      })
      
      output$sub_input <- renderUI({
        plotOutput(ns('plot.subtree'))
      })
    })
    
}
    
## To be copied in the UI
# mod_cluster_ui("cluster_ui_1")
    
## To be copied in the server
# mod_cluster_server("cluster_ui_1")
