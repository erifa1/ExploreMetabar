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
      fluidRow(
        box(title = "Settings:", width = 6, status = "warning", solidHeader = TRUE,
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
          shinyWidgets::materialSwitch(inputId = ns("leaflabels"), label = "Leaf label", status = "info", value = TRUE),
          shinyWidgets::materialSwitch(inputId = ns("branchcolor"), label = "Branch color based on cluster or metadata", status = "info"),
          selectInput(
            ns("clust_fact1"),
            label = "Select metadata column to replace labels:",
            choices = ''
          ),
          actionButton(ns("launch_clust"), "Run Clustering", icon = icon("play-circle"),
                       style="color: #fff; background-color: #3b9ef5; border-color: #1a4469"
          )
        )
      ),
      fluidRow(
        box(
            verbatimTextOutput(ns('nb_clstr')),
            plotOutput(ns('dendro.plot'), height = "800px"),
            width=12, height = "1000px", title='Dendrogram', status = "primary", solidHeader = TRUE)
      ),
      fluidRow(
        box(
          plotOutput(ns('sample.by.clstr')),
          width=12, title = 'Number of samples by cluster', status = "primary", solidHeader = TRUE
        )
      ),
      fluidRow(
        box(width = 12, title = "Heatmap & Multilevel pattern analysis", status = "primary", solidHeader = TRUE,
          h2("Display heatmap of selected cluster"),
          selectInput(
            ns("clust_nb"),
            label = "Select cluster number",
            choices=1
          ),
          numericRangeInput(ns("clstr_minAb"), "Minimum taxa overall raw abundance:", c(1,1), width = NULL, separator = " to "),
          selectInput(
            ns("clstr_rank_glom"),
            label='Select rank to merge taxonomy table',
            choices='',
            selected = 1,
          ),
          plotOutput(ns('subclstr_plot'), height = "600px"),
          h2("Determine taxa specific to each cluster"),
          verbatimTextOutput(ns('indicSpe'))
        )
      )
    )
  )
}
    
#' cluster Server Functions
#'
#' @noRd 
#' @import colourvalues
#' @import dendextend
#' @import phyloseq
#' @import cluster
#' @import indicspecies
#' @import reshape2

mod_cluster_server <- function(input, output, session, r = r){

    ns <- session$ns

    observe({
      req(r$phyloseq_filtered_norm())
      updateSelectInput(session, "clstr_rank_glom",
                        choices = c( phyloseq::rank_names(r$phyloseq_filtered_norm()), "ASV" ),
                        selected = "ASV")
    })

    observe({
      req(r$phyloseq_filtered_norm())
      updateSelectInput(session, "clust_fact1",
                        choices = r$var_list())
    })

    observe({
      req(r$phyloseq_filtered_norm())
      updateNumericRangeInput(session, 'clstr_minAb',
                              "Minimum taxa overall raw abundance:", value=c(0,max(taxa_sums(r$phyloseq_filtered_norm()))))
    })

    observe({
      req(r$phyloseq_filtered_norm())
      if(is.null(phy_tree(r$phyloseq_filtered_norm(), errorIfNULL=FALSE))){
        ch1 = list("bray", "jaccard")
      }else{
        ch1 = list("bray", "jaccard", "unifrac", "wunifrac")
      }
      updateRadioButtons(session, "dist.meth",
                         choices = ch1, inline = TRUE)
    })


    compute.dist <- reactive({
      req(r$phyloseq_filtered_norm())
      dd <- phyloseq::distance(r$phyloseq_filtered_norm(), method = input$dist.meth, type="sample")
      return(dd)
    })

    compute.clust <- reactive({
      req(compute.dist())
      hc <- hclust(compute.dist(), method = input$hclust.meth)
      return(hc)
    })

    compute.k <- eventReactive(input$launch_clust, {
      req(r$phyloseq_filtered_norm(), compute.dist(), compute.clust())
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
      req(compute.clust(), compute.k())
      ct <- cutree(compute.clust(), k=compute.k())
      ct <- as.data.frame(ct)
      colnames(ct) <- c('clstr')
      ct$fact <- r$sdat()[rownames(ct),input$clust_fact1]
      return(ct)
    })

    plot.dendro <- eventReactive(input$launch_clust,{
      dd <- as.dendrogram(compute.clust())
      # browser()
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
      if(is.numeric(ct$fact)){
        p <- ct %>% ggplot(aes(x=as.factor(clstr), y=fact, fill=as.factor(clstr))) + 
          geom_boxplot() + xlab('Cluster number') + ylab(input$clust_fact1)
      } else{
        p <- ct %>%
          group_by(clstr, fact) %>%
          summarise(count = n()) %>% 
          ggplot(aes(x=(fact), y=count, fill=fact)) +
          geom_bar(stat="identity") +
          facet_grid(. ~ clstr) +
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
          xlab("") + ylab("counts")
      }
      return(p)
    })

    output$sample.by.clstr <- renderPlot({
      plot.sample.by.clstr()
    })
    
    
    get_glom_table <- reactive({
      req(input$clstr_rank_glom)
      sub.phy <- r$phyloseq_filtered_norm()
      if(input$clstr_rank_glom != 'ASV'){
        tmp <- fast_tax_glom(sub.phy, input$clstr_rank_glom)
        # FGnames <- tax_table(tmp)[,input$clstr_rank_glom]
        # nnames <- paste(substr(FGnames, 1, 50), taxa_names(tmp), sep="_")
        taxa_names(tmp) <- tax_table(tmp)[,input$clstr_rank_glom]
        sub.phy <- tmp
      }
      return(sub.phy)
    })


    plot.subtree <- reactive({
      req(compute.k(), compute.clust(), get_glom_table())
      k <- compute.k()
      dend <- compute.clust()
      
      dend_list <- get_subdendrograms(as.dendrogram(dend), k, order_clusters_as_data = TRUE)
      dd <- dend_list[[as.numeric(input$clust_nb)]]
      sub.phy <- get_glom_table()
      sub.phy <- phyloseq::prune_samples(labels(dd), sub.phy)

      sub.phy <- phyloseq::prune_taxa(phyloseq::taxa_sums(sub.phy) > input$clstr_minAb[1], sub.phy)
      sub.phy <- phyloseq::prune_taxa(phyloseq::taxa_sums(sub.phy) < input$clstr_minAb[2], sub.phy)


      otable <- phyloseq::otu_table(sub.phy)
      data.com <- reshape2::melt(otable)
      #browser()
      data.com$xlabel <- as.factor(r$sdat()[as.character(data.com$Var2),input$clust_fact1])
      names(data.com) <- c("Tax", "Sample", "Abundance", "xlabel")
      # data.com$Tax = factor(data.com$Tax, levels = sort(unique(as.character(data.com$Tax))))

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
      return(p.heat)
    })


    compute.indicSpe <- reactive({
      req(compute.k(), compute.clust(), compute.cutree(), get_glom_table())
      k <- compute.k()
      dend <- compute.clust()
      clstr <- compute.cutree()
      otable <- t(as.data.frame(otu_table(get_glom_table())))
      grp <- clstr[rownames(otable),'clstr']
      indval <- indicspecies::multipatt(otable, grp, control = permute::how(nperm = 99), duleg = TRUE)
      summary(indval)
    })


    output$indicSpe <- renderPrint({
      req(compute.k(), compute.clust(), compute.cutree())
      compute.indicSpe()
    })


    output$subclstr_plot <- renderPlot({
      p <- plot.subtree()
      p
    })


    output$nb_clstr <- renderPrint({
      req(compute.k())
      k <- compute.k()
      cat('Number of clusters: ', k)
    })

    observe({
      updateSelectInput(inputId = "clust_nb",
                        choices = 1:compute.k())
    })

  
    
}
    
## To be copied in the UI
# mod_cluster_ui("cluster_ui_1")
    
## To be copied in the server
# mod_cluster_server("cluster_ui_1")
