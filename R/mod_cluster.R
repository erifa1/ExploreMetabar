#' cluster UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom bslib layout_sidebar sidebar accordion accordion_panel navset_card_underline nav_panel
#' @importFrom bsicons bs_icon

mod_cluster_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "Configuration",
      open = "desktop",
      width = "350px",
      htmltools::p(
        "Hierarchical clustering of samples. Optimal cluster number is estimated automatically from the chosen method.",
        style = "font-size: 0.9em; color: grey;"
      ),
      tags$hr(),
      accordion(
        id = ns("config_accordion"),
        open = "Clustering settings",
        multiple = TRUE,
        accordion_panel(
          "Clustering settings",
          icon = bs_icon("diagram-3"),
          radioButtons(ns("dist.meth"), "Distance method:", inline = TRUE,
                       choices = '', selected = "bray"),
          radioButtons(ns("hclust.meth"), "Clustering method:", inline = TRUE,
                       choices = c("ward.D", "ward.D2", "single", "complete",
                                   "average", "mcquitty", "median", "centroid"),
                       selected = "ward.D2"),
          radioButtons(ns("k.meth"), "Optimal k method:", inline = TRUE,
                       choices = c("silhouette", "pearson"), selected = "silhouette"),
          selectInput(ns("clust_fact1"), label = "Metadata column for labels:", choices = ''),
          shinyWidgets::materialSwitch(ns("leaflabels"), label = "Show leaf labels",
                                       status = "info", value = TRUE),
          shinyWidgets::materialSwitch(ns("branchcolor"), label = "Color branches by metadata",
                                       status = "info")
        ),
        accordion_panel(
          "Cluster exploration",
          icon = bs_icon("grid"),
          htmltools::p(
            "Available after running clustering.",
            style = "font-size: 0.9em; color: grey;"
          ),
          selectInput(ns("clust_nb"), label = "Cluster number:", choices = 1),
          numericRangeInput(ns("clstr_minAb"), "Taxa raw abundance range:",
                            c(1, 1), width = NULL, separator = " to "),
          selectInput(ns("clstr_rank_glom"), label = "Rank to merge taxonomy:",
                      choices = '', selected = 1)
        )
      ),
      tags$hr(),
      actionButton(
        ns("launch_clust"),
        "Run Clustering",
        icon = bs_icon("play-fill"),
        class = "btn-primary w-100 btn-lg"
      )
    ),
    navset_card_underline(
      title = "Results",
      full_screen = TRUE,
      nav_panel(
        "Dendrogram",
        icon = bs_icon("diagram-3"),
        verbatimTextOutput(ns("nb_clstr")),
        plotOutput(ns("dendro.plot"), height = "700px")
      ),
      nav_panel(
        "Samples by cluster",
        icon = bs_icon("bar-chart"),
        plotOutput(ns("sample.by.clstr"))
      ),
      nav_panel(
        "Cluster heatmap",
        icon = bs_icon("grid-3x3"),
        plotOutput(ns("subclstr_plot"), height = "600px")
      ),
      nav_panel(
        "Indicator species",
        icon = bs_icon("star"),
        verbatimTextOutput(ns("indicSpe"))
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

mod_cluster_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
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
      fact_vals <- r$sdat()[labels(dd), input$clust_fact1]
      if (input$clust_fact1 %in% names(r$numeric_palettes())) {
        # numeric metadata: derive per-leaf colors from a continuous palette
        pal_id <- r$numeric_palettes()[[input$clust_fact1]]
        ramp <- paletteer::paletteer_c(palette = pal_id, n = 256L)
        norm_vals <- (fact_vals - min(fact_vals, na.rm = TRUE)) /
                     diff(range(fact_vals, na.rm = TRUE))
        idx <- pmax(1L, pmin(256L, round(norm_vals * 255) + 1L))
        colours <- as.character(ramp)[idx]
      } else {
        colours <- r$factor_colors()[[input$clust_fact1]][fact_vals]
      }
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
          summarise(count = n(), .groups = "drop") %>%
          ggplot(aes(x=as.factor(clstr), y=count, fill=fact)) +
          geom_bar(stat="identity", position = "stack") +
          xlab("Cluster number") + ylab("counts") +
          scale_fill_manual(values=r$factor_colors()[[input$clust_fact1]])
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
        tmp <- speedyseq::tax_glom(sub.phy, input$clstr_rank_glom)
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
      data.com <- as.data.frame(as.table(as.matrix(otable)))
      colnames(data.com) <- c("Tax", "Sample", "Abundance")
      data.com$xlabel <- as.factor(r$sdat()[as.character(data.com$Sample), input$clust_fact1])
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
      flog.info('done')
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
      updateSelectInput(session, inputId = "clust_nb",
                        choices = 1:compute.k())
    })

  
    
  })
}
    
## To be copied in the UI
# mod_cluster_ui("cluster_ui_1")
    
## To be copied in the server
# mod_cluster_server("cluster_ui_1", r = r)
