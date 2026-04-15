#' mixomics UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_mixomics_ui <- function(id){
  ns <- NS(id)
  tagList(
    card(
      card_header("sPLS-DA settings"),
      selectInput(ns("factor_spls_da"), label = "Factor", choices = ""),
      radioButtons(ns("spls_da_type"), label = "sPLS-DA type", inline = TRUE,
                   choices = list("initial" = "initial", "optimised (can take a long time)" = "optimised"),
                   selected = "initial"),
      uiOutput(ns("ui_nb_comp")),
      uiOutput(ns("ui_nb_feat")),
      uiOutput(ns("ui_optimised_dist")),
      uiOutput(ns("ui_optimised_measure")),
      actionButton(ns("launch_spls_da"), "Launch sPLS-DA", icon = icon("play-circle"),
                   style="color: #fff; background-color: #3b9ef5; border-color: #1a4469")
    ),

    uiOutput(ns("ui_comp_axis")),
    uiOutput(ns("ui_perform_plots")),

    card(
      full_screen = TRUE,
      card_header("Plot of individuals"),
      layout_columns(
        col_widths = c(6, 6),
        checkboxInput(ns("plot_indiv_labels"), label = "Display sample labels", value = FALSE),
        checkboxInput(ns("plot_indiv_ellipses"), label = "Display ellipses", value = TRUE)
      ),
      downloadButton(ns("spls_da_indiv_download"), label = "Download plot"),
      plotOutput(ns('spls_da_indiv'), width = "100%", height = "800px")
    ),

    card(
      full_screen = TRUE,
      card_header("Plot of variables"),
      numericInput(ns("spls_da_var_corr"), label = "Features with correlations below this threshold will not be plotted",
                   min = 0, max = 1, value = 0, step = 0.1),
      downloadButton(ns("spls_da_var_download"), label = "Download plot"),
      plotOutput(ns('spls_da_var'), width = "100%", height = "800px")
    ),

    card(
      full_screen = TRUE,
      card_header("Biplot"),
      numericInput(ns("spls_da_biplot_corr"), label = "Features with correlations below this threshold will not be plotted",
                   min = 0, max = 1, value = 0, step = 0.1),
      layout_columns(
        col_widths = c(6, 6),
        checkboxInput(ns("biplot_labels"), label = "Display sample labels", value = FALSE),
        checkboxInput(ns("biplot_arrows"), label = "Display arrows", value = TRUE)
      ),
      downloadButton(ns("spls_da_biplot_download"), label = "Download plot"),
      plotOutput(ns('spls_da_biplot'), width = "100%", height = "800px")
    ),

    uiOutput(ns("ui_spls_da_loadings")),

    card(
      full_screen = TRUE,
      card_header("Features contribution"),
      uiOutput(ns("ui_comp_select_var")),
      downloadButton(ns("select_var_download"), label = "Download table"),
      DT::dataTableOutput(ns("spls_da_select_var"))
    ),

    card(
      full_screen = TRUE,
      card_header("CIM"),
      downloadButton(ns("spls_da_cim_download"), label = "Download plot"),
      plotOutput(ns('spls_da_cim'), width = "100%", height = "800px")
    ),

    card(
      card_header("Download all sPLS-DA results"),
      downloadButton(ns("download_all"), label = "Download zip")
    )
  )
}
    
#' mixomics Server Functions
#'
#' @param input,output,session,r Internal parameters.
#'
#' @noRd 
#' @import mixOmics
#' @import svglite
#' @import zip
#' 
mod_mixomics_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observe({
      req(r$phyloseq_filtered())
      updateSelectInput(session, "factor_spls_da",
                        choices = r$factor_list())
    })
    
    output$ui_nb_comp <- renderUI({
      req(input$spls_da_type)
      if(input$spls_da_type == "initial"){
        numericInput(ns("nb_comp"),
                     label = "Number of components",
                     min = 2,
                     value = 2
        )
      }
    })
    
    output$ui_nb_feat <- renderUI({
      req(input$spls_da_type, input$nb_comp)
      if(input$spls_da_type == "initial"){
        lapply(1:input$nb_comp, FUN = function(i){
        numericInput(ns(paste("nb_feat", i, sep = "_")),
                   label = paste("Number of features for component", i),
                   value = phyloseq::ntaxa(r$phyloseq_filtered()),
                   min = 1,
                   max = phyloseq::ntaxa(r$phyloseq_filtered()))
        })
      }
    })
    
    # measure and distance used for tuning
    output$ui_optimised_dist <- renderUI({
      req(input$spls_da_type)
      if(input$spls_da_type == "optimised"){
        selectInput(
          ns("tuning_distance"),
          label = "Distance for tuned sPLA-DA",
          choices = c("max.dist", "centroids.dist", "mahalanobis.dist")
        )
      }
    })
    
    output$ui_optimised_measure <- renderUI({
      req(input$spls_da_type)
      if(input$spls_da_type == "optimised"){
        selectInput(
          ns("tuning_measure"),
          label = "Measure for tuned sPLS-DA",
          choices = c("BER", "overall")
        )
      }
    })

    observeEvent(input$launch_spls_da, {
      output$ui_perform_plots <- renderUI({
        req(input$spls_da_type)
        if(input$spls_da_type == "optimised"){
          card(full_screen = TRUE, card_header("sPLS-DA performance evaluation plots"),
              plotOutput(ns('spls_da_ncomp')),
              plotOutput(ns('spls_da_keepX'))
          )
        }
      })
    })
    
    output$ui_comp_axis <- renderUI({
      req(final_ncomp())
      card(card_header("Choose components"),
          numericInput(ns("comp_axis_1"),
                       label = "Component used on the horizontal axis",
                       min = 1,
                       max = final_ncomp(),
                       value = 1
          ),
          numericInput(ns("comp_axis_2"),
                       label = "Component used on the vertical axis",
                       min = 1,
                       max = final_ncomp(),
                       value = 2
          )
      )
    })
    
    output$ui_comp_select_var <- renderUI({
      req(final_ncomp())
      numericInput(ns("comp_select_var"),
                   label = "Component used for feature selection",
                   min = 1,
                   max = final_ncomp(),
                   value = 1)
    })
    
    output$ui_spls_da_loadings <- renderUI({
      req(final_ncomp(), final_keepX())
      card(card_header("Loadings"),
          lapply(1:final_ncomp(), FUN = function(i){
            numericInput(ns(paste("nb_feat_load", i, sep = "_")),
                         label = paste("Number of features to display for component", i),
                         value = final_keepX()[i],
                         min = 1,
                         max = final_keepX()[i])
          }),
          downloadButton(ns("spls_da_loadings_download"), label = "Download plot"),
          lapply(1:final_ncomp(), FUN = function(i){
            plotOutput(ns(paste("spls_da_loadings", i, sep = "_")))
          })
      )
    })
    
    x <- eventReactive(input$launch_spls_da, {
      req(r$phyloseq_filtered(), input$factor_spls_da)
      otu <- t(otu_table(r$phyloseq_filtered()))
      fact <- data.frame(sample_data(r$phyloseq_filtered())[, input$factor_spls_da])[, input$factor_spls_da]
      otu <- otu[!is.na(fact)]
      return(otu)
    })
    
    y <- reactive({
      req(r$phyloseq_filtered(), input$factor_spls_da)
      if(r$tabs$tabselected == "tab_mixomics"){
        fact <- data.frame(sample_data(r$phyloseq_filtered())[, input$factor_spls_da])[, input$factor_spls_da]
        fact <- fact[!is.na(fact)]
        levels(fact) <- sort(unique(fact))
        return(fact)
      }
    })
    
    observeEvent(input$launch_spls_da, {
      req(r$tabs$tabselected, y())
      if(r$tabs$tabselected == "tab_mixomics"){
        if(1 %in% table(y())){
          levels_pb <- paste(names(which(table(y()) == 1)), collapse = ", ")
          shinyalert::shinyalert(title = "Oops", text = paste0("The following levels of the factor ", input$factor_spls_da, " have a single associated sample.\n", levels_pb), type = "error")
          req(FALSE)
        }
      }
    })
    
    spls_da <- eventReactive(input$launch_spls_da, {
      withProgress({
        req(x(), y(), input$spls_da_type)
        if(input$spls_da_type == "initial"){
          list_keepX <- sapply(1:input$nb_comp, FUN = function(i){
            input[[paste("nb_feat", i, sep = "_")]]
          })
          return(mixOmics::splsda(X = x(), Y = y(), ncomp = input$nb_comp, keepX = list_keepX))
        }
      }, message = "Computing sPLS-DA...")
    })
    
    # undergo performance evaluation in order to tune the number of components to use
    perform_spls_da <- eventReactive(input$launch_spls_da, {
      withProgress({
        req(req(x(), y()))
        if(input$spls_da_type == "optimised"){
          spls_da_init <- mixOmics::splsda(X = x(), Y = y(), ncomp = 10)
          mixOmics::perf(spls_da_init, validation = "Mfold", 
                         folds = 10, nrepeat = 10, # use repeated cross-validation
                         progressBar = FALSE, auc = TRUE)
        }
      }, message = "Tuning sPLS-DA, please wait...")
    })
    
    tune_spls_da <- eventReactive(input$launch_spls_da, {
      withProgress({
        req(x(), y(), perform_spls_da(), input$tuning_distance, input$tuning_measure)
        if(input$spls_da_type == "optimised"){
          ncomp <- perform_spls_da()$choice.ncomp[input$tuning_measure, input$tuning_distance]
          list_keepX <- c(1:10,  seq(20, 300, 10))
          tune <- tune.splsda(X = x(), Y = y(),
                              ncomp = ifelse(ncomp == 1, 2, ncomp), # calculate for first components
                              validation = 'Mfold',
                              folds = 5, nrepeat = 10, # use repeated cross-validation
                              dist =  input$tuning_distance,
                              measure =  input$tuning_measure,
                              test.keepX = list_keepX,
                              cpus = 2) # allow for parallelisation to decrease runtime
          return(tune)
        }
      }, message = "Computing optimised sPLS-DA, please wait...")
    })
    
    spls_da_optimal <- eventReactive(input$launch_spls_da, {
      req(x(), y(), tune_spls_da())
      spls_da_opt <- mixOmics::splsda(X = x(), Y = y(),
                                      ncomp = final_ncomp(),
                                      keepX = tune_spls_da()$choice.keepX[1:final_ncomp()]
      )
      return(spls_da_opt)
    })
    
    final_ncomp <- eventReactive(input$launch_spls_da, {
      req(input$spls_da_type)
      if(input$spls_da_type == "initial"){
        ncomp <- input$nb_comp
      }else if(input$spls_da_type == "optimised"){
        ncomp <- ifelse(tune_spls_da()$choice.ncomp$ncomp == 1, 2, tune_spls_da()$choice.ncomp$ncomp)
      }
      return(ncomp)
    })
    
    final_keepX <- eventReactive(input$launch_spls_da, {
      req(input$spls_da_type)
      if(input$spls_da_type == "initial"){
        list_keepX <- sapply(1:input$nb_comp, FUN = function(i){
          input[[paste("nb_feat", i, sep = "_")]]
        })
      }else if(input$spls_da_type == "optimised"){
        list_keepX <- tune_spls_da()$choice.keepX[1:final_ncomp()]
      }
      return(list_keepX)
    })
    
    final_spls_da <- eventReactive(input$launch_spls_da, {
      if(input$spls_da_type == "initial"){
        spls_da <- spls_da()
      }else if(input$spls_da_type == "optimised"){
        spls_da <- spls_da_optimal()
      }
      return(spls_da)
    })
    
    list_colors <- reactive({
      req(r$factor_colors(), input$factor_spls_da)
      df_color <- r$factor_colors()[input$factor_spls_da][[1]]
      df_color <- df_color[levels(y())]
      return(df_color)
    })
    
    observeEvent(input$launch_spls_da, {
      output$spls_da_ncomp <- renderPlot(
        get_spls_da_ncomp()
      )
      
      output$spls_da_keepX <- renderPlot(
        get_spls_da_keepX()
      )
    })
    
    
    get_spls_da_ncomp <- eventReactive(input$launch_spls_da, {
      req(perform_spls_da())
      pl <- plot(perform_spls_da(), sd = TRUE, legend.position = "horizontal")
      return(pl)
    })
    
    get_spls_da_keepX <- eventReactive(input$launch_spls_da, {
      req(tune_spls_da())
      pl <- plot(tune_spls_da())
      return(pl)
    })
    
    get_spls_da_indiv <- reactive({
      withProgress({
        req(final_spls_da(), y(), list_colors())
        plot_ind <- mixOmics::plotIndiv(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2), group = y(),
                                        ind.names = input$plot_indiv_labels, ellipse = input$plot_indiv_ellipses, legend = TRUE,
                                        legend.title = input$factor_spls_da, col.per.group = list_colors()
        )
        return(plot_ind)
      }, message = "Plot of individuals loading, please wait...")
    })
    
    get_spls_da_var <- reactive({
      withProgress({
        req(final_spls_da(), input$spls_da_var_corr)
        plot_var <- mixOmics::plotVar(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2),
                                      cutoff = input$spls_da_var_corr)
        return(plot_var)
      }, message = "Plot of variables loading, please wait...")
    })
    
    get_spls_da_biplot <- reactive({
      withProgress({
        req(final_spls_da(), input$spls_da_biplot_corr, list_colors())
        if(input$biplot_arrows){
          var_arrow <- "black"
        }else{
          var_arrow <- NULL
        }
        plot_biplot <- biplot(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2),
                              var.arrow.col = var_arrow, var.names.col = "black", ind.names = input$biplot_labels,
                              cutoff = input$spls_da_biplot_corr, legend.title = input$factor_spls_da,
                              col.per.group = list_colors())
        return(plot_biplot)
      }, message = "Biplot loading, please wait...")
    })
    
    get_spls_da_cim <- reactive({
      withProgress({
        req(final_spls_da())
        legend_cim <- list(legend = levels(y()), col = list_colors(), title = input$factor_spls_da, cex = 0.7)
        plot_cim <- mixOmics::cim(final_spls_da())
        color <- r$factor_colors()[[input$factor_spls_da]][r$sdat()[plot_cim$row.names, input$factor_spls_da]]
        plot_cim <- mixOmics::cim(final_spls_da(), row.sideColors = color, legend = legend_cim)
        return(plot_cim)
      }, message = "CIM loading, please wait...")
      
    })
    
    get_select_var <- reactive({
      withProgress({
        req(final_spls_da(), input$comp_select_var)
        select_var_tab <- mixOmics::selectVar(final_spls_da(), comp = input$comp_select_var)
        t_table <- as.data.frame(phyloseq::tax_table(r$phyloseq_filtered())) %>%
          rownames_to_column()
        tab_result <- select_var_tab$value %>%
          rownames_to_column() %>%
          left_join(t_table, by = "rowname")
        return(tab_result)
      }, message = "Features contribution loading, please wait...")
    })
    
    output$spls_da_select_var <- DT::renderDataTable({
      req(get_select_var())
      get_select_var()
    }, filter = "top", options = list(scrollX = TRUE))
    
    output$spls_da_indiv <- renderPlot(
      get_spls_da_indiv()
    )
    
    output$spls_da_var <- renderPlot(
      get_spls_da_var()
    )
    
    output$spls_da_biplot <- renderPlot(
      get_spls_da_biplot()
    )
    
    output$spls_da_cim <- renderPlot(
      get_spls_da_cim()
    )
    
    observeEvent(input$launch_spls_da, {
      req(final_ncomp())
      lapply(1:final_ncomp(), FUN = function(i){
        output[[paste("spls_da_loadings", i, sep = "_")]] <- renderPlot(
          withProgress({
            mixOmics::plotLoadings(final_spls_da(), comp = i, contrib = 'max', method = 'mean',
                                 ndisplay = input[[paste("nb_feat_load", i, sep = "_")]], # nb of features to display
                                 legend.color = list_colors(),
                                 size.name = 0.6,
                                 show.ties = FALSE,
                                 layout = c(1, 2))
          }, message = paste0("Contribution on comp", i, " loading, please wait..."))
        )
      })
    })

    output$spls_da_indiv_download <- downloadHandler(
      filename = "splsda_indiv.svg",
      content = function(file){
        grDevices::svg(filename = file, width = 10, height = 10)
        print(get_spls_da_indiv()$graph)
        dev.off()
      }
    )
    
    output$spls_da_var_download <- downloadHandler(
      filename = "splsda_var.svg",
      content = function(file){
        grDevices::svg(filename = file, width = 10, height = 10)
        invisible(mixOmics::plotVar(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2), cutoff = input$spls_da_var_corr))
        dev.off()
      }
    )
    
    output$spls_da_biplot_download <- downloadHandler(
      filename = "splsda_biplot.svg",
      content = function(file){
        req(get_spls_da_biplot())
        grDevices::svg(filename = file, width = 10, height = 10)
        print(get_spls_da_biplot())
        dev.off()
      }
    )
    
    output$spls_da_loadings_download <- downloadHandler(
      filename = "splsda_loadings.zip",
      content = function(file){
        file_list <- c()
        for(i in 1:final_ncomp()){
          grDevices::svg(filename = paste0('splsda_loadings_comp', i, '.svg'), width = 10, height = 10)
          mixOmics::plotLoadings(final_spls_da(), comp = i, contrib = 'max', method = 'mean',
                                 ndisplay = input[[paste("nb_feat_load", i, sep = "_")]],
                                 legend.color = list_colors(),
                                 size.name = 0.6,
                                 show.ties = FALSE)
          dev.off()
          file_list <- c(file_list, paste0('splsda_loadings_comp', i, '.svg'))
        }
        zip::zip(zipfile = file, files = file_list)
      }
    )
    
    output$select_var_download <- downloadHandler(
      filename = "features_contribution.csv",
      content = function(file){
        req(get_select_var())
        write.table(get_select_var(), file, sep = ",", row.names = FALSE)
      }
    )
    
    output$spls_da_cim_download <- downloadHandler(
      filename = "splsda_cim.svg",
      content = function(file){
        req(get_spls_da_cim())
        legend_cim <- list(legend = levels(y()), col = list_colors(), title = input$factor_spls_da, cex = 0.7)
        plot_cim <- mixOmics::cim(final_spls_da())
        color <- r$factor_colors()[[input$factor_spls_da]][r$sdat()[plot_cim$row.names, input$factor_spls_da]]
        grDevices::svg(filename = file, width = 10, height = 10)
        invisible(mixOmics::cim(final_spls_da(), row.sideColors = color, legend = legend_cim))
        dev.off()
      }
    )
    
    output$download_all <- downloadHandler(
      filename = paste0('splsda_results_comp', input$comp_axis_1,'_comp', input$comp_axis_2,'.zip'),
      content = function(file){
        tmpdir <- tempdir()
        plot_indiv <- ggplot2::ggsave(file.path(tmpdir, 'splsda_indiv.svg'), plot = get_spls_da_indiv()$graph, device = 'svg', width = 10, height = 10)

        grDevices::svg(filename = file.path(tmpdir, 'splsda_var.svg'), width = 10, height = 10)
        invisible(mixOmics::plotVar(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2), cutoff = input$spls_da_var_corr))
        dev.off()
        
        biplot <- ggplot2::ggsave(file.path(tmpdir, 'splsda_biplot.svg'), plot = get_spls_da_biplot(), device = 'svg', width = 10, height = 10)
        
        file_list <- c(plot_indiv, file.path(tmpdir, 'splsda_var.svg'), biplot)
        
        t_table <- as.data.frame(phyloseq::tax_table(r$phyloseq_filtered())) %>%
            rownames_to_column()
        
        for(i in c(input$comp_axis_1, input$comp_axis_2)){
          grDevices::svg(filename = file.path(tmpdir, paste0('splsda_loadings_comp', i, '.svg')), width = 10, height = 10)
          mixOmics::plotLoadings(final_spls_da(), comp = i, contrib = 'max', method = 'mean',
                                 ndisplay = input[[paste("nb_feat_load", i, sep = "_")]],
                                 legend.color = list_colors(),
                                 size.name = 0.6,
                                 show.ties = FALSE)
          dev.off()
          file_list <- c(file_list, file.path(tmpdir, paste0('splsda_loadings_comp', i, '.svg')))
          
          
          select_var_tab <- mixOmics::selectVar(final_spls_da(), comp = i)
          tab_result <- select_var_tab$value %>%
            rownames_to_column() %>%
            left_join(t_table, by = "rowname")
          write.table(tab_result, file.path(tmpdir, paste0('splsda_features_contrib_comp', i, '.csv')), sep = ",", row.names = FALSE)
          file_list <- c(file_list, file.path(tmpdir, paste0('splsda_features_contrib_comp', i, '.csv')))
        }
        
        legend_cim <- list(legend = levels(y()), col = list_colors(), title = input$factor_spls_da, cex = 0.7)
        plot_cim <- mixOmics::cim(final_spls_da())
        color <- r$factor_colors()[[input$factor_spls_da]][r$sdat()[plot_cim$row.names, input$factor_spls_da]]
        
        grDevices::svg(filename = file.path(tmpdir, 'splsda_cim.svg'), width = 10, height = 10)
        invisible(mixOmics::cim(final_spls_da(), row.sideColors = color, legend = legend_cim))
        dev.off()
        
        file_list <- c(file_list, file.path(tmpdir, 'splsda_cim.svg'))
        
        zip::zipr(zipfile = file, files = file_list)
      }
    )
  })
}
    
## To be copied in the UI
# mod_mixomics_ui("mixomics_1")
    
## To be copied in the server
# mod_mixomics_server("mixomics_1")
