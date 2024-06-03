#' mixomics UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' @import mixOmics
mod_mixomics_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        box(title = "sPLS-DA settings", width = 12, status = "warning", solidHeader = TRUE,
            selectInput(
              ns("factor_spls_da"),
              label = "Factor",
              choices = ""
            ),
            tags$h5("Optimised sPLS-DA can take time."),
            radioButtons(ns("spls_da_type"),
                         label = "sPLS-DA type",
                         inline = TRUE,
                         choices = c("initial", "optimised"),
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
        box(id = ns("plot_indiv_box"), title = "Plot of individuals", width = 12, status = "primary", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            downloadButton(ns("spls_da_indiv_download"), label = "Download plot"),
            shinycustomloader::withLoader(
              plotOutput(ns('spls_da_indiv'), width = "1000px", height = "1000px"),
              type = "html", loader = "loader4"
            )
        ),
        box(id = ns("plot_var_box"), title = "Plot of variables", width = 12, status = "primary", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            numericInput(ns("spls_da_var_corr"),
                         label = "Features with correlations below this threshold will not be plotted",
                         min = 0,
                         max = 1,
                         value = 0,
                         step = 0.1
            ),
            downloadButton(ns("spls_da_var_download"), label = "Download plot"),
            shinycustomloader::withLoader(
              plotOutput(ns('spls_da_var'), width = "1000px", height = "1000px"),
              type = "html", loader = "loader4"
            )
        ),
        box(title = "Biplot", width = 12, status = "primary", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            numericInput(ns("spls_da_biplot_corr"),
                         label = "Features with correlations below this threshold will not be plotted",
                         min = 0,
                         max = 1,
                         value = 0,
                         step = 0.1
            ),
            downloadButton(ns("spls_da_biplot_download"), label = "Download plot"),
            shinycustomloader::withLoader(
              plotOutput(ns('spls_da_biplot'), width = "1000px", height = "1000px"),
              type = "html", loader = "loader4"
            )
        ),
        uiOutput(ns("ui_spls_da_loadings")),
        box(title = "Selected features", width = 12, status = "primary", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            uiOutput(ns("ui_comp_select_var")),
            shinycustomloader::withLoader(
              DT::dataTableOutput(ns("spls_da_select_var")),
              type = "html", loader = "loader4"
            )
        ),
        box(title = "CIM", width = 12, status = "primary", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            downloadButton(ns("spls_da_cim_download"), label = "Download plot"),
            shinycustomloader::withLoader(
              plotOutput(ns('spls_da_cim')),
              type = "html", loader = "loader4"
            )
        )
      )
    )
  )
}
    
#' mixomics Server Functions
#'
#' @param input,output,session,r Internal parameters.
#'
#' @noRd 
#' @import mixOmics
mod_mixomics_server <- function(input, output, session, r){
    ns <- session$ns

    observe({
      req(r$phyloseq_filtered())
      updateSelectInput(session, "factor_spls_da",
                        choices = r$var_list())
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
                   value = dim(x())[2],
                   min = 1,
                   max = dim(x())[2])
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

    output$ui_perform_plots <- renderUI({
      req(input$spls_da_type)
      if(input$spls_da_type == "optimised"){
        box(title = "sPLS-DA performance evaluation plots", width = 12, status = "primary", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
            shinycustomloader::withLoader(
              plotOutput(ns('spls_da_ncomp')),
              type = "html", loader = "loader4"
            ),
            shinycustomloader::withLoader(
              plotOutput(ns('spls_da_keepX')),
              type = "html", loader = "loader4"
            )
        )
      }
    })
    
    
    output$ui_comp_axis <- renderUI({
      req(final_ncomp())
      box(title = "Choose components", width = 12, status = "warning", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
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
      box(title = "Loadings", width = 12, status = "primary", collapsible = TRUE, collapsed = TRUE, solidHeader = TRUE,
          lapply(1:final_ncomp(), FUN = function(i){
            numericInput(ns(paste("nb_feat_load", i, sep = "_")),
                         label = paste("Number of features to display for component", i),
                         value = final_keepX()[i],
                         min = 1,
                         max = final_keepX()[i])
          }),
          downloadButton(ns("spls_da_loadings_download"), label = "Download plot"),
          lapply(1:final_ncomp(), FUN = function(i){
            shinycustomloader::withLoader(
              plotOutput(ns(paste("spls_da_loadings", i, sep = "_"))),
              type = "html", loader = "loader4"
            )
          })
      )
    })
    
    x <- reactive({
      req(r$phyloseq_filtered())
      return(t(otu_table(r$phyloseq_filtered())))
    })
    
    y <- eventReactive(input$launch_spls_da,{
      fact <- data.frame(sample_data(r$phyloseq_filtered())[, input$factor_spls_da])[, input$factor_spls_da]
      return(fact)
    })
    
    spls_da <- eventReactive(input$launch_spls_da, {
      req(x(), y(), input$spls_da_type)
      if(input$spls_da_type == "initial"){
        list_keepX <- sapply(1:input$nb_comp, FUN = function(i){
          input[[paste("nb_feat", i, sep = "_")]]
        })
        return(mixOmics::splsda(X = x(), Y = y(), ncomp = input$nb_comp, keepX = list_keepX))
      }
    })
    
    # undergo performance evaluation in order to tune the number of components to use
    perform_spls_da <- eventReactive(input$launch_spls_da, {
      req(req(x(), y()))
      if(input$spls_da_type == "optimised"){
      spls_da_init <- mixOmics::splsda(X = x(), Y = y(), ncomp = 10)
      mixOmics::perf(spls_da_init, validation = "Mfold", 
                     folds = 5, nrepeat = 10, # use repeated cross-validation
                     progressBar = FALSE, auc = TRUE) # include AUC values
      }
    })
    
    tune_spls_da <- reactive({
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
    })
    
    spls_da_optimal <- reactive({
      req(x(), y(), tune_spls_da())
      spls_da_opt <- mixOmics::splsda(X = x(), Y = y(),
                                      ncomp = final_ncomp(),
                                      keepX = tune_spls_da()$choice.keepX[1:final_ncomp()]
      )
      return(spls_da_opt)
    })
    
    final_ncomp <- reactive({
      req(input$spls_da_type)
      if(input$spls_da_type == "initial"){
        ncomp <- input$nb_comp
      }else if(input$spls_da_type == "optimised"){
        ncomp <- ifelse(tune_spls_da()$choice.ncomp$ncomp == 1, 2, tune_spls_da()$choice.ncomp$ncomp)
      }
      return(ncomp)
    })
    
    final_keepX <- reactive({
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
    
    final_spls_da <- reactive({
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
    
    # modality <- reactive({
    #   req(input$factor_spls_da)
    #   modal <- get_cat_colors("fact", input$factor_spls_da, r$phyloseq_filtered())
    #   return(modal)
    # })
    
    output$spls_da_ncomp <- renderPlot(
      get_spls_da_ncomp()
    )
    
    output$spls_da_keepX <- renderPlot(
      get_spls_da_keepX()
    )
    
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
      req(final_spls_da(), y(), list_colors())
      plot_ind <- mixOmics::plotIndiv(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2), group = y(),
                          ind.names = FALSE, ellipse = TRUE, legend = TRUE,
                          col.per.group = list_colors()
      )
      return(plot_ind)
    })
    
    get_spls_da_var <- reactive({
      req(final_spls_da(), input$spls_da_var_corr)
      plot_var <- mixOmics::plotVar(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2),
                                    cutoff = input$spls_da_var_corr)
      return(plot_var)
    })
    
    get_spls_da_biplot <- reactive({
      req(final_spls_da(), input$spls_da_biplot_corr, list_colors())
      plot_biplot <- biplot(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2),
                            cutoff = input$spls_da_biplot_corr, legend.title = "Legend",
                            col.per.group = list_colors())
      return(plot_biplot)
    })
    
    get_spls_da_cim <- reactive({
      req(final_spls_da())
      legend_cim <- list(legend = levels(y()), col = list_colors(), title = input$factor_spls_da, cex = 0.7)
      plot_cim <- mixOmics::cim(final_spls_da())
      color <- r$factor_colors()[[input$factor_spls_da]][r$sdat()[plot_cim$row.names, input$factor_spls_da]]
      plot_cim <- mixOmics::cim(final_spls_da(), row.sideColors = color, legend = legend_cim)
      return(plot_cim)
    })
    
    get_select_var <- reactive({
      req(final_spls_da(), input$comp_select_var)
      select_var_tab <- mixOmics::selectVar(final_spls_da(), comp = input$comp_select_var)
      return(select_var_tab)
    })
    
    output$spls_da_select_var <- DT::renderDataTable({
      req(get_select_var())
      t_table <- as.data.frame(phyloseq::tax_table(r$phyloseq_filtered())) %>%
        rownames_to_column()
      tab_result <- get_select_var()$value %>%
        rownames_to_column() %>%
        left_join(t_table, by = "rowname")
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
    
    observe({
      req(final_ncomp())
      lapply(1:final_ncomp(), FUN = function(i){
        output[[paste("spls_da_loadings", i, sep = "_")]] <- renderPlot(
          mixOmics::plotLoadings(final_spls_da(), comp = i, contrib = 'max', method = 'mean',
                                 ndisplay = input[[paste("nb_feat_load", i, sep = "_")]], # nb of features to display
                                 legend.color = list_colors(),
                                 size.name = 0.8,
                                 show.ties = FALSE)
        )
      })
    })

    output$spls_da_indiv_download <- downloadHandler(
      filename = "splsda_indiv.svg",
      content = function(file){
        svg(filename = file)
        print(get_spls_da_indiv()$graph)
        dev.off()
      }
    )
    
    output$spls_da_var_download <- downloadHandler(
      filename = "splsda_var.svg",
      content = function(file){
        svg(filename = file)
        invisible(mixOmics::plotVar(final_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2), cutoff = input$spls_da_var_corr))
        dev.off()
      }
    )
    
    output$spls_da_biplot_download <- downloadHandler(
      filename = "splsda_biplot.svg",
      content = function(file){
        req(get_spls_da_biplot())
        svg(filename = file)
        print(get_spls_da_biplot())
        dev.off()
      }
    )
    
    output$spls_da_loadings_download <- downloadHandler(
      filename = "splsda_loadings.svg",
      content = function(file){
        svg(filename = file)
        invisible(mixOmics::plotLoadings(final_spls_da()))
        dev.off()
      }
    )
    
    output$spls_da_cim_download <- downloadHandler(
      filename = "splsda_cim.svg",
      content = function(file){
        req(get_spls_da_cim())
        svg(filename = file)
        invisible(mixOmics::cim(final_spls_da()))
        dev.off()
      }
    )
}
    
## To be copied in the UI
# mod_mixomics_ui("mixomics_1")
    
## To be copied in the server
# mod_mixomics_server("mixomics_1")
