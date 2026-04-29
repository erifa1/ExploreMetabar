#' mixomics UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import bslib
#' @importFrom bsicons bs_icon
mod_mixomics_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "sPLS-DA Configuration",
      open = "desktop",
      width = "350px",

      htmltools::p(
        "Sparse Partial Least Squares Discriminant Analysis for multivariate feature selection.",
        style = "font-size: 0.9em; color: grey;"
      ),

      tags$hr(),

      selectInput(ns("factor_spls_da"), label = "Factor", choices = ""),

      accordion(
        id = ns("config_accordion"),
        open = "1. Initial sPLS-DA",
        multiple = TRUE,

        accordion_panel(
          "1. Initial sPLS-DA",
          icon = bs_icon("play-circle"),
          htmltools::p("Quick exploration with manual parameters.",
                       style = "font-size: 0.85em; color: grey;"),
          numericInput(ns("nb_comp_initial"),
                       label = "Number of components",
                       min = 2, value = 2),
          uiOutput(ns("ui_nb_feat_initial")),
          div(
            style = "margin: 1rem 0;",
            actionButton(ns("launch_initial"), "Run initial sPLS-DA",
                         icon = bs_icon("play-fill"),
                         class = "btn-primary w-100")
          )
        ),

        accordion_panel(
          "2. Tune parameters",
          icon = bs_icon("gear"),
          htmltools::p("Find optimal ncomp and keepX via cross-validation (slow).",
                       style = "font-size: 0.85em; color: grey;"),
          selectInput(ns("tuning_distance"),
                      label = "Distance",
                      choices = c("max.dist", "centroids.dist", "mahalanobis.dist")),
          selectInput(ns("tuning_measure"),
                      label = "Measure",
                      choices = c("BER", "overall")),
          div(
            style = "margin: 1rem 0;",
            actionButton(ns("launch_tune"), "Tune parameters",
                         icon = bs_icon("gear-fill"),
                         class = "btn-secondary w-100")
          ),
          uiOutput(ns("ui_tune_results"))
        ),

        accordion_panel(
          "3. Final sPLS-DA",
          icon = bs_icon("check-circle"),
          htmltools::p("Run sPLS-DA with tuned parameters (overridable).",
                       style = "font-size: 0.85em; color: grey;"),
          uiOutput(ns("ui_nb_comp_final")),
          uiOutput(ns("ui_nb_feat_final")),
          div(
            style = "margin: 1rem 0;",
            actionButton(ns("launch_final"), "Run final sPLS-DA",
                         icon = bs_icon("play-fill"),
                         class = "btn-success w-100")
          )
        ),

        accordion_panel(
          "Axes & Display",
          icon = bs_icon("sliders"),
          htmltools::p("Available after running any sPLS-DA stage.",
                       style = "font-size: 0.85em; color: grey;"),
          uiOutput(ns("ui_comp_axis")),
          tags$hr(),
          tags$strong("Individuals"),
          checkboxInput(ns("plot_indiv_labels"), label = "Display sample labels", value = FALSE),
          checkboxInput(ns("plot_indiv_ellipses"), label = "Display ellipses", value = TRUE),
          tags$hr(),
          tags$strong("Variables"),
          numericInput(ns("spls_da_var_corr"),
                       label = "Correlation threshold (features below are hidden)",
                       min = 0, max = 1, value = 0, step = 0.1),
          tags$hr(),
          tags$strong("Biplot"),
          numericInput(ns("spls_da_biplot_corr"),
                       label = "Correlation threshold (features below are hidden)",
                       min = 0, max = 1, value = 0, step = 0.1),
          checkboxInput(ns("biplot_labels"), label = "Display sample labels", value = FALSE),
          checkboxInput(ns("biplot_arrows"), label = "Display arrows", value = TRUE)
        )
      )
    ),

    navset_card_underline(
      title = "sPLS-DA Results",
      full_screen = TRUE,

      nav_panel(
        "Individuals",
        icon = bs_icon("people-fill"),
        downloadButton(ns("spls_da_indiv_download"), label = "Download plot"),
        plotOutput(ns('spls_da_indiv'), width = "100%", height = "800px")
      ),

      nav_panel(
        "Variables",
        icon = bs_icon("graph-up-arrow"),
        downloadButton(ns("spls_da_var_download"), label = "Download plot"),
        plotOutput(ns('spls_da_var'), width = "100%", height = "800px")
      ),

      nav_panel(
        "Biplot",
        icon = bs_icon("arrows-angle-expand"),
        downloadButton(ns("spls_da_biplot_download"), label = "Download plot"),
        plotOutput(ns('spls_da_biplot'), width = "100%", height = "800px")
      ),

      nav_panel(
        "Loadings",
        icon = bs_icon("bar-chart-fill"),
        uiOutput(ns("ui_spls_da_loadings"))
      ),

      nav_panel(
        "Features Contribution",
        icon = bs_icon("table"),
        uiOutput(ns("ui_comp_select_var")),
        downloadButton(ns("select_var_download"), label = "Download table"),
        DT::dataTableOutput(ns("spls_da_select_var"))
      ),

      # nav_panel(
      #   "CIM",
      #   icon = bs_icon("grid-3x3-gap"),
      #   downloadButton(ns("spls_da_cim_download"), label = "Download plot"),
      #   plotOutput(ns('spls_da_cim'), width = "100%", height = "800px")
      # ),

      nav_panel(
        "Download All",
        icon = bs_icon("download"),
        downloadButton(ns("download_all"), label = "Download all results (.zip)")
      ),

      nav_panel(
        "Tuning",
        icon = bs_icon("graph-up"),
        uiOutput(ns("ui_perform_plots"))
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
#' @import svglite
#' @import zip
#' @importFrom futile.logger flog.info flog.debug
#'
mod_mixomics_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Active model state — fed by Stage 1 (initial) OR Stage 3 (final)
    active_spls_da <- reactiveVal(NULL)
    active_ncomp <- reactiveVal(NULL)
    active_keepX <- reactiveVal(NULL)

    # Tuning suggestions (used to pre-fill Stage 3 inputs)
    last_tuned <- reactiveVal(NULL)

    observe({
      req(r$phyloseq_filtered())
      updateSelectInput(session, "factor_spls_da",
                        choices = r$factor_list())
    })

    # ===== Common data reactives =====

    x <- reactive({
      req(r$phyloseq_filtered(), input$factor_spls_da)
      otu <- t(otu_table(r$phyloseq_filtered()))
      fact <- data.frame(sample_data(r$phyloseq_filtered())[, input$factor_spls_da])[, input$factor_spls_da]
      otu <- otu[!is.na(fact)]
      flog.debug('x(): factor=%s, %d samples after NA removal, %d taxa',
                 input$factor_spls_da, nrow(otu), ncol(otu))
      return(otu)
    })

    y <- reactive({
      req(r$phyloseq_filtered(), input$factor_spls_da)
      fact <- data.frame(sample_data(r$phyloseq_filtered())[, input$factor_spls_da])[, input$factor_spls_da]
      fact <- fact[!is.na(fact)]
      levels(fact) <- sort(unique(fact))
      flog.debug('y(): factor=%s, n=%d, levels=[%s]',
                 input$factor_spls_da, length(fact),
                 paste(levels(fact), collapse=', '))
      return(fact)
    })

    list_colors <- reactive({
      req(r$factor_colors(), input$factor_spls_da)
      df_color <- r$factor_colors()[input$factor_spls_da][[1]]
      df_color <- df_color[levels(y())]
      return(df_color)
    })

    # Validate factor levels on any launch
    observeEvent(c(input$launch_initial, input$launch_tune, input$launch_final), {
      req(y())
      if(1 %in% table(y())){
        levels_pb <- paste(names(which(table(y()) == 1)), collapse = ", ")
        flog.info('validation failed: single-sample levels detected: %s', levels_pb)
        shinyalert::shinyalert(
          title = "Oops",
          text = paste0("The following levels of the factor ", input$factor_spls_da,
                        " have a single associated sample.\n", levels_pb),
          type = "error"
        )
        req(FALSE)
      }
    }, ignoreInit = TRUE)

    # ===== Stage 1: Initial sPLS-DA =====

    output$ui_nb_feat_initial <- renderUI({
      req(input$nb_comp_initial, r$phyloseq_filtered())
      lapply(1:input$nb_comp_initial, FUN = function(i){
        numericInput(ns(paste("nb_feat_init", i, sep = "_")),
                     label = paste("Number of features for component", i),
                     value = phyloseq::ntaxa(r$phyloseq_filtered()),
                     min = 1,
                     max = phyloseq::ntaxa(r$phyloseq_filtered()))
      })
    })

    splsda_initial <- eventReactive(input$launch_initial, {
      withProgress({
        req(x(), y(), input$nb_comp_initial)
        list_keepX <- sapply(1:input$nb_comp_initial, FUN = function(i){
          input[[paste("nb_feat_init", i, sep = "_")]]
        })
        flog.info('splsda_initial(): ncomp=%d, keepX=[%s]',
                  input$nb_comp_initial, paste(list_keepX, collapse=', '))
        result <- mixOmics::splsda(X = x(), Y = y(),
                                   ncomp = input$nb_comp_initial,
                                   keepX = list_keepX)
        flog.info('splsda_initial() end.')
        list(model = result, ncomp = input$nb_comp_initial, keepX = list_keepX)
      }, message = "Computing initial sPLS-DA...")
    })

    observeEvent(splsda_initial(), {
      res <- splsda_initial()
      active_spls_da(res$model)
      active_ncomp(res$ncomp)
      active_keepX(res$keepX)
      flog.info('active model set from Stage 1 (initial)')
    })

    # ===== Stage 2: Tune parameters =====

    perform_spls_da <- eventReactive(input$launch_tune, {
      withProgress({
        req(x(), y())
        flog.info('perform_spls_da() starting (perf step)...')
        spls_da_init <- mixOmics::splsda(X = x(), Y = y(), ncomp = 10)
        result <- mixOmics::perf(spls_da_init, validation = "Mfold",
                                 folds = 10, nrepeat = 10,
                                 progressBar = FALSE, auc = TRUE)
        flog.info('perform_spls_da() end.')
        result
      }, message = "Tuning ncomp via perf(), please wait...")
    })

    tune_spls_da <- eventReactive(input$launch_tune, {
      withProgress({
        req(x(), y(), perform_spls_da(), input$tuning_distance, input$tuning_measure)
        flog.info('tune_spls_da() starting (tune.splsda step)...')
        ncomp <- perform_spls_da()$choice.ncomp[input$tuning_measure, input$tuning_distance]
        flog.debug('tune_spls_da(): ncomp from perf=%d (dist=%s, measure=%s)',
                   ncomp, input$tuning_distance, input$tuning_measure)
        list_keepX <- c(1:10, seq(20, 300, 10))
        tune <- mixOmics::tune.splsda(X = x(), Y = y(),
                                      ncomp = ifelse(ncomp == 1, 2, ncomp),
                                      validation = 'Mfold',
                                      folds = 5, nrepeat = 10,
                                      dist = input$tuning_distance,
                                      measure = input$tuning_measure,
                                      test.keepX = list_keepX)
        flog.info('tune_spls_da(): optimal ncomp=%d, keepX=[%s]',
                  tune$choice.ncomp$ncomp, paste(tune$choice.keepX, collapse=', '))
        return(tune)
      }, message = "Tuning keepX via tune.splsda(), please wait...")
    })

    observeEvent(tune_spls_da(), {
      tune <- tune_spls_da()
      ncomp <- ifelse(tune$choice.ncomp$ncomp == 1, 2, tune$choice.ncomp$ncomp)
      keepX <- as.integer(tune$choice.keepX[1:ncomp])
      last_tuned(list(ncomp = ncomp, keepX = keepX))
      flog.info('last_tuned set: ncomp=%d, keepX=[%s]', ncomp, paste(keepX, collapse=', '))
    })

    output$ui_tune_results <- renderUI({
      req(last_tuned())
      tagList(
        tags$hr(),
        tags$strong("Suggested parameters:"),
        tags$ul(
          tags$li(paste0("ncomp: ", last_tuned()$ncomp)),
          tags$li(paste0("keepX: ", paste(last_tuned()$keepX, collapse = ", ")))
        ),
        htmltools::p("These have been pre-filled in Stage 3 below.",
                     style = "font-size: 0.85em; color: grey;")
      )
    })

    # Tuning tab content
    output$ui_perform_plots <- renderUI({
      if(input$launch_tune == 0){
        htmltools::p(
          "Run 'Tune parameters' from the sidebar to see the performance plots.",
          style = "color: grey; padding: 2rem;"
        )
      } else {
        tagList(
          tags$h5("Performance evaluation — choosing ncomp"),
          plotOutput(ns('spls_da_ncomp'), height = "400px"),
          tags$hr(),
          tags$h5("Tuning — choosing keepX"),
          plotOutput(ns('spls_da_keepX'), height = "400px")
        )
      }
    })

    output$spls_da_ncomp <- renderPlot({
      req(perform_spls_da())
      plot(perform_spls_da(), sd = TRUE, legend.position = "horizontal")
    })

    output$spls_da_keepX <- renderPlot({
      req(tune_spls_da())
      plot(tune_spls_da())
    })

    # ===== Stage 3: Final sPLS-DA =====

    output$ui_nb_comp_final <- renderUI({
      val <- if(!is.null(last_tuned())) last_tuned()$ncomp else 2
      numericInput(ns("nb_comp_final"),
                   label = "Number of components (override allowed)",
                   min = 2, value = val)
    })

    output$ui_nb_feat_final <- renderUI({
      req(input$nb_comp_final, r$phyloseq_filtered())
      ntx <- phyloseq::ntaxa(r$phyloseq_filtered())
      defaults <- if(!is.null(last_tuned())){
        kx <- last_tuned()$keepX
        sapply(1:input$nb_comp_final, function(i){
          if(i <= length(kx)) kx[i] else ntx
        })
      } else {
        rep(ntx, input$nb_comp_final)
      }
      lapply(1:input$nb_comp_final, FUN = function(i){
        numericInput(ns(paste("nb_feat_final", i, sep = "_")),
                     label = paste("keepX for component", i, "(override allowed)"),
                     value = defaults[i],
                     min = 1,
                     max = ntx)
      })
    })

    splsda_final <- eventReactive(input$launch_final, {
      withProgress({
        req(x(), y(), input$nb_comp_final)
        list_keepX <- sapply(1:input$nb_comp_final, FUN = function(i){
          input[[paste("nb_feat_final", i, sep = "_")]]
        })
        flog.info('splsda_final(): ncomp=%d, keepX=[%s]',
                  input$nb_comp_final, paste(list_keepX, collapse=', '))
        result <- mixOmics::splsda(X = x(), Y = y(),
                                   ncomp = input$nb_comp_final,
                                   keepX = list_keepX)
        flog.info('splsda_final() end.')
        list(model = result, ncomp = input$nb_comp_final, keepX = list_keepX)
      }, message = "Computing final sPLS-DA...")
    })

    observeEvent(splsda_final(), {
      res <- splsda_final()
      active_spls_da(res$model)
      active_ncomp(res$ncomp)
      active_keepX(res$keepX)
      flog.info('active model set from Stage 3 (final)')
    })

    # ===== Common downstream UI (driven by active_*) =====

    output$ui_comp_axis <- renderUI({
      req(active_ncomp())
      tagList(
        numericInput(ns("comp_axis_1"),
                     label = "Component used on the horizontal axis",
                     min = 1, max = active_ncomp(), value = 1),
        numericInput(ns("comp_axis_2"),
                     label = "Component used on the vertical axis",
                     min = 1, max = active_ncomp(), value = 2)
      )
    })

    output$ui_comp_select_var <- renderUI({
      req(active_ncomp())
      numericInput(ns("comp_select_var"),
                   label = "Component used for feature selection",
                   min = 1, max = active_ncomp(), value = 1)
    })

    output$ui_spls_da_loadings <- renderUI({
      req(active_ncomp(), active_keepX())
      tagList(
        lapply(1:active_ncomp(), FUN = function(i){
          numericInput(ns(paste("nb_feat_load", i, sep = "_")),
                       label = paste("Number of features to display for component", i),
                       value = active_keepX()[i],
                       min = 1,
                       max = active_keepX()[i])
        }),
        downloadButton(ns("spls_da_loadings_download"), label = "Download plot"),
        lapply(1:active_ncomp(), FUN = function(i){
          plotOutput(ns(paste("spls_da_loadings", i, sep = "_")))
        })
      )
    })

    # ===== Plots — all read from active_spls_da() =====

    get_spls_da_indiv <- reactive({
      withProgress({
        req(active_spls_da(), y(), list_colors())
        comp1 <- input$comp_axis_1 %||% 1L
        comp2 <- input$comp_axis_2 %||% 2L
        flog.debug('get_spls_da_indiv(): comp=[%s, %s], labels=%s, ellipses=%s',
                   comp1, comp2, input$plot_indiv_labels, input$plot_indiv_ellipses)
        plot_ind <- mixOmics::plotIndiv(active_spls_da(), comp = c(comp1, comp2), group = y(),
                                        ind.names = input$plot_indiv_labels,
                                        ellipse = input$plot_indiv_ellipses,
                                        legend = TRUE,
                                        legend.title = input$factor_spls_da,
                                        col = list_colors())
        return(plot_ind)
      }, message = "Plot of individuals loading, please wait...")
    })

    get_spls_da_var <- reactive({
      withProgress({
        req(active_spls_da(), input$spls_da_var_corr)
        comp1 <- input$comp_axis_1 %||% 1L
        comp2 <- input$comp_axis_2 %||% 2L
        flog.debug('get_spls_da_var(): comp=[%s, %s], cutoff=%g',
                   comp1, comp2, input$spls_da_var_corr)
        plot_var <- mixOmics::plotVar(active_spls_da(), comp = c(comp1, comp2),
                                      cutoff = input$spls_da_var_corr)
        return(plot_var)
      }, message = "Plot of variables loading, please wait...")
    })

    get_spls_da_biplot <- reactive({
      withProgress({
        req(active_spls_da(), input$spls_da_biplot_corr, list_colors())
        comp1 <- input$comp_axis_1 %||% 1L
        comp2 <- input$comp_axis_2 %||% 2L
        flog.debug('get_spls_da_biplot(): comp=[%s, %s], cutoff=%g, arrows=%s, labels=%s',
                   comp1, comp2, input$spls_da_biplot_corr,
                   input$biplot_arrows, input$biplot_labels)
        var_arrow <- if(input$biplot_arrows) "black" else NULL
        plot_biplot <- biplot(active_spls_da(), comp = c(comp1, comp2),
                              var.arrow.col = var_arrow, var.names.col = "black",
                              ind.names = input$biplot_labels,
                              cutoff = input$spls_da_biplot_corr,
                              legend.title = input$factor_spls_da,
                              col = list_colors())
        return(plot_biplot)
      }, message = "Biplot loading, please wait...")
    })

    # get_spls_da_cim <- reactive({
    #   withProgress({
    #     req(active_spls_da())
    #     legend_cim <- list(legend = levels(y()), col = list_colors(),
    #                        title = input$factor_spls_da, cex = 0.7)
    #     plot_cim <- mixOmics::cim(active_spls_da())
    #     color <- r$factor_colors()[[input$factor_spls_da]][r$sdat()[plot_cim$row.names, input$factor_spls_da]]
    #     plot_cim <- mixOmics::cim(active_spls_da(), row.sideColors = color, legend = legend_cim)
    #     return(plot_cim)
    #   }, message = "CIM loading, please wait...")
    # })

    get_select_var <- reactive({
      withProgress({
        req(active_spls_da())
        comp_var <- input$comp_select_var %||% 1L
        flog.debug('get_select_var(): comp=%d', comp_var)
        select_var_tab <- mixOmics::selectVar(active_spls_da(), comp = comp_var)
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

    output$spls_da_indiv <- renderPlot(get_spls_da_indiv())
    output$spls_da_var <- renderPlot(get_spls_da_var())
    output$spls_da_biplot <- renderPlot(get_spls_da_biplot())
    # output$spls_da_cim <- renderPlot(get_spls_da_cim())

    observe({
      req(active_ncomp(), active_spls_da())
      lapply(1:active_ncomp(), FUN = function(i){
        output[[paste("spls_da_loadings", i, sep = "_")]] <- renderPlot(
          withProgress({
            mixOmics::plotLoadings(active_spls_da(), comp = i, contrib = 'max', method = 'mean',
                                   ndisplay = input[[paste("nb_feat_load", i, sep = "_")]],
                                   legend.color = list_colors(),
                                   size.name = 0.6,
                                   show.ties = FALSE,
                                   layout = c(1, 2))
          }, message = paste0("Contribution on comp", i, " loading, please wait..."))
        )
      })
    })

    # ===== Downloads — all read from active_spls_da() =====

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
        invisible(mixOmics::plotVar(active_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2), cutoff = input$spls_da_var_corr))
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
        for(i in 1:active_ncomp()){
          grDevices::svg(filename = paste0('splsda_loadings_comp', i, '.svg'), width = 10, height = 10)
          mixOmics::plotLoadings(active_spls_da(), comp = i, contrib = 'max', method = 'mean',
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

    # output$spls_da_cim_download <- downloadHandler(
    #   filename = "splsda_cim.svg",
    #   content = function(file){
    #     req(get_spls_da_cim())
    #     legend_cim <- list(legend = levels(y()), col = list_colors(), title = input$factor_spls_da, cex = 0.7)
    #     plot_cim <- mixOmics::cim(active_spls_da())
    #     color <- r$factor_colors()[[input$factor_spls_da]][r$sdat()[plot_cim$row.names, input$factor_spls_da]]
    #     grDevices::svg(filename = file, width = 10, height = 10)
    #     invisible(mixOmics::cim(active_spls_da(), row.sideColors = color, legend = legend_cim))
    #     dev.off()
    #   }
    # )

    output$download_all <- downloadHandler(
      filename = paste0('splsda_results_comp', input$comp_axis_1,'_comp', input$comp_axis_2,'.zip'),
      content = function(file){
        tmpdir <- tempdir()
        plot_indiv <- ggplot2::ggsave(file.path(tmpdir, 'splsda_indiv.svg'), plot = get_spls_da_indiv()$graph, device = 'svg', width = 10, height = 10)

        grDevices::svg(filename = file.path(tmpdir, 'splsda_var.svg'), width = 10, height = 10)
        invisible(mixOmics::plotVar(active_spls_da(), comp = c(input$comp_axis_1, input$comp_axis_2), cutoff = input$spls_da_var_corr))
        dev.off()

        biplot <- ggplot2::ggsave(file.path(tmpdir, 'splsda_biplot.svg'), plot = get_spls_da_biplot(), device = 'svg', width = 10, height = 10)

        file_list <- c(plot_indiv, file.path(tmpdir, 'splsda_var.svg'), biplot)

        t_table <- as.data.frame(phyloseq::tax_table(r$phyloseq_filtered())) %>%
            rownames_to_column()

        for(i in c(input$comp_axis_1, input$comp_axis_2)){
          grDevices::svg(filename = file.path(tmpdir, paste0('splsda_loadings_comp', i, '.svg')), width = 10, height = 10)
          mixOmics::plotLoadings(active_spls_da(), comp = i, contrib = 'max', method = 'mean',
                                 ndisplay = input[[paste("nb_feat_load", i, sep = "_")]],
                                 legend.color = list_colors(),
                                 size.name = 0.6,
                                 show.ties = FALSE)
          dev.off()
          file_list <- c(file_list, file.path(tmpdir, paste0('splsda_loadings_comp', i, '.svg')))


          select_var_tab <- mixOmics::selectVar(active_spls_da(), comp = i)
          tab_result <- select_var_tab$value %>%
            rownames_to_column() %>%
            left_join(t_table, by = "rowname")
          write.table(tab_result, file.path(tmpdir, paste0('splsda_features_contrib_comp', i, '.csv')), sep = ",", row.names = FALSE)
          file_list <- c(file_list, file.path(tmpdir, paste0('splsda_features_contrib_comp', i, '.csv')))
        }

        # legend_cim <- list(legend = levels(y()), col = list_colors(), title = input$factor_spls_da, cex = 0.7)
        # plot_cim <- mixOmics::cim(active_spls_da())
        # color <- r$factor_colors()[[input$factor_spls_da]][r$sdat()[plot_cim$row.names, input$factor_spls_da]]

        # grDevices::svg(filename = file.path(tmpdir, 'splsda_cim.svg'), width = 10, height = 10)
        # invisible(mixOmics::cim(active_spls_da(), row.sideColors = color, legend = legend_cim))
        # dev.off()

        # file_list <- c(file_list, file.path(tmpdir, 'splsda_cim.svg'))

        zip::zipr(zipfile = file, files = file_list)
      }
    )
  })
}

## To be copied in the UI
# mod_mixomics_ui("mixomics_1")

## To be copied in the server
# mod_mixomics_server("mixomics_1")
