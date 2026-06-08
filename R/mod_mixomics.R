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
        open = "1. Run sPLS-DA",
        multiple = TRUE,

        accordion_panel(
          "1. Run sPLS-DA",
          icon = bs_icon("play-circle"),
          htmltools::p("Manual parameters — all features by default, or enter your tuned values.",
                       style = "font-size: 0.85em; color: grey;"),
          numericInput(ns("nb_comp_initial"),
                       label = "Number of components",
                       min = 2, value = 2),
          uiOutput(ns("ui_nb_feat_initial")),
          div(
            style = "margin: 1rem 0;",
            input_task_button(ns("launch_initial"), "Run sPLS-DA",
                         icon = bs_icon("play-fill"),
                         class = "btn-primary w-100",
                         label_busy = "Running…")
          )
        ),

        accordion_panel(
          "2. Tune offline",
          icon = bs_icon("gear"),
          htmltools::p(
            "Finding the optimal ncomp and keepX by cross-validation is slow, so it runs in your own R session. Download the kit below, run the script, then enter the suggested values in panel 1 and re-run.",
            style = "font-size: 0.85em; color: grey;"),
          selectInput(ns("tuning_distance"),
                      label = "Distance",
                      choices = c("max.dist", "centroids.dist", "mahalanobis.dist")),
          selectInput(ns("tuning_measure"),
                      label = "Measure",
                      choices = c("BER", "overall")),
          div(
            style = "margin: 1rem 0;",
            downloadButton(ns("download_tuning_kit"),
                           label = "Download tuning kit (.zip)",
                           icon = bs_icon("download"),
                           class = "btn-secondary w-100")
          ),
          htmltools::p(
            "The kit bundles your data object, a commented R script (pre-filled with the factor / distance / measure above) and a README.",
            style = "font-size: 0.8em; color: grey;")
        ),

        accordion_panel(
          "Axes & Display",
          icon = bs_icon("sliders"),
          htmltools::p("Available after running the sPLS-DA.",
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

    # Active model state — fed by the sPLS-DA run (panel 1)
    active_spls_da <- reactiveVal(NULL)
    active_ncomp <- reactiveVal(NULL)
    active_keepX <- reactiveVal(NULL)

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
      validate(need(!input$factor_spls_da %in% names(r$numeric_palettes()),
                    "sPLS-DA requires a categorical factor."))
      df_color <- r$factor_colors()[input$factor_spls_da][[1]]
      df_color <- df_color[levels(y())]
      return(df_color)
    })

    # Validate factor levels on any launch
    observeEvent(input$launch_initial, {
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

    # ===== Step 1: Run sPLS-DA =====

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
      }, message = "Computing sPLS-DA...")
    })

    observeEvent(splsda_initial(), {
      res <- splsda_initial()
      active_spls_da(res$model)
      active_ncomp(res$ncomp)
      active_keepX(res$keepX)
      flog.info('active model set from sPLS-DA run')
    })

    # ===== Step 2: Tune offline =====

    output$download_tuning_kit <- downloadHandler(
      filename = function() paste0("splsda_tuning_kit_", input$factor_spls_da, ".zip"),
      content = function(file) {
        req(x(), y(), input$factor_spls_da)
        flog.info('download_tuning_kit(): factor=%s, dist=%s, measure=%s',
                  input$factor_spls_da, input$tuning_distance, input$tuning_measure)
        tmpdir <- tempdir()

        saveRDS(list(X = x(), Y = y(), factor = input$factor_spls_da),
                file.path(tmpdir, "splsda_input.rds"))

        tmpl <- readLines(app_sys("templates", "tune_splsda_template.R"))
        tmpl <- gsub("{{FACTOR}}",   input$factor_spls_da,  tmpl, fixed = TRUE)
        tmpl <- gsub("{{DISTANCE}}", input$tuning_distance, tmpl, fixed = TRUE)
        tmpl <- gsub("{{MEASURE}}",  input$tuning_measure,  tmpl, fixed = TRUE)
        writeLines(tmpl, file.path(tmpdir, "tune_splsda.R"))

        writeLines(c(
          "ExploreMetabar — sPLS-DA offline tuning kit",
          "",
          "This zip contains:",
          "  - splsda_input.rds : the X / Y data the app prepared for this factor",
          "  - tune_splsda.R    : a commented script that tunes ncomp and keepX",
          "  - README.txt       : this file",
          "",
          "Run tune_splsda.R in your own R session (it can be slow). When it",
          "finishes it prints the suggested 'Number of components' and 'keepX' —",
          "type those into panel 1 (Run sPLS-DA) of the app, then re-run it."
        ), file.path(tmpdir, "README.txt"))

        zip::zipr(zipfile = file,
                  files = file.path(tmpdir, c("splsda_input.rds",
                                              "tune_splsda.R",
                                              "README.txt")))
      }
    )

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
