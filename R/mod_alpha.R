# Internal helpers (pure, Shiny-free) ------------------------------------------

# Tukey box-whisker statistics for a numeric vector. Uses quantile type 7
# (plotly's "linear" method) and whiskers that extend to the most extreme value
# within 1.5*IQR of the quartiles. Returns the five values plotly needs to draw
# a box from precomputed stats (so the hover shows each value once instead of
# duplicating min/max with the fences).
box_whisker_stats <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) == 0) {
    return(list(q1 = NA_real_, median = NA_real_, q3 = NA_real_,
                lowerfence = NA_real_, upperfence = NA_real_))
  }
  q <- stats::quantile(v, probs = c(0.25, 0.5, 0.75), type = 7, names = FALSE)
  iqr <- q[3] - q[1]
  list(
    q1 = q[1], median = q[2], q3 = q[3],
    lowerfence = min(v[v >= q[1] - 1.5 * iqr]),
    upperfence = max(v[v <= q[3] + 1.5 * iqr])
  )
}

# Significance stars from a p-value (vectorised).
p_to_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "ns")))
}

# Build plotly shapes + annotations for pairwise significance brackets.
#   levels  : ordered factor levels (== x-axis category order)
#   sig_df  : Tukey rows (cols `comparison`, `p adj`) already filtered to p<0.05
#   ymax    : max data value (bracket baseline)
#   yrange  : data spread, used to size the vertical stacking step
# Comparisons are matched against level *pairs* rather than split on "-" so that
# level names containing "-" still resolve correctly. Returns shapes (one path
# per bracket), annotations (stars per bracket) and ytop (y headroom).
make_sig_brackets <- function(levels, sig_df, ymax, yrange) {
  empty <- list(shapes = list(), annotations = list(), ytop = ymax)
  if (is.null(sig_df) || nrow(sig_df) == 0) return(empty)

  step <- 0.08 * yrange
  if (!is.finite(step) || step <= 0) step <- 0.08 * abs(ymax)
  if (!is.finite(step) || step <= 0) step <- 1
  tick <- 0.35 * step

  pairs <- list()
  for (k in seq_len(nrow(sig_df))) {
    comp <- as.character(sig_df$comparison[k])
    found <- NULL
    for (i in seq_along(levels)) {
      for (j in seq_along(levels)) {
        if (i == j) next
        if (identical(comp, paste0(levels[i], "-", levels[j]))) { found <- c(i, j); break }
      }
      if (!is.null(found)) break
    }
    if (is.null(found)) next
    pairs[[length(pairs) + 1]] <- list(
      i = min(found) - 1L, j = max(found) - 1L, p = sig_df[["p adj"]][k]
    )
  }
  if (length(pairs) == 0) return(empty)

  shapes <- vector("list", length(pairs))
  annotations <- vector("list", length(pairs))
  for (k in seq_along(pairs)) {
    pr <- pairs[[k]]
    h <- ymax + k * step
    shapes[[k]] <- list(
      type = "path",
      path = sprintf("M %d,%f L %d,%f L %d,%f L %d,%f",
                     pr$i, h - tick, pr$i, h, pr$j, h, pr$j, h - tick),
      xref = "x", yref = "y", line = list(color = "black", width = 1)
    )
    annotations[[k]] <- list(
      x = (pr$i + pr$j) / 2, y = h, text = p_to_stars(pr$p),
      showarrow = FALSE, xref = "x", yref = "y",
      yanchor = "bottom", font = list(size = 14)
    )
  }
  list(shapes = shapes, annotations = annotations,
       ytop = ymax + (length(pairs) + 1) * step)
}

# Module UI

#' @title   mod_alpha_ui and mod_alpha_server
#' @description  A shiny Module.
#'
#' @param id shiny id
#' @param input internal
#' @param output internal
#' @param session internal
#'
#' @rdname mod_alpha
#'
#' @keywords internal
#' @noRd
#' @importFrom shiny NS tagList
#' @importFrom plotly plotlyOutput
#' @importFrom shinyalert shinyalert useShinyalert
#' @import bslib
#' @import bsicons
mod_alpha_ui <- function(id){
  ns <- NS(id)
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      title = "Settings",
      open = "desktop",
      width = "350px",
      htmltools::p(
        "Use the phyloseq object without performing the taxa merging step.",
        style = "font-size: 0.9em; color: grey;"
      ),
      tags$hr(),
      uiOutput(ns('ui_alpha_factor')),
      checkboxInput(ns("checkbox1"), label = "Automatic order factor", value = TRUE),
      input_task_button(ns("launch_alpha"), "Run Alpha Diversity", icon = bs_icon("play-fill"), class = "btn-primary w-100 btn-lg", label_busy = "Computing…")
    ),

    navset_card_underline(
      title = "Alpha Indexes",
      full_screen = TRUE,
      nav_panel(
        title = "Alpha indexes table",
        icon = bs_icon("table"),
        DT::dataTableOutput(ns("alphaout")),
        downloadButton(outputId = ns("alpha_download"), label = "Download Table")
      ),
      nav_panel(
        title = "Alpha indexes by group",
        icon = bs_icon("people"),
        DT::dataTableOutput(ns("alphagrp")),
        downloadButton(outputId = ns("alphagrp_download"), label = "Download group Table")
      ),
      nav_panel(
        title = "Boxplot",
        icon = bs_icon("bar-chart"),
        radioButtons(ns("metrics"), "Choose one index:", inline = TRUE,
                     choices = list("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"),
                     selected = c("Shannon")
        ),
        checkboxInput(ns("show_stats"), "Show statistics on plot (ANOVA + pairwise)", value = TRUE),
        plotly::plotlyOutput(ns("boxplot"))
      ),
      nav_panel(
        title = "Statistics",
        icon = bs_icon("calculator"),
        uiOutput(ns('anovaBox'))
      )
    )
  )
}

# Module Server

#' @rdname mod_alpha
#' @noRd
#' @keywords internal
#' @import phyloseq
#' @importFrom DT renderDataTable
#' @importFrom plotly renderPlotly config layout
#' @importFrom agricolae HSD.test
#' @importFrom gtools mixedsort
#' @import futile.logger
#' @importFrom futile.logger flog.info

mod_alpha_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
  ns <- session$ns

  factor_input <- reactive({ input$Fact1 })
  get_meta_col  <- make_get_meta_col(factor_input, r)
  local_metadata <- make_local_metadata(factor_input, get_meta_col, r)
  isNumFactor   <- make_is_num_factor(get_meta_col, local_metadata)
  local_physeq  <- make_local_physeq(local_metadata, r)


  observe({
    req(r$phyloseq_filtered(), r$var_list())
    shinyWidgets::updatePickerInput(session, "Fact1",
                      choices = r$var_list())
  })

  output$ui_alpha_factor <- renderUI({
        shinyWidgets::pickerInput(
          ns("Fact1"),
          label = "Select one or more factor to test (when multiple selection, qualitative factors are concatenated): ",
          choices = r$var_list(),
          selected = r$var_list()[1],
          multiple = FALSE
        )
  })


  alpha1 <- eventReactive(input$launch_alpha, {
    withProgress(message = 'Computing alpha diversity tables', min=0, max=10, value = 0,{
      flog.info('computing alpha1...')
      req(local_physeq())

      data <- local_physeq()
      setProgress(value = 5, detail = 'estimate richness')
      alphatab <- estimate_richness(data, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
                                                       "InvSimpson") )
      row.names(alphatab) = sample_names(data)

      LL=list()
      LL$alphatab = as.data.frame(alphatab)
      LL$data = data
      flog.info('computing alpha1 done.')
      setProgress(value = 10, detail = 'done')
      return(LL)
    })
  })


  output$alphaout <- DT::renderDataTable({
    LL = alpha1()
    round_df(LL$alphatab, 4)
  }, options = list(
    pageLength = 5,
    scrollX = TRUE
  ))


  alphagrp_table <- eventReactive(input$launch_alpha, {
    req(alpha1(), local_metadata(), get_meta_col())
    withProgress(message = 'Group table', min=0, max=10, value = 0,{
    alpha.table <- alpha1()$alphatab

    alpha.table =  tibble::rownames_to_column(alpha.table)
    alpha.table <- dplyr::left_join(local_metadata(), alpha.table, by = c( "sample.id" = "rowname"))

    alpha.table[,'rowname'] <- NULL

    if(!isNumFactor()){
      alpha.table <- alpha.table %>%
        group_by_at(get_meta_col()) %>%
        summarise(
          tibble(
            across(where(is.numeric), ~mean(.x), .names = "mean_{.col}"),
            across(where(is.numeric), ~median(.x), .names = "median_{.col}")
          )
        )
    }
    setProgress(value = 10, detail = 'done')
    return(alpha.table)
    })
  })

  output$alphagrp <- DT::renderDataTable({
    round_df(as.data.frame(alphagrp_table()), 4)
  }, options = list(
    pageLength = 5,
    scrollX = TRUE
  ))

  output$alphagrp_download <- downloadHandler(
    filename = "alphagrp_index.csv",
    content = function(file) {
      write.table(alphagrp_table(), file, sep="\t", col.names=NA)}
  )

  output$alpha_download <- downloadHandler(
    filename = "alpha_index.csv",
    content = function(file) {
      LL = alpha1()
      write.table(LL$alphatab, file, sep="\t", col.names=NA)}
  )


  boxtab <- eventReactive(input$launch_alpha, {
    req(r$sdat(), input$Fact1)
    withProgress(message = 'Boxplot table', min=0, max=10, value = 0,{
    flog.info('boxtab function...')
    LL = alpha1()
    alphatab =  tibble::rownames_to_column(LL$alphatab)

    boxtab <- dplyr::left_join(r$sdat(), alphatab, by = c('sample.id' = "rowname"))

    if(! is.numeric(boxtab[, input$Fact1])){
      if(input$checkbox1){
        boxtab[[input$Fact1]] <- factor(boxtab[[input$Fact1]], levels = gtools::mixedsort(levels(as.factor(boxtab[[input$Fact1]]))))
      }
    }

    if( !any(names(boxtab)=="sample.id") ) {
      boxtab <- dplyr::rename(boxtab, sample.id = rowname)
    }

    boxtab$Depth <- sample_sums(local_physeq())
    setProgress(value = 10, detail = 'done')
    flog.info('boxtab done.')
    return(boxtab)
    })
  }
)


  get_box_plot <- reactive({
    req(boxtab(), get_meta_col())
    flog.info('renderPloty...')
    withProgress(message = 'Rendering plot...', min=0, max=10, value = 0,{
      dt <- boxtab()
      metric <- input$metrics
      mcol <- get_meta_col()
      show_stats <- isTRUE(input$show_stats)

      # ---- numeric factor: scatter (no boxplot / no pairwise) ----
      if (is.numeric(dt[, mcol])) {
        p <- plot_ly(dt, x = as.formula(glue("~{mcol}")), y = as.formula(glue("~{metric}")),
                     color = as.formula(glue("~{mcol}")), type = 'scatter', mode = 'markers')
        title_txt <- metric
        if (show_stats) {
          sub <- tryCatch({
            fs <- alpha_lm()$fstatistic
            sprintf("Linear model: F(%d,%d) = %.2f, p = %.3g",
                    as.integer(fs[2]), as.integer(fs[3]), fs[1],
                    stats::pf(fs[1], fs[2], fs[3], lower.tail = FALSE))
          }, error = function(e) NULL)
          if (!is.null(sub)) title_txt <- paste0(metric, "<br><sub>", sub, "</sub>")
        }
        return(
          p %>% layout(title = list(text = title_txt), margin = list(t = 60),
                       yaxis = list(title = glue('{metric}'))) %>%
            config(toImageButtonOptions = list(format = "svg"))
        )
      }

      # ---- categorical factor: boxplot from precomputed quartiles ----
      # Drawing the box from precomputed stats (no raw y) makes the hover show
      # each value once, fixing the min/max == fence duplication that plotly.js
      # 2.25.2 cannot avoid for data-driven boxes.
      dt[[mcol]] <- factor(dt[[mcol]])
      lvls <- levels(dt[[mcol]])

      stats_df <- do.call(rbind, lapply(lvls, function(g) {
        s <- box_whisker_stats(dt[[metric]][dt[[mcol]] == g])
        data.frame(grp = g, q1 = s$q1, median = s$median, q3 = s$q3,
                   lowerfence = s$lowerfence, upperfence = s$upperfence,
                   stringsAsFactors = FALSE)
      }))
      stats_df$grp <- factor(stats_df$grp, levels = lvls)

      p <- plot_ly(stats_df, x = ~grp, q1 = ~q1, median = ~median, q3 = ~q3,
                   lowerfence = ~lowerfence, upperfence = ~upperfence,
                   type = "box", color = ~grp, colors = r$factor_colors()[[input$Fact1]])

      # Precomputed boxes draw no points, so overlay the outliers explicitly.
      idx <- match(as.character(dt[[mcol]]), as.character(stats_df$grp))
      out_idx <- dt[[metric]] < stats_df$lowerfence[idx] | dt[[metric]] > stats_df$upperfence[idx]
      out_idx[is.na(out_idx)] <- FALSE
      if (any(out_idx)) {
        od <- dt[out_idx, , drop = FALSE]
        p <- p %>% plotly::add_markers(
          x = as.character(od[[mcol]]), y = od[[metric]],
          marker = list(color = "black", size = 6),
          text = paste0(od$sample.id, ": ", round(od[[metric]], 3)),
          hoverinfo = "text", showlegend = FALSE, inherit = FALSE
        )
      }

      # ---- optional ANOVA subtitle + pairwise significance brackets ----
      yvals <- dt[[metric]][is.finite(dt[[metric]])]
      ymax <- max(yvals); ymin <- min(yvals); yrange <- ymax - ymin
      title_txt <- metric
      brackets <- list(shapes = list(), annotations = list(), ytop = ymax)
      if (show_stats) {
        res <- tryCatch({
          rc <- reacalpha()
          aov_df <- as.data.frame(rc$aov1[[1]])
          rn <- trimws(rownames(aov_df))
          fi <- match(mcol, rn); ri <- match("Residuals", rn)
          sub <- sprintf("ANOVA (%s): F(%d,%d) = %.2f, p = %.3g", mcol,
                         as.integer(aov_df[["Df"]][fi]), as.integer(aov_df[["Df"]][ri]),
                         aov_df[["F value"]][fi], aov_df[["Pr(>F)"]][fi])
          sig <- rc$groups1[rc$groups1[["p adj"]] < 0.05, , drop = FALSE]
          list(sub = sub, br = make_sig_brackets(lvls, sig, ymax, yrange))
        }, error = function(e) {
          flog.info(paste("alpha stats overlay skipped:", conditionMessage(e))); NULL
        })
        if (!is.null(res)) {
          title_txt <- paste0(metric, "<br><sub>", res$sub, "</sub>")
          brackets <- res$br
        }
      }

      yaxis_cfg <- list(title = glue('{metric}'))
      if (length(brackets$shapes) > 0) {
        yaxis_cfg$range <- c(ymin - 0.05 * yrange, brackets$ytop + 0.05 * yrange)
      }
      p %>% layout(
        title = list(text = title_txt), margin = list(t = 60),
        xaxis = list(categoryorder = "array", categoryarray = lvls),
        yaxis = yaxis_cfg,
        shapes = brackets$shapes, annotations = brackets$annotations
      ) %>% config(toImageButtonOptions = list(format = "svg"))
    })
  })


  output$boxplot <- renderPlotly({
    get_box_plot()
  })

  output$alphalrm <- renderPrint(
    print(alpha_lm())
  )

  alpha_lm <- reactive({
    dt <- boxtab()
    form1 = glue::glue("{input$metrics} ~ Depth + {get_meta_col()}")
    fit = lm(as.formula(form1), data=dt)
    return(summary(fit))
  })

  output$anovaBox <- renderUI({
    if(isNumFactor()){
      tagList(
        h3("Linear regression model"),
        verbatimTextOutput(ns("alphalrm"))
      )
    } else {
      tagList(
        h3("ANOVA Results"),
        DT::dataTableOutput(ns("testalpha")),
        tags$hr(),
        h3("Post-hoc Tukey HSD Test"),
        downloadButton(outputId = ns("boxtab_download"), label = "Download Tukey Table"),
        DT::dataTableOutput(ns("boxstats"))
      )
    }
  })


  reacalpha <- reactive({
    req(input$metrics, get_meta_col(), boxtab())

    flog.info('Alpha tests...')
    withProgress(message = 'Statistics...', min=0, max=10, value = 0,{

    form1 = glue::glue("{input$metrics} ~ Depth + {get_meta_col()}")
    anova_res1 <- aov( as.formula(form1), boxtab())

    tukey_hsd <- TukeyHSD(anova_res1, get_meta_col())

    LL = list()
    LL$form1 = form1
    LL$aov1 = summary(anova_res1)
    LL$groups1 <- tukey_hsd[[get_meta_col()]]

    LL$groups1 <- LL$groups1 %>% as.data.frame() %>% rownames_to_column('comparison')

    setProgress(value = 10, detail = 'done')

    })
    flog.info('Done...')
    return(LL)
 })


 output$testalpha <- DT::renderDataTable({
   req(input$metrics)
   tt <- reacalpha()
   anova_df <- as.data.frame(tt$aov1[[1]])
   anova_df <- tibble::rownames_to_column(anova_df, "Term")
   datatable(round_df(anova_df, 4), options = list(
     pageLength = 10,
     scrollX = TRUE,
     dom = 't'
   ))
 }, server = FALSE)


 output$boxstats <- DT::renderDataTable({
   req(reacalpha())
   LL = reacalpha()
   round_df(as.data.frame(LL$groups1), 4)
 }, filter="top", options = list(
   pageLength = 5,
   scrollX = TRUE
 ))


 output$boxtab_download <- downloadHandler(
   filename = "alpha_boxplot_stats.csv",
   content = function(file) {
     req(reacalpha())
     LL = reacalpha()
     write.table(LL$groups1, file, sep="\t", col.names=NA)}
 )
  })
}

## To be copied in the UI
# mod_alpha_ui("alpha_ui_1")

## To be copied in the server
# mod_alpha_server("alpha_ui_1", r = r)