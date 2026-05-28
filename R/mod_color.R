#' color UI Function
#'
#' @description Module that owns the per-variable / per-rank color maps used
#' across the app. Reads from the active phyloseq object via the shared `r`
#' object and writes back `r$factor_colors()`, `r$numeric_palettes()` and
#' `r$taxa_colors()`.
#'
#' UI is rendered inside the "Colors" nav panel of [mod_data_loading_ui()].
#'
#' @param id Internal parameter for `{shiny}`.
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import bslib
#' @importFrom bsicons bs_icon
mod_color_ui <- function(id) {
  ns <- NS(id)
  tagList(
    accordion(
      open = FALSE,
      accordion_panel(
        "Color file format",
        icon = bs_icon("info-circle"),
        tagList(
          tags$p("Long-format CSV with three columns: ",
                 tags$code("variable"), ", ",
                 tags$code("modality"), ", ",
                 tags$code("color"), "."),
          tags$ul(
            tags$li(tags$b("variable"),
                    " — metadata variable name (factor file) or ",
                    "taxonomic rank name (taxa file). The special value ",
                    tags$code("__numeric__"),
                    " marks a row that assigns a continuous palette to a ",
                    "numeric metadata variable."),
            tags$li(tags$b("modality"),
                    " — the modality of a categorical variable, the taxon ",
                    "name for a rank, or (when ",
                    tags$code("variable = __numeric__"),
                    ") the name of the numeric variable."),
            tags$li(tags$b("color"),
                    " — 7-character hex code ",
                    tags$code("#RRGGBB"),
                    " for categorical / taxa rows, ",
                    tags$b("or"),
                    " a ",
                    tags$code("package::palette"),
                    " string for ",
                    tags$code("__numeric__"),
                    " rows (only ", tags$code("viridis"), " and ",
                    tags$code("ggthemes"),
                    " continuous palettes are accepted).")
          ),
          tags$p("Example:"),
          tags$pre(paste(
            "variable,modality,color",
            "treatment,control,#1f77b4",
            "treatment,case,#ff7f0e",
            "Phylum,Firmicutes,#2ca02c",
            "Phylum,Bacteroidetes,#d62728",
            "__numeric__,age,viridis::viridis",
            sep = "\n"
          )),
          tags$p(tags$b("Validation:"),
                 " rows with unknown variables, unknown modalities or ",
                 "invalid color codes are dropped silently and reported in ",
                 "the log below each file input. Missing modalities are ",
                 "auto-filled from a deterministic palette."),
          tags$p(tags$b("Round-trip:"),
                 " use ", tags$em("Download CSV file"),
                 " buttons to obtain a template seeded with the current ",
                 "color assignments, then edit and re-upload.")
        )
      )
    ),
    layout_columns(
      col_widths = c(6, 6),
      div(
        uiOutput(ns("var_color")),
        downloadButton(ns("var_color_download"),
                       label = "Download CSV file"),
        fileInput(ns("file_var_color"),
                  label = "var-color file",
                  placeholder = "var_colors.csv",
                  accept = c(".csv", ".tsv", "text/csv")),
        verbatimTextOutput(ns("fact_color_file_log"))
      ),
      div(
        uiOutput(ns("taxa_color")),
        downloadButton(ns("taxa_color_download"),
                       label = "Download CSV file"),
        fileInput(ns("file_taxa_color"),
                  label = "taxa-color file",
                  placeholder = "taxa_colors.csv",
                  accept = c(".csv", ".tsv", "text/csv")),
        verbatimTextOutput(ns("taxa_color_file_log"))
      )
    )
  )
}


#' color Server Function
#'
#' Owns `r$factor_colors()`, `r$numeric_palettes()` and `r$taxa_colors()`.
#' Pure logic lives in [color_helpers] — this server is glue.
#'
#' @noRd
#' @import phyloseq
mod_color_server <- function(id, r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- ground-truth reactive sources -----------------------------------

    rank_list <- reactive({
      req(r$phyloseq_filtered())
      phy <- r$phyloseq_filtered()
      validate(need(!is.null(phyloseq::access(phy, "tax_table")),
                    "Tax table not available yet."))
      phyloseq::rank_names(phy)
    })

    factor_modalities <- reactive({
      req(r$phyloseq_filtered(), r$factor_list())
      extract_modalities(r$phyloseq_filtered(), r$factor_list(), type = "fact")
    })

    taxa_modalities <- reactive({
      req(r$phyloseq_filtered(), rank_list())
      extract_modalities(r$phyloseq_filtered(), rank_list(), type = "taxa")
    })

    numeric_vars <- reactive({
      req(r$sdat(), r$var_list(), r$factor_list())
      setdiff(r$var_list(), r$factor_list())
    })

    # ---- automatic assignments -------------------------------------------

    auto_factor_colors <- reactive({
      req(factor_modalities())
      assign_qualitative_colors(factor_modalities(),
                                pool = qualitative_palette_pool(),
                                type = "fact")
    })

    auto_taxa_colors <- reactive({
      req(taxa_modalities())
      assign_qualitative_colors(taxa_modalities(),
                                pool = qualitative_palette_pool(),
                                type = "taxa")
    })

    auto_numeric_palettes <- reactive({
      assign_numeric_palettes(numeric_vars(),
                              pool = continuous_palette_pool())
    })

    # ---- optional CSV override -------------------------------------------

    parsed_factor_csv <- reactive({
      req(input$file_var_color, factor_modalities())
      df <- read_color_csv(input$file_var_color$datapath)
      validate(need(!is.null(df),
                    "Could not parse CSV: expected columns variable, modality, color."))
      apply_color_csv(df, factor_modalities(), accept_numeric = TRUE)
    })

    parsed_taxa_csv <- reactive({
      req(input$file_taxa_color, taxa_modalities())
      df <- read_color_csv(input$file_taxa_color$datapath)
      validate(need(!is.null(df),
                    "Could not parse CSV: expected columns variable, modality, color."))
      apply_color_csv(df, taxa_modalities(), accept_numeric = FALSE)
    })

    # ---- merged outputs (auto + override) --------------------------------

    factor_colors_merged <- reactive({
      req(auto_factor_colors())
      if (is.null(input$file_var_color)) return(auto_factor_colors())
      loaded <- parsed_factor_csv()$colors
      if (length(loaded) == 0L) return(auto_factor_colors())
      filled <- fill_missing_colors(loaded, factor_modalities())
      # variables not covered by the upload keep their auto assignment
      not_covered <- setdiff(names(auto_factor_colors()), names(filled))
      c(filled, auto_factor_colors()[not_covered])
    })

    numeric_palettes_merged <- reactive({
      req(auto_numeric_palettes())
      if (is.null(input$file_var_color)) return(auto_numeric_palettes())
      loaded <- parsed_factor_csv()$numerics
      if (length(loaded) == 0L) return(auto_numeric_palettes())
      not_covered <- setdiff(names(auto_numeric_palettes()), names(loaded))
      c(loaded, auto_numeric_palettes()[not_covered])
    })

    taxa_colors_merged <- reactive({
      req(auto_taxa_colors())
      if (is.null(input$file_taxa_color)) return(auto_taxa_colors())
      loaded <- parsed_taxa_csv()$colors
      if (length(loaded) == 0L) return(auto_taxa_colors())
      filled <- fill_missing_colors(loaded, taxa_modalities())
      not_covered <- setdiff(names(auto_taxa_colors()), names(filled))
      c(filled, auto_taxa_colors()[not_covered])
    })

    # ---- expose contract to other modules --------------------------------

    r$factor_colors    <- factor_colors_merged
    r$numeric_palettes <- numeric_palettes_merged
    r$taxa_colors      <- taxa_colors_merged

    # ---- selection widgets ------------------------------------------------

    output$var_color <- renderUI({
      req(r$var_list(), r$sdat())
      vars <- r$var_list()
      shinyWidgets::pickerInput(
        ns("list_var_color"),
        label = "Variables whose associated colors are exported in the CSV file",
        choices = vars,
        selected = vars[1],
        multiple = TRUE,
        options = shinyWidgets::pickerOptions(
          actionsBox = TRUE,
          liveSearch = TRUE,
          showContent = FALSE
        ),
        choicesOpt = list(
          content = unlist(lapply(vars, function(x) {
            htmltools::doRenderTags(
              tags$div(
                shiny::splitLayout(
                  cellWidths = 200,
                  tags$div(style = htmltools::css(fontWeight = "bold"), x),
                  tags$div(style = htmltools::css(color = "grey"),
                           class(r$sdat()[, x]))
                )
              )
            )
          }))
        )
      )
    })

    output$taxa_color <- renderUI({
      req(rank_list())
      shinyWidgets::pickerInput(
        ns("list_taxa_color"),
        label = "Taxonomic ranks whose associated colors are exported in the CSV file",
        choices = rank_list(),
        selected = rank_list()[1],
        multiple = TRUE,
        options = shinyWidgets::pickerOptions(
          actionsBox = TRUE,
          liveSearch = TRUE,
          showContent = FALSE
        )
      )
    })

    # ---- import logs ------------------------------------------------------

    output$fact_color_file_log <- renderPrint({
      if (is.null(input$file_var_color)) {
        cat("No colors have been loaded.")
        return(invisible())
      }
      parsed <- parsed_factor_csv()
      loaded_vars <- unique(c(names(parsed$colors), names(parsed$numerics)))
      if (length(loaded_vars) == 0L) {
        cat("No valid rows were loaded from the CSV file.\n")
      } else {
        cat("Loaded color assignments for: ",
            paste(loaded_vars, collapse = ", "), "\n", sep = "")
      }
      filled <- attr(fill_missing_colors(parsed$colors, factor_modalities()),
                     "missing")
      if (length(filled) > 0L) {
        any_miss <- sapply(filled, length) > 0
        if (any(any_miss)) {
          cat("\nMissing modalities (auto-filled from the fallback palette):\n")
          print(filled[any_miss])
        }
      }
      if (nrow(parsed$dropped) > 0L) {
        cat("\nDropped rows:\n")
        print(parsed$dropped)
      }
    })

    output$taxa_color_file_log <- renderPrint({
      if (is.null(input$file_taxa_color)) {
        cat("No colors have been loaded.")
        return(invisible())
      }
      parsed <- parsed_taxa_csv()
      loaded_ranks <- names(parsed$colors)
      if (length(loaded_ranks) == 0L) {
        cat("No valid rows were loaded from the CSV file.\n")
      } else {
        cat("Loaded color assignments for ranks: ",
            paste(loaded_ranks, collapse = ", "), "\n", sep = "")
      }
      filled <- attr(fill_missing_colors(parsed$colors, taxa_modalities()),
                     "missing")
      if (length(filled) > 0L) {
        nb <- sapply(filled, length)
        if (any(nb > 0L)) {
          cat("\nNumber of missing taxa auto-filled per rank:\n")
          print(nb[nb > 0L])
        }
      }
      if (nrow(parsed$dropped) > 0L) {
        cat("\nDropped rows:\n")
        print(parsed$dropped)
      }
    })

    # ---- downloads --------------------------------------------------------

    output$var_color_download <- downloadHandler(
      filename = "var_colors.csv",
      content = function(file) {
        req(r$factor_colors(), r$numeric_palettes())
        selected <- input$list_var_color
        if (length(selected) == 0L) {
          showNotification("No variable selected.", type = "warning",
                           duration = 5)
          return()
        }
        if (length(selected) >= 50L) {
          showNotification("Too many variables to export.", type = "warning",
                           duration = 10)
          return()
        }
        df <- build_factor_color_df(
          factor_colors    = r$factor_colors(),
          numeric_palettes = r$numeric_palettes(),
          selected_vars    = selected
        )
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    output$taxa_color_download <- downloadHandler(
      filename = "taxa_colors.csv",
      content = function(file) {
        req(r$taxa_colors())
        selected <- input$list_taxa_color
        if (length(selected) == 0L) {
          showNotification("No rank selected.", type = "warning", duration = 5)
          return()
        }
        df <- build_taxa_color_df(
          taxa_colors    = r$taxa_colors(),
          selected_ranks = selected
        )
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

  })
}

## To be copied in the UI
# mod_color_ui("color_ui_1")

## To be copied in the server
# mod_color_server("color_ui_1", r = r)
