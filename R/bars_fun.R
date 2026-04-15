
#' aggregate_top_taxa from microbiome package
#'
#' @param x phyloseq object
#' @param top Keep the top-n taxa, and merge the rest under the category
#'   'Other'. Can also be a character vector listing the groups to combine.
#' @param level Summarization level (from 'rank_names(pseq)')
#'
#' @importFrom microbiome aggregate_taxa
#' @importFrom microbiome top_taxa
#'
#' @export
aggregate_top_taxa <- function(x, top, level) {
  x <- aggregate_taxa(x, level)
  tops <- top_taxa(x, top)
  tax <- tax_table(x)
  inds <- which(!rownames(tax) %in% tops)
  tax[inds, level] <- "Other"
  tax_table(x) <- tax
  tt <- tax_table(x)[, level]
  tax_table(x) <- tax_table(tt)
  aggregate_taxa(x, level)
}


#' Barplots plotly
#'
#' @param data a phyloseq object
#' @param rank Taxonomy rank to merge features (among rank_names(data), or
#'   'ASV' for no glom)
#' @param top Number of top taxa to plot
#' @param Ord1 Variable used to order samples (X axis) or split the barplot
#' @param sample_labels If TRUE, x axis labels are sample IDs
#' @param split If TRUE, make a facet_wrap-like plot grouped by Ord1
#' @param split_sid_order If TRUE, keep sample IDs order from metadata in split
#' @param relative Plot relative (TRUE) or raw abundance (FALSE)
#' @param autoorder Automatic ordering of x axis labels with mixedsort
#' @param ylab Y axis title
#' @param outfile Output html file (NULL to skip saving)
#' @param pal A color palette
#' @param verbose Print log messages
#'
#' @return An interactive plotly barplot
#'
#' @importFrom plotly plot_ly subplot layout
#' @importFrom tidyr pivot_longer
#' @importFrom gtools mixedsort
#' @importFrom dplyr filter mutate arrange group_by group_map pull across where
#' @importFrom rlang .data `:=`
#' @import futile.logger
#'
#' @export
bars_fun <- function(data, rank = "Genus", top = 10, Ord1 = NULL,
                     sample_labels = FALSE, split = FALSE,
                     split_sid_order = FALSE, relative = TRUE,
                     autoorder = TRUE, ylab = "Abundance",
                     outfile = "plot_compo.html", verbose = TRUE,
                     pal = NULL) {

  if (verbose) invisible(flog.threshold(INFO)) else invisible(flog.threshold(ERROR))

  if (!Ord1 %in% sample_variables(data)) {
    stop("Wrong value in Ord1, use: ", toString(sample_variables(data)))
  }

  # --- 1. Preprocess phyloseq → sdata + otable ---
  flog.info("Preprocess...")
  psobj.top <- aggregate_top_taxa(data, rank, top = top)

  sdata <- as.data.frame(sample_data(psobj.top), stringsAsFactors = TRUE)
  sdata$sample.id <- sample_names(psobj.top)

  otable <- as.data.frame(otu_table(psobj.top))
  rownames(otable) <- tax_table(psobj.top)[, rank]

  # --- 2. Normalize if relative ---
  if (relative) {
    otable <- apply(otable, 2, function(x) x / sum(x))
  }

  # --- 3. Melt to long format (single pass) ---
  flog.info("Melting table...")
  n_meta <- ncol(sdata)
  dat <- cbind.data.frame(sdata, as.data.frame(t(otable)))

  meltdat <- dat |>
    tidyr::pivot_longer(
      cols = -seq_len(n_meta),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::filter(!is.na(.data[[Ord1]])) |>
    dplyr::mutate(variable = factor(variable))

  # Put "Other" first in the legend
  lvls <- levels(meltdat$variable)
  meltdat$variable <- factor(meltdat$variable, levels = c("Other", lvls[lvls != "Other"]))

  # --- 4. Order samples ---
  if (autoorder) {
    flog.info("Ordering samples...")
    sorted_levels <- gtools::mixedsort(unique(meltdat[[Ord1]]))
    meltdat <- meltdat |>
      dplyr::mutate(!!Ord1 := factor(.data[[Ord1]], levels = sorted_levels)) |>
      dplyr::arrange(.data[[Ord1]], sample.id)
  }

  ordered_ids <- unique(meltdat$sample.id)
  ord_levels <- levels(factor(meltdat[[Ord1]]))

  # --- 5. Build x-axis labels ---
  lab_col <- if (sample_labels) "sample.id" else Ord1
  xform <- list(
    categoryorder = "array",
    categoryarray = ordered_ids,
    title = "Samples",
    tickmode = "array",
    tickvals = seq(0, nrow(sdata)),
    ticktext = sdata[ordered_ids, lab_col],
    tickangle = -90
  )

  # --- 6. Group color bar (subplot below main plot) ---
  df_groups <- data.frame(
    x = ordered_ids,
    g = factor(sdata[ordered_ids, Ord1], levels = ord_levels),
    y = 1
  )

  has_groups <- length(unique(df_groups$g)) > 1

  subp1 <- plotly::plot_ly(
    df_groups, type = "bar", x = ~x, y = ~y,
    color = ~g, legendgroup = ~g,
    showlegend = FALSE, colors = pal
  ) |>
    plotly::layout(
      xaxis = list(zeroline = FALSE, showline = FALSE, showgrid = FALSE),
      yaxis = list(showticklabels = FALSE, title = "", showgrid = FALSE)
    )

  # --- 7. Build main barplot or split ---
  plot_title <- if (relative) "Relative abundance" else "Raw abundance"

  if (!split) {
    flog.info("Plotting %s...", plot_title)
    p1 <- plotly::plot_ly(
      meltdat, x = ~sample.id, y = ~value,
      type = "bar", name = ~variable, color = ~variable, colors = pal
    ) |>
      plotly::layout(
        title = plot_title,
        yaxis = list(title = ylab),
        xaxis = xform, barmode = "stack"
      )

    if (has_groups) {
      p1 <- plotly::subplot(p1, subp1, nrows = 2, shareX = TRUE, heights = c(0.95, 0.05)) |>
        plotly::layout(xaxis = xform)
    }
  } else {
    flog.info("Splitted plot...")
    if (split_sid_order) {
      meltdat$sample.id <- factor(meltdat$sample.id, levels = unique(meltdat$sample.id))
    }

    p1 <- meltdat |>
      dplyr::arrange(.data[[Ord1]]) |>
      dplyr::group_by(.data[[Ord1]]) |>
      dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor)) |>
      dplyr::group_map(
        ~ plotly::plot_ly(
          data = .x, x = ~sample.id, y = ~value, type = "bar",
          name = ~variable, color = ~variable, legendgroup = ~variable,
          showlegend = (dplyr::pull(.y, 1) == ord_levels[1]),
          colors = pal
        ),
        .keep = TRUE
      ) |>
      plotly::subplot(nrows = 1, shareX = TRUE, shareY = TRUE, titleX = FALSE) |>
      plotly::layout(
        title = "",
        xaxis = list(title = glue::glue("{Ord1} =\n{ord_levels[1]}")),
        yaxis = list(title = ylab),
        barmode = "stack"
      )

    for (i in 2:length(ord_levels)) {
      p1$x$layoutAttrs[[1]][[paste0("xaxis", i)]] <- NULL
      p1$x$layoutAttrs[[1]][[paste0("xaxis", i)]]$title <- glue::glue("{Ord1} =\n{ord_levels[i]}")
    }
  }

  # --- 8. Save and return ---
  if (!is.null(outfile)) htmlwidgets::saveWidget(p1, outfile)
  flog.info("Finish.")
  p1
}
