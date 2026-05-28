#' Color helpers
#'
#' Pure (Shiny-free) helpers used by [mod_color_server()] to build the
#' `r$factor_colors()`, `r$numeric_palettes()` and `r$taxa_colors()` reactives.
#'
#' The producer logic used to live in `R/color_functions.R` as two parallel
#' branches (`fact` vs `taxa`). It was collapsed into a single parametrised
#' family. All randomness was removed so colors are stable across reactive
#' invalidations.
#'
#' @keywords internal
#' @name color_helpers
NULL


# ---- palette pools ----------------------------------------------------------

#' Qualitative palette pool used for factor / taxon coloring.
#'
#' @return data.frame with columns `package`, `palette`, `length` (subset of
#'   `paletteer::palettes_d_names`).
#' @keywords internal
qualitative_palette_pool <- function() {
  pal <- paletteer::palettes_d_names
  pal <- pal[pal$type == "qualitative", ]
  pal <- pal[pal$package %in% c("ggsci", "ggthemes", "jcolors", "pals",
                                "Polychrome", "RColorBrewer"), ]
  # Drop monochrome / grayscale / problematic palettes.
  pal <- pal[pal$package != "ggthemes" |
               !pal$palette %in% c("fivethirtyeight", "Seattle_Grays",
                                   "Classic_Gray_5", "excel_Grayscale",
                                   "stata_mono", "stata_economist"), ]
  pal <- pal[pal$package != "Polychrome" |
               !pal$palette %in% c("glasbey", "kelly"), ]
  pal <- pal[pal$package != "pals" | pal$palette != "kelly", ]
  pal
}


#' Continuous palette pool used for numeric metadata variables.
#'
#' @return data.frame with columns `package`, `palette` (subset of
#'   `paletteer::palettes_c_names`).
#' @keywords internal
continuous_palette_pool <- function() {
  pal <- paletteer::palettes_c_names
  pal[pal$package %in% c("viridis", "ggthemes"), ]
}


#' All `package::palette` identifiers accepted in `__numeric__` CSV rows.
#' @keywords internal
continuous_palette_ids <- function() {
  pool <- continuous_palette_pool()
  paste(pool$package, pool$palette, sep = "::")
}


# ---- modality extraction ----------------------------------------------------

#' Extract modalities for each metadata variable or taxonomic rank.
#'
#' @param phy a `phyloseq` object.
#' @param vars character vector of variable names (factors) or rank names.
#' @param type one of `"fact"` or `"taxa"`.
#' @return named list (one entry per `vars`) of character vectors.
#' @keywords internal
extract_modalities <- function(phy, vars, type = c("fact", "taxa")) {
  type <- match.arg(type)
  src <- if (type == "fact") {
    as.data.frame(phyloseq::sample_data(phy))
  } else {
    as.data.frame(phyloseq::tax_table(phy))
  }
  out <- lapply(vars, function(v) {
    if (!v %in% colnames(src)) return(character(0))
    as.character(unique(src[[v]]))
  })
  names(out) <- vars
  out
}


# ---- automatic color assignment --------------------------------------------

#' Build a deterministic large-cardinality qualitative palette.
#'
#' Used when the requested cardinality exceeds the largest single qualitative
#' palette in the pool (typically `ggsci::default_igv`, 51 colors). The
#' previous implementation fell back to `grDevices::hcl.colors("Spectral")`,
#' which is a *sequential* ramp — taxa at high cardinality ended up colored
#' as variations of red sliding to purple. This helper instead concatenates
#' several of the largest qualitative palettes (deduplicated by hex) so the
#' first ~100 colors are visually distinct, then cycles to cover any
#' remainder.
#'
#' @param n integer, number of colors required.
#' @param palette_ids character vector of `package::palette` identifiers,
#'   in concatenation order.
#' @return character vector of length `n`.
#' @keywords internal
high_cardinality_palette <- function(
    n,
    palette_ids = c("ggsci::default_igv", "Polychrome::palette36",
                    "pals::glasbey", "ggsci::default_ucscgb")) {
  hex <- unique(unlist(lapply(palette_ids, function(id) {
    as.character(paletteer::paletteer_d(palette = id))
  }), use.names = FALSE))
  if (n <= length(hex)) hex[seq_len(n)] else rep_len(hex, n)
}


#' Deterministically assign hex colors to a set of grouped modalities.
#'
#' Picks the smallest palette in `pool` that covers each variable's
#' cardinality, then removes that palette from the pool so subsequent
#' variables get a distinct palette. Variables exceeding the largest
#' available palette fall back to [high_cardinality_palette()] — a
#' concatenation of several large qualitative palettes that stays visually
#' distinct well above the 51-color ceiling of the pool.
#'
#' `NA` modalities are dropped from the auto-assignment and mapped to
#' `"#000000"` under the literal name `"NA"`, matching the historical
#' behavior used by the consumer modules.
#'
#' @param modalities_list named list of character vectors (one per variable).
#' @param pool palette pool from [qualitative_palette_pool()].
#' @param type `"fact"` (pick first palette of the smallest-fit tie) or
#'   `"taxa"` (pick last) — preserves the legacy bias.
#' @return named list of named character vectors.
#' @keywords internal
assign_qualitative_colors <- function(modalities_list,
                                      pool = qualitative_palette_pool(),
                                      type = c("fact", "taxa")) {
  type <- match.arg(type)
  out <- list()
  for (var in names(modalities_list)) {
    mods <- modalities_list[[var]]
    has_na <- any(is.na(mods))
    mods <- mods[!is.na(mods)]
    n <- length(mods)
    na_pair <- if (has_na) c(`NA` = "#000000") else NULL
    if (n == 0) {
      out[[var]] <- na_pair
      next
    }
    fit <- pool[pool$length >= n, , drop = FALSE]
    if (nrow(fit) > 0) {
      fit <- fit[fit$length == min(fit$length), , drop = FALSE]
      idx <- if (type == "fact") 1L else nrow(fit)
      pkg <- fit$package[idx]
      plt <- fit$palette[idx]
      hex <- as.character(paletteer::paletteer_d(
        palette = paste(pkg, plt, sep = "::"), n = n
      ))[seq_len(n)]
      pool <- pool[!(pool$package == pkg & pool$palette == plt), , drop = FALSE]
    } else {
      hex <- high_cardinality_palette(n)
    }
    names(hex) <- as.character(mods)
    out[[var]] <- c(hex, na_pair)
  }
  out
}


#' Assign a continuous palette identifier to each numeric variable.
#'
#' Cycles through `pool` if there are more numeric variables than palettes.
#'
#' @param num_vars character vector of numeric metadata variable names.
#' @param pool palette pool from [continuous_palette_pool()].
#' @return named list of `"package::palette"` strings.
#' @keywords internal
assign_numeric_palettes <- function(num_vars,
                                    pool = continuous_palette_pool()) {
  if (length(num_vars) == 0) return(list())
  idx <- ((seq_along(num_vars) - 1L) %% nrow(pool)) + 1L
  pals <- paste(pool$package[idx], pool$palette[idx], sep = "::")
  names(pals) <- num_vars
  as.list(pals)
}


# ---- CSV import / export ----------------------------------------------------

#' Read a long-format color CSV.
#'
#' Expected columns: `variable`, `modality`, `color`. Extra columns are
#' ignored. Returns a character-only data.frame.
#'
#' @param path file path.
#' @return data.frame with columns `variable`, `modality`, `color`, all
#'   character. NULL if the file cannot be parsed.
#' @keywords internal
read_color_csv <- function(path) {
  df <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)
  required <- c("variable", "modality", "color")
  if (!all(required %in% colnames(df))) return(NULL)
  df <- df[, required, drop = FALSE]
  df$variable <- as.character(df$variable)
  df$modality <- as.character(df$modality)
  df$color    <- as.character(df$color)
  df
}


#' Apply a parsed color CSV against the ground-truth modality list.
#'
#' Validates each row and silently drops rows that don't match the
#' ground-truth (with a structured report so the UI can surface what was
#' rejected). Rows whose `variable == "__numeric__"` are routed to the
#' numeric palette map; `modality` then carries the numeric variable name.
#'
#' @param df data.frame returned by [read_color_csv()].
#' @param modalities_list ground-truth named list (var -> character vector
#'   of valid modalities).
#' @param accept_numeric whether to honor `__numeric__` rows. Set to `FALSE`
#'   when parsing the taxa CSV.
#' @return list with elements:
#'   * `colors` — named list of named hex vectors (factors only).
#'   * `numerics` — named list of `package::palette` strings.
#'   * `dropped` — data.frame of dropped rows with a `reason` column.
#' @keywords internal
apply_color_csv <- function(df, modalities_list, accept_numeric = TRUE) {
  dropped <- data.frame(variable = character(0), modality = character(0),
                        color = character(0), reason = character(0),
                        stringsAsFactors = FALSE)
  colors <- list()
  numerics <- list()
  if (is.null(df) || nrow(df) == 0) {
    return(list(colors = colors, numerics = numerics, dropped = dropped))
  }
  is_num <- df$variable == "__numeric__"
  if (accept_numeric && any(is_num)) {
    valid_pals <- continuous_palette_ids()
    num_df <- df[is_num, , drop = FALSE]
    for (i in seq_len(nrow(num_df))) {
      v <- num_df$modality[i]
      p <- num_df$color[i]
      if (p %in% valid_pals) {
        numerics[[v]] <- p
      } else {
        dropped <- rbind(dropped, cbind(num_df[i, , drop = FALSE],
                                         reason = "unknown continuous palette"))
      }
    }
  } else if (any(is_num)) {
    num_df <- df[is_num, , drop = FALSE]
    dropped <- rbind(dropped, cbind(num_df,
                                     reason = "numeric row not accepted here"))
  }
  fact_df <- df[!is_num, , drop = FALSE]
  for (var in unique(fact_df$variable)) {
    sub <- fact_df[fact_df$variable == var, , drop = FALSE]
    if (!var %in% names(modalities_list)) {
      dropped <- rbind(dropped,
                       cbind(sub, reason = "unknown variable / rank"))
      next
    }
    known <- as.character(modalities_list[[var]])
    valid_color <- plotfunctions::isColor(sub$color)
    valid_mod   <- sub$modality %in% known
    bad <- !(valid_color & valid_mod)
    if (any(bad)) {
      reason <- ifelse(!valid_color, "invalid color",
                ifelse(!valid_mod, "unknown modality", "ok"))[bad]
      dropped <- rbind(dropped, cbind(sub[bad, , drop = FALSE],
                                       reason = reason))
    }
    sub <- sub[!bad, , drop = FALSE]
    if (nrow(sub) > 0) {
      hex <- sub$color
      names(hex) <- sub$modality
      colors[[var]] <- hex
    }
  }
  list(colors = colors, numerics = numerics, dropped = dropped)
}


#' Fill any modality missing from `loaded` against the ground-truth list.
#'
#' Deterministic — uses [high_cardinality_palette()] so re-runs produce
#' identical, visually-distinct fills even for many missing modalities.
#'
#' @param loaded named list of named hex vectors (subset).
#' @param full_modalities ground-truth named list (var -> all modalities).
#' @return named list of named hex vectors with `attr(., "missing")`
#'   giving the names that were filled in for each variable.
#' @keywords internal
fill_missing_colors <- function(loaded, full_modalities) {
  out <- list()
  missing_report <- list()
  for (var in names(loaded)) {
    have <- loaded[[var]]
    full <- as.character(full_modalities[[var]])
    full <- full[!is.na(full)]
    miss <- setdiff(full, names(have))
    if (length(miss) > 0) {
      fill <- high_cardinality_palette(length(miss))
      names(fill) <- miss
      out[[var]] <- c(have, fill)
      missing_report[[var]] <- miss
    } else {
      out[[var]] <- have
      missing_report[[var]] <- character(0)
    }
  }
  attr(out, "missing") <- missing_report
  out
}


#' Build a long-format data.frame for CSV export of factor + numeric colors.
#'
#' @param factor_colors named list of named hex vectors.
#' @param numeric_palettes named list of `package::palette` strings.
#' @param selected_vars optional character vector restricting export.
#' @return data.frame with columns `variable`, `modality`, `color`.
#' @keywords internal
build_factor_color_df <- function(factor_colors, numeric_palettes = list(),
                                  selected_vars = NULL) {
  vars <- if (is.null(selected_vars)) {
    c(names(factor_colors), names(numeric_palettes))
  } else selected_vars
  rows <- list()
  for (v in vars) {
    if (v %in% names(factor_colors)) {
      hex <- factor_colors[[v]]
      rows[[length(rows) + 1L]] <- data.frame(
        variable = v, modality = names(hex), color = unname(hex),
        stringsAsFactors = FALSE
      )
    } else if (v %in% names(numeric_palettes)) {
      rows[[length(rows) + 1L]] <- data.frame(
        variable = "__numeric__", modality = v,
        color = numeric_palettes[[v]], stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(variable = character(0), modality = character(0),
                      color = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}


#' Build a long-format data.frame for CSV export of taxa colors.
#' @keywords internal
build_taxa_color_df <- function(taxa_colors, selected_ranks = NULL) {
  ranks <- if (is.null(selected_ranks)) names(taxa_colors) else selected_ranks
  rows <- list()
  for (r in ranks) {
    hex <- taxa_colors[[r]]
    if (is.null(hex)) next
    rows[[length(rows) + 1L]] <- data.frame(
      variable = r, modality = names(hex), color = unname(hex),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(data.frame(variable = character(0), modality = character(0),
                      color = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}
