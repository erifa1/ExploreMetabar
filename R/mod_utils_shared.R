#' Shared reactive utilities for ExploreMetabar modules
#'
#' Factory functions that create commonly-used reactives for modules that
#' work with phyloseq metadata factors. Each function returns a reactive
#' expression that can be used inside a module server.
#'
#' @noRd
#' @import shiny
#' @import dplyr
#' @import tidyr

#' Create a reactive that returns the metadata column name(s) for the selected factor(s)
#'
#' Handles single and multi-factor selection. When multiple factors are selected,
#' they are collapsed into a single column name with '_' separator.
#'
#' @param factor_input A reactive expression returning the selected factor(s)
#' @param r The shared reactive values object
#' @return A reactive returning the metadata column name
#' @noRd
make_get_meta_col <- function(factor_input, r) {
  reactive({
    req(factor_input(), r$sdat())
    metadata <- r$sdat()
    facts <- factor_input()
    if (length(facts) == 1) {
      meta.col <- facts
    } else if (length(facts) > 1) {
      validate(
        need(!any(sapply(metadata[, facts], is.numeric)),
             message = "You can't select multiple numeric factors")
      )
      meta.col <- paste0(facts, collapse = '_')
    }
    return(meta.col)
  })
}


#' Create a reactive that checks if the selected factor is numeric
#'
#' @param get_meta_col A reactive returning the metadata column name
#' @param local_metadata A reactive returning the local metadata data.frame
#' @return A reactive returning TRUE if the factor is numeric
#' @noRd
make_is_num_factor <- function(get_meta_col, local_metadata) {
  reactive({
    req(get_meta_col(), local_metadata())
    metadata <- local_metadata()
    is.numeric(metadata[, get_meta_col()])
  })
}


#' Create a reactive that returns the local metadata with united factor columns
#'
#' When multiple factors are selected, they are united into a single column.
#' Categorical factors are converted to factors.
#'
#' @param factor_input A reactive expression returning the selected factor(s)
#' @param get_meta_col A reactive returning the metadata column name
#' @param r The shared reactive values object
#' @param keep_all_cols If TRUE, keep all columns (used by mod_beta). If FALSE,
#'   select only sample.id and the factor column.
#' @return A reactive returning the metadata data.frame
#' @noRd
make_local_metadata <- function(factor_input, get_meta_col, r, keep_all_cols = FALSE) {
  reactive({
    req(factor_input(), r$sdat())
    metadata <- r$sdat()
    facts <- factor_input()
    if (!all(sapply(metadata[, facts], is.numeric))) {
      metadata <- tidyr::unite(metadata, !!get_meta_col(), facts, na.rm = TRUE)
      metadata[, get_meta_col()] <- as.factor(metadata[, get_meta_col()])
      if (!keep_all_cols) {
        metadata <- dplyr::select(metadata, "sample.id", get_meta_col())
      }
    } else {
      if (!keep_all_cols) {
        metadata <- dplyr::select(metadata, "sample.id", facts)
      }
    }
    return(metadata)
  })
}


#' Create a reactive that returns a phyloseq object with updated sample_data
#'
#' @param local_metadata A reactive returning the local metadata data.frame
#' @param r The shared reactive values object
#' @return A reactive returning the phyloseq object
#' @noRd
make_local_physeq <- function(local_metadata, r) {
  reactive({
    req(r$phyloseq_filtered(), local_metadata())
    phy <- r$phyloseq_filtered()
    phyloseq::sample_data(phy) <- phyloseq::sample_data(local_metadata())
    return(phy)
  })
}


#' Round numeric columns of a data.frame
#'
#' Helper used at the rendering layer of DT::renderDataTable to limit the
#' number of displayed decimals without altering the upstream reactive data
#' (which other reactives may depend on with full precision).
#'
#' Non-numeric columns (character, factor, logical, integer kept as-is when
#' not rounded meaningfully) are left untouched. If the input is not a
#' data.frame it is returned unchanged.
#'
#' @param x A data.frame (or coercible object). Non-data.frames are returned as-is.
#' @param digits Integer number of decimal places (default: 4).
#' @return A data.frame with numeric columns rounded to `digits` decimals.
#' @noRd
round_df <- function(x, digits = 4) {
  if (is.data.frame(x)) {
    x[] <- lapply(x, function(col) if (is.numeric(col)) round(col, digits) else col)
  }
  x
}
