# Column detection helpers
#' @keywords internal
#' @noRd
is_repaired_numeric_name <- function(x) {
  grepl("^\\d+(\\.\\d+)?(\\.\\.\\.[0-9]+)?$", x)
}

#' @keywords internal
#' @noRd
is_summary_col_name <- function(x) {
  x <- tolower(trimws(as.character(x)))
  stringr::str_detect(x, "^(average|stdev)(\\.\\.\\.[0-9]+)?$") |
    is_repaired_numeric_name(x)
}

#' @keywords internal
#' @noRd
is_blank_label <- function(x) {
  stringr::str_detect(tolower(as.character(x)), "(blank|blk)")
}

#' @keywords internal
#' @noRd
is_qc_label <- function(x) {
  stringr::str_detect(
    tolower(as.character(x)),
    "(quality\\s*control|pooled|(^|[_ -])qc([_ -]|$)|^qc[0-9])"
  )
}

#' @keywords internal
#' @noRd
is_standard_label <- function(x) {
  stringr::str_detect(
    tolower(as.character(x)),
    "(standard|std($|[_ -]|[0-9]))"
  )
}

#' @keywords internal
#' @noRd
get_column_names <- function(x) {
  if (is.character(x)) {
    return(x)
  }
  names(as_feature_table(x))
}

#' @keywords internal
#' @noRd
get_cols_by_role <- function(x, roles) {
  sample_meta <- get_sample_meta(x)
  if (is.null(sample_meta)) {
    return(character(0))
  }
  sample_meta$column[sample_meta$role %in% roles]
}

#' @keywords internal
#' @noRd
get_sample_cols <- function(df) {
  cols <- get_cols_by_role(df, "sample")
  if (length(cols) > 0) {
    return(cols)
  }

  cand <- setdiff(get_column_names(df), c(META_COLS, INTERNAL_COLS))
  cand <- cand[!is_summary_col_name(cand)]
  cand <- setdiff(cand, c(get_blank_cols(cand), get_qc_cols(cand), get_std_cols(cand)))
  cand
}

#' @keywords internal
#' @noRd
get_blank_cols <- function(x) {
  cols <- get_cols_by_role(x, "blank")
  if (length(cols) > 0) {
    return(cols)
  }
  cols <- get_column_names(x)
  cols[is_blank_label(cols)]
}

#' @keywords internal
#' @noRd
get_qc_cols <- function(x) {
  cols <- get_cols_by_role(x, "qc")
  if (length(cols) > 0) {
    return(cols)
  }
  cols <- get_column_names(x)
  cols[is_qc_label(cols)]
}

#' @keywords internal
#' @noRd
get_std_cols <- function(x) {
  cols <- get_cols_by_role(x, "standard")
  if (length(cols) > 0) {
    return(cols)
  }
  cols <- get_column_names(x)
  cols[is_standard_label(cols)]
}

#' @keywords internal
#' @noRd
get_quant_cols <- function(x, include_qc = TRUE, include_blank = TRUE, include_std = TRUE) {
  sample_meta <- get_sample_meta(x)
  if (!is.null(sample_meta)) {
    roles <- c("sample")
    if (isTRUE(include_qc)) roles <- c(roles, "qc")
    if (isTRUE(include_blank)) roles <- c(roles, "blank")
    if (isTRUE(include_std)) roles <- c(roles, "standard")
    cols <- sample_meta$column[sample_meta$role %in% roles]
    return(cols)
  }

  cols <- get_sample_cols(x)
  if (isTRUE(include_qc)) cols <- unique(c(cols, get_qc_cols(x)))
  if (isTRUE(include_blank)) cols <- unique(c(cols, get_blank_cols(x)))
  if (isTRUE(include_std)) cols <- unique(c(cols, get_std_cols(x)))
  all_cols <- get_column_names(x)
  all_cols[all_cols %in% cols]
}

#' @keywords internal
#' @noRd
get_signal_cols <- function(x, include_qc = FALSE) {
  sample_meta <- get_sample_meta(x)
  if (!is.null(sample_meta)) {
    roles <- "sample"
    if (isTRUE(include_qc)) {
      roles <- c(roles, "qc")
    }
    return(sample_meta$column[sample_meta$role %in% roles])
  }

  cols <- get_sample_cols(x)
  if (isTRUE(include_qc)) {
    cols <- unique(c(cols, get_qc_cols(x)))
    all_cols <- get_column_names(x)
    cols <- all_cols[all_cols %in% cols]
  }
  cols
}

#' @keywords internal
#' @noRd
derive_compound_id <- function(df) {
  dplyr::if_else(
    df$`Annotation tag (VS1.0)` == "4",
    stringr::str_c(df$`Metabolite name`, df$`Alignment ID`, sep = "_"),
    df$`Metabolite name`
  )
}

#' @keywords internal
#' @noRd
resolve_reference_row <- function(df, ref) {
  if (is.null(ref)) return(NA_integer_)
  if (is.numeric(ref)) {
    idx <- as.integer(ref[1])
    if (is.na(idx) || idx < 1 || idx > nrow(df)) {
      stop("normalize_ref row index is out of range.")
    }
    return(idx)
  }
  ref_chr <- as.character(ref[1])
  idx <- which(as.character(df$`Alignment ID`) == ref_chr)
  if (length(idx) == 0) {
    idx <- which(as.character(df$`Metabolite name`) == ref_chr)
  }
  if (length(idx) == 0) {
    stop("normalize_ref not found in Alignment ID or Metabolite name.")
  }
  if (length(idx) > 1) {
    warning("normalize_ref matched multiple rows; using the first match.")
  }
  idx[1]
}
