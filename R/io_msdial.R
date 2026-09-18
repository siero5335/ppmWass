#' Read MS-DIAL Export
#'
#' @param file Path to MS-DIAL export file.
#' @param skip Number of metadata rows before the tabular header (default 4).
#' @param show_col_types Passed to readr::read_tsv().
#' @return An `msdial_import` object with the main table and preserved header metadata.
#' @export
read_ms_dial <- function(file, skip = 4, show_col_types = FALSE) {
  raw_lines <- readr::read_lines(file, n_max = skip + 1)
  if (length(raw_lines) < skip + 1) {
    stop("MS-DIAL export must contain metadata rows and a tabular header.")
  }

  data <- readr::read_tsv(
    file,
    skip = skip,
    show_col_types = show_col_types,
    name_repair = "unique_quiet"
  )

  preview <- utils::read.delim(
    file = file,
    sep = "\t",
    header = FALSE,
    nrows = skip + 1,
    fill = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = character()
  )

  sample_meta <- build_msdial_sample_meta(data, preview, skip = skip)

  structure(
    list(
      data = data,
      sample_meta = sample_meta,
      raw_header = raw_lines[seq_len(skip)],
      source = file
    ),
    class = "msdial_import"
  )
}

#' @keywords internal
#' @noRd
is_msdial_import <- function(x) {
  inherits(x, "msdial_import")
}

#' @keywords internal
#' @noRd
as_feature_table <- function(x) {
  if (is_msdial_import(x)) {
    return(x$data)
  }
  x
}

#' @keywords internal
#' @noRd
get_sample_meta <- function(x) {
  if (is_msdial_import(x)) {
    return(x$sample_meta)
  }
  NULL
}

#' @keywords internal
#' @noRd
build_msdial_sample_meta <- function(data, preview, skip = 4) {
  data_cols <- names(data)
  n_cols <- length(data_cols)
  preview_tbl <- tibble::as_tibble(preview, .name_repair = "minimal")
  preview_ncol <- ncol(preview_tbl)
  offset <- max(0L, preview_ncol - n_cols)

  header_row <- as.character(preview_tbl[skip + 1, seq_len(min(preview_ncol, n_cols)), drop = TRUE])
  if (length(header_row) < n_cols) {
    header_row <- c(header_row, rep(NA_character_, n_cols - length(header_row)))
  }
  column_raw <- normalize_msdial_meta_value(header_row[seq_len(n_cols)])

  metadata_rows <- seq_len(min(skip, nrow(preview_tbl) - 1L))
  labels <- if (length(metadata_rows) > 0) {
    normalize_msdial_header_key(preview_tbl[[1]][metadata_rows])
  } else {
    character(0)
  }

  find_row <- function(expected, fallback) {
    idx <- which(labels == expected)
    if (length(idx) > 0) {
      return(metadata_rows[idx[1]])
    }
    if (fallback <= length(metadata_rows)) {
      return(metadata_rows[fallback])
    }
    NA_integer_
  }

  extract_values <- function(row_idx) {
    if (is.na(row_idx)) {
      return(rep(NA_character_, n_cols))
    }
    idx <- seq_len(n_cols) + offset
    idx <- idx[idx <= preview_ncol]
    values <- as.character(preview_tbl[row_idx, idx, drop = TRUE])
    if (length(values) < n_cols) {
      values <- c(values, rep(NA_character_, n_cols - length(values)))
    }
    normalize_msdial_meta_value(values[seq_len(n_cols)])
  }

  class_vals <- extract_values(find_row("class", 1))
  file_type_vals <- extract_values(find_row("file_type", 2))
  injection_vals <- extract_values(find_row("injection_order", 3))
  batch_vals <- extract_values(find_row("batch_id", 4))

  tibble::tibble(
    column = data_cols,
    column_raw = column_raw,
    index = seq_along(data_cols),
    class = class_vals,
    file_type = file_type_vals,
    injection_order = injection_vals,
    batch_id = batch_vals
  ) -> sample_meta

  sample_meta$role <- vapply(
    seq_len(nrow(sample_meta)),
    function(i) infer_msdial_role(
      column = sample_meta$column[[i]],
      column_raw = sample_meta$column_raw[[i]],
      class = sample_meta$class[[i]],
      file_type = sample_meta$file_type[[i]],
      injection_order = sample_meta$injection_order[[i]],
      batch_id = sample_meta$batch_id[[i]]
    ),
    character(1)
  )

  sample_meta
}

#' @keywords internal
#' @noRd
normalize_msdial_header_key <- function(x) {
  out <- tolower(trimws(as.character(x)))
  out <- gsub("[^a-z0-9]+", "_", out)
  out <- gsub("^_+|_+$", "", out)
  out
}

#' @keywords internal
#' @noRd
normalize_msdial_meta_value <- function(x) {
  missing <- is.na(x)
  out <- trimws(as.character(x))
  out[missing] <- NA_character_
  out[out == ""] <- NA_character_
  out
}

#' @keywords internal
#' @noRd
infer_msdial_role <- function(column, column_raw, class, file_type, injection_order, batch_id) {
  label <- tolower(paste(stats::na.omit(c(class, file_type, column_raw, column)), collapse = " "))
  has_metadata <- any(!is.na(c(class, file_type, injection_order, batch_id)))

  if (column %in% c(META_COLS, INTERNAL_COLS)) {
    return("feature_meta")
  }
  if (any(is_summary_col_name(c(column, column_raw)))) {
    return("summary")
  }
  if (is_qc_label(label)) {
    return("qc")
  }
  if (is_blank_label(label)) {
    return("blank")
  }
  if (is_standard_label(label)) {
    return("standard")
  }
  if (has_metadata) {
    return("sample")
  }
  "other"
}
