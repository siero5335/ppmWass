# Internal helpers shared by publication drivers.

#' Build a deterministic grouping key without dropping missing-value groups
#'
#' Base `interaction()` returns `NA` whenever any grouping column is missing;
#' `split()` then silently drops those rows. Publication parameter tables use
#' intentional `NA` values for inactive settings (for example Sinkhorn epsilon
#' under exact OT), so those rows must remain a real group.
#' @keywords internal
#' @noRd
na_safe_group_key <- function(data, columns = names(data)) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }
  if (!length(columns)) return(rep("<all>", nrow(data)))
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop(
      "Missing grouping columns: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  encoded <- lapply(data[columns], function(value) {
    value <- as.character(value)
    is_missing <- is.na(value)
    value[is_missing] <- ""
    ifelse(
      is_missing,
      "N",
      paste0("V", nchar(value, type = "bytes"), ":", value)
    )
  })
  do.call(paste, c(encoded, list(sep = "\u001f")))
}

#' Summarize a publication metric with an optional complete-replicate contract
#' @keywords internal
#' @noRd
summarize_publication_metric <- function(value, require_complete = FALSE) {
  finite <- is.finite(value)
  n_finite <- sum(finite)
  complete <- n_finite == length(value)
  eligible <- n_finite > 0L && (!isTRUE(require_complete) || complete)
  list(
    n_finite = n_finite,
    mean = if (eligible) mean(value[finite]) else NA_real_,
    sd = if (eligible && n_finite > 1L) stats::sd(value[finite]) else NA_real_,
    complete = complete
  )
}
