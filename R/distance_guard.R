#' @keywords internal
#' @noRd
batch_distance_methods <- function() {
  c("cosine", "weighted_cosine", "hellinger")
}

#' @keywords internal
#' @noRd
guard_batch_compatible <- function(params, check_method = TRUE) {
  reasons <- character(0)

  if (isTRUE(params$use_typical_loss)) {
    reasons <- c(reasons, "use_typical_loss = TRUE not supported by sparse-bin batch path")
  }
  if (isTRUE(params$use_split_loss)) {
    reasons <- c(reasons, "use_split_loss = TRUE not supported by sparse-bin batch path")
  }
  if (isTRUE(params$use_mref_confidence)) {
    reasons <- c(reasons, "use_mref_confidence = TRUE not supported by sparse-bin batch path")
  }

  method <- params$distance_method %||% "cosine"
  if (isTRUE(check_method) && !(method %in% batch_distance_methods())) {
    reasons <- c(reasons, sprintf("distance_method '%s' has no sparse-bin batch implementation", method))
  }

  list(ok = length(reasons) == 0L, reason = reasons)
}

#' @keywords internal
#' @noRd
report_guard <- function(g, prefix = "[v2 sparse-bin] ") {
  for (r in g$reason) message(prefix, "skipping sparse-bin backend: ", r)
  invisible(NULL)
}
