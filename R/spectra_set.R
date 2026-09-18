#' Build a SpectraSet from raw spectra lists
#'
#' `SpectraSet` is a lightweight container for sparse-bin distance backends.
#' Sparse matrices are built lazily on first use and cached per side/weighting
#' scheme, so callers should reuse the same object across multiple methods.
#'
#' The sparse-bin representation is an approximation of the legacy
#' pair-aligned ppm matching path, not a bit-for-bit replacement.
#'
#' @param frag_list List of two-column fragment spectra matrices.
#' @param loss_list Optional list of two-column loss spectra matrices.
#' @param ids Optional spectrum IDs. Defaults to `names(frag_list)`.
#' @param mz_min,mz_max Mass range for the ppm-spaced bin grid.
#' @param bin_ppm Bin spacing in ppm.
#' @param smear_ppm Half-width of the peak smearing kernel in ppm.
#' @param smear_kind Smearing kernel, either `"tri"` or `"rect"`.
#' @return A `SpectraSet` S3 object.
#' @export
as_spectra_set <- function(frag_list,
                           loss_list = NULL,
                           ids = NULL,
                           mz_min = 35,
                           mz_max = 650,
                           bin_ppm = 2,
                           smear_ppm = 15,
                           smear_kind = c("tri", "rect")) {
  smear_kind <- match.arg(smear_kind)
  if (!is.finite(bin_ppm) || bin_ppm <= 0) stop("bin_ppm must be positive.")
  if (!is.finite(smear_ppm) || smear_ppm <= 0) stop("smear_ppm must be positive.")
  if (!is.finite(mz_min) || !is.finite(mz_max) || mz_min <= 0 || mz_max <= mz_min) {
    stop("mz_min and mz_max must be finite and satisfy 0 < mz_min < mz_max.")
  }
  if (smear_ppm < bin_ppm) {
    warning("smear_ppm < bin_ppm; consider lowering bin_ppm for stable smearing.")
  }

  n <- length(frag_list)
  if (is.null(ids)) {
    ids <- names(frag_list)
    if (is.null(ids)) ids <- as.character(seq_len(n))
  }
  if (length(ids) != n) stop("ids must have the same length as frag_list.")
  if (!is.null(loss_list) && length(loss_list) != n) {
    stop("loss_list must be NULL or have the same length as frag_list.")
  }

  sort_spec <- function(s) {
    if (is.null(s) || nrow(s) == 0) return(matrix(numeric(0), ncol = 2))
    s[order(s[, 1]), , drop = FALSE]
  }

  bin_factor <- 1 + bin_ppm * 1e-6
  n_bins <- ceiling(log(mz_max / mz_min) / log(bin_factor)) + 1L
  mz_bins <- mz_min * bin_factor^seq.int(0, n_bins - 1L)
  bin_edges <- mz_bins / sqrt(bin_factor)

  obj <- list(
    n = n,
    ids = as.character(ids),
    frag_pp = lapply(frag_list, sort_spec),
    loss_pp = if (!is.null(loss_list)) lapply(loss_list, sort_spec) else NULL,
    mz_bins = mz_bins,
    bin_edges = bin_edges,
    bin_ppm = bin_ppm,
    smear_ppm = smear_ppm,
    smear_bins = max(1L, as.integer(round(smear_ppm / bin_ppm))),
    smear_kind = smear_kind,
    n_bins = length(mz_bins),
    cache = new.env(parent = emptyenv())
  )
  class(obj) <- c("SpectraSet", "list")
  obj
}

#' @export
print.SpectraSet <- function(x, ...) {
  cat("<SpectraSet>\n")
  cat("  n spectra  :", x$n, "\n")
  cat("  bin grid   :", x$n_bins, "bins,", x$bin_ppm, "ppm spacing\n")
  cat("  smear      :", x$smear_kind, "kernel,", x$smear_ppm, "ppm half-width\n")
  cat("  has loss   :", !is.null(x$loss_pp), "\n")
  cached <- ls(x$cache)
  cat("  cached     :", if (length(cached)) paste(cached, collapse = ", ") else "(none)", "\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
ensure_sparse_matrix <- function(set, side = c("frag", "loss"),
                                 weight = c("raw", "wcos"),
                                 params = list(mass_power = 3, intensity_power = 0.5)) {
  side <- match.arg(side)
  weight <- match.arg(weight)
  key <- paste0(
    side, "_M_", weight,
    if (weight == "wcos") sprintf("_mp%g_ip%g", params$mass_power, params$intensity_power) else ""
  )

  if (exists(key, envir = set$cache, inherits = FALSE)) {
    return(get(key, envir = set$cache, inherits = FALSE))
  }

  spec_list <- if (side == "frag") set$frag_pp else set$loss_pp
  if (is.null(spec_list)) stop("SpectraSet has no '", side, "' side.")

  M <- build_smeared_sparse(
    spec_list = spec_list,
    bin_edges = set$bin_edges,
    mz_bins = set$mz_bins,
    smear_bins = set$smear_bins,
    smear_kind = set$smear_kind,
    weight = weight,
    params = params
  )

  assign(key, M, envir = set$cache)
  M
}

#' @keywords internal
#' @noRd
build_smeared_sparse <- function(spec_list, bin_edges, mz_bins,
                                 smear_bins, smear_kind, weight, params) {
  n_spec <- length(spec_list)
  n_bins <- length(mz_bins)

  if (smear_kind == "rect") {
    kernel <- rep(1, 2L * smear_bins + 1L)
  } else {
    kernel <- 1 - abs(seq.int(-smear_bins, smear_bins)) / (smear_bins + 1)
  }
  kernel <- kernel / sum(kernel)

  total_peaks <- sum(vapply(spec_list, function(s) if (is.null(s)) 0L else nrow(s), integer(1)))
  if (total_peaks == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(n_spec, n_bins), repr = "C"
    ))
  }

  est_nnz <- total_peaks * length(kernel)
  i_idx <- integer(est_nnz)
  j_idx <- integer(est_nnz)
  v_val <- numeric(est_nnz)
  pos <- 0L

  for (k in seq_len(n_spec)) {
    s <- spec_list[[k]]
    if (is.null(s) || nrow(s) == 0) next

    mz_p <- s[, 1]
    in_p <- s[, 2]
    intensity_used <- switch(weight,
      raw = in_p,
      wcos = (mz_p ^ params$mass_power) * (pmax(in_p, 0) ^ params$intensity_power)
    )

    bin_idx <- findInterval(mz_p, bin_edges, all.inside = TRUE)
    for (pp in seq_along(mz_p)) {
      ctr <- bin_idx[pp]
      lo <- max(1L, ctr - smear_bins)
      hi <- min(n_bins, ctr + smear_bins)
      ker_lo <- lo - (ctr - smear_bins) + 1L
      ker_hi <- ker_lo + (hi - lo)
      span <- hi - lo + 1L

      idx_range <- (pos + 1L):(pos + span)
      i_idx[idx_range] <- k
      j_idx[idx_range] <- lo:hi
      v_val[idx_range] <- intensity_used[pp] * kernel[ker_lo:ker_hi]
      pos <- pos + span
    }
  }

  Matrix::sparseMatrix(
    i = i_idx[seq_len(pos)],
    j = j_idx[seq_len(pos)],
    x = v_val[seq_len(pos)],
    dims = c(n_spec, n_bins),
    repr = "C"
  )
}

#' @keywords internal
#' @noRd
sparse_row_norms <- function(M) {
  sqrt(Matrix::rowSums(M * M))
}

#' @keywords internal
#' @noRd
validate_compatible_grids <- function(A, B) {
  if (!inherits(A, "SpectraSet") || !inherits(B, "SpectraSet")) {
    stop("Both inputs must be SpectraSet objects.")
  }
  if (A$n_bins != B$n_bins) stop("Bin grids differ: ", A$n_bins, " vs ", B$n_bins, " bins.")
  same_grid <- abs(A$mz_bins[1] - B$mz_bins[1]) <= 1e-6 &&
    abs(A$mz_bins[length(A$mz_bins)] - B$mz_bins[length(B$mz_bins)]) <= 1e-6 &&
    abs(A$bin_ppm - B$bin_ppm) <= 1e-9 &&
    abs(A$smear_ppm - B$smear_ppm) <= 1e-9 &&
    identical(A$smear_kind, B$smear_kind)
  if (!same_grid) stop("Bin grid parameters differ between SpectraSets.")
  invisible(TRUE)
}
