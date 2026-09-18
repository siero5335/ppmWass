#' @keywords internal
#' @noRd
cosine_distance_matrix <- function(Q_set, L_set, side = c("frag", "loss")) {
  side <- match.arg(side)
  validate_compatible_grids(Q_set, L_set)

  Q <- ensure_sparse_matrix(Q_set, side, weight = "raw")
  L <- ensure_sparse_matrix(L_set, side, weight = "raw")
  block_distance(Q, L, sparse_row_norms(Q), sparse_row_norms(L), "cosine",
                 row_ids = Q_set$ids, col_ids = L_set$ids)
}

#' @keywords internal
#' @noRd
weighted_cosine_distance_matrix <- function(Q_set, L_set,
                                            side = c("frag", "loss"),
                                            mass_power = 3,
                                            intensity_power = 0.5) {
  side <- match.arg(side)
  validate_compatible_grids(Q_set, L_set)

  pp <- list(mass_power = mass_power, intensity_power = intensity_power)
  Q <- ensure_sparse_matrix(Q_set, side, weight = "wcos", params = pp)
  L <- ensure_sparse_matrix(L_set, side, weight = "wcos", params = pp)
  block_distance(Q, L, sparse_row_norms(Q), sparse_row_norms(L), "weighted_cosine",
                 row_ids = Q_set$ids, col_ids = L_set$ids)
}

#' @keywords internal
#' @noRd
hellinger_distance_matrix <- function(Q_set, L_set, side = c("frag", "loss")) {
  side <- match.arg(side)
  validate_compatible_grids(Q_set, L_set)

  Q <- hellinger_normalise(ensure_sparse_matrix(Q_set, side, weight = "raw"))
  L <- hellinger_normalise(ensure_sparse_matrix(L_set, side, weight = "raw"))
  block_distance(Q, L, NULL, NULL, "hellinger", row_ids = Q_set$ids, col_ids = L_set$ids)
}

#' @keywords internal
#' @noRd
combine_frag_loss <- function(d_frag, d_loss = NULL, w_frag = 0.7, w_loss = 0.3) {
  if (is.null(d_loss)) return(d_frag)
  stopifnot(identical(dim(d_frag), dim(d_loss)))
  out <- sqrt(w_frag * d_frag^2 + w_loss * d_loss^2)
  rownames(out) <- rownames(d_frag)
  colnames(out) <- colnames(d_frag)
  out
}

#' Compute Dense Sparse-Bin Distance Matrix
#'
#' This experimental API materializes a dense `n_query x n_library` matrix and
#' is intended for small/medium GC benchmarks. For LC-scale retrieval, prefer
#' [compute_distance_topk_v2()], which keeps only top-k hits per query.
#'
#' If `params$distance_method` is not batchable, this function falls back to the
#' pair-loop path using the spectra stored in `Q_set` and `L_set`. If options
#' that require additional legacy-side inputs are enabled (`use_typical_loss`,
#' `use_split_loss`, `use_mref_confidence`), this function stops rather than
#' silently dropping those inputs; use the compatibility shim instead.
#'
#' @param Q_set,L_set `SpectraSet` objects sharing the same bin grid.
#' @param params Parameter list.
#' @param use_loss Whether to include loss spectra. `NULL` uses loss when both
#'   sets contain loss spectra.
#' @param max_dense_cells Maximum number of matrix cells allowed.
#' @param progress If `TRUE`, emit progress messages.
#' @return Dense numeric distance matrix.
#' @export
compute_distance_matrix_search_v2 <- function(Q_set, L_set, params,
                                              use_loss = NULL,
                                              max_dense_cells = 1e7,
                                              progress = TRUE) {
  check_common_param_typos(params)
  validate_compatible_grids(Q_set, L_set)

  cells <- as.numeric(Q_set$n) * as.numeric(L_set$n)
  if (cells > max_dense_cells) {
    stop(sprintf(
      paste0("Refusing to materialise dense %d x %d matrix (%g cells > %g). ",
             "Use compute_distance_topk_v2() for LC-scale retrieval, or increase ",
             "max_dense_cells explicitly."),
      Q_set$n, L_set$n, cells, max_dense_cells
    ))
  }

  feature_guard <- guard_batch_compatible(params, check_method = FALSE)
  if (!feature_guard$ok) {
    stop(
      "compute_distance_matrix_search_v2() cannot safely reproduce these options without ",
      "legacy auxiliary inputs: ", paste(feature_guard$reason, collapse = "; "),
      ". Use compute_distance_matrix_search() with backend = 'pair_loop' or the sparse-bin shim fallback."
    )
  }

  method <- params$distance_method %||% "cosine"
  has_loss <- !is.null(Q_set$loss_pp) && !is.null(L_set$loss_pp)
  if (is.null(use_loss)) use_loss <- has_loss

  if (!(method %in% batch_distance_methods())) {
    if (isTRUE(progress)) {
      message("[v2/Tier B] method = ", method, " is not batchable; using pair-loop fallback.")
    }
    return(pair_loop_distance_matrix(Q_set, L_set, params, use_loss = use_loss))
  }

  if (isTRUE(progress)) {
    message(sprintf("[v2/Tier B] dense matrix: %d x %d, method = %s", Q_set$n, L_set$n, method))
  }

  d_frag <- batch_one_side(Q_set, L_set, "frag", method, params)
  if (use_loss && has_loss) {
    d_loss <- batch_one_side(Q_set, L_set, "loss", method, params)
    return(combine_frag_loss(
      d_frag, d_loss,
      w_frag = params$w_frag %||% 0.7,
      w_loss = params$w_loss %||% 0.3
    ))
  }
  d_frag
}

#' @keywords internal
#' @noRd
batch_one_side <- function(Q_set, L_set, side, method, params) {
  switch(method,
    cosine = cosine_distance_matrix(Q_set, L_set, side),
    weighted_cosine = weighted_cosine_distance_matrix(
      Q_set, L_set, side,
      mass_power = params$mass_power %||% 3,
      intensity_power = params$intensity_power %||% 0.5
    ),
    hellinger = hellinger_distance_matrix(Q_set, L_set, side),
    stop("Unsupported batch method: ", method)
  )
}

#' @keywords internal
#' @noRd
pair_loop_distance_matrix <- function(Q_set, L_set, params, use_loss = TRUE) {
  nq <- Q_set$n
  nl <- L_set$n
  D <- matrix(0, nq, nl, dimnames = list(Q_set$ids, L_set$ids))
  empty <- matrix(numeric(0), ncol = 2)

  q_loss <- if (use_loss) Q_set$loss_pp else NULL
  l_loss <- if (use_loss) L_set$loss_pp else NULL

  for (i in seq_len(nq)) {
    for (j in seq_len(nl)) {
      D[i, j] <- combined_distance(
        Q_set$frag_pp[[i]], L_set$frag_pp[[j]],
        if (use_loss) q_loss[[i]] else empty,
        if (use_loss) l_loss[[j]] else empty,
        params
      )
    }
  }
  D
}

#' @keywords internal
#' @noRd
hellinger_normalise <- function(M) {
  rs <- Matrix::rowSums(M)
  rs[rs < .Machine$double.eps] <- 1
  Mn <- Matrix::Diagonal(x = 1 / rs) %*% M
  Mn@x <- sqrt(Mn@x)
  Mn
}

#' @keywords internal
#' @noRd
block_distance <- function(Q, L, qn = NULL, ln = NULL, method,
                           row_ids = NULL, col_ids = NULL) {
  if (method == "hellinger") {
    BC <- as.matrix(Matrix::tcrossprod(Q, L))
    if (!is.matrix(BC)) BC <- matrix(BC, nrow = nrow(Q), ncol = nrow(L))
    BC <- pmin(pmax(BC, 0), 1)
    inner <- 1 - BC
    inner[inner < 0] <- 0
    d <- sqrt(inner)
  } else {
    qn[qn < .Machine$double.eps] <- 1
    ln[ln < .Machine$double.eps] <- 1
    sim <- as.matrix(Matrix::tcrossprod(Q, L))
    if (!is.matrix(sim)) sim <- matrix(sim, nrow = nrow(Q), ncol = nrow(L))
    sim <- sim / outer(qn, ln)
    sim[!is.finite(sim)] <- 0
    sim <- pmin(pmax(sim, 0), 1)
    d <- 1 - sim
  }
  if (!is.matrix(d)) d <- matrix(d, nrow = nrow(Q), ncol = nrow(L))
  if (!is.null(row_ids)) rownames(d) <- row_ids
  if (!is.null(col_ids)) colnames(d) <- col_ids
  d
}
