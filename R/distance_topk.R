#' Compute Top-K Sparse-Bin Library Candidates
#'
#' LC-scale retrieval API. Queries are processed in blocks; each block may
#' materialize a temporary `block_size x n_library` matrix, but only the top-k
#' candidates per query are retained.
#'
#' @param Q_set,L_set `SpectraSet` objects sharing the same bin grid.
#' @param params Parameter list. `params$distance_method` is the cheap stage-1
#'   method and must be one of `"cosine"`, `"weighted_cosine"`, or `"hellinger"`.
#' @param top_k Number of stage-1 candidates to keep per query.
#' @param expensive_method Optional method for rescoring stage-1 candidates with
#'   the legacy pair distance function, e.g. `"ppm_wasserstein"`.
#' @param expensive_top_k Number of candidates to keep after rescoring.
#' @param block_size Query block size.
#' @param use_loss Whether to include loss spectra. `NULL` uses loss when both
#'   sets contain loss spectra.
#' @param return_prefilter If `TRUE`, include stage-1 distances aligned to the
#'   returned candidates.
#' @param progress If `TRUE`, emit block progress messages.
#' @return A `TopKResult` object.
#' @export
compute_distance_topk_v2 <- function(Q_set, L_set, params,
                                     top_k = 200L,
                                     expensive_method = NULL,
                                     expensive_top_k = NULL,
                                     block_size = 1024L,
                                     use_loss = NULL,
                                     return_prefilter = FALSE,
                                     progress = TRUE) {
  check_common_param_typos(params)
  validate_compatible_grids(Q_set, L_set)
  if (Q_set$n < 1L) stop("Q_set must contain at least one spectrum.")
  if (L_set$n < 1L) stop("L_set must contain at least one spectrum.")

  guard <- guard_batch_compatible(params, check_method = TRUE)
  if (!guard$ok) {
    report_guard(guard, prefix = "[v2/Tier A] ")
    stop("compute_distance_topk_v2() requires a batchable stage-1 method.")
  }

  cheap_method <- params$distance_method %||% "cosine"
  has_loss <- !is.null(Q_set$loss_pp) && !is.null(L_set$loss_pp)
  if (is.null(use_loss)) use_loss <- has_loss

  top_k <- as.integer(min(top_k, L_set$n))
  if (top_k < 1L) stop("top_k must be >= 1.")
  expensive_top_k <- if (is.null(expensive_top_k)) top_k else as.integer(min(expensive_top_k, top_k))
  if (expensive_top_k < 1L) stop("expensive_top_k must be >= 1.")

  nq <- Q_set$n
  block_size <- as.integer(min(block_size, nq))
  if (block_size < 1L) stop("block_size must be >= 1.")

  pp <- list(
    mass_power = params$mass_power %||% 3,
    intensity_power = params$intensity_power %||% 0.5
  )
  weight <- switch(cheap_method,
    cosine = "raw",
    weighted_cosine = "wcos",
    hellinger = "raw",
    stop("Unsupported stage-1 method: ", cheap_method)
  )

  L_frag <- ensure_sparse_matrix(L_set, "frag", weight, pp)
  L_loss <- if (use_loss) ensure_sparse_matrix(L_set, "loss", weight, pp) else NULL
  if (cheap_method == "hellinger") {
    L_frag <- hellinger_normalise(L_frag)
    if (!is.null(L_loss)) L_loss <- hellinger_normalise(L_loss)
  }
  L_norm_frag <- if (cheap_method == "hellinger") NULL else sparse_row_norms(L_frag)
  L_norm_loss <- if (!is.null(L_loss) && cheap_method != "hellinger") sparse_row_norms(L_loss) else NULL

  topk_idx <- matrix(NA_integer_, nq, top_k)
  topk_dist <- matrix(NA_real_, nq, top_k)
  rownames(topk_idx) <- rownames(topk_dist) <- Q_set$ids

  block_starts <- seq.int(1L, nq, by = block_size)
  for (bs in block_starts) {
    be <- min(bs + block_size - 1L, nq)
    if (isTRUE(progress)) {
      message(sprintf("[v2/Tier A] block %d-%d / %d (stage1=%s)", bs, be, nq, cheap_method))
    }

    Q_frag <- ensure_block_matrix(Q_set, "frag", weight, pp, bs, be, cheap_method)
    d_frag <- block_distance(
      Q_frag, L_frag,
      if (cheap_method == "hellinger") NULL else sparse_row_norms(Q_frag),
      L_norm_frag,
      cheap_method
    )

    if (use_loss && !is.null(L_loss)) {
      Q_loss <- ensure_block_matrix(Q_set, "loss", weight, pp, bs, be, cheap_method)
      d_loss <- block_distance(
        Q_loss, L_loss,
        if (cheap_method == "hellinger") NULL else sparse_row_norms(Q_loss),
        L_norm_loss,
        cheap_method
      )
      d_block <- sqrt((params$w_frag %||% 0.7) * d_frag^2 + (params$w_loss %||% 0.3) * d_loss^2)
    } else {
      d_block <- d_frag
    }

    block_candidates <- select_topk_candidates(
      d_block,
      top_k = top_k,
      query_ids = Q_set$ids[bs:be],
      lib_ids = L_set$ids,
      exclude_self = FALSE
    )
    topk_idx[bs:be, ] <- block_candidates$candidate_idx
    topk_dist[bs:be, ] <- block_candidates$first_pass_distance
  }

  topk_dist_pre <- if (isTRUE(return_prefilter)) topk_dist else NULL
  refined <- FALSE

  if (!is.null(expensive_method)) {
    if (isTRUE(progress)) {
      message(sprintf("[v2/Tier A] refining top-%d candidates with %s", top_k, expensive_method))
    }

    params_exp <- params
    params_exp$distance_method <- expensive_method
    params_exp$tol_ppm <- params_exp$tol_ppm %||% 20
    params_exp$wasserstein_align <- params_exp$wasserstein_align %||% FALSE
    params_exp$w_frag <- params_exp$w_frag %||% 0.7
    params_exp$w_loss <- params_exp$w_loss %||% 0.3
    refined_dist <- refine_topk_candidates(
      query_frag_list = Q_set$frag_pp,
      query_loss_list = if (use_loss) Q_set$loss_pp else NULL,
      lib_frag_list = L_set$frag_pp,
      lib_loss_list = if (use_loss) L_set$loss_pp else NULL,
      candidate_idx = topk_idx,
      params = params_exp,
      use_parallel = params$use_parallel %||% FALSE,
      n_cores = params$n_cores %||% 1L,
      progress = progress,
      progress_prefix = "[v2/Tier A]",
      progress_step = block_size
    )

    refined_idx <- matrix(NA_integer_, nq, expensive_top_k)
    refined_top_dist <- matrix(NA_real_, nq, expensive_top_k)
    refined_pre_dist <- if (isTRUE(return_prefilter)) matrix(NA_real_, nq, expensive_top_k) else NULL

    for (i in seq_len(nq)) {
      ord <- top_k_order(refined_dist[i, ], expensive_top_k)
      refined_idx[i, ] <- topk_idx[i, ord]
      refined_top_dist[i, ] <- refined_dist[i, ord]
      if (isTRUE(return_prefilter)) refined_pre_dist[i, ] <- topk_dist_pre[i, ord]
    }

    topk_idx <- refined_idx
    topk_dist <- refined_top_dist
    topk_dist_pre <- refined_pre_dist
    rownames(topk_idx) <- rownames(topk_dist) <- Q_set$ids
    if (!is.null(topk_dist_pre)) rownames(topk_dist_pre) <- Q_set$ids
    refined <- TRUE
  }

  out <- list(
    query_ids = Q_set$ids,
    library_ids = L_set$ids,
    topk_idx = topk_idx,
    topk_dist = topk_dist,
    topk_dist_pre = topk_dist_pre,
    method = if (refined) expensive_method else cheap_method,
    prefilter = if (refined) cheap_method else NULL,
    refined = refined,
    block_size = block_size,
    call = match.call()
  )
  class(out) <- c("TopKResult", "list")
  out
}

#' @export
print.TopKResult <- function(x, ...) {
  cat("<TopKResult>\n")
  cat("  queries  :", length(x$query_ids), "\n")
  cat("  library  :", length(x$library_ids), "\n")
  cat("  top_k    :", ncol(x$topk_idx), "\n")
  cat("  method   :", x$method,
      if (!is.null(x$prefilter)) sprintf(" (prefilter: %s)", x$prefilter) else "", "\n")
  cat("  refined  :", x$refined, "\n")
  cat("  block_sz :", x$block_size, "\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
ensure_block_matrix <- function(set, side, weight, params, bs, be, method) {
  M <- ensure_sparse_matrix(set, side, weight, params)
  M_blk <- M[bs:be, , drop = FALSE]
  if (method == "hellinger") M_blk <- hellinger_normalise(M_blk)
  M_blk
}

#' Top-K Retrieval Accuracy
#'
#' Computes simple ID-based hit rates from a `TopKResult`.
#'
#' @param result A `TopKResult`.
#' @param top_k Integer vector of cutoffs.
#' @param ground_truth Optional named list or character vector of acceptable
#'   library IDs for each query. If `NULL`, query IDs themselves are treated as
#'   the expected IDs.
#' @param exclude_self If `TRUE`, candidates with the same ID as the query are
#'   excluded before scoring. This avoids trivial self-hits in diagnostics.
#' @return Named numeric vector of hit rates.
#' @export
topk_accuracy <- function(result, top_k = c(1, 3, 5, 10),
                          ground_truth = NULL,
                          exclude_self = FALSE) {
  qid <- result$query_ids
  lid <- result$library_ids
  nq <- length(qid)

  truth_for <- function(i) {
    if (is.null(ground_truth)) return(qid[i])
    if (is.list(ground_truth)) {
      if (!is.null(names(ground_truth))) return(ground_truth[[qid[i]]])
      return(ground_truth[[i]])
    }
    if (!is.null(names(ground_truth))) return(ground_truth[qid[i]])
    ground_truth[i]
  }

  out <- numeric(length(top_k))
  names(out) <- paste0("top_", top_k)
  for (ki in seq_along(top_k)) {
    K <- min(top_k[ki], ncol(result$topk_idx))
    hit <- 0L
    scored <- 0L
    for (i in seq_len(nq)) {
      cand_idx <- result$topk_idx[i, ]
      cand_idx <- cand_idx[!is.na(cand_idx)]
      if (isTRUE(exclude_self)) {
        cand_idx <- cand_idx[lid[cand_idx] != qid[i]]
      }
      cand_idx <- cand_idx[seq_len(min(K, length(cand_idx)))]
      expected <- truth_for(i)
      expected <- expected[!is.na(expected)]
      if (length(expected) == 0L || length(cand_idx) == 0L) next
      scored <- scored + 1L
      if (any(lid[cand_idx] %in% expected)) hit <- hit + 1L
    }
    out[ki] <- if (scored > 0L) hit / scored else NA_real_
  }
  out
}
