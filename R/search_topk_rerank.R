#' Search a Library with Top-K Candidate Reranking
#'
#' `search_topk_rerank()` is a high-level candidate-set reranking API. It first
#' computes a lightweight query-by-library first-pass distance matrix, keeps the
#' nearest `first_pass_top_k` finite candidates for each query, and reranks only
#' those candidates with a higher-accuracy method such as `"ppm_wasserstein"`.
#'
#' This is not exact all-vs-all reranking: library entries that do not enter the
#' first-pass candidate set are not evaluated by `rerank_method`. In v1 the
#' first-pass stage still materializes a full query-by-library matrix via
#' [compute_distance_matrix_search()]. Future LC-HRMS extensions may add
#' block-wise first-pass search, sparse-bin candidates, or external ANN
#' candidates.
#'
#' @param query_frag_list,lib_frag_list Named lists of two-column fragment
#'   spectra matrices.
#' @param query_loss_list,lib_loss_list Optional named lists of two-column
#'   neutral-loss spectra matrices. Both must be `NULL` for fragment-only
#'   search, or both must be non-`NULL` for fragment-plus-loss search.
#' @param params Parameter list. `params$distance_method` is overridden for the
#'   first-pass and rerank stages.
#' @param first_pass_method Distance method used to generate candidates.
#' @param rerank_method Distance method used to rerank candidates.
#' @param first_pass_top_k Number of first-pass candidates per query. Values
#'   larger than the library size are clamped.
#' @param final_top_k Number of reranked hits returned per query.
#' @param query_ids,lib_ids IDs aligned to the query and library lists.
#' @param exclude_self Whether to remove candidates with identical query/library
#'   IDs. `NULL` auto-detects this when the query and library ID vectors are
#'   identical.
#' @param first_pass_block_size Reserved for future block-wise first-pass
#'   search. Non-`NULL` values currently warn and fall back to full-matrix
#'   first-pass computation.
#' @param return_first_pass If `TRUE`, return the full query-by-library
#'   first-pass distance matrix and retain first-pass rank and distance columns
#'   in `results`. If `FALSE`, the matrix is returned as `NULL` and those
#'   columns are omitted.
#' @param use_parallel,n_cores Parallel settings passed to reranking and the
#'   pair-loop first-pass backend.
#' @param progress If `TRUE`, emit progress messages.
#' @param query_loss_typ_list,lib_loss_typ_list Optional query and library
#'   typical-loss spectra passed through to compatible distance paths.
#' @param query_mref_conf,lib_mref_conf Optional query and library Mref
#'   confidence vectors passed through to compatible distance paths.
#' @param ... Reserved optional arguments. Unsupported names are warned and
#'   ignored in v1.
#' @return A list with long-form `results`, per-query `candidate_summary`, an
#'   aligned `first_pass_distance_matrix` (or `NULL` when
#'   `return_first_pass = FALSE`), stage-specific parameter lists, and timing
#'   information. Matrix rows and columns follow `query_ids` and `lib_ids`,
#'   respectively. When `exclude_self = TRUE`, self-excluded entries are not
#'   counted in `candidate_summary$n_finite_first_pass`.
#' @examples
#' params <- eihrms_default_params()
#' params$distance_method <- "hellinger"
#' frag <- list(
#'   a = matrix(c(100, 1), ncol = 2, byrow = TRUE),
#'   b = matrix(c(101, 1), ncol = 2, byrow = TRUE)
#' )
#' search_topk_rerank(
#'   query_frag_list = frag,
#'   lib_frag_list = frag,
#'   params = params,
#'   first_pass_method = "hellinger",
#'   rerank_method = "hellinger",
#'   first_pass_top_k = 2L,
#'   final_top_k = 1L,
#'   exclude_self = TRUE,
#'   progress = FALSE
#' )
#'
#' # SI-style reranking simulations can sweep first_pass_top_k over
#' # c(5, 10, 20, 50, 100), often with final_top_k = 1L.
#' @export
search_topk_rerank <- function(
    query_frag_list,
    query_loss_list = NULL,
    lib_frag_list,
    lib_loss_list = NULL,
    params,
    first_pass_method = "entropy_weighted",
    rerank_method = "ppm_wasserstein",
    first_pass_top_k = 50L,
    final_top_k = 10L,
    query_ids = names(query_frag_list),
    lib_ids = names(lib_frag_list),
    exclude_self = NULL,
    first_pass_block_size = NULL,
    return_first_pass = TRUE,
    use_parallel = params$use_parallel %||% FALSE,
    n_cores = params$n_cores %||% 1L,
    progress = TRUE,
    query_loss_typ_list = NULL,
    lib_loss_typ_list = NULL,
    query_mref_conf = NULL,
    lib_mref_conf = NULL,
    ...) {
  total_start <- proc.time()[["elapsed"]]
  dots <- list(...)
  if (length(dots) > 0L) {
    warning(
      "Unsupported optional arguments ignored in search_topk_rerank(): ",
      paste(names(dots), collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.list(params)) stop("params must be a list.")
  check_common_param_typos(params)
  validate_distance_method_name(first_pass_method, "first_pass_method")
  validate_distance_method_name(rerank_method, "rerank_method")

  nq <- length(query_frag_list)
  nl <- length(lib_frag_list)
  if (nq < 1L) stop("query_frag_list must contain at least one spectrum.")
  if (nl < 1L) stop("lib_frag_list must contain at least one spectrum.")

  if (is.null(query_ids)) query_ids <- as.character(seq_len(nq))
  if (is.null(lib_ids)) lib_ids <- as.character(seq_len(nl))
  query_ids <- as.character(query_ids)
  lib_ids <- as.character(lib_ids)
  if (length(query_ids) != nq) stop("query_ids length must match query_frag_list length.")
  if (length(lib_ids) != nl) stop("lib_ids length must match lib_frag_list length.")
  names(query_frag_list) <- query_ids
  names(lib_frag_list) <- lib_ids

  first_pass_top_k <- as.integer(first_pass_top_k)
  final_top_k <- as.integer(final_top_k)
  if (length(first_pass_top_k) != 1L || is.na(first_pass_top_k) || first_pass_top_k < 1L) {
    stop("first_pass_top_k must be >= 1.")
  }
  if (length(final_top_k) != 1L || is.na(final_top_k) || final_top_k < 1L) {
    stop("final_top_k must be >= 1.")
  }
  if (final_top_k > first_pass_top_k) {
    stop("final_top_k must be <= first_pass_top_k.")
  }
  if (!is.logical(return_first_pass) || length(return_first_pass) != 1L ||
      is.na(return_first_pass)) {
    stop("return_first_pass must be TRUE or FALSE.")
  }
  if (!is.null(first_pass_block_size)) {
    warning(
      "first_pass_block_size is reserved for future block-wise first-pass search; ",
      "the current implementation computes the full first-pass distance matrix.",
      call. = FALSE
    )
  }

  if (first_pass_top_k > nl) {
    message("first_pass_top_k > n_library; clamping to n_library.")
    first_pass_top_k <- nl
  }

  loss_null <- c(is.null(query_loss_list), is.null(lib_loss_list))
  if (any(loss_null) && !all(loss_null)) {
    stop("query_loss_list and lib_loss_list must both be NULL or both be non-NULL.")
  }

  fragment_only <- all(loss_null)
  empty_loss <- matrix(numeric(0), ncol = 2)
  if (fragment_only) {
    query_loss_list <- rep(list(empty_loss), nq)
    lib_loss_list <- rep(list(empty_loss), nl)
  } else {
    if (length(query_loss_list) != nq) stop("query_loss_list length must match query_frag_list length.")
    if (length(lib_loss_list) != nl) stop("lib_loss_list length must match lib_frag_list length.")
    if (isTRUE(params$use_split_loss)) {
      stop(
        "search_topk_rerank() v1 does not support split-loss reranking. ",
        "Set params$use_split_loss = FALSE or use a lower-level split-loss distance API."
      )
    }
  }
  names(query_loss_list) <- query_ids
  names(lib_loss_list) <- lib_ids

  exclude_self <- resolve_exclude_self(query_ids, lib_ids, exclude_self)

  normalize_stage_params <- function(params_in, method) {
    out <- params_in
    out$distance_method <- method
    out$tol_ppm <- out$tol_ppm %||% 20
    out$wasserstein_align <- out$wasserstein_align %||% FALSE
    out$wasserstein_transition_mult <- out$wasserstein_transition_mult %||% 3
    out$ot_method <- out$ot_method %||% "exact"
    out$sinkhorn_epsilon <- out$sinkhorn_epsilon %||% 0.05
    out$sinkhorn_niter <- out$sinkhorn_niter %||% 100L
    out$mass_power <- out$mass_power %||% 3
    out$intensity_power <- out$intensity_power %||% 0.5
    out$use_parallel <- isTRUE(use_parallel)
    out$n_cores <- as.integer(n_cores %||% 1L)
    if (fragment_only) {
      out$w_frag <- 1
      out$w_loss <- 0
      out$use_typical_loss <- FALSE
      out$use_split_loss <- FALSE
      out$use_mref_confidence <- FALSE
    } else {
      out$w_frag <- out$w_frag %||% 0.7
      out$w_loss <- out$w_loss %||% 0.3
    }
    out
  }

  params_first <- normalize_stage_params(params, first_pass_method)
  params_rerank <- normalize_stage_params(params, rerank_method)

  first_start <- proc.time()[["elapsed"]]
  D_first <- compute_distance_matrix_search(
    query_frag_list = query_frag_list,
    query_loss_list = query_loss_list,
    lib_frag_list = lib_frag_list,
    lib_loss_list = lib_loss_list,
    params = params_first,
    progress = progress,
    query_loss_typ_list = if (fragment_only) NULL else query_loss_typ_list,
    lib_loss_typ_list = if (fragment_only) NULL else lib_loss_typ_list,
    query_mref_conf = if (fragment_only) NULL else query_mref_conf,
    lib_mref_conf = if (fragment_only) NULL else lib_mref_conf
  )
  if (!is.matrix(D_first) || !is.numeric(D_first) ||
      !identical(dim(D_first), c(nq, nl))) {
    stop(
      "First-pass search must return a numeric query-by-library matrix with ",
      nq, " row(s) and ", nl, " column(s)."
    )
  }
  if ((!is.null(rownames(D_first)) &&
       !identical(as.character(rownames(D_first)), query_ids)) ||
      (!is.null(colnames(D_first)) &&
       !identical(as.character(colnames(D_first)), lib_ids))) {
    stop("First-pass search returned distance-matrix IDs out of input order.")
  }
  dimnames(D_first) <- list(query_ids, lib_ids)
  first_pass_sec <- proc.time()[["elapsed"]] - first_start

  candidate_start <- proc.time()[["elapsed"]]
  candidates <- select_topk_candidates(
    D_first,
    top_k = first_pass_top_k,
    query_ids = query_ids,
    lib_ids = lib_ids,
    exclude_self = exclude_self
  )
  candidate_extraction_sec <- proc.time()[["elapsed"]] - candidate_start

  rerank_start <- proc.time()[["elapsed"]]
  if (identical(first_pass_method, rerank_method)) {
    rerank_distance <- candidates$first_pass_distance
  } else {
    rerank_distance <- refine_topk_candidates(
      query_frag_list = query_frag_list,
      query_loss_list = if (fragment_only) NULL else query_loss_list,
      lib_frag_list = lib_frag_list,
      lib_loss_list = if (fragment_only) NULL else lib_loss_list,
      candidate_idx = candidates$candidate_idx,
      params = params_rerank,
      use_parallel = use_parallel,
      n_cores = n_cores,
      progress = progress,
      progress_prefix = "[top-k rerank]",
      query_loss_typ_list = if (fragment_only) NULL else query_loss_typ_list,
      lib_loss_typ_list = if (fragment_only) NULL else lib_loss_typ_list,
      query_mref_conf = if (fragment_only) NULL else query_mref_conf,
      lib_mref_conf = if (fragment_only) NULL else lib_mref_conf
    )
  }
  rerank_sec <- proc.time()[["elapsed"]] - rerank_start

  format_start <- proc.time()[["elapsed"]]
  results <- format_topk_long(
    candidate_idx = candidates$candidate_idx,
    rerank_distance = rerank_distance,
    first_pass_rank = candidates$first_pass_rank,
    first_pass_distance = candidates$first_pass_distance,
    query_ids = query_ids,
    lib_ids = lib_ids,
    final_top_k = final_top_k,
    first_pass_method = first_pass_method,
    rerank_method = rerank_method,
    return_first_pass = return_first_pass
  )

  returned <- tabulate(match(results$query_id, query_ids), nbins = nq)
  candidate_summary <- data.frame(
    query_id = query_ids,
    n_candidates = rowSums(!is.na(candidates$candidate_idx)),
    n_finite_first_pass = candidates$n_finite_first_pass,
    n_returned = returned,
    stringsAsFactors = FALSE
  )

  short <- candidate_summary$n_returned < final_top_k
  if (any(short)) {
    warning(
      "Fewer than final_top_k finite candidates returned for ",
      sum(short), " query(s).",
      call. = FALSE
    )
  }
  result_formatting_sec <- proc.time()[["elapsed"]] - format_start

  list(
    results = results,
    candidate_summary = candidate_summary,
    first_pass_distance_matrix = if (return_first_pass) D_first else NULL,
    params_first_pass = params_first,
    params_rerank = params_rerank,
    timing = list(
      first_pass_sec = first_pass_sec,
      topk_extraction_sec = candidate_extraction_sec,
      rerank_sec = rerank_sec,
      result_formatting_sec = result_formatting_sec,
      total_sec = proc.time()[["elapsed"]] - total_start
    )
  )
}
