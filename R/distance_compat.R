#' Compute Query-vs-Library Distance Matrix
#'
#' Computes a rectangular query-by-library distance matrix. By default this uses
#' the legacy exact pair-loop backend. Set `params$backend = "sparse_bins"` to
#' opt into the experimental sparse-bin backend for batchable methods.
#'
#' The sparse-bin backend is approximate and currently supports only `cosine`,
#' `weighted_cosine`, and `hellinger` without typical-loss, split-loss, or
#' Mref-confidence options. Unsupported sparse-bin configurations fall back to
#' the legacy backend with a message instead of silently ignoring inputs.
#'
#' @param query_frag_list Named list of query fragment spectra.
#' @param query_loss_list Named list of query loss spectra.
#' @param lib_frag_list Named list of library fragment spectra.
#' @param lib_loss_list Named list of library loss spectra.
#' @param params Parameter list from [eihrms_default_params()].
#' @param progress If `TRUE`, emit progress messages.
#' @param query_loss_typ_list,lib_loss_typ_list Optional typical-loss spectra.
#' @param query_mref_conf,lib_mref_conf Optional Mref confidence vectors.
#' @param query_loss_anchor_list,query_loss_pair_list Optional query split-loss spectra.
#' @param lib_loss_anchor_list,lib_loss_pair_list Optional library split-loss spectra.
#' @param query_loss_anchor_typ_list,query_loss_pair_typ_list Optional query typical split-loss spectra.
#' @param lib_loss_anchor_typ_list,lib_loss_pair_typ_list Optional library typical split-loss spectra.
#' @return A numeric matrix of dimension `nquery x nlibrary`.
#' @export
compute_distance_matrix_search <- function(query_frag_list, query_loss_list,
                                           lib_frag_list, lib_loss_list,
                                           params, progress = TRUE,
                                           query_loss_typ_list = NULL,
                                           lib_loss_typ_list = NULL,
                                           query_mref_conf = NULL,
                                           lib_mref_conf = NULL,
                                           query_loss_anchor_list = NULL,
                                           query_loss_pair_list = NULL,
                                           lib_loss_anchor_list = NULL,
                                           lib_loss_pair_list = NULL,
                                           query_loss_anchor_typ_list = NULL,
                                           query_loss_pair_typ_list = NULL,
                                           lib_loss_anchor_typ_list = NULL,
                                           lib_loss_pair_typ_list = NULL) {
  check_common_param_typos(params)
  backend <- params$backend %||% "pair_loop"

  call_pairloop <- function() {
    compute_distance_matrix_search_pairloop(
      query_frag_list = query_frag_list,
      query_loss_list = query_loss_list,
      lib_frag_list = lib_frag_list,
      lib_loss_list = lib_loss_list,
      params = params,
      progress = progress,
      query_loss_typ_list = query_loss_typ_list,
      lib_loss_typ_list = lib_loss_typ_list,
      query_mref_conf = query_mref_conf,
      lib_mref_conf = lib_mref_conf,
      query_loss_anchor_list = query_loss_anchor_list,
      query_loss_pair_list = query_loss_pair_list,
      lib_loss_anchor_list = lib_loss_anchor_list,
      lib_loss_pair_list = lib_loss_pair_list,
      query_loss_anchor_typ_list = query_loss_anchor_typ_list,
      query_loss_pair_typ_list = query_loss_pair_typ_list,
      lib_loss_anchor_typ_list = lib_loss_anchor_typ_list,
      lib_loss_pair_typ_list = lib_loss_pair_typ_list
    )
  }

  if (identical(backend, "pair_loop")) return(call_pairloop())
  if (!identical(backend, "sparse_bins")) {
    stop("Unknown backend: ", backend, ". Use 'pair_loop' or 'sparse_bins'.")
  }

  guard <- guard_batch_compatible(params, check_method = TRUE)
  if (!guard$ok) {
    if (isTRUE(progress)) {
      report_guard(guard, prefix = "[v2 shim] ")
      message("[v2 shim] falling back to pair_loop backend.")
    }
    return(call_pairloop())
  }

  bin_ppm <- params$bin_ppm %||% 2
  smear_ppm <- params$smear_ppm %||% (params$tol_ppm %||% 15)
  smear_kind <- params$smear_kind %||% "tri"
  mz_min <- params$min_mz %||% 35
  mz_max <- params$max_mz %||% 650

  Q <- as_spectra_set(
    query_frag_list, query_loss_list,
    ids = names(query_frag_list),
    mz_min = mz_min, mz_max = mz_max,
    bin_ppm = bin_ppm, smear_ppm = smear_ppm,
    smear_kind = smear_kind
  )
  L <- as_spectra_set(
    lib_frag_list, lib_loss_list,
    ids = names(lib_frag_list),
    mz_min = mz_min, mz_max = mz_max,
    bin_ppm = bin_ppm, smear_ppm = smear_ppm,
    smear_kind = smear_kind
  )

  compute_distance_matrix_search_v2(
    Q, L, params,
    max_dense_cells = params$max_dense_cells %||% 1e7,
    progress = progress
  )
}
