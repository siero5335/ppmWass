#' Compute Similarity Matrices
#'
#' @param frag_list List of fragment spectra.
#' @param loss_list List of loss spectra.
#' @param ri Named vector of retention indices.
#' @param params Parameter list.
#' @param loss_typ_list Optional typical-loss spectra aligned to `frag_list`.
#' @param mref_conf Optional named numeric vector of Mref confidence values.
#' @param loss_anchor_list Optional anchored loss spectra for split-loss mode.
#' @param loss_pair_list Optional pairwise loss spectra for split-loss mode.
#' @param loss_anchor_typ_list Optional typical-loss projections for anchored losses.
#' @param loss_pair_typ_list Optional typical-loss projections for pairwise losses.
#' @return A list with distance and similarity matrices.
#' @export
compute_similarity_matrices <- function(frag_list, loss_list, ri, params,
                                 loss_typ_list = NULL, mref_conf = NULL,
                                 loss_anchor_list = NULL, loss_pair_list = NULL,
                                 loss_anchor_typ_list = NULL, loss_pair_typ_list = NULL) {
  dist_components <- NULL
  if (isTRUE(params$return_distance_components)) {
    dist_res <- compute_distance_matrices(
      frag_list, loss_list, params,
      loss_typ_list = loss_typ_list,
      mref_conf = mref_conf,
      loss_anchor_list = loss_anchor_list,
      loss_pair_list = loss_pair_list,
      loss_anchor_typ_list = loss_anchor_typ_list,
      loss_pair_typ_list = loss_pair_typ_list
    )
    dist_raw <- dist_res$dist
    dist_components <- dist_res$components
  } else {
    dist_raw <- compute_distance_matrix(
      frag_list, loss_list, params,
      loss_typ_list = loss_typ_list,
      mref_conf = mref_conf,
      loss_anchor_list = loss_anchor_list,
      loss_pair_list = loss_pair_list,
      loss_anchor_typ_list = loss_anchor_typ_list,
      loss_pair_typ_list = loss_pair_typ_list
    )
  }
  # Clamp distances to [0, 1] before converting to similarity
  dist_raw <- pmin(pmax(dist_raw, 0), 1)
  sim_raw <- 1 - dist_raw
  diag(sim_raw) <- 1

  ids <- names(frag_list)
  n <- length(ids)
  ri_vec <- ri[ids]

  dRI <- abs(outer(ri_vec, ri_vec, "-"))
  dRI_na <- is.na(dRI)

  if (params$ri_mode_cluster == "none") {
    sim_cluster <- sim_raw
  } else if (params$ri_mode_cluster == "soft") {
    pen <- exp(-(dRI / params$ri_sigma)^2)
    pen[dRI_na] <- 1
    sim_cluster <- sim_raw * pen
  } else if (params$ri_mode_cluster == "hard") {
    sim_cluster <- sim_raw
    sim_cluster[dRI > params$ri_window_hard & !dRI_na] <- 0
  } else {
    stop("Unknown ri_mode_cluster: ", params$ri_mode_cluster)
  }
  diag(sim_cluster) <- 1

  if (params$ri_mode_analog == "none") {
    sim_analog_base <- sim_raw
  } else if (params$ri_mode_analog == "soft") {
    pen <- exp(-(dRI / params$ri_sigma)^2)
    pen[dRI_na] <- 1
    sim_analog_base <- sim_raw * pen
  } else if (params$ri_mode_analog == "hard") {
    sim_analog_base <- sim_raw
    sim_analog_base[dRI > params$ri_window_hard & !dRI_na] <- 0
  } else {
    stop("Unknown ri_mode_analog: ", params$ri_mode_analog)
  }

  sim_analog <- sim_analog_base
  diag(sim_analog) <- 0
  k <- params$analog_k
  keep <- matrix(FALSE, n, n, dimnames = dimnames(sim_analog))
  for (i in seq_len(n)) {
    jj <- order(sim_analog[i, ], decreasing = TRUE)[1:min(k, n - 1)]
    keep[i, jj] <- TRUE
  }
  keep <- keep | t(keep)
  sim_analog[!keep] <- 0
  diag(sim_analog) <- 1

  list(
    dist_raw = dist_raw,
    dist_components = dist_components,
    sim_raw = sim_raw,
    sim_cluster = sim_cluster,
    sim_analog = sim_analog,
    dRI = dRI
  )
}


#' Library Search via Spectral Similarity
#'
#' Computes a query-vs-library distance matrix and returns top-N hits for each
#' query spectrum ranked by similarity.
#'
#' @param query_spectra Output of \code{build_spectra()} or
#'   \code{build_spectra_from_msp()}.
#' @param library_spectra Output of \code{build_spectra_from_msp()} or
#'   \code{build_spectra()}.
#' @param params Parameter list from \code{eihrms_default_params()}.
#' @param top_n Integer; number of top hits to return per query (default 10).
#' @param ri_tolerance Numeric or NULL; if set, library entries with
#'   \code{|RI_query - RI_lib| > ri_tolerance} are excluded before ranking.
#' @param as_list Logical; if TRUE return a named list keyed by query ID,
#'   otherwise return a single data.frame in long format (default FALSE).
#' @return A data.frame (long format) or a named list of data.frames.
#' @export
search_library <- function(query_spectra, library_spectra, params,
                           top_n = 10, ri_tolerance = NULL,
                           as_list = FALSE) {
  # --- Extract spectra lists ---
  q_frag <- query_spectra$frag_list
  q_loss <- query_spectra$loss_list
  l_frag <- library_spectra$frag_list
  l_loss <- library_spectra$loss_list

  if (length(q_frag) == 0 || length(l_frag) == 0) {
    message("Empty query or library spectra.")
    empty <- data.frame(
      query_id = character(), rank = integer(), library_id = character(),
      similarity = numeric(), delta_ri = numeric(),
      stringsAsFactors = FALSE
    )
    if (as_list) return(list())
    return(empty)
  }

  # --- Compute rectangular distance matrix ---
  dist_mat <- compute_distance_matrix_search(
    query_frag_list = q_frag,
    query_loss_list = q_loss,
    lib_frag_list = l_frag,
    lib_loss_list = l_loss,
    params = params,
    query_loss_typ_list = query_spectra$loss_typ_list,
    lib_loss_typ_list = library_spectra$loss_typ_list,
    query_mref_conf = query_spectra$mref_confidence,
    lib_mref_conf = library_spectra$mref_confidence,
    query_loss_anchor_list = query_spectra$loss_anchor_list,
    query_loss_pair_list = query_spectra$loss_pair_list,
    lib_loss_anchor_list = library_spectra$loss_anchor_list,
    lib_loss_pair_list = library_spectra$loss_pair_list,
    query_loss_anchor_typ_list = query_spectra$loss_anchor_typ_list,
    query_loss_pair_typ_list = query_spectra$loss_pair_typ_list,
    lib_loss_anchor_typ_list = library_spectra$loss_anchor_typ_list,
    lib_loss_pair_typ_list = library_spectra$loss_pair_typ_list
  )

  # --- Convert to similarity ---
  sim_mat <- 1 - pmin(pmax(dist_mat, 0), 1)
  sim_mat[!is.finite(dist_mat)] <- NA_real_

  # --- RI filter ---
  q_ri <- query_spectra$ri
  l_ri <- library_spectra$ri
  dRI <- NULL
  if (!is.null(q_ri) && !is.null(l_ri)) {
    q_ri_vec <- q_ri[rownames(sim_mat)]
    l_ri_vec <- l_ri[colnames(sim_mat)]
    dRI <- abs(outer(q_ri_vec, l_ri_vec, "-"))
    dimnames(dRI) <- dimnames(sim_mat)

    if (!is.null(ri_tolerance)) {
      mask <- !is.na(dRI) & dRI <= ri_tolerance
      sim_mat[!mask] <- NA
    }
  }

  # --- Extract top-N per query ---
  q_ids <- rownames(sim_mat)
  l_ids <- colnames(sim_mat)

  results_list <- lapply(seq_along(q_ids), function(i) {
    sims <- sim_mat[i, ]
    valid <- !is.na(sims)
    if (sum(valid) == 0) return(NULL)

    sims_valid <- sims[valid]
    ord <- order(sims_valid, decreasing = TRUE)
    n_take <- min(top_n, length(ord))
    top_idx <- ord[seq_len(n_take)]
    top_ids <- l_ids[valid][top_idx]

    df_hit <- data.frame(
      query_id = q_ids[i],
      rank = seq_len(n_take),
      library_id = top_ids,
      similarity = unname(sims_valid[top_idx]),
      stringsAsFactors = FALSE
    )

    # Add delta RI if available
    if (!is.null(dRI)) {
      ri_vals <- dRI[i, top_ids]
      df_hit$delta_ri <- unname(ri_vals)
    }

    df_hit
  })

  names(results_list) <- q_ids
  results_list <- results_list[!vapply(results_list, is.null, logical(1))]

  # --- Merge library metadata if available ---
  lib_df_spec <- library_spectra$df_spec
  has_meta <- !is.null(library_spectra$msp_metadata) &&
    !is.null(lib_df_spec) && "id" %in% names(lib_df_spec)
  if (has_meta) {
    # df_spec contains the retained records after RI filtering, deduplication,
    # and empty-spectrum removal. Raw msp_metadata still contains every input
    # record, so it cannot be paired with df_spec by position (or by name when
    # an earlier duplicate was filtered out).
    meta <- lib_df_spec
    if (!"name" %in% names(meta)) meta$name <- meta$id
    if ("RI" %in% names(meta)) meta$ri <- meta$RI
    meta_cols <- intersect(
      c("name", "cas", "inchikey", "formula", "mw", "ri"),
      colnames(meta)
    )
    results_list <- lapply(results_list, function(df_hit) {
      index <- match(df_hit$library_id, meta$id)
      data.frame(df_hit, meta[index, meta_cols, drop = FALSE],
                 row.names = NULL, check.names = FALSE)
    })
  }

  if (as_list) return(results_list)

  # Combine to long data.frame
  result_df <- do.call(rbind, results_list)
  if (is.null(result_df)) {
    result_df <- data.frame(
      query_id = character(), rank = integer(), library_id = character(),
      similarity = numeric(), stringsAsFactors = FALSE
    )
  }
  rownames(result_df) <- NULL

  # Reorder columns: query_id, rank, library_id, similarity, delta_ri, then metadata
  core_cols <- c("query_id", "rank", "library_id", "similarity")
  if ("delta_ri" %in% colnames(result_df)) core_cols <- c(core_cols, "delta_ri")
  other_cols <- setdiff(colnames(result_df), core_cols)
  result_df <- result_df[, c(core_cols, other_cols), drop = FALSE]

  result_df
}
