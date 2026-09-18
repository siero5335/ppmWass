#' @keywords internal
#' @noRd
resolve_library_inchikeys <- function(spectra_obj) {
  if (!is.null(spectra_obj$df_spec) && "inchikey" %in% names(spectra_obj$df_spec)) {
    return(spectra_obj$df_spec$inchikey)
  }
  if (!is.null(spectra_obj$msp_metadata) && "inchikey" %in% names(spectra_obj$msp_metadata)) {
    return(spectra_obj$msp_metadata$inchikey)
  }
  NULL
}

#' @keywords internal
#' @noRd
evaluate_library_search_retrieval <- function(dist_mat,
                                              query_inchikeys,
                                              library_inchikeys,
                                              top_k = c(1, 5, 10, 20),
                                              inchikey_chars = 14) {
  q_prefix <- normalize_inchikey_prefix(query_inchikeys, inchikey_chars = inchikey_chars)
  l_prefix <- normalize_inchikey_prefix(library_inchikeys, inchikey_chars = inchikey_chars)

  valid_query <- !is.na(q_prefix)
  has_match <- vapply(seq_along(q_prefix), function(i) {
    valid_query[i] && any(l_prefix == q_prefix[i], na.rm = TRUE)
  }, logical(1))
  query_idx <- which(has_match)

  top_k <- sort(unique(as.integer(top_k)))
  metric_names <- paste0("P@", top_k)

  if (length(query_idx) == 0) {
    out <- as.list(stats::setNames(rep(NA_real_, length(metric_names)), metric_names))
    out$map <- NA_real_
    out$mrr <- NA_real_
    out$n_queries <- 0L
    return(out)
  }

  p_at_k <- numeric(length(top_k))
  ap <- numeric(length(query_idx))
  rr <- numeric(length(query_idx))

  for (ii in seq_along(query_idx)) {
    i <- query_idx[ii]
    candidates <- which(is.finite(dist_mat[i, ]))
    if (!length(candidates)) next
    ranked <- candidates[order(dist_mat[i, candidates])]
    relevant <- l_prefix[ranked] == q_prefix[i]
    relevant[is.na(relevant)] <- FALSE

    if (!any(relevant)) next

    cum_rel <- cumsum(relevant)
    precision_at_rank <- cum_rel / seq_along(relevant)
    rr[ii] <- 1 / which(relevant)[1]
    # Excluded relevant entries remain misses in the fixed-library estimand.
    n_relevant <- sum(l_prefix == q_prefix[i], na.rm = TRUE)
    ap[ii] <- sum(precision_at_rank[relevant]) / n_relevant

    for (k_idx in seq_along(top_k)) {
      k <- min(top_k[k_idx], length(relevant))
      p_at_k[k_idx] <- p_at_k[k_idx] + sum(relevant[seq_len(k)]) / k
    }
  }

  out <- as.list(p_at_k / length(query_idx))
  names(out) <- metric_names
  out$map <- mean(ap)
  out$mrr <- mean(rr)
  out$n_queries <- length(query_idx)
  out
}

#' Benchmark Query-vs-Library Search Across Methods
#'
#' Runs `compute_distance_matrix_search()` for one or more methods and evaluates
#' retrieval quality against library entries sharing the same InChIKey prefix.
#' Nonfinite distances and RI-excluded entries are not ranked. Queries with a
#' relevant entry in the original library remain in the evaluation denominator;
#' if no candidate survives, their metrics are zero. Average precision uses the
#' original number of relevant library entries, and Precision@K divides by the
#' smaller of K and the number of surviving candidates. With an RI tolerance,
#' pairs with missing RI values are excluded, as in `search_library()`.
#'
#' @param query_spectra Output from `build_spectra()` or `build_spectra_from_msp()`.
#' @param library_spectra Output from `build_spectra()` or `build_spectra_from_msp()`.
#' @param methods Character vector of distance methods to compare.
#' @param params Base parameter list.
#' @param top_k Integer vector of cutoffs for Precision@K.
#' @param inchikey_chars Number of InChIKey prefix characters used for matching.
#' @param ri_tolerance Optional RI tolerance passed to candidate filtering after distance calculation.
#' @param keep_search_results If `TRUE`, also return top-ranked hit tables from `search_library()`.
#' @param search_top_n Number of hits to export per query when `keep_search_results = TRUE`.
#' @param progress If `TRUE`, print progress messages.
#' @return A list with per-method results and a `summary` data frame.
#' @export
benchmark_library_search <- function(query_spectra,
                                     library_spectra,
                                     methods = c("cosine", "entropy", "wasserstein"),
                                     params = eihrms_default_params(),
                                     top_k = c(1, 5, 10, 20),
                                     inchikey_chars = 14,
                                     ri_tolerance = NULL,
                                     keep_search_results = TRUE,
                                     search_top_n = 10,
                                     progress = TRUE) {
  query_inchikeys <- resolve_library_inchikeys(query_spectra)
  library_inchikeys <- resolve_library_inchikeys(library_spectra)
  if (is.null(query_inchikeys) || is.null(library_inchikeys)) {
    stop("Both query_spectra and library_spectra must include InChIKey metadata for retrieval benchmarking.")
  }

  results <- list()
  summary_rows <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    method <- methods[i]
    params_i <- params
    params_i$distance_method <- method
    params_i <- validate_params(params_i)

    if (isTRUE(progress)) {
      message("Library benchmark: ", method, " (", i, "/", length(methods), ")")
    }

    dist_mat <- compute_distance_matrix_search(
      query_frag_list = query_spectra$frag_list,
      query_loss_list = query_spectra$loss_list,
      lib_frag_list = library_spectra$frag_list,
      lib_loss_list = library_spectra$loss_list,
      params = params_i,
      progress = FALSE,
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

    if (!is.null(ri_tolerance) && !is.null(query_spectra$ri) && !is.null(library_spectra$ri)) {
      q_ri <- query_spectra$ri[rownames(dist_mat)]
      l_ri <- library_spectra$ri[colnames(dist_mat)]
      dRI <- abs(outer(q_ri, l_ri, "-"))
      eligible <- !is.na(dRI) & dRI <= ri_tolerance
      dist_mat[!eligible] <- Inf
    }

    metrics <- evaluate_library_search_retrieval(
      dist_mat = dist_mat,
      query_inchikeys = query_inchikeys,
      library_inchikeys = library_inchikeys,
      top_k = top_k,
      inchikey_chars = inchikey_chars
    )

    summary_rows[[i]] <- data.frame(
      method = method,
      n_queries = metrics$n_queries,
      map = metrics$map,
      mrr = metrics$mrr,
      matrix(unlist(metrics[paste0("P@", sort(unique(as.integer(top_k))))]),
             nrow = 1,
             dimnames = list(NULL, paste0("p_at_", sort(unique(as.integer(top_k))))))
    )

    results[[method]] <- list(
      method = method,
      dist_matrix = dist_mat,
      metrics = metrics
    )

    if (isTRUE(keep_search_results)) {
      results[[method]]$hits <- search_library(
        query_spectra = query_spectra,
        library_spectra = library_spectra,
        params = params_i,
        top_n = search_top_n,
        ri_tolerance = ri_tolerance,
        as_list = FALSE
      )
    }
  }

  results$summary <- do.call(rbind, summary_rows)
  rownames(results$summary) <- NULL
  results
}
