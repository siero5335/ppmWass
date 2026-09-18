#' Compute Distance Matrix
#'
#' @param frag_list List of fragment spectra.
#' @param loss_list List of loss spectra.
#' @param params Parameter list.
#' @param progress If TRUE, print progress messages.
#' @param loss_typ_list Optional typical-loss spectra aligned to `frag_list`.
#' @param mref_conf Optional named numeric vector of Mref confidence values.
#' @param loss_anchor_list Optional anchored loss spectra for split-loss mode.
#' @param loss_pair_list Optional pairwise loss spectra for split-loss mode.
#' @param loss_anchor_typ_list Optional typical-loss projections for anchored losses.
#' @param loss_pair_typ_list Optional typical-loss projections for pairwise losses.
#' @return A square query-by-library distance matrix. Symmetric methods are
#'   computed by one triangle and mirroring; directional methods are evaluated
#'   in explicit row-query to column-library orientation.
#' @export
compute_distance_matrix <- function(frag_list, loss_list, params, progress = TRUE,
                                    loss_typ_list = NULL, mref_conf = NULL,
                                    loss_anchor_list = NULL, loss_pair_list = NULL,
                                    loss_anchor_typ_list = NULL, loss_pair_typ_list = NULL) {
  ids <- names(frag_list)
  n <- length(ids)

  if (isTRUE(progress)) {
    message("Computing distance matrix for ", n, " compounds...")
    message("Method: ", params$distance_method)
    message("Weights: frag = ", params$w_frag, ", loss = ", params$w_loss)
    if (isTRUE(params$use_split_loss)) {
      message("Split loss: anchored weight = ", if (is.null(params$loss_anchor_weight)) 0.6 else params$loss_anchor_weight)
    }
  }

  dist_mat <- matrix(0, n, n, dimnames = list(ids, ids))
  if (n < 2) return(dist_mat)

  if (!is.null(loss_list) && !is.null(names(loss_list))) loss_list <- loss_list[ids]
  if (!is.null(loss_typ_list) && !is.null(names(loss_typ_list))) loss_typ_list <- loss_typ_list[ids]

  if (!is.null(loss_anchor_list) && !is.null(names(loss_anchor_list))) loss_anchor_list <- loss_anchor_list[ids]
  if (!is.null(loss_pair_list) && !is.null(names(loss_pair_list))) loss_pair_list <- loss_pair_list[ids]
  if (!is.null(loss_anchor_typ_list) && !is.null(names(loss_anchor_typ_list))) loss_anchor_typ_list <- loss_anchor_typ_list[ids]
  if (!is.null(loss_pair_typ_list) && !is.null(names(loss_pair_typ_list))) loss_pair_typ_list <- loss_pair_typ_list[ids]

  use_typical_loss_combined <- isTRUE(params$use_typical_loss) && !is.null(loss_typ_list)
  use_split_loss <- isTRUE(params$use_split_loss) &&
    !is.null(loss_anchor_list) && !is.null(loss_pair_list)
  use_typical_loss_split <- isTRUE(params$use_typical_loss) &&
    !is.null(loss_anchor_typ_list) && !is.null(loss_pair_typ_list)

  use_conf <- isTRUE(params$use_mref_confidence) && !is.null(mref_conf)

  # Stein-style composite is directional because its ratio term depends on the
  # number of positive peaks in the first (query) argument. Never mirror a
  # directional score: evaluate every row-query/column-library cell explicitly.
  if (!distance_method_is_symmetric(params$distance_method, params$ot_method)) {
    if (isTRUE(progress)) {
      message("Directional method: computing the full query-by-library matrix")
    }
    directional_dist <- compute_distance_matrix_search(
      query_frag_list = frag_list,
      query_loss_list = loss_list,
      lib_frag_list = frag_list,
      lib_loss_list = loss_list,
      params = params,
      progress = progress,
      query_loss_typ_list = loss_typ_list,
      lib_loss_typ_list = loss_typ_list,
      query_mref_conf = mref_conf,
      lib_mref_conf = mref_conf,
      query_loss_anchor_list = loss_anchor_list,
      query_loss_pair_list = loss_pair_list,
      lib_loss_anchor_list = loss_anchor_list,
      lib_loss_pair_list = loss_pair_list,
      query_loss_anchor_typ_list = loss_anchor_typ_list,
      query_loss_pair_typ_list = loss_pair_typ_list,
      lib_loss_anchor_typ_list = loss_anchor_typ_list,
      lib_loss_pair_typ_list = loss_pair_typ_list
    )
    # Preserve the square distance-matrix contract even when an empty channel
    # makes the underlying pairwise implementation return a nonzero self-score.
    diag(directional_dist) <- 0
    return(directional_dist)
  }

  if (can_use_prepared_pairloop(params) &&
      !use_typical_loss_combined && !use_split_loss &&
      !use_typical_loss_split && !use_conf) {
    return(compute_distance_matrix_prepared(
      frag_list, loss_list, params, progress = progress
    ))
  }

  conf_vec <- NULL
  if (use_conf) {
    conf_vec <- mref_conf
    if (!is.null(names(conf_vec))) {
      conf_vec <- conf_vec[ids]
    }
    conf_vec <- as.numeric(conf_vec)
  }

  compute_row <- function(i) {
    d <- numeric(n - i)

    confA <- if (use_conf) conf_vec[i] else NULL

    if (use_split_loss) {
      lossA_anchor <- loss_anchor_list[[i]]
      lossA_pair <- loss_pair_list[[i]]
      lossA_anchor_typ <- if (use_typical_loss_split) loss_anchor_typ_list[[i]] else NULL
      lossA_pair_typ <- if (use_typical_loss_split) loss_pair_typ_list[[i]] else NULL

      for (j in (i + 1):n) {
        confB <- if (use_conf) conf_vec[j] else NULL
        d[j - i] <- combined_distance_split(
          frag_list[[i]], frag_list[[j]],
          lossA_anchor, loss_anchor_list[[j]],
          lossA_pair, loss_pair_list[[j]],
          params,
          lossA_anchor_typ = lossA_anchor_typ,
          lossB_anchor_typ = if (use_typical_loss_split) loss_anchor_typ_list[[j]] else NULL,
          lossA_pair_typ = lossA_pair_typ,
          lossB_pair_typ = if (use_typical_loss_split) loss_pair_typ_list[[j]] else NULL,
          confA = confA, confB = confB
        )
      }
    } else {
      lossA_typ <- if (use_typical_loss_combined) loss_typ_list[[i]] else NULL

      for (j in (i + 1):n) {
        confB <- if (use_conf) conf_vec[j] else NULL
        lossB_typ <- if (use_typical_loss_combined) loss_typ_list[[j]] else NULL

        d[j - i] <- combined_distance(
          frag_list[[i]], frag_list[[j]],
          loss_list[[i]], loss_list[[j]], params,
          lossA_typ = lossA_typ, lossB_typ = lossB_typ,
          confA = confA, confB = confB
        )
      }
    }
    list(i = i, d = d)
  }

  if (!params$use_parallel || params$n_cores <= 1) {
    for (i in 1:(n - 1)) {
      res <- compute_row(i)
      dist_mat[i, (i + 1):n] <- res$d
      dist_mat[(i + 1):n, i] <- res$d
      if (isTRUE(progress) && i %% 10 == 0) message("Progress: ", i, " / ", n - 1)
    }
    return(dist_mat)
  }

  cores <- resolve_parallel_cores(params$n_cores)
  idx <- 1:(n - 1)

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c(
        "frag_list", "loss_list",
        "loss_typ_list",
        "loss_anchor_list", "loss_pair_list",
        "loss_anchor_typ_list", "loss_pair_typ_list",
        "conf_vec", "params",
        "combined_distance", "combined_distance_split",
        "compute_distance", "hellinger_distance", "cosine_distance",
        "entropy_distance", "entropy_unweighted_distance",
        "weighted_cosine_distance", "composite_distance",
        "wasserstein_distance", "ppm_wasserstein_distance",
        "PPMWASS_OT_MARGINAL_TOLERANCE", "assess_ot_plan",
        "append_ot_diagnostic", "empty_ot_attempts", "ot_attempt_row",
        "align_spectra", "is_near_zero"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      if (params$distance_method == "ppm_wasserstein" &&
          requireNamespace("approxOT", quietly = TRUE)) {
        library(approxOT)
      }
      if (params$distance_method %in% c("wasserstein", "ppm_wasserstein") &&
          requireNamespace("transport", quietly = TRUE)) {
        library(transport)
      }
      if (params$distance_method %in% c("entropy", "entropy_weighted", "entropy_unweighted") &&
          requireNamespace("msentropy", quietly = TRUE)) {
        library(msentropy)
      }
      NULL
    })
    results <- parallel::parLapply(cl, idx, compute_row)
  } else {
    results <- parallel::mclapply(idx, compute_row, mc.cores = cores)
  }

  for (res in results) {
    i <- res$i
    dist_mat[i, (i + 1):n] <- res$d
    dist_mat[(i + 1):n, i] <- res$d
  }

  dist_mat
}

#' Legacy Query-vs-Library Distance Matrix
#'
#' Pair-loop implementation retained for exact compatibility and as the
#' fallback backend behind \code{compute_distance_matrix_search()}.
#' @keywords internal
#' @noRd
compute_distance_matrix_search_pairloop <- function(query_frag_list, query_loss_list,
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
  q_ids <- names(query_frag_list)
  l_ids <- names(lib_frag_list)
  nq <- length(q_ids)
  nl <- length(l_ids)

  if (isTRUE(progress)) {
    message("Library search: ", nq, " queries x ", nl, " library entries")
    message("Method: ", params$distance_method)
  }

  dist_mat <- matrix(0, nq, nl, dimnames = list(q_ids, l_ids))
  if (nq == 0 || nl == 0) return(dist_mat)

  use_typical_loss_combined <- isTRUE(params$use_typical_loss) &&
    !is.null(query_loss_typ_list) && !is.null(lib_loss_typ_list)
  use_split_loss <- isTRUE(params$use_split_loss) &&
    !is.null(query_loss_anchor_list) && !is.null(lib_loss_anchor_list) &&
    !is.null(query_loss_pair_list) && !is.null(lib_loss_pair_list)
  use_typical_loss_split <- isTRUE(params$use_typical_loss) &&
    !is.null(query_loss_anchor_typ_list) && !is.null(lib_loss_anchor_typ_list) &&
    !is.null(query_loss_pair_typ_list) && !is.null(lib_loss_pair_typ_list)
  use_conf <- isTRUE(params$use_mref_confidence) &&
    !is.null(query_mref_conf) && !is.null(lib_mref_conf)

  if (can_use_prepared_pairloop(params) &&
      !use_typical_loss_combined && !use_split_loss &&
      !use_typical_loss_split && !use_conf) {
    return(compute_distance_matrix_search_prepared(
      query_frag_list, query_loss_list,
      lib_frag_list, lib_loss_list,
      params, progress = progress
    ))
  }

  q_conf_vec <- NULL
  l_conf_vec <- NULL
  if (use_conf) {
    q_conf_vec <- query_mref_conf
    l_conf_vec <- lib_mref_conf
    if (!is.null(names(q_conf_vec))) q_conf_vec <- q_conf_vec[q_ids]
    if (!is.null(names(l_conf_vec))) l_conf_vec <- l_conf_vec[l_ids]
    q_conf_vec <- as.numeric(q_conf_vec)
    l_conf_vec <- as.numeric(l_conf_vec)
  }

  compute_query_row <- function(i) {
    d <- numeric(nl)
    confA <- if (use_conf) q_conf_vec[i] else NULL

    if (use_split_loss) {
      qA_anchor <- query_loss_anchor_list[[i]]
      qA_pair <- query_loss_pair_list[[i]]
      qA_anchor_typ <- if (use_typical_loss_split) query_loss_anchor_typ_list[[i]] else NULL
      qA_pair_typ <- if (use_typical_loss_split) query_loss_pair_typ_list[[i]] else NULL

      for (j in seq_len(nl)) {
        confB <- if (use_conf) l_conf_vec[j] else NULL
        d[j] <- combined_distance_split(
          query_frag_list[[i]], lib_frag_list[[j]],
          qA_anchor, lib_loss_anchor_list[[j]],
          qA_pair, lib_loss_pair_list[[j]],
          params,
          lossA_anchor_typ = qA_anchor_typ,
          lossB_anchor_typ = if (use_typical_loss_split) lib_loss_anchor_typ_list[[j]] else NULL,
          lossA_pair_typ = qA_pair_typ,
          lossB_pair_typ = if (use_typical_loss_split) lib_loss_pair_typ_list[[j]] else NULL,
          confA = confA, confB = confB
        )
      }
    } else {
      qA_typ <- if (use_typical_loss_combined) query_loss_typ_list[[i]] else NULL

      for (j in seq_len(nl)) {
        confB <- if (use_conf) l_conf_vec[j] else NULL
        lB_typ <- if (use_typical_loss_combined) lib_loss_typ_list[[j]] else NULL

        d[j] <- combined_distance(
          query_frag_list[[i]], lib_frag_list[[j]],
          query_loss_list[[i]], lib_loss_list[[j]], params,
          lossA_typ = qA_typ, lossB_typ = lB_typ,
          confA = confA, confB = confB
        )
      }
    }
    list(i = i, d = d)
  }

  if (!params$use_parallel || params$n_cores <= 1) {
    for (i in seq_len(nq)) {
      res <- compute_query_row(i)
      dist_mat[i, ] <- res$d
      if (isTRUE(progress) && i %% 10 == 0) message("Progress: ", i, " / ", nq)
    }
    return(dist_mat)
  }

  cores <- resolve_parallel_cores(params$n_cores)
  idx <- seq_len(nq)

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c(
        "query_frag_list", "query_loss_list",
        "lib_frag_list", "lib_loss_list",
        "query_loss_typ_list", "lib_loss_typ_list",
        "query_loss_anchor_list", "query_loss_pair_list",
        "lib_loss_anchor_list", "lib_loss_pair_list",
        "query_loss_anchor_typ_list", "query_loss_pair_typ_list",
        "lib_loss_anchor_typ_list", "lib_loss_pair_typ_list",
        "q_conf_vec", "l_conf_vec", "nl", "params",
        "combined_distance", "combined_distance_split",
        "compute_distance", "hellinger_distance", "cosine_distance",
        "entropy_distance", "entropy_unweighted_distance",
        "weighted_cosine_distance", "composite_distance",
        "wasserstein_distance", "ppm_wasserstein_distance",
        "PPMWASS_OT_MARGINAL_TOLERANCE", "assess_ot_plan",
        "append_ot_diagnostic", "empty_ot_attempts", "ot_attempt_row",
        "align_spectra", "is_near_zero"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      if (params$distance_method == "ppm_wasserstein" &&
          requireNamespace("approxOT", quietly = TRUE)) {
        library(approxOT)
      }
      if (params$distance_method %in% c("wasserstein", "ppm_wasserstein") &&
          requireNamespace("transport", quietly = TRUE)) {
        library(transport)
      }
      if (params$distance_method %in% c("entropy", "entropy_weighted", "entropy_unweighted") &&
          requireNamespace("msentropy", quietly = TRUE)) {
        library(msentropy)
      }
      NULL
    })
    results <- parallel::parLapply(cl, idx, compute_query_row)
  } else {
    results <- parallel::mclapply(idx, compute_query_row, mc.cores = cores)
  }

  for (res in results) {
    dist_mat[res$i, ] <- res$d
  }

  dist_mat
}

#' Compute Distance Matrices with Component Breakdown
#'
#' This function computes the main distance matrix **and** (optionally) the distances/weights
#' for each sub-component (fragment vs loss, and for split-loss mode: anchored vs pairwise).
#' It is intended for diagnostics, interpretation, and plotting.
#'
#' @param frag_list Named list of fragment spectra.
#' @param loss_list Named list of loss spectra (combined/legacy).
#' @param params Parameter list.
#' @param progress If TRUE, print progress messages.
#' @param loss_typ_list Optional list of typical-loss spectra for combined loss.
#' @param mref_conf Optional named numeric vector of Mref confidence (0..1).
#' @param loss_anchor_list Optional list of anchored-loss spectra (A).
#' @param loss_pair_list Optional list of pairwise-loss spectra (B).
#' @param loss_anchor_typ_list Optional list of typical-loss spectra for anchored losses.
#' @param loss_pair_typ_list Optional list of typical-loss spectra for pairwise losses.
#' @return A list with:
#'   - dist: distance matrix (same as compute_distance_matrix())
#'   - components: named list of matrices (per-component distances/weights)
#' @export
compute_distance_matrices <- function(frag_list, loss_list, params, progress = TRUE,
                                      loss_typ_list = NULL, mref_conf = NULL,
                                      loss_anchor_list = NULL, loss_pair_list = NULL,
                                      loss_anchor_typ_list = NULL, loss_pair_typ_list = NULL) {
  ids <- names(frag_list)
  n <- length(ids)
  directional <- !distance_method_is_symmetric(
    params$distance_method, params$ot_method
  )

  if (isTRUE(progress)) {
    message("Computing distance matrices (with components) for ", n, " compounds...")
    message("Method: ", params$distance_method)
    message("Weights: frag = ", params$w_frag, ", loss = ", params$w_loss)
    if (isTRUE(params$use_split_loss)) {
      message("Split loss: anchored weight = ", if (is.null(params$loss_anchor_weight)) 0.6 else params$loss_anchor_weight)
    }
  }

  make_mat <- function(fill = NA_real_) matrix(fill, n, n, dimnames = list(ids, ids))

  dist_mat <- make_mat(0)
  if (n < 2) {
    return(list(dist = dist_mat, components = list()))
  }

  if (!is.null(loss_list) && !is.null(names(loss_list))) loss_list <- loss_list[ids]
  if (!is.null(loss_typ_list) && !is.null(names(loss_typ_list))) loss_typ_list <- loss_typ_list[ids]

  if (!is.null(loss_anchor_list) && !is.null(names(loss_anchor_list))) loss_anchor_list <- loss_anchor_list[ids]
  if (!is.null(loss_pair_list) && !is.null(names(loss_pair_list))) loss_pair_list <- loss_pair_list[ids]
  if (!is.null(loss_anchor_typ_list) && !is.null(names(loss_anchor_typ_list))) loss_anchor_typ_list <- loss_anchor_typ_list[ids]
  if (!is.null(loss_pair_typ_list) && !is.null(names(loss_pair_typ_list))) loss_pair_typ_list <- loss_pair_typ_list[ids]

  use_typical_loss_combined <- isTRUE(params$use_typical_loss) && !is.null(loss_typ_list)
  use_split_loss <- isTRUE(params$use_split_loss) &&
    !is.null(loss_anchor_list) && !is.null(loss_pair_list)
  use_typical_loss_split <- isTRUE(params$use_typical_loss) &&
    !is.null(loss_anchor_typ_list) && !is.null(loss_pair_typ_list)

  use_conf <- isTRUE(params$use_mref_confidence) && !is.null(mref_conf)
  conf_vec <- NULL
  if (use_conf) {
    conf_vec <- mref_conf
    if (!is.null(names(conf_vec))) conf_vec <- conf_vec[ids]
    conf_vec <- as.numeric(conf_vec)
  }

  comps <- list(
    d_frag = make_mat(0),
    d_loss = make_mat(0)
  )

  if (!use_split_loss) {
    comps$d_loss_raw <- make_mat(0)
    if (use_typical_loss_combined) {
      comps$d_loss_typ <- make_mat(NA_real_)
      comps$w_loss_typical <- make_mat(0)
      comps$w_loss_raw <- make_mat(1)
    }
  } else {
    comps$d_anchor <- make_mat(0)
    comps$d_anchor_raw <- make_mat(0)
    comps$d_pair <- make_mat(0)
    comps$d_pair_raw <- make_mat(0)
    comps$w_anchor_final <- make_mat(NA_real_)
    comps$w_pair_final <- make_mat(NA_real_)

    if (use_typical_loss_split) {
      comps$d_anchor_typ <- make_mat(NA_real_)
      comps$w_anchor_typical <- make_mat(0)
      comps$w_anchor_raw <- make_mat(1)
      comps$d_pair_typ <- make_mat(NA_real_)
      comps$w_pair_typical <- make_mat(0)
      comps$w_pair_raw <- make_mat(1)
    }
  }

  if (use_conf) {
    comps$conf_pair <- make_mat(NA_real_)
  }

  component_names <- names(comps)

  compute_row <- function(i) {
    js <- if (directional) setdiff(seq_len(n), i) else seq.int(i + 1L, n)
    out <- list(i = i, j = js)
    for (nm in component_names) {
      out[[nm]] <- numeric(length(js))
    }
    d_total <- numeric(length(js))

    confA <- if (use_conf) conf_vec[i] else NULL

    if (use_split_loss) {
      lossA_anchor <- loss_anchor_list[[i]]
      lossA_pair <- loss_pair_list[[i]]
      lossA_anchor_typ <- if (use_typical_loss_split) loss_anchor_typ_list[[i]] else NULL
      lossA_pair_typ <- if (use_typical_loss_split) loss_pair_typ_list[[i]] else NULL

      for (k in seq_along(js)) {
        j <- js[[k]]
        confB <- if (use_conf) conf_vec[j] else NULL

        det <- combined_distance_split_details(
          frag_list[[i]], frag_list[[j]],
          lossA_anchor, loss_anchor_list[[j]],
          lossA_pair, loss_pair_list[[j]],
          params,
          lossA_anchor_typ = lossA_anchor_typ,
          lossB_anchor_typ = if (use_typical_loss_split) loss_anchor_typ_list[[j]] else NULL,
          lossA_pair_typ = lossA_pair_typ,
          lossB_pair_typ = if (use_typical_loss_split) loss_pair_typ_list[[j]] else NULL,
          confA = confA, confB = confB
        )

        d_total[k] <- det$d_total
        out[["d_frag"]][k] <- det$d_frag
        out[["d_loss"]][k] <- det$d_loss
        out[["d_anchor"]][k] <- det$d_anchor
        out[["d_anchor_raw"]][k] <- det$d_anchor_raw
        out[["d_pair"]][k] <- det$d_pair
        out[["d_pair_raw"]][k] <- det$d_pair_raw
        out[["w_anchor_final"]][k] <- det$w_anchor_final
        out[["w_pair_final"]][k] <- det$w_pair_final

        if (use_typical_loss_split) {
          out[["d_anchor_typ"]][k] <- det$d_anchor_typ
          out[["w_anchor_typical"]][k] <- det$w_anchor_typical
          out[["w_anchor_raw"]][k] <- det$w_anchor_raw
          out[["d_pair_typ"]][k] <- det$d_pair_typ
          out[["w_pair_typical"]][k] <- det$w_pair_typical
          out[["w_pair_raw"]][k] <- det$w_pair_raw
        }

        if (use_conf) {
          out[["conf_pair"]][k] <- det$conf_pair
        }
      }
    } else {
      lossA_typ <- if (use_typical_loss_combined) loss_typ_list[[i]] else NULL

      for (k in seq_along(js)) {
        j <- js[[k]]
        confB <- if (use_conf) conf_vec[j] else NULL
        lossB_typ <- if (use_typical_loss_combined) loss_typ_list[[j]] else NULL

        det <- combined_distance_details(
          frag_list[[i]], frag_list[[j]],
          loss_list[[i]], loss_list[[j]], params,
          lossA_typ = lossA_typ, lossB_typ = lossB_typ,
          confA = confA, confB = confB
        )

        d_total[k] <- det$d_total
        out[["d_frag"]][k] <- det$d_frag
        out[["d_loss"]][k] <- det$d_loss
        out[["d_loss_raw"]][k] <- det$d_loss_raw
        if (use_typical_loss_combined) {
          out[["d_loss_typ"]][k] <- det$d_loss_typ
          out[["w_loss_typical"]][k] <- det$w_loss_typical
          out[["w_loss_raw"]][k] <- det$w_loss_raw
        }
        if (use_conf) {
          a <- if (is.finite(confA)) confA else 0
          b <- if (is.finite(confB)) confB else 0
          a <- min(max(a, 0), 1)
          b <- min(max(b, 0), 1)
          out[["conf_pair"]][k] <- sqrt(a * b)
        }
      }
    }

    out$d_total <- d_total
    out
  }

  fill_from_row <- function(dist_mat, comps, res) {
    i <- res$i
    js <- res$j
    dist_mat[i, js] <- res$d_total
    if (!directional) dist_mat[js, i] <- res$d_total
    for (nm in component_names) {
      comps[[nm]][i, js] <- res[[nm]]
      if (!directional) comps[[nm]][js, i] <- res[[nm]]
    }
    list(dist_mat = dist_mat, comps = comps)
  }

  row_indices <- if (directional) seq_len(n) else seq_len(n - 1L)
  if (!params$use_parallel || params$n_cores <= 1) {
    for (i in row_indices) {
      res <- compute_row(i)
      filled <- fill_from_row(dist_mat, comps, res)
      dist_mat <- filled$dist_mat
      comps <- filled$comps
      if (isTRUE(progress) && i %% 10 == 0) {
        message("Progress: ", i, " / ", length(row_indices))
      }
    }
  } else {
    cores <- resolve_parallel_cores(params$n_cores)
    idx <- row_indices

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(
        cl,
        varlist = c(
          "frag_list", "loss_list",
          "loss_typ_list",
          "loss_anchor_list", "loss_pair_list",
          "loss_anchor_typ_list", "loss_pair_typ_list",
          "conf_vec", "params",
          "combined_distance_details", "combined_distance_split_details",
          "compute_distance", "hellinger_distance", "cosine_distance",
          "entropy_distance", "entropy_unweighted_distance",
          "weighted_cosine_distance", "composite_distance",
          "wasserstein_distance", "ppm_wasserstein_distance",
          "PPMWASS_OT_MARGINAL_TOLERANCE", "assess_ot_plan",
          "append_ot_diagnostic", "empty_ot_attempts", "ot_attempt_row",
          "align_spectra", "is_near_zero"
        ),
        envir = environment()
      )
      parallel::clusterEvalQ(cl, {
        if (params$distance_method == "ppm_wasserstein" &&
            requireNamespace("approxOT", quietly = TRUE)) {
          library(approxOT)
        }
        if (params$distance_method %in% c("wasserstein", "ppm_wasserstein") &&
            requireNamespace("transport", quietly = TRUE)) {
          library(transport)
        }
        if (params$distance_method %in% c("entropy", "entropy_weighted", "entropy_unweighted") &&
            requireNamespace("msentropy", quietly = TRUE)) {
          library(msentropy)
        }
        NULL
      })
      results <- parallel::parLapply(cl, idx, compute_row)
    } else {
      results <- parallel::mclapply(idx, compute_row, mc.cores = cores)
    }

    for (res in results) {
      filled <- fill_from_row(dist_mat, comps, res)
      dist_mat <- filled$dist_mat
      comps <- filled$comps
    }
  }

  diag(dist_mat) <- 0
  for (nm in names(comps)) {
    if (startsWith(nm, "d_")) diag(comps[[nm]]) <- 0
  }

  ratio_mat <- function(num, den) {
    out <- make_mat(0)
    ok <- is.finite(num) & is.finite(den) & den > 0
    out[ok] <- num[ok] / den[ok]
    out
  }

  D2 <- params$w_frag * (comps$d_frag^2) + params$w_loss * (comps$d_loss^2)
  comps$c_frag_total <- ratio_mat(params$w_frag * (comps$d_frag^2), D2)
  comps$c_loss_total <- ratio_mat(params$w_loss * (comps$d_loss^2), D2)

  if (!use_split_loss) {
    w_raw <- if (!is.null(comps$w_loss_raw)) comps$w_loss_raw else make_mat(1)
    w_typ <- if (!is.null(comps$w_loss_typical)) comps$w_loss_typical else make_mat(0)

    d_raw2 <- comps$d_loss_raw^2
    d_typ2 <- if (!is.null(comps$d_loss_typ)) comps$d_loss_typ^2 else make_mat(0)

    d_raw2[!is.finite(d_raw2)] <- 0
    d_typ2[!is.finite(d_typ2)] <- 0
    w_raw[!is.finite(w_raw)] <- 0
    w_typ[!is.finite(w_typ)] <- 0

    loss2 <- comps$d_loss^2
    comps$c_loss_raw_in_loss <- ratio_mat(w_raw * d_raw2, loss2)
    comps$c_loss_typ_in_loss <- ratio_mat(w_typ * d_typ2, loss2)

    comps$c_loss_raw_total <- ratio_mat(params$w_loss * w_raw * d_raw2, D2)
    comps$c_loss_typ_total <- ratio_mat(params$w_loss * w_typ * d_typ2, D2)
  } else {
    wA <- comps$w_anchor_final
    wB <- comps$w_pair_final
    wA[!is.finite(wA)] <- 0
    wB[!is.finite(wB)] <- 0

    dA2 <- comps$d_anchor^2
    dB2 <- comps$d_pair^2
    dA2[!is.finite(dA2)] <- 0
    dB2[!is.finite(dB2)] <- 0

    loss2 <- comps$d_loss^2
    comps$c_anchor_in_loss <- ratio_mat(wA * dA2, loss2)
    comps$c_pair_in_loss <- ratio_mat(wB * dB2, loss2)

    comps$c_anchor_total <- ratio_mat(params$w_loss * wA * dA2, D2)
    comps$c_pair_total <- ratio_mat(params$w_loss * wB * dB2, D2)

    wA_raw <- if (!is.null(comps$w_anchor_raw)) comps$w_anchor_raw else make_mat(1)
    wA_typ <- if (!is.null(comps$w_anchor_typical)) comps$w_anchor_typical else make_mat(0)
    dA_raw2 <- comps$d_anchor_raw^2
    dA_typ2 <- if (!is.null(comps$d_anchor_typ)) comps$d_anchor_typ^2 else make_mat(0)

    wA_raw[!is.finite(wA_raw)] <- 0
    wA_typ[!is.finite(wA_typ)] <- 0
    dA_raw2[!is.finite(dA_raw2)] <- 0
    dA_typ2[!is.finite(dA_typ2)] <- 0

    anchor2 <- comps$d_anchor^2
    comps$c_anchor_raw_in_anchor <- ratio_mat(wA_raw * dA_raw2, anchor2)
    comps$c_anchor_typ_in_anchor <- ratio_mat(wA_typ * dA_typ2, anchor2)

    comps$c_anchor_raw_total <- ratio_mat(params$w_loss * wA * wA_raw * dA_raw2, D2)
    comps$c_anchor_typ_total <- ratio_mat(params$w_loss * wA * wA_typ * dA_typ2, D2)

    wB_raw <- if (!is.null(comps$w_pair_raw)) comps$w_pair_raw else make_mat(1)
    wB_typ <- if (!is.null(comps$w_pair_typical)) comps$w_pair_typical else make_mat(0)
    dB_raw2 <- comps$d_pair_raw^2
    dB_typ2 <- if (!is.null(comps$d_pair_typ)) comps$d_pair_typ^2 else make_mat(0)

    wB_raw[!is.finite(wB_raw)] <- 0
    wB_typ[!is.finite(wB_typ)] <- 0
    dB_raw2[!is.finite(dB_raw2)] <- 0
    dB_typ2[!is.finite(dB_typ2)] <- 0

    pair2 <- comps$d_pair^2
    comps$c_pair_raw_in_pair <- ratio_mat(wB_raw * dB_raw2, pair2)
    comps$c_pair_typ_in_pair <- ratio_mat(wB_typ * dB_typ2, pair2)

    comps$c_pair_raw_total <- ratio_mat(params$w_loss * wB * wB_raw * dB_raw2, D2)
    comps$c_pair_typ_total <- ratio_mat(params$w_loss * wB * wB_typ * dB_typ2, D2)
  }

  list(dist = dist_mat, components = comps)
}
