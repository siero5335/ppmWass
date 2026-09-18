#' @keywords internal
#' @noRd
valid_distance_methods <- function() {
  c(
    "hellinger", "wasserstein", "ppm_wasserstein", "cosine", "entropy",
    "entropy_weighted", "entropy_unweighted", "weighted_cosine", "composite"
  )
}

#' @keywords internal
#' @noRd
validate_distance_method_name <- function(method, arg = "distance_method") {
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% valid_distance_methods()) {
    stop(arg, " must be one of: ", paste(valid_distance_methods(), collapse = ", "))
  }
  invisible(method)
}

#' @keywords internal
#' @noRd
top_k_order <- function(x, k) {
  ok <- which(is.finite(x))
  if (length(ok) == 0L) return(rep(NA_integer_, k))
  if (length(ok) <= k) return(c(ok[order(x[ok], ok)], rep(NA_integer_, k - length(ok))))

  vals <- x[ok]
  kth <- sort(vals, partial = k)[k]
  below <- ok[vals < kth]
  equal <- ok[vals == kth]
  below <- below[order(x[below], below)]
  equal <- equal[order(equal)]
  c(below, utils::head(equal, k - length(below)))
}

#' @keywords internal
#' @noRd
resolve_exclude_self <- function(query_ids, lib_ids, exclude_self = NULL) {
  if (is.null(exclude_self)) {
    return(identical(as.character(query_ids), as.character(lib_ids)))
  }
  if (!is.logical(exclude_self) || length(exclude_self) != 1L || is.na(exclude_self)) {
    stop("exclude_self must be TRUE, FALSE, or NULL.")
  }
  exclude_self
}

#' @keywords internal
#' @noRd
select_topk_candidates <- function(D_first, top_k,
                                   query_ids = rownames(D_first),
                                   lib_ids = colnames(D_first),
                                   exclude_self = FALSE) {
  if (!is.matrix(D_first)) stop("D_first must be a matrix.")
  nq <- nrow(D_first)
  nl <- ncol(D_first)
  top_k <- as.integer(top_k)
  if (length(top_k) != 1L || is.na(top_k) || top_k < 1L) {
    stop("top_k must be >= 1.")
  }
  top_k <- min(top_k, nl)

  if (is.null(query_ids)) query_ids <- as.character(seq_len(nq))
  if (is.null(lib_ids)) lib_ids <- as.character(seq_len(nl))
  if (length(query_ids) != nq) stop("query_ids length must match nrow(D_first).")
  if (length(lib_ids) != nl) stop("lib_ids length must match ncol(D_first).")

  candidate_idx <- matrix(NA_integer_, nq, top_k)
  first_pass_distance <- matrix(NA_real_, nq, top_k)
  first_pass_rank <- matrix(NA_integer_, nq, top_k)
  n_finite_first_pass <- integer(nq)

  for (i in seq_len(nq)) {
    qrow <- D_first[i, ]
    if (isTRUE(exclude_self)) {
      qrow[lib_ids == query_ids[i]] <- Inf
    }
    n_finite_first_pass[i] <- sum(is.finite(qrow))
    ord <- top_k_order(qrow, top_k)
    candidate_idx[i, ] <- ord
    ok <- !is.na(ord)
    first_pass_distance[i, ok] <- qrow[ord[ok]]
    first_pass_rank[i, ok] <- seq_len(sum(ok))
  }

  rownames(candidate_idx) <- rownames(first_pass_distance) <- rownames(first_pass_rank) <- query_ids

  list(
    candidate_idx = candidate_idx,
    first_pass_distance = first_pass_distance,
    first_pass_rank = first_pass_rank,
    n_finite_first_pass = n_finite_first_pass
  )
}

#' @keywords internal
#' @noRd
refine_topk_candidates <- function(query_frag_list, query_loss_list,
                                   lib_frag_list, lib_loss_list,
                                   candidate_idx, params,
                                   use_parallel = params$use_parallel %||% FALSE,
                                   n_cores = params$n_cores %||% 1L,
                                   progress = TRUE,
                                   progress_prefix = "[top-k]",
                                   progress_step = NULL,
                                   query_loss_typ_list = NULL,
                                   lib_loss_typ_list = NULL,
                                   query_mref_conf = NULL,
                                   lib_mref_conf = NULL) {
  if (!is.matrix(candidate_idx)) stop("candidate_idx must be a matrix.")
  nq <- nrow(candidate_idx)
  top_k <- ncol(candidate_idx)
  empty <- matrix(numeric(0), ncol = 2)
  use_loss <- !is.null(query_loss_list) && !is.null(lib_loss_list)
  use_typical_loss <- isTRUE(params$use_typical_loss) &&
    !is.null(query_loss_typ_list) && !is.null(lib_loss_typ_list)
  use_conf <- isTRUE(params$use_mref_confidence) &&
    !is.null(query_mref_conf) && !is.null(lib_mref_conf)
  q_conf_vec <- if (use_conf) as.numeric(query_mref_conf) else NULL
  l_conf_vec <- if (use_conf) as.numeric(lib_mref_conf) else NULL

  refine_one <- function(i) {
    d_ref <- rep(NA_real_, top_k)
    lossA <- if (use_loss) query_loss_list[[i]] else empty
    lossA_typ <- if (use_typical_loss) query_loss_typ_list[[i]] else NULL
    confA <- if (use_conf) q_conf_vec[i] else NULL

    for (m in seq_len(top_k)) {
      j <- candidate_idx[i, m]
      if (is.na(j)) next
      d <- combined_distance(
        query_frag_list[[i]], lib_frag_list[[j]],
        lossA,
        if (use_loss) lib_loss_list[[j]] else empty,
        params,
        lossA_typ = lossA_typ,
        lossB_typ = if (use_typical_loss) lib_loss_typ_list[[j]] else NULL,
        confA = confA,
        confB = if (use_conf) l_conf_vec[j] else NULL
      )
      if (is.finite(d)) d_ref[m] <- d
    }
    d_ref
  }

  if (isTRUE(use_parallel) && (n_cores %||% 1L) > 1L) {
    cores <- resolve_parallel_cores(n_cores)
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, library(ppmWass))
      parallel::clusterExport(
        cl,
        varlist = c(
          "candidate_idx", "top_k", "query_frag_list", "query_loss_list",
          "lib_frag_list", "lib_loss_list", "query_loss_typ_list",
          "lib_loss_typ_list", "q_conf_vec", "l_conf_vec", "use_loss",
          "use_typical_loss", "use_conf", "params", "empty",
          "combined_distance"
        ),
        envir = environment()
      )
      rows <- parallel::parLapply(cl, seq_len(nq), refine_one)
    } else {
      rows <- parallel::mclapply(seq_len(nq), refine_one, mc.cores = cores)
    }
  } else {
    if (is.null(progress_step)) progress_step <- max(1L, min(1024L, nq))
    rows <- vector("list", nq)
    for (i in seq_len(nq)) {
      if (isTRUE(progress) && (i == 1L || i %% progress_step == 0L || i == nq)) {
        message(sprintf("%s refine query %d / %d", progress_prefix, i, nq))
      }
      rows[[i]] <- refine_one(i)
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- rownames(candidate_idx)
  out
}

#' @keywords internal
#' @noRd
format_topk_long <- function(candidate_idx, rerank_distance,
                             first_pass_rank = NULL,
                             first_pass_distance = NULL,
                             query_ids,
                             lib_ids,
                             final_top_k,
                             first_pass_method,
                             rerank_method,
                             return_first_pass = TRUE) {
  if (!is.matrix(candidate_idx) || !is.matrix(rerank_distance)) {
    stop("candidate_idx and rerank_distance must be matrices.")
  }
  final_top_k <- as.integer(final_top_k)
  rows <- vector("list", nrow(candidate_idx))

  for (i in seq_len(nrow(candidate_idx))) {
    cand <- candidate_idx[i, ]
    dist <- rerank_distance[i, ]
    ok <- which(!is.na(cand) & is.finite(dist))
    if (length(ok) == 0L) {
      rows[[i]] <- data.frame(
        query_id = character(0),
        lib_id = character(0),
        rerank_rank = integer(0),
        rerank_distance = numeric(0),
        first_pass_rank = integer(0),
        first_pass_distance = numeric(0),
        first_pass_method = character(0),
        rerank_method = character(0),
        stringsAsFactors = FALSE
      )
      next
    }

    ok <- ok[order(dist[ok], cand[ok])]
    ok <- utils::head(ok, final_top_k)
    out <- data.frame(
      query_id = rep(query_ids[i], length(ok)),
      lib_id = lib_ids[cand[ok]],
      rerank_rank = seq_along(ok),
      rerank_distance = dist[ok],
      first_pass_rank = if (is.null(first_pass_rank)) ok else first_pass_rank[i, ok],
      first_pass_distance = if (is.null(first_pass_distance)) NA_real_ else first_pass_distance[i, ok],
      first_pass_method = rep(first_pass_method, length(ok)),
      rerank_method = rep(rerank_method, length(ok)),
      stringsAsFactors = FALSE
    )
    rows[[i]] <- out
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (!isTRUE(return_first_pass)) {
    out$first_pass_rank <- NULL
    out$first_pass_distance <- NULL
  }
  out
}
