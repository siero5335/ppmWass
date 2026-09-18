#' @keywords internal
#' @noRd
prepare_spectrum_alignment <- function(spec) {
  if (is.null(spec) || nrow(spec) == 0) {
    return(list(mz = numeric(0), int = numeric(0), cs = 0))
  }
  ord <- order(spec[, 1])
  int <- spec[ord, 2]
  list(
    mz = spec[ord, 1],
    int = int,
    cs = c(0, cumsum(int))
  )
}

#' @keywords internal
#' @noRd
prepare_spectra_alignment_list <- function(spec_list) {
  lapply(spec_list, prepare_spectrum_alignment)
}

#' @keywords internal
#' @noRd
align_spectra_prepared <- function(specA, specB, ppm = 20) {
  mzA <- specA$mz
  mzB <- specB$mz
  if (length(mzA) == 0 || length(mzB) == 0) {
    return(list(p = numeric(0), q = numeric(0), mz = numeric(0)))
  }

  all_mz <- sort(unique(c(mzA, mzB)))
  n <- length(all_mz)
  merged_mz <- vector("list", n)
  k <- 0
  i <- 1
  while (i <= n) {
    current_mz <- all_mz[i]
    tol <- current_mz * ppm * 1e-6
    j <- i
    while (j <= n && abs(all_mz[j] - current_mz) <= tol) {
      j <- j + 1
    }
    k <- k + 1
    merged_mz[[k]] <- mean(all_mz[i:(j - 1)])
    i <- j
  }
  merged_mz <- unlist(merged_mz[seq_len(k)], use.names = FALSE)

  tol_vec <- merged_mz * ppm * 1e-6
  lower <- merged_mz - tol_vec
  upper <- merged_mz + tol_vec

  leftA <- findInterval(lower, mzA, left.open = TRUE) + 1
  rightA <- findInterval(upper, mzA)
  leftB <- findInterval(lower, mzB, left.open = TRUE) + 1
  rightB <- findInterval(upper, mzB)

  p <- numeric(length(merged_mz))
  q <- numeric(length(merged_mz))

  okA <- leftA <= rightA
  okB <- leftB <= rightB
  if (any(okA)) p[okA] <- specA$cs[rightA[okA] + 1] - specA$cs[leftA[okA]]
  if (any(okB)) q[okB] <- specB$cs[rightB[okB] + 1] - specB$cs[leftB[okB]]

  if (sum(p) > 0) p <- p / sum(p)
  if (sum(q) > 0) q <- q / sum(q)

  list(p = p, q = q, mz = merged_mz)
}

#' @keywords internal
#' @noRd
compute_distance_prepared <- function(specA, specB, method = "hellinger", ppm = 20,
                                      mass_power = 3, intensity_power = 0.5) {
  aligned <- align_spectra_prepared(specA, specB, ppm)
  p <- aligned$p
  q <- aligned$q
  if (length(p) == 0) return(1)

  switch(method,
    hellinger = sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2),
    cosine = {
      dot <- sum(p * q)
      norm_p <- sqrt(sum(p^2))
      norm_q <- sqrt(sum(q^2))
      if (is_near_zero(norm_p) || is_near_zero(norm_q)) return(1)
      1 - dot / (norm_p * norm_q)
    },
    weighted_cosine = {
      mz <- aligned$mz
      wp <- (mz^mass_power) * (p^intensity_power)
      wq <- (mz^mass_power) * (q^intensity_power)
      dot <- sum(wp * wq)
      norm_p <- sqrt(sum(wp^2))
      norm_q <- sqrt(sum(wq^2))
      if (is_near_zero(norm_p) || is_near_zero(norm_q)) return(1)
      similarity <- dot / (norm_p * norm_q)
      similarity <- min(max(similarity, 0), 1)
      1 - similarity
    },
    composite = {
      mz <- aligned$mz
      wp <- (mz^mass_power) * (p^intensity_power)
      wq <- (mz^mass_power) * (q^intensity_power)
      dot <- sum(wp * wq)
      norm_p <- sqrt(sum(wp^2))
      norm_q <- sqrt(sum(wq^2))
      if (is_near_zero(norm_p) || is_near_zero(norm_q)) return(1)
      sw <- dot / (norm_p * norm_q)
      both_nonzero <- which(p > 0 & q > 0)
      nxy <- length(both_nonzero)
      nx <- sum(p > 0)
      if (nxy < 2 || nx == 0) {
        sw <- min(max(sw, 0), 1)
        return(1 - sw)
      }
      p_shared <- p[both_nonzero]
      q_shared <- q[both_nonzero]
      ratio_sum <- 0
      valid_pairs <- 0
      for (k in seq_len(nxy - 1)) {
        ratio_p <- p_shared[k] / p_shared[k + 1]
        ratio_q <- q_shared[k] / q_shared[k + 1]
        if (!is.finite(ratio_p) || !is.finite(ratio_q) ||
            is_near_zero(ratio_p) || is_near_zero(ratio_q)) next
        ratio_sum <- ratio_sum + min(ratio_p, ratio_q) / max(ratio_p, ratio_q)
        valid_pairs <- valid_pairs + 1
      }
      if (valid_pairs == 0) {
        sw <- min(max(sw, 0), 1)
        return(1 - sw)
      }
      sr <- ratio_sum / valid_pairs
      similarity <- (nx * sw + nxy * sr) / (nx + nxy)
      similarity <- min(max(similarity, 0), 1)
      1 - similarity
    },
    stop("Unsupported prepared distance method: ", method)
  )
}

#' @keywords internal
#' @noRd
combined_distance_prepared <- function(fragA, fragB, lossA, lossB, params) {
  mass_power <- if (is.null(params$mass_power)) 3 else params$mass_power
  intensity_power <- if (is.null(params$intensity_power)) 0.5 else params$intensity_power
  method <- params$distance_method
  ppm <- params$tol_ppm

  d_frag <- compute_distance_prepared(
    fragA, fragB, method, ppm,
    mass_power = mass_power,
    intensity_power = intensity_power
  )
  d_loss <- compute_distance_prepared(
    lossA, lossB, method, ppm,
    mass_power = mass_power,
    intensity_power = intensity_power
  )

  sqrt(params$w_frag * d_frag^2 + params$w_loss * d_loss^2)
}

#' @keywords internal
#' @noRd
can_use_prepared_pairloop <- function(params) {
  (params$distance_method %in% c("hellinger", "cosine", "weighted_cosine", "composite")) &&
    !isTRUE(params$use_typical_loss) &&
    !isTRUE(params$use_split_loss) &&
    !isTRUE(params$use_mref_confidence)
}

#' Whether a distance method is symmetric in its two spectrum arguments
#'
#' The Stein-style composite implementation is intentionally directional: its
#' ratio term uses the number of positive peaks in the first (query) spectrum.
#' Finite-iteration approximate OT is evaluated directionally because the
#' current approxOT iteration order is not guaranteed to be invariant to
#' exchanging its two marginals. Exact ppm-Wasserstein and the other currently
#' implemented publication methods are symmetric.
#'
#' @param method Distance-method identifier.
#' @param ot_method Optimal-transport backend when `method` is
#'   `"ppm_wasserstein"`.
#' @return A single logical value.
#' @keywords internal
#' @noRd
distance_method_is_symmetric <- function(method, ot_method = "exact") {
  if (is.null(ot_method)) ot_method <- "exact"
  if (identical(method, "ppm_wasserstein")) {
    return(ot_method %in% c("exact", "exact_sparse"))
  }
  !method %in% c("composite")
}

#' @keywords internal
#' @noRd
compute_distance_matrix_prepared <- function(frag_list, loss_list, params,
                                             progress = TRUE) {
  ids <- names(frag_list)
  n <- length(ids)
  dist_mat <- matrix(0, n, n, dimnames = list(ids, ids))
  if (n < 2) return(dist_mat)

  if (!distance_method_is_symmetric(params$distance_method, params$ot_method)) {
    directional_dist <- compute_distance_matrix_search_prepared(
      query_frag_list = frag_list,
      query_loss_list = loss_list,
      lib_frag_list = frag_list,
      lib_loss_list = loss_list,
      params = params,
      progress = progress
    )
    diag(directional_dist) <- 0
    return(directional_dist)
  }

  frag_pp <- prepare_spectra_alignment_list(frag_list)
  loss_pp <- prepare_spectra_alignment_list(loss_list)

  compute_row <- function(i) {
    d <- numeric(n - i)
    for (j in (i + 1):n) {
      d[j - i] <- combined_distance_prepared(
        frag_pp[[i]], frag_pp[[j]],
        loss_pp[[i]], loss_pp[[j]],
        params
      )
    }
    list(i = i, d = d)
  }

  if (!params$use_parallel || params$n_cores <= 1) {
    for (i in seq_len(n - 1)) {
      res <- compute_row(i)
      dist_mat[i, (i + 1):n] <- res$d
      dist_mat[(i + 1):n, i] <- res$d
      if (isTRUE(progress) && i %% 10 == 0) message("Progress: ", i, " / ", n - 1)
    }
    return(dist_mat)
  }

  cores <- resolve_parallel_cores(params$n_cores)
  idx <- seq_len(n - 1)

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c(
        "frag_pp", "loss_pp", "params", "n",
        "combined_distance_prepared", "compute_distance_prepared",
        "align_spectra_prepared", "is_near_zero"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      if (requireNamespace("ppmWass", quietly = TRUE)) {
        suppressPackageStartupMessages(library(ppmWass))
      }
      NULL
    })
    results <- parallel::parLapply(cl, idx, compute_row)
  } else {
    results <- parallel::mclapply(idx, compute_row, mc.cores = cores)
  }

  for (res in results) {
    dist_mat[res$i, (res$i + 1):n] <- res$d
    dist_mat[(res$i + 1):n, res$i] <- res$d
  }

  dist_mat
}

#' @keywords internal
#' @noRd
compute_distance_matrix_search_prepared <- function(query_frag_list, query_loss_list,
                                                    lib_frag_list, lib_loss_list,
                                                    params, progress = TRUE) {
  q_ids <- names(query_frag_list)
  l_ids <- names(lib_frag_list)
  nq <- length(q_ids)
  nl <- length(l_ids)
  dist_mat <- matrix(0, nq, nl, dimnames = list(q_ids, l_ids))
  if (nq == 0 || nl == 0) return(dist_mat)

  q_frag_pp <- prepare_spectra_alignment_list(query_frag_list)
  q_loss_pp <- prepare_spectra_alignment_list(query_loss_list)
  l_frag_pp <- prepare_spectra_alignment_list(lib_frag_list)
  l_loss_pp <- prepare_spectra_alignment_list(lib_loss_list)

  compute_query_row <- function(i) {
    d <- numeric(nl)
    for (j in seq_len(nl)) {
      d[j] <- combined_distance_prepared(
        q_frag_pp[[i]], l_frag_pp[[j]],
        q_loss_pp[[i]], l_loss_pp[[j]],
        params
      )
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
        "q_frag_pp", "q_loss_pp", "l_frag_pp", "l_loss_pp",
        "params", "nl",
        "combined_distance_prepared", "compute_distance_prepared",
        "align_spectra_prepared", "is_near_zero"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      if (requireNamespace("ppmWass", quietly = TRUE)) {
        suppressPackageStartupMessages(library(ppmWass))
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
