#' @keywords internal
#' @noRd
hellinger_distance <- function(specA, specB, ppm = 20) {
  if (nrow(specA) == 0 || nrow(specB) == 0) return(1)

  aligned <- align_spectra(specA, specB, ppm)
  p <- aligned$p
  q <- aligned$q

  if (length(p) == 0) return(1)

  sqrt(sum((sqrt(p) - sqrt(q))^2)) / sqrt(2)
}

#' @keywords internal
#' @noRd
cosine_distance <- function(specA, specB, ppm = 20) {
  if (nrow(specA) == 0 || nrow(specB) == 0) return(1)

  aligned <- align_spectra(specA, specB, ppm)
  p <- aligned$p
  q <- aligned$q

  if (length(p) == 0) return(1)

  dot <- sum(p * q)
  norm_p <- sqrt(sum(p^2))
  norm_q <- sqrt(sum(q^2))

  if (is_near_zero(norm_p) || is_near_zero(norm_q)) return(1)

  1 - dot / (norm_p * norm_q)
}

#' @keywords internal
#' @noRd
weighted_cosine_distance <- function(specA, specB, ppm = 20, mass_power = 3, intensity_power = 0.5) {
  if (nrow(specA) == 0 || nrow(specB) == 0) return(1)

  aligned <- align_spectra(specA, specB, ppm)
  p <- aligned$p
  q <- aligned$q
  mz <- aligned$mz

  if (length(p) == 0) return(1)

  wp <- (mz^mass_power) * (p^intensity_power)
  wq <- (mz^mass_power) * (q^intensity_power)

  dot <- sum(wp * wq)
  norm_p <- sqrt(sum(wp^2))
  norm_q <- sqrt(sum(wq^2))

  if (is_near_zero(norm_p) || is_near_zero(norm_q)) return(1)

  similarity <- dot / (norm_p * norm_q)
  similarity <- min(max(similarity, 0), 1)
  1 - similarity
}

#' @keywords internal
#' @noRd
composite_distance <- function(specA, specB, ppm = 20, mass_power = 3, intensity_power = 0.5) {
  if (nrow(specA) == 0 || nrow(specB) == 0) return(1)

  aligned <- align_spectra(specA, specB, ppm)
  p <- aligned$p
  q <- aligned$q
  mz <- aligned$mz

  if (length(p) == 0) return(1)

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
}

#' @keywords internal
#' @noRd
wasserstein_distance <- function(specA, specB, ppm = 20, align = FALSE) {
  if (!requireNamespace("transport", quietly = TRUE)) {
    warning("transport package not available, using Hellinger instead")
    return(hellinger_distance(specA, specB, ppm))
  }

  if (nrow(specA) == 0 || nrow(specB) == 0) return(1)

  if (align) {
    aligned <- align_spectra(specA, specB, ppm)
    if (length(aligned$mz) == 0) return(1)

    mz <- aligned$mz
    p <- aligned$p
    q <- aligned$q

    idx_p <- which(p > 0)
    idx_q <- which(q > 0)
    if (length(idx_p) == 0 || length(idx_q) == 0) return(1)

    mz_p <- mz[idx_p]
    mz_q <- mz[idx_q]
    w_p <- p[idx_p]
    w_q <- q[idx_q]
  } else {
    mz_p <- specA[, 1]
    mz_q <- specB[, 1]
    w_p <- specA[, 2]
    w_q <- specB[, 2]
  }

  if (sum(w_p) <= 0 || sum(w_q) <= 0) return(1)

  w_p <- w_p / sum(w_p)
  w_q <- w_q / sum(w_q)

  mz_all <- c(mz_p, mz_q)
  mz_min <- min(mz_all)
  mz_max <- max(mz_all)
  mz_range <- mz_max - mz_min

  if (is_near_zero(mz_range)) {
    return(0)
  }

  mz_p_scaled <- (mz_p - mz_min) / mz_range
  mz_q_scaled <- (mz_q - mz_min) / mz_range

  w_dist <- transport::wasserstein1d(mz_p_scaled, mz_q_scaled, p = 1,
                                     wa = w_p, wb = w_q)

  min(max(w_dist, 0), 1)
}

# PPM-aware Wasserstein distance for HRMS spectra.
#
# This uses a ppm-aware ground cost that saturates at 1 beyond a few multiples
# of a ppm base-cost scale. The base scale is not a zero-cost tolerance: every
# nonzero displacement has nonzero cost. approxOT supplies the approximate
# Sinkhorn/Greenkhorn paths; transport supplies exact OT and the deterministic
# fallback. Each returned plan is accepted only after explicit marginal checks.
PPMWASS_OT_MARGINAL_TOLERANCE <- 1e-8

#' Assess whether a sparse OT plan satisfies both requested marginals
#'
#' Finite objective values alone do not establish solver convergence. This
#' helper validates plan structure, mass finiteness/non-negativity, and the
#' maximum absolute row/column marginal residual.
#' @keywords internal
#' @noRd
assess_ot_plan <- function(plan, cost_matrix, target_row, target_col,
                           tolerance = PPMWASS_OT_MARGINAL_TOLERANCE) {
  empty <- list(
    total_cost = NA_real_, accepted = FALSE,
    validation_error = "invalid_plan_structure",
    row_residual_l1 = NA_real_, row_residual_linf = NA_real_,
    col_residual_l1 = NA_real_, col_residual_linf = NA_real_,
    nonfinite_mass = NA_integer_, negative_mass = NA_integer_,
    invalid_index = NA_integer_
  )
  if (is.null(plan) || !all(c("from", "to", "mass") %in% names(plan))) {
    return(empty)
  }

  mass <- suppressWarnings(as.numeric(plan$mass))
  if (!length(mass) || length(plan$from) != length(mass) ||
      length(plan$to) != length(mass)) {
    return(empty)
  }

  # approxOT returns a dense column-major `transport.plan`. Validate the index
  # orientation at column boundaries, then use matrix row/column sums. This
  # avoids an R-level loop and repeated grouping in the publication all-pairs
  # path. Exact `transport` plans are sparse and take the general branch below.
  nr <- length(target_row)
  nc <- length(target_col)
  dense_starts <- if (nr > 0L && length(mass) == nr * nc) {
    seq.int(1L, length(mass), by = nr)
  } else integer(0)
  dense_ordered <- inherits(plan, "transport.plan") &&
    is.integer(plan$from) && is.integer(plan$to) &&
    length(dense_starts) == nc &&
    identical(plan$from[seq_len(nr)], seq_len(nr)) &&
    all(plan$from[dense_starts] == 1L) &&
    all(plan$from[dense_starts + nr - 1L] == nr) &&
    identical(plan$to[dense_starts], seq_len(nc)) &&
    identical(plan$to[dense_starts + nr - 1L], seq_len(nc))

  if (dense_ordered) {
    from <- plan$from
    to <- plan$to
    invalid_index <- 0L
  } else {
    numeric_index <- function(x) {
      if (is.factor(x)) suppressWarnings(as.numeric(as.character(x))) else
        suppressWarnings(as.numeric(x))
    }
    from_raw <- numeric_index(plan$from)
    to_raw <- numeric_index(plan$to)
    from <- suppressWarnings(as.integer(from_raw))
    to <- suppressWarnings(as.integer(to_raw))
    invalid_index <- sum(
      !is.finite(from_raw) | !is.finite(to_raw) |
        from_raw != from | to_raw != to |
        from < 1L | from > nr | to < 1L | to > nc
    )
  }
  nonfinite_mass <- sum(!is.finite(mass))
  negative_mass <- sum(is.finite(mass) & mass < -tolerance)
  empty$invalid_index <- as.integer(invalid_index)
  empty$nonfinite_mass <- as.integer(nonfinite_mass)
  empty$negative_mass <- as.integer(negative_mass)
  if (invalid_index > 0L) {
    empty$validation_error <- "invalid_plan_index"
    return(empty)
  }
  if (nonfinite_mass > 0L) {
    empty$validation_error <- "nonfinite_plan_mass"
    return(empty)
  }
  if (negative_mass > 0L) {
    empty$validation_error <- "negative_plan_mass"
    return(empty)
  }

  if (dense_ordered) {
    mass_matrix <- matrix(mass, nrow = nr, ncol = nc)
    row_mass <- rowSums(mass_matrix)
    col_mass <- colSums(mass_matrix)
  } else {
    row_mass <- numeric(nr)
    col_mass <- numeric(nc)
    row_grouped <- rowsum(mass, from, reorder = FALSE)
    col_grouped <- rowsum(mass, to, reorder = FALSE)
    row_mass[as.integer(rownames(row_grouped))] <- row_grouped[, 1L]
    col_mass[as.integer(rownames(col_grouped))] <- col_grouped[, 1L]
  }
  row_delta <- row_mass - target_row
  col_delta <- col_mass - target_col
  empty$row_residual_l1 <- sum(abs(row_delta))
  empty$row_residual_linf <- max(abs(row_delta))
  empty$col_residual_l1 <- sum(abs(col_delta))
  empty$col_residual_linf <- max(abs(col_delta))

  total_cost <- if (dense_ordered) {
    sum(mass * as.numeric(cost_matrix))
  } else {
    tryCatch(
      sum(mass * cost_matrix[cbind(from, to)]),
      error = function(e) NA_real_
    )
  }
  empty$total_cost <- total_cost
  if (!is.finite(total_cost)) {
    empty$validation_error <- "nonfinite_total_cost"
    return(empty)
  }
  max_residual <- max(empty$row_residual_linf, empty$col_residual_linf)
  if (!is.finite(max_residual) || max_residual > tolerance) {
    empty$validation_error <- "marginal_residual_exceeds_tolerance"
    return(empty)
  }
  empty$accepted <- TRUE
  empty$validation_error <- NA_character_
  empty
}

#' @keywords internal
#' @noRd
append_ot_diagnostic <- function(collector, record) {
  if (is.null(collector)) return(invisible(record))
  if (is.environment(collector)) {
    records <- collector$records
    if (is.null(records)) records <- list()
    records[[length(records) + 1L]] <- record
    collector$records <- records
  } else if (is.function(collector)) {
    collector(record)
  } else {
    stop("OT diagnostics collector must be an environment or a function.",
         call. = FALSE)
  }
  invisible(record)
}

#' @keywords internal
#' @noRd
empty_ot_attempts <- function() {
  data.frame(
    step = character(), backend = character(), niter = integer(),
    epsilon = numeric(), total_cost = numeric(), accepted = logical(),
    validation_error = character(), solver_error = character(),
    row_residual_l1 = numeric(), row_residual_linf = numeric(),
    col_residual_l1 = numeric(), col_residual_linf = numeric(),
    nonfinite_mass = integer(), negative_mass = integer(),
    invalid_index = integer(), stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
ot_attempt_row <- function(step, backend, niter, epsilon, assessment,
                           solver_error = NA_character_) {
  data.frame(
    step = step, backend = backend,
    niter = if (is.null(niter)) NA_integer_ else as.integer(niter),
    epsilon = if (is.null(epsilon)) NA_real_ else as.numeric(epsilon),
    total_cost = assessment$total_cost,
    accepted = isTRUE(assessment$accepted),
    validation_error = assessment$validation_error,
    solver_error = solver_error,
    row_residual_l1 = assessment$row_residual_l1,
    row_residual_linf = assessment$row_residual_linf,
    col_residual_l1 = assessment$col_residual_l1,
    col_residual_linf = assessment$col_residual_linf,
    nonfinite_mass = assessment$nonfinite_mass,
    negative_mass = assessment$negative_mass,
    invalid_index = assessment$invalid_index,
    stringsAsFactors = FALSE
  )
}

#' PPM-aware Wasserstein Distance for HRMS Spectra
#'
#' Computes a Wasserstein-like distance using a ppm-aware ground cost function.
#' The cost saturates at one after `transition_mult` multiples of the `ppm` base
#' scale. Exact unregularized OT is the default. Finite-iteration approximate
#' plans are retained as explicitly requested alternative backends and are
#' retried deterministically before exact fallback is considered.
#'
#' @param specA Two-column matrix (mz, intensity) for spectrum A.
#' @param specB Two-column matrix (mz, intensity) for spectrum B.
#' @param ppm PPM base-cost scale (legacy argument name; default: 15).
#' @param transition_mult Multiples of the base-cost scale before saturation.
#' @param align Whether to pre-align spectra.
#' @param ot_method One of `"sinkhorn"`, `"greenkhorn"`, `"exact"` (default),
#'   or `"exact_sparse"`. The sparse backend computes the same unregularized
#'   objective using connected components of edges below the saturated cost.
#'   It can be slower on dense inputs and does not automatically change backends.
#' @param sinkhorn_epsilon Regularisation parameter for approximate OT.
#' @param sinkhorn_niter Requested maximum approximate-solver iterations.
#' @param .diagnostics Optional environment or callback receiving diagnostics.
#' @param .context Optional caller-supplied pair context stored in diagnostics.
#' @param .marginal_tolerance Maximum accepted row/column marginal residual.
#' @param .approx_requested_only Internal diagnostic mode that executes exactly
#'   one approximate-OT call at the requested epsilon and iteration count. It
#'   never retries and never invokes exact OT. A finite raw objective is
#'   returned even when its marginal residual exceeds `.marginal_tolerance`;
#'   the failed residual check remains available in `.diagnostics`.
#' @param .approx_solver,.exact_solver Optional injected solvers for testing.
#' @return Distance value between 0 and 1.
#' @keywords internal
#' @noRd
ppm_wasserstein_distance <- function(specA, specB, ppm = 15,
                                     transition_mult = 3,
                                     align = FALSE,
                                     ot_method = "exact",
                                     sinkhorn_epsilon = 0.05,
                                     sinkhorn_niter = 100L,
                                     .diagnostics = NULL,
                                     .context = NULL,
                                     .marginal_tolerance = PPMWASS_OT_MARGINAL_TOLERANCE,
                                     .approx_requested_only = FALSE,
                                     .approx_solver = NULL,
                                     .exact_solver = NULL) {
  valid_ot_methods <- c("exact", "sinkhorn", "greenkhorn", "exact_sparse")
  if (length(ot_method) != 1L || is.na(ot_method) ||
      !ot_method %in% valid_ot_methods) {
    stop(
      "ot_method must be one of: ", paste(valid_ot_methods, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(.marginal_tolerance) != 1L ||
      !is.finite(.marginal_tolerance) || .marginal_tolerance <= 0) {
    stop(".marginal_tolerance must be one positive finite value.", call. = FALSE)
  }
  if (length(.approx_requested_only) != 1L ||
      is.na(.approx_requested_only) || !is.logical(.approx_requested_only)) {
    stop(".approx_requested_only must be one nonmissing logical value.",
         call. = FALSE)
  }
  if (isTRUE(.approx_requested_only) &&
      !ot_method %in% c("sinkhorn", "greenkhorn")) {
    stop(
      ".approx_requested_only is available only for sinkhorn or greenkhorn.",
      call. = FALSE
    )
  }
  has_approxOT <- !is.null(.approx_solver) ||
    requireNamespace("approxOT", quietly = TRUE)
  has_transport <- !is.null(.exact_solver) ||
    requireNamespace("transport", quietly = TRUE)

  attempts <- list()
  sparse_details <- NULL
  cost_matrix <- NULL
  finalize_diagnostic <- function(selected_path, selected_total,
                                  status = "ok") {
    if (is.null(.diagnostics)) return(invisible(NULL))
    attempt_df <- if (length(attempts)) do.call(rbind, attempts) else
      empty_ot_attempts()
    solver_errors <- attempt_df$solver_error[
      !is.na(attempt_df$solver_error) & nzchar(attempt_df$solver_error)
    ]
    validation_errors <- attempt_df$validation_error[
      !is.na(attempt_df$validation_error) & nzchar(attempt_df$validation_error)
    ]
    record <- list(
      context = if (is.null(.context)) list() else .context,
      requested_method = ot_method,
      requested_epsilon = if (ot_method %in% c("sinkhorn", "greenkhorn")) {
        as.numeric(sinkhorn_epsilon)
      } else NA_real_,
      requested_niter = if (ot_method %in% c("sinkhorn", "greenkhorn")) {
        as.integer(sinkhorn_niter)
      } else NA_integer_,
      requested_only = isTRUE(.approx_requested_only),
      selected_path = selected_path,
      fallback_used = !selected_path %in% c(
        "analytic_constant_cost", "empty_spectrum", "empty_aligned_spectrum",
        "empty_positive_marginal", "invalid_or_nonpositive_marginal",
        "approx_primary", "approx_requested_only",
        "approx_requested_only_invalid",
        "analytic_constant_cost_with_raw_requested_diagnostic",
        "exact_primary", "exact_sparse_primary"
      ),
      status = status,
      selected_total = selected_total,
      marginal_tolerance = .marginal_tolerance,
      query_peaks = nrow(specA),
      library_peaks = nrow(specB),
      cost_rows = if (is.null(cost_matrix)) NA_integer_ else nrow(cost_matrix),
      cost_cols = if (is.null(cost_matrix)) NA_integer_ else ncol(cost_matrix),
      cost_nonfinite = if (is.null(cost_matrix)) NA_integer_ else
        sum(!is.finite(cost_matrix)),
      original_solver_error = if (length(solver_errors)) solver_errors[[1L]] else
        NA_character_,
      original_validation_error = if (length(validation_errors)) {
        validation_errors[[1L]]
      } else NA_character_,
      attempts = attempt_df
    )
    if (!is.null(sparse_details)) record$sparse_details <- sparse_details
    append_ot_diagnostic(.diagnostics, record)
    invisible(record)
  }

  if (ot_method %in% c("exact", "exact_sparse") && !has_transport) {
    stop(
      "ot_method='exact' requires the transport package; no approximate fallback is allowed.",
      call. = FALSE
    )
  }
  if (isTRUE(.approx_requested_only) && !has_approxOT) {
    stop(
      "Requested-only approximate OT requires the approxOT package; no other backend is allowed.",
      call. = FALSE
    )
  }
  if (!has_approxOT && !has_transport) {
    stop(
      "ppm_wasserstein requires transport for exact OT or approxOT for an approximate backend.",
      call. = FALSE
    )
  }

  if (nrow(specA) == 0 || nrow(specB) == 0) {
    finalize_diagnostic("empty_spectrum", 1)
    return(1)
  }

  # --- spectrum preparation (shared with previous implementation) ---
  if (align) {
    aligned <- align_spectra(specA, specB, ppm)
    if (length(aligned$mz) == 0) {
      finalize_diagnostic("empty_aligned_spectrum", 1)
      return(1)
    }
    mz <- aligned$mz
    p <- aligned$p
    q <- aligned$q
    idx_p <- which(p > 0)
    idx_q <- which(q > 0)
    if (length(idx_p) == 0 || length(idx_q) == 0) {
      finalize_diagnostic("empty_positive_marginal", 1)
      return(1)
    }
    mz_p <- mz[idx_p]; mz_q <- mz[idx_q]
    w_p <- p[idx_p];   w_q <- q[idx_q]
  } else {
    mz_p <- specA[, 1]; mz_q <- specB[, 1]
    w_p  <- specA[, 2]; w_q  <- specB[, 2]
  }

  total_p <- sum(w_p)
  total_q <- sum(w_q)
  valid_input <- all(is.finite(mz_p)) && all(is.finite(mz_q)) &&
    all(mz_p > 0) && all(mz_q > 0) &&
    all(is.finite(w_p)) && all(is.finite(w_q)) &&
    all(w_p >= 0) && all(w_q >= 0) &&
    is.finite(total_p) && is.finite(total_q) && total_p > 0 && total_q > 0
  if (!valid_input) {
    invalid_value <- if (isTRUE(.approx_requested_only)) NA_real_ else 1
    finalize_diagnostic("invalid_or_nonpositive_marginal", invalid_value,
                        status = "invalid_input")
    return(invalid_value)
  }

  w_p <- w_p / total_p
  w_q <- w_q / total_q

  if (identical(ot_method, "exact_sparse")) {
    result <- sparse_exact_ot(
      mz_p, mz_q, w_p, w_q, ppm * transition_mult,
      tolerance = .marginal_tolerance, solver = .exact_solver,
      diagnostics = !is.null(.diagnostics)
    )
    attempts <- result$attempts
    sparse_details <- result$details
    if (!is.finite(result$total)) {
      finalize_diagnostic("exact_sparse_failed", NA_real_, status = "error")
      stop("Sparse exact OT failed plan validation; no alternate backend was used.",
           call. = FALSE)
    }
    value <- min(max(result$total, 0), 1)
    finalize_diagnostic("exact_sparse_primary", value)
    return(value)
  }

  # --- ppm-aware ground cost matrix ---
  mz_mean  <- outer(mz_p, mz_q, function(a, b) (a + b) / 2)
  mz_diff  <- abs(outer(mz_p, mz_q, "-"))
  delta_ppm <- mz_diff / mz_mean * 1e6
  cost_matrix <- pmin(delta_ppm / (ppm * transition_mult), 1)
  if (!length(cost_matrix) || any(!is.finite(cost_matrix))) {
    finalize_diagnostic("invalid_cost_matrix", NA_real_, status = "error")
    stop("PPM Wasserstein ground-cost matrix contains nonfinite values.",
         call. = FALSE)
  }

  # A constant ground-cost matrix has an analytic OT value equal to that
  # constant for every valid coupling. The approxOT Sinkhorn backend can return
  # NaN for the common all-ones case; use the exact identity rather than
  # replacing the resulting distance later at matrix level.
  if (!isTRUE(.approx_requested_only) &&
      length(cost_matrix) && all(is.finite(cost_matrix)) &&
      all(cost_matrix == cost_matrix[[1]])) {
    value <- min(max(cost_matrix[[1]], 0), 1)
    finalize_diagnostic("analytic_constant_cost", value)
    return(value)
  }

  # --- optimal transport computation ---
  approx_solver <- if (!is.null(.approx_solver)) .approx_solver else
    function(...) approxOT::transport_plan_given_C(...)
  exact_solver <- if (!is.null(.exact_solver)) .exact_solver else
    function(...) transport::transport(...)

  run_approx <- function(niter, method, step) {
    solver_error <- NA_character_
    plan <- tryCatch(
      approx_solver(
        mass_x = w_p, mass_y = w_q, p = 1, cost = cost_matrix,
        method = method, epsilon = sinkhorn_epsilon,
        niter = as.integer(niter)
      ),
      error = function(e) {
        solver_error <<- conditionMessage(e)
        NULL
      }
    )
    assessment <- assess_ot_plan(
      plan, cost_matrix, w_p, w_q, tolerance = .marginal_tolerance
    )
    if (!is.na(solver_error)) assessment$validation_error <- "solver_error"
    if (!is.null(.diagnostics)) {
      attempts[[length(attempts) + 1L]] <<- ot_attempt_row(
        step, paste0("approxOT_", method), niter, sinkhorn_epsilon,
        assessment, solver_error
      )
    }
    assessment
  }

  run_exact <- function(step) {
    solver_error <- NA_character_
    plan <- tryCatch(
      exact_solver(a = w_p, b = w_q, costm = cost_matrix),
      error = function(e) {
        solver_error <<- conditionMessage(e)
        NULL
      }
    )
    assessment <- assess_ot_plan(
      plan, cost_matrix, w_p, w_q, tolerance = .marginal_tolerance
    )
    if (!is.na(solver_error)) assessment$validation_error <- "solver_error"
    if (!is.null(.diagnostics)) {
      attempts[[length(attempts) + 1L]] <<- ot_attempt_row(
        step, "transport_exact", NULL, NULL, assessment, solver_error
      )
    }
    assessment
  }

  run_approx_sequence <- function(method, primary_step = "approx_primary",
                                  retry_prefix = "approx_retry_") {
    retry_iterations <- unique(c(
      as.integer(sinkhorn_niter),
      max(as.integer(sinkhorn_niter), 500L),
      max(as.integer(sinkhorn_niter), 2000L)
    ))
    for (attempt_index in seq_along(retry_iterations)) {
      retry_niter <- retry_iterations[[attempt_index]]
      step <- if (attempt_index == 1L) primary_step else
        paste0(retry_prefix, retry_niter)
      assessment <- run_approx(retry_niter, method, step)
      if (isTRUE(assessment$accepted)) {
        return(list(total = assessment$total_cost, path = step))
      }
    }
    list(total = NA_real_, path = NA_character_)
  }

  # "exact" always uses transport::transport() when available. Approximate
  # methods retry at the requested epsilon and then use exact OT only if every
  # finite/nonfinite plan fails the explicit marginal-residual criterion.
  use_approxOT_path <- has_approxOT && ot_method %in% c("sinkhorn", "greenkhorn")
  selected <- list(total = NA_real_, path = NA_character_, status = "ok")
  if (isTRUE(.approx_requested_only)) {
    raw <- run_approx(
      as.integer(sinkhorn_niter), ot_method, "approx_requested_only"
    )
    constant_ground_cost <- length(cost_matrix) &&
      all(cost_matrix == cost_matrix[[1]])
    if (constant_ground_cost) {
      selected <- list(
        total = as.numeric(cost_matrix[[1]]),
        path = "analytic_constant_cost_with_raw_requested_diagnostic",
        status = if (isTRUE(raw$accepted)) {
          "analytic_constant_cost_raw_diagnostic_ok"
        } else {
          "analytic_constant_cost_raw_diagnostic_invalid"
        }
      )
    } else {
      raw_usable <- is.finite(raw$total_cost) &&
        (is.na(raw$validation_error) ||
         identical(raw$validation_error, "marginal_residual_exceeds_tolerance"))
      if (raw_usable) {
        selected <- list(
          total = raw$total_cost,
          path = "approx_requested_only",
          status = "ok"
        )
      }
    }
  } else if (use_approxOT_path) {
    selected <- run_approx_sequence(ot_method)
    if (!is.finite(selected$total) && has_transport) {
      exact <- run_exact("exact_fallback")
      if (isTRUE(exact$accepted)) {
        selected <- list(total = exact$total_cost, path = "exact_fallback")
      }
    }
  } else if (has_transport) {
    step <- if (identical(ot_method, "exact")) "exact_primary" else
      "exact_no_approx_backend"
    exact <- run_exact(step)
    if (isTRUE(exact$accepted)) {
      selected <- list(total = exact$total_cost, path = step)
    }
  } else {
    # Never change the requested estimand silently. In particular, replacing
    # exact unregularized OT by finite-iteration Sinkhorn can introduce both
    # regularization bias and query/library order dependence.
    stop(
      "ot_method='exact' requires the transport package; no approximate fallback is allowed.",
      call. = FALSE
    )
  }

  if (!is.finite(selected$total)) {
    if (isTRUE(.approx_requested_only)) {
      finalize_diagnostic(
        "approx_requested_only_invalid", NA_real_, status = "raw_invalid"
      )
      return(NA_real_)
    }
    finalize_diagnostic("failed_after_fallback", NA_real_, status = "error")
    stop(
      paste0(
        "Optimal-transport computation remained invalid after deterministic ",
        "fallback (nonfinite value, malformed plan, or marginal residual > ",
        format(.marginal_tolerance, scientific = TRUE), ")."
      ),
      call. = FALSE
    )
  }

  value <- min(max(selected$total, 0), 1)
  selected_status <- if (is.null(selected$status)) "ok" else selected$status
  finalize_diagnostic(selected$path, value, status = selected_status)
  value
}

#' Execute one raw requested finite-iteration ppm-Wasserstein calculation
#'
#' This publication-diagnostic helper deliberately bypasses the ordinary
#' retry/exact-fallback policy. It invokes approxOT exactly once using the
#' supplied method, epsilon, and iteration count. Use `.diagnostics` to audit
#' the selected path, plan validity, and marginal residuals. For a constant
#' ground-cost matrix, the raw call is still executed and audited, while the
#' returned distance uses the coupling-independent analytic constant; this is
#' an input identity, not a retry, fallback backend, or imputation.
#' @keywords internal
#' @noRd
ppm_wasserstein_requested_approx <- function(
    specA, specB, ppm = 15, transition_mult = 3, align = FALSE,
    ot_method = "sinkhorn", sinkhorn_epsilon = 0.05,
    sinkhorn_niter = 100L, .diagnostics = NULL, .context = NULL,
    .marginal_tolerance = PPMWASS_OT_MARGINAL_TOLERANCE,
    .approx_solver = NULL) {
  ppm_wasserstein_distance(
    specA = specA, specB = specB, ppm = ppm,
    transition_mult = transition_mult, align = align,
    ot_method = ot_method, sinkhorn_epsilon = sinkhorn_epsilon,
    sinkhorn_niter = sinkhorn_niter, .diagnostics = .diagnostics,
    .context = .context, .marginal_tolerance = .marginal_tolerance,
    .approx_requested_only = TRUE, .approx_solver = .approx_solver
  )
}

#' @keywords internal
#' @noRd
entropy_distance <- function(specA, specB, ppm = 20) {
  if (!requireNamespace("msentropy", quietly = TRUE)) {
    warning("msentropy package not available, using Hellinger instead")
    return(hellinger_distance(specA, specB, ppm))
  }

  if (nrow(specA) == 0 || nrow(specB) == 0) return(1)

  sim <- msentropy::msentropy_similarity(
    specA, specB,
    ms2_tolerance_in_da = -1,
    ms2_tolerance_in_ppm = ppm,
    clean_spectra = FALSE,
    weighted = TRUE
  )

  1 - sim
}

#' @keywords internal
#' @noRd
entropy_unweighted_distance <- function(specA, specB, ppm = 20) {
  if (!requireNamespace("msentropy", quietly = TRUE)) {
    warning("msentropy package not available, using Hellinger instead")
    return(hellinger_distance(specA, specB, ppm))
  }

  if (nrow(specA) == 0 || nrow(specB) == 0) return(1)

  sim <- msentropy::msentropy_similarity(
    specA, specB,
    ms2_tolerance_in_da = -1,
    ms2_tolerance_in_ppm = ppm,
    clean_spectra = FALSE,
    weighted = FALSE
  )

  1 - sim
}

#' @keywords internal
#' @noRd
compute_distance <- function(specA, specB, method = "hellinger", ppm = 20,
                             align_wasserstein = FALSE, mass_power = 3,
                             intensity_power = 0.5,
                             wasserstein_transition_mult = 3,
                             ot_method = "exact",
                             sinkhorn_epsilon = 0.05,
                             sinkhorn_niter = 100L,
                             ot_diagnostics = NULL,
                             ot_context = NULL,
                             ot_marginal_tolerance = PPMWASS_OT_MARGINAL_TOLERANCE) {
  switch(method,
         "hellinger" = hellinger_distance(specA, specB, ppm),
         "wasserstein" = wasserstein_distance(specA, specB, ppm,
                                              align = align_wasserstein),
         "ppm_wasserstein" = ppm_wasserstein_distance(
           specA, specB, ppm,
           transition_mult = wasserstein_transition_mult,
           align = align_wasserstein,
           ot_method = ot_method,
           sinkhorn_epsilon = sinkhorn_epsilon,
           sinkhorn_niter = sinkhorn_niter,
           .diagnostics = ot_diagnostics,
           .context = ot_context,
           .marginal_tolerance = ot_marginal_tolerance
         ),
         "cosine" = cosine_distance(specA, specB, ppm),
         "entropy" = entropy_distance(specA, specB, ppm),
         "entropy_weighted" = entropy_distance(specA, specB, ppm),
         "entropy_unweighted" = entropy_unweighted_distance(specA, specB, ppm),
         "weighted_cosine" = weighted_cosine_distance(specA, specB, ppm, mass_power, intensity_power),
         "composite" = composite_distance(specA, specB, ppm, mass_power, intensity_power),
         hellinger_distance(specA, specB, ppm)
  )
}
