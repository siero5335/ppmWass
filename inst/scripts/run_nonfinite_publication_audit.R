#!/usr/bin/env Rscript

# Publication audit of raw approximate-OT failures and all input conditions
# that bypass or prevent the raw Sinkhorn call.  The selected distance exactly
# follows the package path used by the main retrieval run:
#
#   empty/nonpositive marginal -> 1
#   constant cost              -> analytic constant
#   requested Sinkhorn         -> 500 -> 2,000 iterations -> exact OT
#   each selected channel      -> clamp to [0, 1]
#
# Alternative epsilon, strict-residual, and Greenkhorn runs are diagnostic
# only; they never change the publication distance.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  substring(hit[[1L]], nchar(prefix) + 1L)
}
split_arg <- function(name, default) {
  value <- get_arg(name, default)
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

repo_dir <- normalizePath(get_arg("repo", getwd()), mustWork = TRUE)
main_dir <- normalizePath(get_arg("main-dir"), mustWork = TRUE)
output_dir <- normalizePath(
  get_arg("output-dir", file.path(dirname(main_dir), "..", "audit")),
  mustWork = FALSE
)
datasets <- split_arg("datasets", "RECETOX,HREI-MSDB")
n_cores <- as.integer(get_arg("n-cores", "8"))
ppm <- as.numeric(get_arg("base-cost-ppm", "15"))
transition_mult <- as.numeric(get_arg("transition-mult", "3"))
epsilon <- as.numeric(get_arg("sinkhorn-epsilon", "0.05"))
primary_niter <- as.integer(get_arg("sinkhorn-niter", "100"))
strict_niter <- as.integer(get_arg("strict-niter", "5000"))
strict_residual_tolerance <- as.numeric(
  get_arg("strict-residual-tolerance", "1e-10")
)
alternate_backend <- tolower(get_arg("alternate-backend", "greenkhorn"))
main_reproduction_tolerance <- as.numeric(
  get_arg("main-reproduction-tolerance", "1e-12")
)
cluster_boot <- as.integer(get_arg("cluster-boot", "20000"))
seed <- as.integer(get_arg("seed", "20260718"))
w_frag <- as.numeric(get_arg("w-frag", "0.70"))
w_loss <- as.numeric(get_arg("w-loss", "0.30"))

if (!length(datasets) || any(!datasets %in% c("RECETOX", "HREI-MSDB"))) {
  stop("datasets must be RECETOX and/or HREI-MSDB.")
}
numeric_positive <- c(
  n_cores = n_cores, ppm = ppm, transition_mult = transition_mult,
  epsilon = epsilon, primary_niter = primary_niter,
  strict_niter = strict_niter, strict_residual_tolerance = strict_residual_tolerance,
  cluster_boot = cluster_boot
)
if (any(!is.finite(numeric_positive)) || any(numeric_positive <= 0)) {
  stop("Core counts, solver settings, and tolerances must be finite and positive.")
}
if (!is.finite(seed) || !is.finite(main_reproduction_tolerance) ||
    main_reproduction_tolerance < 0 || !is.finite(w_frag) || !is.finite(w_loss) ||
    w_frag < 0 || w_loss < 0 || w_frag + w_loss <= 0) {
  stop("Invalid seed, main reproduction tolerance, or channel weights.")
}
if (!alternate_backend %in% c("greenkhorn", "none")) {
  stop("alternate-backend must be greenkhorn or none.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
if (!requireNamespace("approxOT", quietly = TRUE)) stop("approxOT is required")
if (!requireNamespace("transport", quietly = TRUE)) stop("transport is required")
pkgload::load_all(repo_dir, quiet = TRUE)

stats_helper <- file.path(
  repo_dir, "inst", "scripts", "lib", "publication_retrieval_statistics.R"
)
if (!file.exists(stats_helper)) stop("Missing statistics helper: ", stats_helper)
source(stats_helper, local = TRUE)

per_query_file <- file.path(main_dir, "combined_retrieval_per_query_tol15.csv")
if (!file.exists(per_query_file)) {
  stop(
    "Required combined per-query results are missing: ", per_query_file,
    ". Run the main finalizer before the nonfinite audit."
  )
}

dataset_dir_name <- function(dataset) gsub("[^A-Za-z0-9]+", "_", dataset)
clamp_channel <- function(x) {
  if (!length(x) || !is.finite(x[[1L]])) return(NA_real_)
  min(max(as.numeric(x[[1L]]), 0), 1)
}
safe_sum <- function(x) {
  if (is.null(x) || !length(x)) return(0)
  suppressWarnings(sum(x))
}
safe_min <- function(x) {
  if (is.null(x) || !length(x) || !any(is.finite(x))) return(NA_real_)
  min(x[is.finite(x)])
}
safe_max <- function(x) {
  if (is.null(x) || !length(x) || !any(is.finite(x))) return(NA_real_)
  max(x[is.finite(x)])
}

coerce_spectrum <- function(x, label) {
  if (is.null(x)) {
    return(list(
      ok = TRUE, empty = TRUE, matrix = matrix(numeric(), nrow = 0L, ncol = 2L),
      reason = paste0(label, "_is_null")
    ))
  }
  converted <- tryCatch(as.matrix(x), error = function(e) e)
  if (inherits(converted, "error")) {
    return(list(
      ok = FALSE, empty = FALSE, matrix = matrix(numeric(), nrow = 0L, ncol = 2L),
      reason = paste0(label, "_matrix_conversion_error: ", conditionMessage(converted))
    ))
  }
  if (is.null(dim(converted)) || length(dim(converted)) != 2L || ncol(converted) < 2L) {
    return(list(
      ok = FALSE, empty = FALSE, matrix = matrix(numeric(), nrow = 0L, ncol = 2L),
      reason = paste0(label, "_must_be_a_two_column_matrix")
    ))
  }
  if (!nrow(converted)) {
    return(list(
      ok = TRUE, empty = TRUE, matrix = matrix(numeric(), nrow = 0L, ncol = 2L),
      reason = paste0(label, "_has_zero_rows")
    ))
  }
  first_two <- converted[, 1:2, drop = FALSE]
  numeric_values <- suppressWarnings(as.numeric(first_two))
  numeric_matrix <- matrix(
    numeric_values, nrow = nrow(first_two), ncol = 2L,
    dimnames = list(rownames(first_two), colnames(first_two))
  )
  list(ok = TRUE, empty = FALSE, matrix = numeric_matrix, reason = NA_character_)
}

prepare_channel <- function(a, b) {
  aa <- coerce_spectrum(a, "query")
  bb <- coerce_spectrum(b, "library")
  base <- list(
    valid = FALSE, solver_ready = FALSE, package_return_one = FALSE,
    preparation_status = "malformed_input", preparation_reason = NA_character_,
    a = aa$matrix, b = bb$matrix, wa_raw = numeric(), wb_raw = numeric(),
    wa = numeric(), wb = numeric(), cost = matrix(numeric(), 0L, 0L),
    constant = FALSE, finite_marginals = FALSE,
    nonnegative_marginals = FALSE, positive_marginal_totals = FALSE
  )
  # Match ppm_wasserstein_distance(): an empty argument returns 1 before any
  # peak-column or marginal inspection of the other argument.
  if (aa$empty || bb$empty) {
    base$preparation_status <- "empty_channel"
    base$preparation_reason <- paste(
      na.omit(c(if (aa$empty) aa$reason, if (bb$empty) bb$reason)),
      collapse = "; "
    )
    base$package_return_one <- TRUE
    return(base)
  }
  if (!aa$ok || !bb$ok) {
    base$preparation_reason <- paste(
      na.omit(c(if (!aa$ok) aa$reason, if (!bb$ok) bb$reason)),
      collapse = "; "
    )
    return(base)
  }

  base$a <- aa$matrix
  base$b <- bb$matrix
  base$wa_raw <- aa$matrix[, 2L]
  base$wb_raw <- bb$matrix[, 2L]
  finite_values <- all(is.finite(aa$matrix[, 1:2, drop = FALSE])) &&
    all(is.finite(bb$matrix[, 1:2, drop = FALSE]))
  if (!finite_values) {
    base$preparation_reason <- "nonfinite_mz_or_intensity"
    return(base)
  }
  if (any(aa$matrix[, 1L] <= 0) || any(bb$matrix[, 1L] <= 0)) {
    base$preparation_reason <- "nonpositive_mz"
    return(base)
  }

  sum_a <- sum(base$wa_raw)
  sum_b <- sum(base$wb_raw)
  base$finite_marginals <- all(is.finite(c(base$wa_raw, base$wb_raw, sum_a, sum_b)))
  base$nonnegative_marginals <- all(base$wa_raw >= 0) && all(base$wb_raw >= 0)
  base$positive_marginal_totals <- is.finite(sum_a) && is.finite(sum_b) &&
    sum_a > 0 && sum_b > 0

  if (!is.finite(sum_a) || !is.finite(sum_b)) {
    base$preparation_status <- "invalid_marginal"
    base$preparation_reason <- "nonfinite_marginal_total"
    return(base)
  }
  if (sum_a <= 0 || sum_b <= 0) {
    base$preparation_status <- "zero_or_invalid_marginal"
    base$preparation_reason <- "nonpositive_marginal_total"
    base$package_return_one <- TRUE
    return(base)
  }
  if (!base$nonnegative_marginals) {
    base$preparation_status <- "invalid_marginal"
    base$preparation_reason <- "negative_marginal_entry"
    return(base)
  }

  base$wa <- base$wa_raw / sum_a
  base$wb <- base$wb_raw / sum_b
  mean_mz <- outer(aa$matrix[, 1L], bb$matrix[, 1L], function(x, y) (x + y) / 2)
  delta <- abs(outer(aa$matrix[, 1L], bb$matrix[, 1L], "-")) / mean_mz * 1e6
  base$cost <- pmin(delta / (ppm * transition_mult), 1)
  if (!length(base$cost) || any(!is.finite(base$cost))) {
    base$preparation_reason <- "nonfinite_or_empty_ground_cost"
    return(base)
  }
  base$constant <- all(base$cost == base$cost[[1L]])
  base$valid <- TRUE
  base$solver_ready <- TRUE
  base$preparation_status <- "valid"
  base$preparation_reason <- NA_character_
  base
}

empty_residual <- function() list(
  row_l1 = NA_real_, row_linf = NA_real_, col_l1 = NA_real_,
  col_linf = NA_real_, nonfinite_mass = NA_integer_, negative_mass = NA_integer_
)

plan_residuals <- function(plan, wa, wb) {
  out <- empty_residual()
  if (is.null(plan) || is.null(plan$mass) || is.null(plan$from) ||
      is.null(plan$to)) return(out)
  mass <- as.numeric(plan$mass)
  out$nonfinite_mass <- sum(!is.finite(mass))
  out$negative_mass <- sum(is.finite(mass) & mass < 0)
  if (any(!is.finite(mass))) return(out)
  row_mass <- numeric(length(wa))
  col_mass <- numeric(length(wb))
  row_sum <- tapply(mass, plan$from, sum)
  col_sum <- tapply(mass, plan$to, sum)
  row_mass[as.integer(names(row_sum))] <- row_sum
  col_mass[as.integer(names(col_sum))] <- col_sum
  row_delta <- row_mass - wa
  col_delta <- col_mass - wb
  out$row_l1 <- sum(abs(row_delta))
  out$row_linf <- max(abs(row_delta))
  out$col_l1 <- sum(abs(col_delta))
  out$col_linf <- max(abs(col_delta))
  out
}

not_attempted <- function(backend, reason, available = FALSE) list(
  backend = backend, available = available, attempted = FALSE,
  status = "not_attempted", total = NA_real_, residual = empty_residual(),
  error = reason
)

attempt_approx <- function(prepared, niter, eps, method = "sinkhorn",
                           available = TRUE) {
  if (!available) {
    return(not_attempted(method, paste0(method, " backend unavailable"), FALSE))
  }
  if (!isTRUE(prepared$solver_ready)) {
    return(not_attempted(
      method, paste0("input_not_solver_ready: ", prepared$preparation_status), TRUE
    ))
  }
  error_text <- NA_character_
  plan <- tryCatch(
    approxOT::transport_plan_given_C(
      mass_x = prepared$wa, mass_y = prepared$wb, p = 1,
      cost = prepared$cost, method = method, epsilon = eps,
      niter = as.integer(niter)
    ),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )
  total <- if (is.null(plan)) NA_real_ else tryCatch(
    sum(plan$mass * prepared$cost[cbind(plan$from, plan$to)]),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NA_real_
    }
  )
  list(
    backend = method, available = TRUE, attempted = TRUE,
    status = if (!is.na(error_text)) "error" else if (is.finite(total)) "finite" else "nonfinite",
    total = total, residual = plan_residuals(plan, prepared$wa, prepared$wb),
    error = error_text
  )
}

attempt_exact <- function(prepared) {
  if (!isTRUE(prepared$solver_ready)) {
    return(not_attempted(
      "transport_network_simplex",
      paste0("input_not_solver_ready: ", prepared$preparation_status), TRUE
    ))
  }
  error_text <- NA_character_
  plan <- tryCatch(
    transport::transport(a = prepared$wa, b = prepared$wb, costm = prepared$cost),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )
  total <- if (is.null(plan)) NA_real_ else tryCatch(
    sum(plan$mass * prepared$cost[cbind(plan$from, plan$to)]),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NA_real_
    }
  )
  list(
    backend = "transport_network_simplex", available = TRUE, attempted = TRUE,
    status = if (!is.na(error_text)) "error" else if (is.finite(total)) "finite" else "nonfinite",
    total = total, residual = plan_residuals(plan, prepared$wa, prepared$wb),
    error = error_text
  )
}

residual_passes <- function(attempt, tolerance) {
  residual <- attempt$residual
  values <- c(
    residual$row_l1, residual$row_linf, residual$col_l1, residual$col_linf
  )
  isTRUE(attempt$attempted) && is.finite(attempt$total) &&
    all(is.finite(values)) && all(values <= tolerance) &&
    identical(as.integer(residual$nonfinite_mass), 0L) &&
    identical(as.integer(residual$negative_mass), 0L)
}

classify_failure <- function(prepared, primary, retry500, retry2000,
                             eps002, eps010, strict, alternate, exact) {
  if (identical(prepared$preparation_status, "empty_channel")) {
    return(c(
      cause = "empty_channel",
      detail = "Package bypasses OT and returns channel distance 1."
    ))
  }
  if (identical(prepared$preparation_status, "malformed_input")) {
    return(c(cause = "malformed_input", detail = prepared$preparation_reason))
  }
  if (prepared$preparation_status %in% c("invalid_marginal", "zero_or_invalid_marginal")) {
    return(c(cause = "zero_or_invalid_marginal", detail = prepared$preparation_reason))
  }
  if (isTRUE(prepared$constant) && !is.finite(primary$total)) {
    return(c(
      cause = "backend_specific_constant_cost_failure",
      detail = "Raw Sinkhorn was nonfinite on a constant cost; analytic OT identity recovers it."
    ))
  }
  if (!is.na(primary$error) && nzchar(primary$error)) {
    recovered <- c(alternate$total, exact$total)
    return(c(
      cause = if (any(is.finite(recovered))) "sinkhorn_backend_error" else "unresolved_backend_error",
      detail = primary$error
    ))
  }
  if (is.finite(retry500$total) || is.finite(retry2000$total) || is.finite(strict$total)) {
    return(c(
      cause = "solver_nonconvergence_at_primary_iterations",
      detail = "Same-epsilon Sinkhorn became finite with additional iterations."
    ))
  }
  if (is.finite(eps002$total) || is.finite(eps010$total)) {
    return(c(
      cause = "epsilon_sensitive_numerical_failure",
      detail = "Alternative epsilon produced a finite approximate-OT total."
    ))
  }
  if (is.finite(alternate$total)) {
    return(c(
      cause = "sinkhorn_backend_specific_failure",
      detail = paste0(alternate$backend, " produced a finite approximate-OT total.")
    ))
  }
  if (is.finite(exact$total)) {
    nonfinite_mass <- primary$residual$nonfinite_mass
    return(c(
      cause = if (is.finite(nonfinite_mass) && nonfinite_mass > 0)
        "possible_overflow_underflow_nonfinite_plan_mass" else
        "approximate_solver_failure_exact_recovers",
      detail = "Exact network-simplex OT was finite after approximate diagnostics failed."
    ))
  }
  c(
    cause = "unresolved_numerical_failure",
    detail = "No selected or diagnostic solver produced a finite result."
  )
}

select_main_channel <- function(prepared, primary, retry500, retry2000, exact) {
  if (isTRUE(prepared$package_return_one)) {
    step <- if (identical(prepared$preparation_status, "empty_channel"))
      "package_empty_channel_return_1" else
      "package_nonpositive_marginal_return_1"
    return(list(raw = 1, clamped = 1, step = step))
  }
  if (!isTRUE(prepared$valid)) {
    return(list(raw = NA_real_, clamped = NA_real_, step = "unresolved_invalid_input"))
  }
  if (isTRUE(prepared$constant)) {
    raw <- prepared$cost[[1L]]
    return(list(raw = raw, clamped = clamp_channel(raw), step = "analytic_constant_cost"))
  }
  attempts <- list(
    primary = primary, same_epsilon_retry_500 = retry500,
    same_epsilon_retry_2000 = retry2000, exact_transport = exact
  )
  for (step in names(attempts)) {
    value <- attempts[[step]]$total
    if (is.finite(value)) {
      return(list(raw = value, clamped = clamp_channel(value), step = step))
    }
  }
  list(raw = NA_real_, clamped = NA_real_, step = "unresolved")
}

diagnose_channel <- function(dataset, query_id, library_id, i, j, channel,
                             prepared, primary) {
  can_run <- isTRUE(prepared$solver_ready)
  retry500 <- if (can_run) {
    attempt_approx(prepared, max(primary_niter, 500L), epsilon, "sinkhorn")
  } else not_attempted("sinkhorn", "input_not_solver_ready", TRUE)
  retry2000 <- if (can_run) {
    attempt_approx(prepared, max(primary_niter, 2000L), epsilon, "sinkhorn")
  } else not_attempted("sinkhorn", "input_not_solver_ready", TRUE)
  eps002 <- if (can_run) {
    attempt_approx(prepared, primary_niter, 0.02, "sinkhorn")
  } else not_attempted("sinkhorn", "input_not_solver_ready", TRUE)
  eps010 <- if (can_run) {
    attempt_approx(prepared, primary_niter, 0.10, "sinkhorn")
  } else not_attempted("sinkhorn", "input_not_solver_ready", TRUE)
  strict <- if (can_run) {
    attempt_approx(prepared, max(primary_niter, strict_niter), epsilon, "sinkhorn")
  } else not_attempted("sinkhorn", "input_not_solver_ready", TRUE)
  alternate_available <- !identical(alternate_backend, "none")
  alternate <- if (can_run && alternate_available) {
    attempt_approx(
      prepared, primary_niter, epsilon, alternate_backend,
      available = alternate_available
    )
  } else {
    not_attempted(
      alternate_backend,
      if (can_run) "alternate_backend_disabled" else "input_not_solver_ready",
      alternate_available
    )
  }
  exact <- if (can_run) attempt_exact(prepared) else
    not_attempted("transport_network_simplex", "input_not_solver_ready", TRUE)
  selected <- select_main_channel(prepared, primary, retry500, retry2000, exact)
  classified <- classify_failure(
    prepared, primary, retry500, retry2000, eps002, eps010,
    strict, alternate, exact
  )
  r <- primary$residual
  sr <- strict$residual
  ar <- alternate$residual
  er <- exact$residual
  audit_trigger <- if (!isTRUE(prepared$valid)) {
    prepared$preparation_status
  } else {
    "primary_nonfinite"
  }

  data.frame(
    dataset = dataset,
    query_id = query_id,
    library_id = library_id,
    query_index = i,
    library_index = j,
    channel = channel,
    audit_trigger = audit_trigger,
    failure_cause = unname(classified[["cause"]]),
    failure_cause_detail = unname(classified[["detail"]]),
    preparation_status = prepared$preparation_status,
    preparation_reason = prepared$preparation_reason,
    query_peaks = nrow(prepared$a),
    library_peaks = nrow(prepared$b),
    query_intensity_sum_raw = safe_sum(prepared$wa_raw),
    library_intensity_sum_raw = safe_sum(prepared$wb_raw),
    query_intensity_sum_normalized = if (length(prepared$wa)) safe_sum(prepared$wa) else NA_real_,
    library_intensity_sum_normalized = if (length(prepared$wb)) safe_sum(prepared$wb) else NA_real_,
    finite_positive_marginals = prepared$valid,
    finite_marginals = prepared$finite_marginals,
    nonnegative_marginals = prepared$nonnegative_marginals,
    positive_marginal_totals = prepared$positive_marginal_totals,
    package_return_one = prepared$package_return_one,
    solver_ready = prepared$solver_ready,
    cost_rows = nrow(prepared$cost),
    cost_cols = ncol(prepared$cost),
    cost_min = safe_min(prepared$cost),
    cost_max = safe_max(prepared$cost),
    cost_nonfinite = sum(!is.finite(prepared$cost)),
    constant_cost = prepared$constant,
    primary_attempted = primary$attempted,
    primary_status = primary$status,
    primary_total = primary$total,
    primary_row_residual_l1 = r$row_l1,
    primary_row_residual_linf = r$row_linf,
    primary_col_residual_l1 = r$col_l1,
    primary_col_residual_linf = r$col_linf,
    primary_nonfinite_mass = r$nonfinite_mass,
    primary_negative_mass = r$negative_mass,
    primary_error = primary$error,
    retry_500_attempted = retry500$attempted,
    retry_500_status = retry500$status,
    retry_500_total = retry500$total,
    retry_500_error = retry500$error,
    retry_2000_attempted = retry2000$attempted,
    retry_2000_status = retry2000$status,
    retry_2000_total = retry2000$total,
    retry_2000_error = retry2000$error,
    epsilon_002_attempted = eps002$attempted,
    epsilon_002_status = eps002$status,
    epsilon_002_total = eps002$total,
    epsilon_002_error = eps002$error,
    epsilon_010_attempted = eps010$attempted,
    epsilon_010_status = eps010$status,
    epsilon_010_total = eps010$total,
    epsilon_010_error = eps010$error,
    strict_tolerance_argument_available = FALSE,
    strict_tolerance_equivalent =
      "fixed-iteration Sinkhorn followed by explicit marginal-residual threshold",
    strict_iterations = max(primary_niter, strict_niter),
    strict_residual_tolerance = strict_residual_tolerance,
    strict_attempted = strict$attempted,
    strict_status = strict$status,
    strict_total = strict$total,
    strict_residual_pass = residual_passes(strict, strict_residual_tolerance),
    strict_row_residual_l1 = sr$row_l1,
    strict_row_residual_linf = sr$row_linf,
    strict_col_residual_l1 = sr$col_l1,
    strict_col_residual_linf = sr$col_linf,
    strict_error = strict$error,
    alternate_backend = alternate_backend,
    alternate_backend_available = alternate$available,
    alternate_backend_attempted = alternate$attempted,
    alternate_backend_status = alternate$status,
    alternate_backend_total = alternate$total,
    alternate_backend_row_residual_l1 = ar$row_l1,
    alternate_backend_row_residual_linf = ar$row_linf,
    alternate_backend_col_residual_l1 = ar$col_l1,
    alternate_backend_col_residual_linf = ar$col_linf,
    alternate_backend_error = alternate$error,
    exact_attempted = exact$attempted,
    exact_status = exact$status,
    exact_total = exact$total,
    exact_row_residual_l1 = er$row_l1,
    exact_row_residual_linf = er$row_linf,
    exact_col_residual_l1 = er$col_l1,
    exact_col_residual_linf = er$col_linf,
    exact_error = exact$error,
    selected_total_raw = selected$raw,
    selected_total = selected$clamped,
    selected_step = selected$step,
    main_combined_distance = NA_real_,
    recomputed_combined_distance = NA_real_,
    main_reproduction_abs_difference = NA_real_,
    stringsAsFactors = FALSE
  )
}

diagnostic_prototype <- function() {
  prepared <- prepare_channel(
    matrix(numeric(), nrow = 0L, ncol = 2L),
    matrix(numeric(), nrow = 0L, ncol = 2L)
  )
  primary <- not_attempted("sinkhorn", "prototype", TRUE)
  diagnose_channel("", "", "", 0L, 0L, "", prepared, primary)[0, , drop = FALSE]
}

select_unflagged_channel <- function(prepared, primary) {
  if (isTRUE(prepared$constant)) {
    value <- prepared$cost[[1L]]
    return(list(raw = value, clamped = clamp_channel(value), step = "analytic_constant_cost"))
  }
  if (is.finite(primary$total)) {
    return(list(
      raw = primary$total, clamped = clamp_channel(primary$total), step = "primary"
    ))
  }
  stop("Internal error: an unflagged channel has no finite main-equivalent value.")
}

scan_dataset <- function(dataset) {
  ds_dir <- file.path(main_dir, dataset_dir_name(dataset))
  spectra_path <- file.path(ds_dir, "spectra_filtered_tol15.rds")
  matrix_path <- file.path(ds_dir, "dist_ppm_wasserstein_tol15.rds")
  if (!file.exists(spectra_path) || !file.exists(matrix_path)) {
    stop("Missing main nonfinite-audit input for ", dataset, ": ", ds_dir)
  }
  spectra <- readRDS(spectra_path)
  main_dm <- as.matrix(readRDS(matrix_path))
  ids <- as.character(spectra$df_spec$id)
  if (anyDuplicated(ids) || length(ids) != length(spectra$frag_list) ||
      length(ids) != length(spectra$loss_list)) {
    stop(dataset, " spectra IDs are duplicated or channel lengths are inconsistent.")
  }
  names(spectra$frag_list) <- ids
  names(spectra$loss_list) <- ids
  if (!all(dim(main_dm) == c(length(ids), length(ids)))) {
    stop(dataset, " main distance matrix dimensions do not match spectra.")
  }
  if (!is.null(rownames(main_dm)) && !is.null(colnames(main_dm))) {
    if (!setequal(ids, rownames(main_dm)) || !setequal(ids, colnames(main_dm))) {
      stop(dataset, " main distance matrix IDs do not match spectra.")
    }
    main_dm <- main_dm[ids, ids, drop = FALSE]
  }
  n <- length(ids)
  total_pairs <- n * (n - 1L) / 2L
  message("[", dataset, "] raw Sinkhorn scan over ", total_pairs,
          " unique unordered pairs")

  scan_i <- function(i) {
    if (i >= n) return(NULL)
    rows <- list()
    for (j in (i + 1L):n) {
      channels <- list(
        fragment = prepare_channel(spectra$frag_list[[i]], spectra$frag_list[[j]]),
        pooled_derived = prepare_channel(spectra$loss_list[[i]], spectra$loss_list[[j]])
      )
      primary <- lapply(channels, function(prepared) {
        if (isTRUE(prepared$solver_ready)) {
          attempt_approx(prepared, primary_niter, epsilon, "sinkhorn")
        } else {
          not_attempted(
            "sinkhorn", paste0("input_not_solver_ready: ", prepared$preparation_status),
            TRUE
          )
        }
      })
      triggered <- vapply(names(channels), function(channel) {
        !isTRUE(channels[[channel]]$valid) || !is.finite(primary[[channel]]$total)
      }, logical(1))
      if (!any(triggered)) next

      channel_rows <- list()
      selected <- list()
      for (channel in names(channels)) {
        if (triggered[[channel]]) {
          row <- diagnose_channel(
            dataset, ids[[i]], ids[[j]], i, j, channel,
            channels[[channel]], primary[[channel]]
          )
          channel_rows[[channel]] <- row
          selected[[channel]] <- list(
            raw = row$selected_total_raw[[1L]],
            clamped = row$selected_total[[1L]],
            step = row$selected_step[[1L]]
          )
        } else {
          selected[[channel]] <- select_unflagged_channel(
            channels[[channel]], primary[[channel]]
          )
        }
      }

      selected_values <- vapply(selected, function(x) x$clamped, numeric(1))
      recomputed <- if (all(is.finite(selected_values))) {
        sqrt(
          w_frag * selected_values[["fragment"]]^2 +
            w_loss * selected_values[["pooled_derived"]]^2
        )
      } else {
        NA_real_
      }
      main_value <- main_dm[i, j]
      difference <- if (is.finite(main_value) && is.finite(recomputed)) {
        abs(main_value - recomputed)
      } else {
        NA_real_
      }
      for (channel in names(channel_rows)) {
        channel_rows[[channel]]$main_combined_distance <- main_value
        channel_rows[[channel]]$recomputed_combined_distance <- recomputed
        channel_rows[[channel]]$main_reproduction_abs_difference <- difference
        rows[[length(rows) + 1L]] <- channel_rows[[channel]]
      }
    }
    if (!length(rows)) NULL else do.call(rbind, rows)
  }

  by_i <- if (.Platform$OS.type != "windows" && n_cores > 1L && n > 1L) {
    parallel::mclapply(
      seq_len(n - 1L), scan_i, mc.cores = n_cores, mc.preschedule = TRUE
    )
  } else if (n > 1L) {
    lapply(seq_len(n - 1L), scan_i)
  } else {
    list()
  }
  by_i <- by_i[vapply(by_i, function(x) !is.null(x) && nrow(x), logical(1))]
  diagnostics <- if (length(by_i)) do.call(rbind, by_i) else diagnostic_prototype()
  list(
    diagnostics = diagnostics,
    scan = data.frame(
      dataset = dataset,
      spectra = n,
      total_unordered_pairs_scanned = total_pairs,
      total_symmetric_offdiagonal_cells = 2L * total_pairs,
      stringsAsFactors = FALSE
    )
  )
}

scan_results <- lapply(datasets, scan_dataset)
diagnostics <- do.call(rbind, lapply(scan_results, `[[`, "diagnostics"))
scan_manifest <- do.call(rbind, lapply(scan_results, `[[`, "scan"))
rownames(diagnostics) <- NULL

utils::write.csv(
  diagnostics, file.path(output_dir, "nonfinite_pair_diagnostics.csv"),
  row.names = FALSE, na = ""
)

max_finite_or <- function(x, default = NA_real_) {
  x <- x[is.finite(x)]
  if (!length(x)) default else max(x)
}

summary_rows <- lapply(datasets, function(dataset) {
  x <- diagnostics[diagnostics$dataset == dataset, , drop = FALSE]
  scan <- scan_manifest[scan_manifest$dataset == dataset, , drop = FALSE]
  keys <- if (nrow(x)) unique(paste(x$query_index, x$library_index, sep = "|")) else character()
  compared_keys <- if (nrow(x)) unique(paste(
    x$query_index[is.finite(x$main_reproduction_abs_difference)],
    x$library_index[is.finite(x$main_reproduction_abs_difference)], sep = "|"
  )) else character()
  violations <- if (nrow(x)) {
    x$main_reproduction_abs_difference > main_reproduction_tolerance
  } else logical()
  data.frame(
    dataset = dataset,
    spectra = scan$spectra[[1L]],
    total_unordered_pairs_scanned = scan$total_unordered_pairs_scanned[[1L]],
    total_symmetric_offdiagonal_cells = scan$total_symmetric_offdiagonal_cells[[1L]],
    diagnostic_channels = nrow(x),
    primary_nonfinite_channels = sum(x$audit_trigger == "primary_nonfinite"),
    invalid_or_empty_channels = sum(x$audit_trigger != "primary_nonfinite"),
    unique_unordered_pairs = length(keys),
    symmetric_matrix_cells = 2L * length(keys),
    analytic_constant_cost = sum(x$selected_step == "analytic_constant_cost"),
    package_empty_channel_return_1 = sum(x$selected_step == "package_empty_channel_return_1"),
    package_nonpositive_marginal_return_1 =
      sum(x$selected_step == "package_nonpositive_marginal_return_1"),
    same_epsilon_retry_500 = sum(x$selected_step == "same_epsilon_retry_500"),
    same_epsilon_retry_2000 = sum(x$selected_step == "same_epsilon_retry_2000"),
    exact_transport = sum(x$selected_step == "exact_transport"),
    strict_diagnostic_finite = sum(is.finite(x$strict_total)),
    strict_diagnostic_residual_pass = sum(x$strict_residual_pass, na.rm = TRUE),
    alternate_backend_finite = sum(is.finite(x$alternate_backend_total)),
    unresolved_channels = sum(!is.finite(x$selected_total)),
    reproduction_compared_pairs = length(compared_keys),
    main_reproduction_tolerance = main_reproduction_tolerance,
    reproduction_violating_channels = sum(violations, na.rm = TRUE),
    max_abs_main_recomputed_difference = max_finite_or(
      x$main_reproduction_abs_difference, if (nrow(x)) NA_real_ else 0
    ),
    stringsAsFactors = FALSE
  )
})
fallback_summary <- do.call(rbind, summary_rows)
utils::write.csv(
  fallback_summary, file.path(output_dir, "nonfinite_fallback_summary.csv"),
  row.names = FALSE, na = ""
)

cause_rows <- lapply(datasets, function(dataset) {
  x <- diagnostics[diagnostics$dataset == dataset, , drop = FALSE]
  if (!nrow(x)) {
    return(data.frame(
      dataset = dataset, failure_cause = "none", n_channels = 0L,
      n_unique_pairs = 0L, stringsAsFactors = FALSE
    ))
  }
  causes <- unique(x$failure_cause)
  do.call(rbind, lapply(causes, function(cause) {
    y <- x[x$failure_cause == cause, , drop = FALSE]
    data.frame(
      dataset = dataset,
      failure_cause = cause,
      n_channels = nrow(y),
      n_unique_pairs = length(unique(paste(y$query_index, y$library_index, sep = "|"))),
      stringsAsFactors = FALSE
    )
  }))
})
failure_cause_summary <- do.call(rbind, cause_rows)
utils::write.csv(
  failure_cause_summary,
  file.path(output_dir, "nonfinite_failure_cause_summary.csv"),
  row.names = FALSE, na = ""
)

affected <- if (nrow(diagnostics)) {
  endpoints <- rbind(
    diagnostics[, c("dataset", "query_id", "failure_cause"), drop = FALSE],
    setNames(
      diagnostics[, c("dataset", "library_id", "failure_cause"), drop = FALSE],
      c("dataset", "query_id", "failure_cause")
    )
  )
  endpoints <- unique(endpoints)
  grouped <- split(endpoints, interaction(
    endpoints$dataset, endpoints$query_id, drop = TRUE, lex.order = TRUE
  ))
  do.call(rbind, lapply(grouped, function(x) data.frame(
    dataset = x$dataset[[1L]],
    query_id = x$query_id[[1L]],
    failure_causes = paste(sort(unique(x$failure_cause)), collapse = ";"),
    stringsAsFactors = FALSE
  )))
} else {
  data.frame(
    dataset = character(), query_id = character(), failure_causes = character(),
    stringsAsFactors = FALSE
  )
}
rownames(affected) <- NULL
utils::write.csv(
  affected, file.path(output_dir, "nonfinite_affected_query_ids.csv"),
  row.names = FALSE, na = ""
)

per_query <- utils::read.csv(
  per_query_file, check.names = FALSE, stringsAsFactors = FALSE
)
required_per_query <- c(
  "dataset", "method", "query_id", "inchikey", "inchikey_prefix",
  "top1_fractional", "rr_fractional", "p_at_1_fractional"
)
missing_per_query <- setdiff(required_per_query, names(per_query))
if (length(missing_per_query)) {
  stop(
    "Combined per-query results lack required columns: ",
    paste(missing_per_query, collapse = ", ")
  )
}
if (!all(datasets %in% unique(per_query$dataset))) {
  stop(
    "Combined per-query results do not contain requested dataset(s): ",
    paste(setdiff(datasets, unique(per_query$dataset)), collapse = ", ")
  )
}

metric_specs <- data.frame(
  metric = c("top1", "mrr", "p_at_1"),
  value_col = c("top1_fractional", "rr_fractional", "p_at_1_fractional"),
  cluster_col = c("inchikey", "inchikey", "inchikey_prefix"),
  cluster_unit = c("full_inchikey", "full_inchikey", "inchikey_prefix"),
  stringsAsFactors = FALSE
)
influence_rows <- list()
influence_index <- 0L
for (di in seq_along(datasets)) {
  dataset <- datasets[[di]]
  dataset_data <- per_query[per_query$dataset == dataset, , drop = FALSE]
  methods <- sort(unique(as.character(dataset_data$method)))
  affected_ids <- affected$query_id[affected$dataset == dataset]
  for (mi in seq_along(methods)) {
    method <- methods[[mi]]
    x <- dataset_data[dataset_data$method == method, , drop = FALSE]
    if (anyDuplicated(x$query_id)) {
      stop("Duplicate query IDs in combined per-query results for ", dataset, "/", method)
    }
    keep <- !x$query_id %in% affected_ids
    for (si in seq_len(nrow(metric_specs))) {
      spec <- metric_specs[si, ]
      full <- cluster_bootstrap_mean(
        x, spec$value_col, spec$cluster_col,
        R = cluster_boot,
        seed = seed + di * 100000L + mi * 1000L + si
      )
      excluded <- cluster_bootstrap_mean(
        x[keep, , drop = FALSE], spec$value_col, spec$cluster_col,
        R = cluster_boot,
        seed = seed + di * 100000L + mi * 1000L + si
      )
      eligible <- is.finite(x[[spec$value_col]])
      influence_index <- influence_index + 1L
      influence_rows[[influence_index]] <- data.frame(
        dataset = dataset,
        method = method,
        metric = spec$metric,
        tie_policy = "fractional_expected",
        cluster_unit = spec$cluster_unit,
        cluster_boot_R = cluster_boot,
        exclusion_scope = "affected queries only; candidate library fixed",
        difference_direction = "affected-query-excluded minus full",
        full_estimate = full[["estimate"]],
        full_cluster_ci_low = full[["ci_low"]],
        full_cluster_ci_high = full[["ci_high"]],
        n_full = as.integer(full[["n_queries"]]),
        n_full_clusters = as.integer(full[["n_clusters"]]),
        fallback_query_excluded_estimate = excluded[["estimate"]],
        fallback_query_excluded_cluster_ci_low = excluded[["ci_low"]],
        fallback_query_excluded_cluster_ci_high = excluded[["ci_high"]],
        n_after_exclusion = as.integer(excluded[["n_queries"]]),
        n_after_exclusion_clusters = as.integer(excluded[["n_clusters"]]),
        difference = excluded[["estimate"]] - full[["estimate"]],
        n_affected_query_ids = length(intersect(unique(x$query_id), affected_ids)),
        n_affected_eligible_queries = sum(eligible & !keep),
        stringsAsFactors = FALSE
      )
    }
  }
}
influence <- do.call(rbind, influence_rows)
utils::write.csv(
  influence, file.path(output_dir, "nonfinite_affected_query_sensitivity.csv"),
  row.names = FALSE, na = ""
)

git_commit <- tryCatch(
  trimws(system2(
    "git", c("-C", repo_dir, "rev-parse", "HEAD"),
    stdout = TRUE, stderr = FALSE
  )),
  error = function(e) NA_character_
)
if (length(git_commit) != 1L || !nzchar(git_commit)) git_commit <- NA_character_
git_status <- tryCatch(
  system2(
    "git", c("-C", repo_dir, "status", "--porcelain"),
    stdout = TRUE, stderr = FALSE
  ),
  error = function(e) character()
)
exact_command <- paste(
  c(file.path(R.home("bin"), "Rscript"), commandArgs()), collapse = " "
)

parameters <- data.frame(
  parameter = c(
    "timestamp_utc", "package_commit", "package_tree_dirty", "command",
    "repo_dir", "main_dir", "datasets", "ppmWass_base_cost_ppm",
    "transition_multiplier", "saturation_width_ppm", "sinkhorn_epsilon",
    "primary_iterations", "retry_iterations", "diagnostic_epsilons",
    "strict_tolerance_argument_available", "strict_tolerance_equivalent",
    "strict_iterations", "strict_residual_tolerance", "alternate_backend",
    "automatic_fallback_order", "channel_clamp", "w_frag", "w_loss",
    "main_reproduction_tolerance", "cluster_boot_R", "seed", "n_cores",
    "approxOT_version", "transport_version"
  ),
  value = c(
    format(Sys.time(), tz = "UTC"), git_commit, length(git_status) > 0L,
    exact_command, repo_dir, main_dir, paste(datasets, collapse = ","), ppm,
    transition_mult, ppm * transition_mult, epsilon, primary_niter,
    "500,2000", "0.02,0.10 (diagnostic only)", FALSE,
    "same-epsilon fixed iterations followed by explicit marginal-residual threshold",
    max(primary_niter, strict_niter), strict_residual_tolerance,
    alternate_backend,
    paste(
      "empty/nonpositive marginal return 1; analytic constant cost;",
      "same-epsilon 500; same-epsilon 2000; exact transport; otherwise stop"
    ),
    "each selected channel clamped to [0,1] before weighted RMS",
    w_frag, w_loss, main_reproduction_tolerance, cluster_boot, seed, n_cores,
    as.character(utils::packageVersion("approxOT")),
    as.character(utils::packageVersion("transport"))
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  parameters, file.path(output_dir, "nonfinite_audit_parameters.csv"),
  row.names = FALSE, na = ""
)
writeLines(
  capture.output(sessionInfo()), file.path(output_dir, "nonfinite_sessionInfo.txt"),
  useBytes = TRUE
)

print(fallback_summary)
if (any(fallback_summary$unresolved_channels > 0L)) {
  stop("Unresolved invalid/nonfinite channels remain; see diagnostics.")
}
if (any(fallback_summary$reproduction_violating_channels > 0L)) {
  stop(
    "Audit fallback reconstruction differs from the main matrix by more than ",
    format(main_reproduction_tolerance, scientific = TRUE),
    "; see nonfinite_fallback_summary.csv."
  )
}
message("Nonfinite audit complete: ", output_dir)
