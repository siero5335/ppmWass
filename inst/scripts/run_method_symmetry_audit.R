#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]])
}

repo_dir <- normalizePath(get_arg("repo", getwd()), mustWork = TRUE)
output_file <- get_arg(
  "output",
  file.path(dirname(repo_dir), "audit", "method_symmetry_check.csv")
)
sinkhorn_diagnostic_file <- file.path(
  dirname(output_file), "sinkhorn_directionality_diagnostic.csv"
)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
pkgload::load_all(repo_dir, quiet = TRUE)

compute_one <- getFromNamespace("compute_distance", "ppmWass")
compute_requested_approx <- getFromNamespace(
  "ppm_wasserstein_requested_approx", "ppmWass"
)
is_symmetric <- getFromNamespace("distance_method_is_symmetric", "ppmWass")

normalize_spectrum <- function(mz, intensity) {
  intensity <- intensity / sum(intensity)
  cbind(mz = mz, intensity = intensity)
}

# Unequal peak counts are deliberate: the Stein-style composite denominator
# depends on the first (query) spectrum and should therefore be directional.
pairs <- list(
  unequal_3v2 = list(
    a = normalize_spectrum(c(50, 75, 100), c(0.60, 0.30, 0.10)),
    b = normalize_spectrum(c(50.0001, 100.0001), c(0.70, 0.30))
  ),
  unequal_4v3 = list(
    a = normalize_spectrum(c(60, 85, 110, 150), c(0.45, 0.25, 0.20, 0.10)),
    b = normalize_spectrum(c(60.0002, 110.0001, 150.0003), c(0.50, 0.35, 0.15))
  ),
  unequal_5v2 = list(
    a = normalize_spectrum(c(55, 70, 95, 125, 180), c(0.35, 0.25, 0.20, 0.12, 0.08)),
    b = normalize_spectrum(c(55.0002, 125.0002), c(0.65, 0.35))
  ),
  sinkhorn_stress_fragment_4v5 = list(
    a = normalize_spectrum(
      c(91.05425, 105.069885, 119.08556, 161.13245),
      c(0.0640527962151771, 0.180230435344229,
        0.12268667663643, 0.0738341491041677)
    ),
    b = normalize_spectrum(
      c(81.069855, 91.05417, 105.06977, 119.08542, 161.13237),
      c(0.0405897497350223, 0.0732282724366475,
        0.164381906032416, 0.176824758042847,
        0.0993508930752914)
    )
  ),
  sinkhorn_stress_derived_6v6 = list(
    a = normalize_spectrum(
      c(12, 14.016, 26.016, 28.031, 42.047, 100.089),
      c(0.0354009246790235, 0.058915903693447,
        0.0315229488467013, 0.0333187289949651,
        0.0233039583929867, 0.0246389415012392)
    ),
    b = normalize_spectrum(
      c(12, 14.016, 26.016, 28.031, 42.047, 86.073),
      c(0.0354873249116053, 0.0628615896782458,
        0.028761943408224, 0.0325367815459892,
        0.0317338080483333, 0.0251867765806948)
    )
  )
)

methods <- c(
  "ppm_wasserstein", "wasserstein", "composite", "entropy",
  "entropy_weighted",
  "entropy_unweighted", "cosine", "weighted_cosine", "hellinger"
)
# The publication-grade ppm-Wasserstein backend is exact OT. Approximate
# finite-iteration direction dependence is checked separately below and is
# never accepted as evidence that the mathematical distance is symmetric.
tolerance <- 1e-8
rows <- list()
k <- 0L

for (method in methods) {
  expected <- is_symmetric(method, "exact")
  for (pair_name in names(pairs)) {
    k <- k + 1L
    pair <- pairs[[pair_name]]
    d_ab <- compute_one(
      pair$a, pair$b, method = method, ppm = 15,
      align_wasserstein = FALSE, mass_power = 3, intensity_power = 0.5,
      wasserstein_transition_mult = 3, ot_method = "exact",
      sinkhorn_epsilon = 0.05, sinkhorn_niter = 100L
    )
    d_ba <- compute_one(
      pair$b, pair$a, method = method, ppm = 15,
      align_wasserstein = FALSE, mass_power = 3, intensity_power = 0.5,
      wasserstein_transition_mult = 3, ot_method = "exact",
      sinkhorn_epsilon = 0.05, sinkhorn_niter = 100L
    )
    difference <- abs(d_ab - d_ba)
    rows[[k]] <- data.frame(
      method = method,
      spectrum_pair = pair_name,
      d_ab = d_ab,
      d_ba = d_ba,
      absolute_difference = difference,
      expected_symmetric = expected,
      test_pass = if (expected) difference <= tolerance else TRUE,
      comparison_tolerance = tolerance,
      stringsAsFactors = FALSE
    )
  }
}

result <- do.call(rbind, rows)
directional_methods <- unique(result$method[!result$expected_symmetric])
directionality_demonstrated <- vapply(directional_methods, function(method) {
  any(result$absolute_difference[result$method == method] > tolerance)
}, logical(1))
if (length(directional_methods)) {
  result$directionality_demonstrated_for_method <- NA
  for (method in directional_methods) {
    result$directionality_demonstrated_for_method[result$method == method] <-
      directionality_demonstrated[[method]]
  }
} else {
  result$directionality_demonstrated_for_method <- NA
}
utils::write.csv(result, output_file, row.names = FALSE)
print(result)
if (!all(result$test_pass) || any(!directionality_demonstrated)) {
  stop("At least one method symmetry check failed; see ", output_file)
}

stress_names <- grep("^sinkhorn_stress_", names(pairs), value = TRUE)
sinkhorn_rows <- lapply(stress_names, function(pair_name) {
  pair <- pairs[[pair_name]]
  evaluate_exact <- function(a, b) {
    compute_one(
      a, b, method = "ppm_wasserstein", ppm = 15,
      align_wasserstein = FALSE, mass_power = 3, intensity_power = 0.5,
      wasserstein_transition_mult = 3, ot_method = "exact",
      sinkhorn_epsilon = 0.05, sinkhorn_niter = 100L
    )
  }
  evaluate_requested_sinkhorn <- function(a, b, orientation) {
    collector <- new.env(parent = emptyenv())
    value <- compute_requested_approx(
      a, b, ppm = 15, transition_mult = 3, align = FALSE,
      ot_method = "sinkhorn", sinkhorn_epsilon = 0.05,
      sinkhorn_niter = 100L, .diagnostics = collector,
      .context = list(spectrum_pair = pair_name, orientation = orientation)
    )
    records <- collector$records
    if (length(records) != 1L) {
      stop("Requested-only Sinkhorn audit expected exactly one diagnostic record.")
    }
    record <- records[[1L]]
    attempts <- record$attempts
    setting_match <- nrow(attempts) == 1L &&
      identical(attempts$step[[1]], "approx_requested_only") &&
      identical(attempts$backend[[1]], "approxOT_sinkhorn") &&
      identical(attempts$niter[[1]], 100L) &&
      isTRUE(all.equal(attempts$epsilon[[1]], 0.05, tolerance = 0))
    exact_attempts <- if (nrow(attempts)) {
      sum(grepl("exact", attempts$step, fixed = TRUE) |
          grepl("transport_exact", attempts$backend, fixed = TRUE))
    } else 0L
    retry_attempts <- if (nrow(attempts)) {
      sum(grepl("retry", attempts$step, fixed = TRUE))
    } else 0L
    list(
      value = value, record = record, attempt = attempts[1L, , drop = FALSE],
      setting_match = setting_match,
      provenance_gate_pass = setting_match &&
        identical(record$selected_path, "approx_requested_only") &&
        !isTRUE(record$fallback_used) && retry_attempts == 0L &&
        exact_attempts == 0L,
      retry_attempts = retry_attempts, exact_attempts = exact_attempts
    )
  }
  sinkhorn_ab <- evaluate_requested_sinkhorn(pair$a, pair$b, "AB")
  sinkhorn_ba <- evaluate_requested_sinkhorn(pair$b, pair$a, "BA")
  exact_ab <- evaluate_exact(pair$a, pair$b)
  exact_ba <- evaluate_exact(pair$b, pair$a)
  data.frame(
    spectrum_pair = pair_name,
    approximate_backend = "approxOT_sinkhorn",
    approximate_ot_estimand = "entropic_finite_iteration_raw_requested_only",
    solver_execution_policy = "raw_requested_only_no_retry_no_exact_fallback",
    sinkhorn_epsilon = 0.05,
    sinkhorn_niter = 100L,
    sinkhorn_d_ab = sinkhorn_ab$value,
    sinkhorn_d_ba = sinkhorn_ba$value,
    sinkhorn_absolute_direction_gap = abs(sinkhorn_ab$value - sinkhorn_ba$value),
    sinkhorn_selected_path_ab = sinkhorn_ab$record$selected_path,
    sinkhorn_selected_path_ba = sinkhorn_ba$record$selected_path,
    sinkhorn_solver_attempts_ab = nrow(sinkhorn_ab$record$attempts),
    sinkhorn_solver_attempts_ba = nrow(sinkhorn_ba$record$attempts),
    sinkhorn_requested_setting_match_ab = sinkhorn_ab$setting_match,
    sinkhorn_requested_setting_match_ba = sinkhorn_ba$setting_match,
    sinkhorn_retry_attempts_ab = sinkhorn_ab$retry_attempts,
    sinkhorn_retry_attempts_ba = sinkhorn_ba$retry_attempts,
    sinkhorn_exact_attempts_ab = sinkhorn_ab$exact_attempts,
    sinkhorn_exact_attempts_ba = sinkhorn_ba$exact_attempts,
    sinkhorn_validation_error_ab = sinkhorn_ab$attempt$validation_error,
    sinkhorn_validation_error_ba = sinkhorn_ba$attempt$validation_error,
    sinkhorn_row_residual_linf_ab = sinkhorn_ab$attempt$row_residual_linf,
    sinkhorn_col_residual_linf_ab = sinkhorn_ab$attempt$col_residual_linf,
    sinkhorn_row_residual_linf_ba = sinkhorn_ba$attempt$row_residual_linf,
    sinkhorn_col_residual_linf_ba = sinkhorn_ba$attempt$col_residual_linf,
    sinkhorn_nonfinite_plan_mass_ab = sinkhorn_ab$attempt$nonfinite_mass,
    sinkhorn_nonfinite_plan_mass_ba = sinkhorn_ba$attempt$nonfinite_mass,
    requested_only_provenance_gate_pass =
      sinkhorn_ab$provenance_gate_pass && sinkhorn_ba$provenance_gate_pass,
    exact_d_ab = exact_ab,
    exact_d_ba = exact_ba,
    exact_absolute_direction_gap = abs(exact_ab - exact_ba),
    finite_iteration_directionality_detected =
      is.finite(sinkhorn_ab$value) && is.finite(sinkhorn_ba$value) &&
      abs(sinkhorn_ab$value - sinkhorn_ba$value) > 1e-6,
    exact_symmetry_pass = abs(exact_ab - exact_ba) <= 1e-12,
    publication_primary_backend = "transport_exact_unregularized",
    stringsAsFactors = FALSE
  )
})
sinkhorn_diagnostic <- do.call(rbind, sinkhorn_rows)
utils::write.csv(
  sinkhorn_diagnostic, sinkhorn_diagnostic_file, row.names = FALSE
)
print(sinkhorn_diagnostic)
if (!all(sinkhorn_diagnostic$requested_only_provenance_gate_pass) ||
    !all(sinkhorn_diagnostic$finite_iteration_directionality_detected) ||
    !all(sinkhorn_diagnostic$exact_symmetry_pass)) {
  stop(
    "Sinkhorn directionality or exact-OT symmetry stress contract failed; see ",
    sinkhorn_diagnostic_file
  )
}
git_commit <- tryCatch(
  trimws(system2("git", c("-C", repo_dir, "rev-parse", "HEAD"),
                  stdout = TRUE, stderr = FALSE)),
  error = function(e) NA_character_
)
git_status <- tryCatch(
  system2("git", c("-C", repo_dir, "status", "--porcelain"),
          stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
utils::write.csv(
  data.frame(
    parameter = c("timestamp_utc", "package_commit", "package_tree_dirty",
                  "command", "symmetry_tolerance", "methods",
                  "publication_primary_ot_backend",
                  "approximate_backend_diagnostic"),
    value = c(format(Sys.time(), tz = "UTC"), git_commit,
              length(git_status) > 0L, paste(commandArgs(), collapse = " "),
              tolerance, paste(methods, collapse = ","),
              "transport_exact_unregularized",
              paste(
                "raw approxOT Sinkhorn epsilon=0.05 niter=100; exactly one",
                "requested-setting call per orientation; no retry/exact fallback"
              )),
    stringsAsFactors = FALSE
  ),
  file.path(dirname(output_file), "method_symmetry_check_parameters.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()),
           file.path(dirname(output_file), "method_symmetry_check_sessionInfo.txt"))
message("Wrote method symmetry audit: ", output_file)
