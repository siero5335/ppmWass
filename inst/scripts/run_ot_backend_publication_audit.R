#!/usr/bin/env Rscript

# Publication audit for the exact-primary ppm-Wasserstein backend and the
# finite-iteration Sinkhorn backend retained for diagnostics.  The main
# publication matrix must be finite, symmetric, zero-diagonal, and traceable to
# an exact unregularized OT run.  In contrast, nonfinite values and direction
# gaps from raw approxOT Sinkhorn (epsilon 0.05, 100 iterations) are measured
# and reported but deliberately do not fail the audit.

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
  get_arg("output-dir", file.path(dirname(main_dir), "audit", "ot_backend")),
  mustWork = FALSE
)
datasets <- split_arg("datasets", "RECETOX,HREI-MSDB")
file_tag <- get_arg("file-tag", "tol15")
n_cores_requested <- as.integer(get_arg("n-cores", "8"))
seed <- as.integer(get_arg("seed", "20260718"))
exact_sample_size <- as.integer(get_arg("exact-sample-size", "1000"))
ppm <- as.numeric(get_arg("base-cost-ppm", "15"))
transition_mult <- as.numeric(get_arg("transition-mult", "3"))
sinkhorn_epsilon <- as.numeric(get_arg("sinkhorn-epsilon", "0.05"))
sinkhorn_niter <- as.integer(get_arg("sinkhorn-niter", "100"))
w_frag <- as.numeric(get_arg("w-frag", "0.70"))
w_loss <- as.numeric(get_arg("w-loss", "0.30"))
exact_tolerance <- as.numeric(get_arg("exact-tolerance", "1e-12"))
main_tolerance <- as.numeric(get_arg("main-tolerance", "1e-12"))
top_n <- as.integer(get_arg("top-n", "50"))
# Zero means the complete unordered-pair population.  A positive limit exists
# only for a labelled development smoke test; publication runs must leave it 0.
max_pairs_per_dataset <- as.integer(
  get_arg("max-pairs-per-dataset", "0")
)
parameter_file_arg <- get_arg("parameter-file", "")

if (!length(datasets) || any(!datasets %in% c("RECETOX", "HREI-MSDB"))) {
  stop("--datasets must contain RECETOX and/or HREI-MSDB.", call. = FALSE)
}
if (!grepl("^[A-Za-z0-9._-]+$", file_tag)) {
  stop("--file-tag contains unsupported characters.", call. = FALSE)
}
positive_values <- c(
  n_cores = n_cores_requested,
  exact_sample_size = exact_sample_size,
  ppm = ppm,
  transition_mult = transition_mult,
  sinkhorn_epsilon = sinkhorn_epsilon,
  sinkhorn_niter = sinkhorn_niter,
  exact_tolerance = exact_tolerance,
  main_tolerance = main_tolerance,
  top_n = top_n
)
if (any(!is.finite(positive_values)) || any(positive_values <= 0)) {
  stop("Core counts, sample sizes, solver settings, and tolerances must be positive.",
       call. = FALSE)
}
if (!is.finite(seed) || !is.finite(max_pairs_per_dataset) ||
    max_pairs_per_dataset < 0L || any(!is.finite(c(w_frag, w_loss))) ||
    w_frag < 0 || w_loss < 0 || w_frag + w_loss <= 0) {
  stop("Invalid seed, pair limit, or channel weights.", call. = FALSE)
}
if (!identical(sinkhorn_niter, 100L)) {
  stop("This publication diagnostic is defined for --sinkhorn-niter=100.",
       call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
if (!requireNamespace("transport", quietly = TRUE)) stop("transport is required")
if (!requireNamespace("approxOT", quietly = TRUE)) stop("approxOT is required")
pkgload::load_all(repo_dir, quiet = TRUE)

compute_one <- getFromNamespace("compute_distance", "ppmWass")
assess_plan <- getFromNamespace("assess_ot_plan", "ppmWass")

dataset_dir_name <- function(dataset) gsub("[^A-Za-z0-9]+", "_", dataset)
tagged_name <- function(stem, extension) paste0(stem, "_", file_tag, extension)

git_one_line <- function(arguments) {
  value <- tryCatch(
    suppressWarnings(system2("git", c("-C", repo_dir, arguments),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (!length(value)) NA_character_ else trimws(value[[1L]])
}

git_commit <- git_one_line(c("rev-parse", "HEAD"))
git_status <- tryCatch(
  system2("git", c("-C", repo_dir, "status", "--porcelain"),
          stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)

resolve_parameter_file <- function() {
  if (nzchar(parameter_file_arg)) {
    return(normalizePath(parameter_file_arg, mustWork = TRUE))
  }
  preferred <- file.path(main_dir, tagged_name("main_run_parameters", ".csv"))
  if (file.exists(preferred)) return(normalizePath(preferred, mustWork = TRUE))
  legacy <- file.path(main_dir, "main_run_parameters.csv")
  if (file.exists(legacy)) return(normalizePath(legacy, mustWork = TRUE))
  candidates <- Sys.glob(file.path(main_dir, "main_run_parameters*.csv"))
  if (length(candidates) != 1L) {
    stop(
      "Could not resolve one main-run parameter file under ", main_dir,
      "; supply --parameter-file=...",
      call. = FALSE
    )
  }
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

parameter_file <- resolve_parameter_file()
main_parameters <- utils::read.csv(
  parameter_file, stringsAsFactors = FALSE, check.names = FALSE
)
if (!all(c("parameter", "value") %in% names(main_parameters))) {
  stop("Main-run parameter file must contain parameter and value columns.",
       call. = FALSE)
}

parameter_value <- function(name) {
  values <- unique(as.character(
    main_parameters$value[main_parameters$parameter == name]
  ))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 1L) values[[1L]] else NA_character_
}

is_false_text <- function(x) {
  length(x) == 1L && !is.na(x) && tolower(x) %in% c("false", "0", "no")
}

is_true_text <- function(x) {
  length(x) == 1L && !is.na(x) && tolower(x) %in% c("true", "1", "yes")
}

provenance_values <- list(
  ot_method = parameter_value("ot_method"),
  ot_estimand = parameter_value("ot_estimand"),
  solver_backend = parameter_value("solver_backend"),
  transport_available = parameter_value("transport_available"),
  package_commit = parameter_value("package_commit"),
  package_tree_dirty = parameter_value("package_tree_dirty"),
  distance_generation_commit = parameter_value("distance_generation_commit")
)

provenance_exact_pass <- identical(tolower(provenance_values$ot_method), "exact") &&
  identical(
    provenance_values$ot_estimand,
    "unregularized_exact_transport_cost"
  ) &&
  !is.na(provenance_values$solver_backend) &&
  grepl("transport", provenance_values$solver_backend, ignore.case = TRUE) &&
  grepl("exact", provenance_values$solver_backend, ignore.case = TRUE) &&
  is_true_text(provenance_values$transport_available) &&
  !is.na(provenance_values$package_commit) &&
  grepl("^[0-9a-fA-F]{7,40}$", provenance_values$package_commit) &&
  is_false_text(provenance_values$package_tree_dirty) &&
  !is.na(provenance_values$distance_generation_commit) &&
  grepl("^[0-9a-fA-F]{7,40}$", provenance_values$distance_generation_commit)

loaded <- list()
health_rows <- vector("list", length(datasets))

for (di in seq_along(datasets)) {
  dataset <- datasets[[di]]
  dataset_dir <- file.path(main_dir, dataset_dir_name(dataset))
  spectra_file <- file.path(
    dataset_dir, tagged_name("spectra_filtered", ".rds")
  )
  matrix_file <- file.path(
    dataset_dir, tagged_name("dist_ppm_wasserstein", ".rds")
  )
  spectra_exists <- file.exists(spectra_file)
  matrix_exists <- file.exists(matrix_file)
  spectra <- if (spectra_exists) readRDS(spectra_file) else NULL
  distance_matrix <- if (matrix_exists) readRDS(matrix_file) else NULL

  spectra_structure_pass <- is.list(spectra) &&
    !is.null(spectra$df_spec) && "id" %in% names(spectra$df_spec) &&
    is.list(spectra$frag_list) && is.list(spectra$loss_list) &&
    length(spectra$frag_list) == nrow(spectra$df_spec) &&
    length(spectra$loss_list) == nrow(spectra$df_spec)
  expected_ids <- if (spectra_structure_pass) {
    as.character(spectra$df_spec$id)
  } else character()
  matrix_structure_pass <- is.matrix(distance_matrix) &&
    is.numeric(distance_matrix) && length(expected_ids) > 0L &&
    identical(dim(distance_matrix), c(length(expected_ids), length(expected_ids))) &&
    identical(rownames(distance_matrix), expected_ids) &&
    identical(colnames(distance_matrix), expected_ids)
  n_nonfinite <- if (is.matrix(distance_matrix) && is.numeric(distance_matrix)) {
    sum(!is.finite(distance_matrix))
  } else NA_integer_
  finite_pass <- isTRUE(matrix_structure_pass) && identical(n_nonfinite, 0L)
  symmetry_max_gap <- if (finite_pass) {
    max(abs(distance_matrix - t(distance_matrix)))
  } else NA_real_
  diagonal_max_abs <- if (finite_pass) max(abs(diag(distance_matrix))) else NA_real_
  range_min <- if (finite_pass) min(distance_matrix) else NA_real_
  range_max <- if (finite_pass) max(distance_matrix) else NA_real_
  symmetric_pass <- finite_pass && is.finite(symmetry_max_gap) &&
    symmetry_max_gap <= exact_tolerance
  diagonal_zero_pass <- finite_pass && is.finite(diagonal_max_abs) &&
    diagonal_max_abs <= exact_tolerance
  range_pass <- finite_pass && range_min >= -exact_tolerance &&
    range_max <= 1 + exact_tolerance
  health_pass <- spectra_exists && matrix_exists && spectra_structure_pass &&
    matrix_structure_pass && finite_pass && symmetric_pass &&
    diagonal_zero_pass && range_pass && provenance_exact_pass

  health_rows[[di]] <- data.frame(
    dataset = dataset,
    spectra_file = normalizePath(spectra_file, mustWork = FALSE),
    spectra_md5 = if (spectra_exists) unname(tools::md5sum(spectra_file)) else NA_character_,
    matrix_file = normalizePath(matrix_file, mustWork = FALSE),
    matrix_md5 = if (matrix_exists) unname(tools::md5sum(matrix_file)) else NA_character_,
    parameter_file = parameter_file,
    parameter_file_md5 = unname(tools::md5sum(parameter_file)),
    n_spectra = length(expected_ids),
    matrix_rows = if (is.matrix(distance_matrix)) nrow(distance_matrix) else NA_integer_,
    matrix_cols = if (is.matrix(distance_matrix)) ncol(distance_matrix) else NA_integer_,
    n_nonfinite = n_nonfinite,
    symmetry_max_absolute_gap = symmetry_max_gap,
    diagonal_max_absolute_value = diagonal_max_abs,
    matrix_min = range_min,
    matrix_max = range_max,
    spectra_structure_pass = spectra_structure_pass,
    matrix_structure_and_id_pass = matrix_structure_pass,
    finite_pass = finite_pass,
    symmetric_pass = symmetric_pass,
    diagonal_zero_pass = diagonal_zero_pass,
    range_0_1_pass = range_pass,
    provenance_exact_primary_pass = provenance_exact_pass,
    provenance_ot_method = provenance_values$ot_method,
    provenance_ot_estimand = provenance_values$ot_estimand,
    provenance_solver_backend = provenance_values$solver_backend,
    provenance_transport_available = provenance_values$transport_available,
    provenance_package_commit = provenance_values$package_commit,
    provenance_package_tree_dirty = provenance_values$package_tree_dirty,
    provenance_distance_generation_commit =
      provenance_values$distance_generation_commit,
    health_pass = health_pass,
    stringsAsFactors = FALSE
  )
  if (spectra_structure_pass && matrix_structure_pass) {
    names(spectra$frag_list) <- expected_ids
    names(spectra$loss_list) <- expected_ids
    loaded[[dataset]] <- list(
      spectra = spectra, matrix = distance_matrix,
      spectra_file = spectra_file, matrix_file = matrix_file
    )
  }
}

health <- do.call(rbind, health_rows)
health_file <- file.path(output_dir, "ot_exact_primary_main_health.csv")
utils::write.csv(health, health_file, row.names = FALSE)

available_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (!is.finite(available_cores) || available_cores < 1L) available_cores <- 1L
n_cores_effective <- max(1L, min(n_cores_requested, available_cores))

exact_command <- paste(
  c(file.path(R.home("bin"), "Rscript"), commandArgs()), collapse = " "
)
scan_scope <- if (max_pairs_per_dataset == 0L) {
  "complete_all_unordered_pairs"
} else {
  paste0("development_smoke_first_", max_pairs_per_dataset, "_unordered_pairs")
}
parameter_output <- data.frame(
  parameter = c(
    "timestamp_utc", "command", "repository", "package_commit_at_audit",
    "package_tree_dirty_at_audit", "main_dir", "main_parameter_file",
    "output_dir", "datasets", "file_tag", "scan_scope",
    "max_pairs_per_dataset", "seed", "exact_sample_size_per_dataset",
    "n_cores_requested", "n_cores_effective", "base_cost_ppm",
    "transition_multiplier", "cost_saturation_ppm", "publication_ot_method",
    "publication_ot_estimand", "approximate_diagnostic_backend",
    "sinkhorn_epsilon", "sinkhorn_iterations", "fragment_weight",
    "pooled_derived_weight", "exact_symmetry_tolerance",
    "main_matrix_reproduction_tolerance", "top_extremes_per_dataset_channel",
    "exact_primary_failure_policy", "approximate_diagnostic_failure_policy"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE), exact_command, repo_dir,
    git_commit, length(git_status) > 0L, main_dir, parameter_file, output_dir,
    paste(datasets, collapse = ","), file_tag, scan_scope,
    max_pairs_per_dataset, seed, exact_sample_size, n_cores_requested,
    n_cores_effective, ppm, transition_mult, ppm * transition_mult,
    "transport_exact_unregularized", "unregularized_exact_transport_cost",
    "raw_approxOT_sinkhorn_AB_and_BA", sinkhorn_epsilon, sinkhorn_niter,
    w_frag, w_loss, exact_tolerance, main_tolerance, top_n,
    paste(
      "abort if main health fails or any exact sample is nonfinite,",
      "directionally inconsistent, or differs from the main matrix"
    ),
    "record direction gaps/nonfinite values; never fail on approximate behavior"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  parameter_output,
  file.path(output_dir, "ot_backend_publication_audit_parameters.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "ot_backend_publication_audit_sessionInfo.txt"),
  useBytes = TRUE
)

if (!all(health$health_pass)) {
  stop(
    "Exact-primary main matrix/provenance health failed; see ", health_file,
    call. = FALSE
  )
}

coerce_channel <- function(x) {
  if (is.null(x)) {
    return(list(status = "empty_null", matrix = matrix(numeric(), 0L, 2L),
                empty = TRUE, valid = TRUE))
  }
  converted <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (is.null(converted) || length(dim(converted)) != 2L || ncol(converted) < 2L) {
    return(list(status = "malformed", matrix = matrix(numeric(), 0L, 2L),
                empty = FALSE, valid = FALSE))
  }
  if (!nrow(converted)) {
    return(list(status = "empty_zero_rows", matrix = matrix(numeric(), 0L, 2L),
                empty = TRUE, valid = TRUE))
  }
  values <- suppressWarnings(as.numeric(converted[, 1:2, drop = FALSE]))
  matrix_value <- matrix(values, nrow = nrow(converted), ncol = 2L)
  list(status = "present", matrix = matrix_value, empty = FALSE, valid = TRUE)
}

prepare_channel_pair <- function(a, b) {
  aa <- coerce_channel(a)
  bb <- coerce_channel(b)
  base <- list(
    status = "invalid_input", solver_ready = FALSE,
    empty_a = aa$empty, empty_b = bb$empty,
    n_peaks_a = nrow(aa$matrix), n_peaks_b = nrow(bb$matrix),
    constant_cost = NA, constant_cost_value = NA_real_,
    wa = numeric(), wb = numeric(), cost = matrix(numeric(), 0L, 0L)
  )
  # This mirrors the package's early empty-channel return while preserving the
  # distinction from malformed/nonpositive inputs in the audit output.
  if (aa$empty || bb$empty) {
    base$status <- "empty_channel"
    return(base)
  }
  if (!aa$valid || !bb$valid) return(base)
  mz_a <- aa$matrix[, 1L]
  mz_b <- bb$matrix[, 1L]
  wa <- aa$matrix[, 2L]
  wb <- bb$matrix[, 2L]
  total_a <- sum(wa)
  total_b <- sum(wb)
  valid <- all(is.finite(c(mz_a, mz_b, wa, wb, total_a, total_b))) &&
    all(mz_a > 0) && all(mz_b > 0) && all(wa >= 0) && all(wb >= 0) &&
    total_a > 0 && total_b > 0
  if (!valid) return(base)
  wa <- wa / total_a
  wb <- wb / total_b
  mean_mz <- outer(mz_a, mz_b, function(x, y) (x + y) / 2)
  delta_ppm <- abs(outer(mz_a, mz_b, "-")) / mean_mz * 1e6
  cost <- pmin(delta_ppm / (ppm * transition_mult), 1)
  if (!length(cost) || any(!is.finite(cost))) return(base)
  constant_cost <- all(cost == cost[[1L]])
  base$status <- if (constant_cost) {
    "valid_constant_cost"
  } else {
    "valid_nonconstant_cost"
  }
  base$solver_ready <- TRUE
  base$constant_cost <- constant_cost
  base$constant_cost_value <- if (constant_cost) cost[[1L]] else NA_real_
  base$wa <- wa
  base$wb <- wb
  base$cost <- cost
  base
}

empty_solver_result <- function(error = NA_character_) {
  list(
    total = NA_real_, accepted = FALSE, validation_error = "solver_not_run",
    solver_error = error, row_residual_linf = NA_real_,
    col_residual_linf = NA_real_, nonfinite_mass = NA_integer_
  )
}

run_raw_sinkhorn <- function(wa, wb, cost) {
  solver_error <- NA_character_
  plan <- tryCatch(
    approxOT::transport_plan_given_C(
      mass_x = wa, mass_y = wb, p = 1, cost = cost,
      method = "sinkhorn", epsilon = sinkhorn_epsilon,
      niter = sinkhorn_niter
    ),
    error = function(e) {
      solver_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(plan)) return(empty_solver_result(solver_error))
  assessed <- assess_plan(plan, cost, wa, wb, tolerance = 1e-8)
  list(
    total = assessed$total_cost,
    accepted = isTRUE(assessed$accepted),
    validation_error = assessed$validation_error,
    solver_error = solver_error,
    row_residual_linf = assessed$row_residual_linf,
    col_residual_linf = assessed$col_residual_linf,
    nonfinite_mass = assessed$nonfinite_mass
  )
}

scan_pair_channel <- function(dataset, pair_index, i, j, ids, channel,
                              spectrum_a, spectrum_b) {
  prepared <- prepare_channel_pair(spectrum_a, spectrum_b)
  if (prepared$solver_ready) {
    ab <- run_raw_sinkhorn(prepared$wa, prepared$wb, prepared$cost)
    ba <- run_raw_sinkhorn(prepared$wb, prepared$wa, t(prepared$cost))
  } else {
    ab <- empty_solver_result()
    ba <- empty_solver_result()
  }
  finite_ab <- is.finite(ab$total)
  finite_ba <- is.finite(ba$total)
  data.frame(
    dataset = dataset,
    pair_index = pair_index,
    spectrum_index_a = i,
    spectrum_index_b = j,
    spectrum_id_a = ids[[i]],
    spectrum_id_b = ids[[j]],
    channel = channel,
    n_peaks_a = prepared$n_peaks_a,
    n_peaks_b = prepared$n_peaks_b,
    input_status = prepared$status,
    empty_a = prepared$empty_a,
    empty_b = prepared$empty_b,
    empty_either = prepared$empty_a || prepared$empty_b,
    invalid_input = identical(prepared$status, "invalid_input"),
    constant_cost = prepared$constant_cost,
    constant_cost_value = prepared$constant_cost_value,
    solver_called = prepared$solver_ready,
    sinkhorn_ab = ab$total,
    sinkhorn_ba = ba$total,
    finite_ab = finite_ab,
    finite_ba = finite_ba,
    any_nonfinite = prepared$solver_ready && (!finite_ab || !finite_ba),
    absolute_direction_gap = if (finite_ab && finite_ba) {
      abs(ab$total - ba$total)
    } else NA_real_,
    ab_plan_accepted = ab$accepted,
    ba_plan_accepted = ba$accepted,
    ab_validation_error = ab$validation_error,
    ba_validation_error = ba$validation_error,
    ab_solver_error = ab$solver_error,
    ba_solver_error = ba$solver_error,
    ab_row_residual_linf = ab$row_residual_linf,
    ab_col_residual_linf = ab$col_residual_linf,
    ba_row_residual_linf = ba$row_residual_linf,
    ba_col_residual_linf = ba$col_residual_linf,
    ab_nonfinite_plan_mass = ab$nonfinite_mass,
    ba_nonfinite_plan_mass = ba$nonfinite_mass,
    stringsAsFactors = FALSE
  )
}

pair_rank <- function(i, j, n) {
  as.integer((i - 1L) * (2L * n - i) / 2 + (j - i))
}

parallel_map <- function(x, fun, cores) {
  if (cores <= 1L || length(x) <= 1L) return(lapply(x, fun))
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(
      x, fun, mc.cores = min(cores, length(x)), mc.preschedule = TRUE
    ))
  }
  # PSOCK workers receive the closure (including the selected spectra) and load
  # the two solver namespaces explicitly.  This branch is not used on the
  # publication macOS host but keeps --n-cores meaningful on Windows.
  cluster <- parallel::makeCluster(min(cores, length(x)))
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterEvalQ(cluster, {
    library(approxOT)
    library(transport)
    NULL
  })
  parallel::parLapply(cluster, x, fun)
}

summary_columns <- c(
  "dataset", "channel", "solver_called", "empty_either", "invalid_input",
  "constant_cost", "finite_ab", "finite_ba", "any_nonfinite",
  "ab_plan_accepted", "ba_plan_accepted", "absolute_direction_gap"
)

all_pair_file <- file.path(
  output_dir, "ot_sinkhorn100_all_unordered_pairs.csv.gz"
)
all_pair_temp <- paste0(all_pair_file, ".tmp-", Sys.getpid())
if (file.exists(all_pair_temp)) unlink(all_pair_temp)
on.exit(if (file.exists(all_pair_temp)) unlink(all_pair_temp), add = TRUE)
all_pair_connection <- gzfile(all_pair_temp, open = "wt", compression = 6)
header_written <- FALSE
summary_parts <- list()
nonfinite_parts <- list()
top_parts <- list()
part_index <- 0L
nonfinite_index <- 0L
top_index <- 0L

for (dataset in datasets) {
  spectra <- loaded[[dataset]]$spectra
  ids <- as.character(spectra$df_spec$id)
  n <- length(ids)
  i_values <- if (n >= 2L) seq_len(n - 1L) else integer()
  total_pairs <- n * (n - 1L) / 2
  if (max_pairs_per_dataset > 0L && total_pairs > max_pairs_per_dataset) {
    keep_pair <- function(i, j) pair_rank(i, j, n) <= max_pairs_per_dataset
  } else {
    keep_pair <- function(i, j) TRUE
  }

  scan_i <- function(i) {
    rows <- vector("list", 2L * (n - i))
    row_index <- 0L
    for (j in seq.int(i + 1L, n)) {
      if (!keep_pair(i, j)) next
      index <- pair_rank(i, j, n)
      row_index <- row_index + 1L
      rows[[row_index]] <- scan_pair_channel(
        dataset, index, i, j, ids, "fragment",
        spectra$frag_list[[i]], spectra$frag_list[[j]]
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- scan_pair_channel(
        dataset, index, i, j, ids, "pooled_derived",
        spectra$loss_list[[i]], spectra$loss_list[[j]]
      )
    }
    if (!row_index) return(NULL)
    do.call(rbind, rows[seq_len(row_index)])
  }

  block_size <- max(8L, n_cores_effective * 2L)
  blocks <- split(i_values, ceiling(seq_along(i_values) / block_size))
  for (block_number in seq_along(blocks)) {
    block_result <- parallel_map(blocks[[block_number]], scan_i, n_cores_effective)
    block_result <- block_result[!vapply(block_result, is.null, logical(1))]
    if (!length(block_result)) next
    chunk <- do.call(rbind, block_result)
    utils::write.table(
      chunk, all_pair_connection, sep = ",", row.names = FALSE,
      col.names = !header_written, append = header_written, quote = TRUE,
      na = "", qmethod = "double"
    )
    header_written <- TRUE

    part_index <- part_index + 1L
    summary_parts[[part_index]] <- chunk[, summary_columns, drop = FALSE]
    nonfinite_rows <- chunk[
      chunk$any_nonfinite |
        (!is.na(chunk$ab_nonfinite_plan_mass) & chunk$ab_nonfinite_plan_mass > 0L) |
        (!is.na(chunk$ba_nonfinite_plan_mass) & chunk$ba_nonfinite_plan_mass > 0L) |
        (!is.na(chunk$ab_solver_error) & nzchar(chunk$ab_solver_error)) |
        (!is.na(chunk$ba_solver_error) & nzchar(chunk$ba_solver_error)),
      , drop = FALSE
    ]
    if (nrow(nonfinite_rows)) {
      nonfinite_index <- nonfinite_index + 1L
      nonfinite_parts[[nonfinite_index]] <- nonfinite_rows
    }
    finite_gap_rows <- chunk[is.finite(chunk$absolute_direction_gap), , drop = FALSE]
    if (nrow(finite_gap_rows)) {
      split_top <- split(finite_gap_rows, finite_gap_rows$channel)
      for (candidate in split_top) {
        candidate <- candidate[
          order(candidate$absolute_direction_gap, decreasing = TRUE),
          , drop = FALSE
        ]
        top_index <- top_index + 1L
        top_parts[[top_index]] <- head(candidate, top_n)
      }
    }
    message(
      "[", dataset, "] Sinkhorn scan block ", block_number, "/",
      length(blocks), ": ", nrow(chunk), " channel-pair rows"
    )
  }
}
close(all_pair_connection)
if (!header_written) {
  stop("No unordered spectrum pairs were available for the Sinkhorn scan.",
       call. = FALSE)
}
if (file.exists(all_pair_file) && !file.remove(all_pair_file)) {
  stop("Could not replace prior all-pairs output: ", all_pair_file,
       call. = FALSE)
}
if (!file.rename(all_pair_temp, all_pair_file)) {
  stop("Could not finalize compressed all-pairs output.", call. = FALSE)
}

summary_data <- do.call(rbind, summary_parts)

summarize_group <- function(x, dataset_label, channel_label) {
  gaps <- x$absolute_direction_gap[is.finite(x$absolute_direction_gap)]
  quantile_or_na <- function(prob) {
    if (!length(gaps)) NA_real_ else
      unname(stats::quantile(gaps, prob, names = FALSE, type = 8))
  }
  data.frame(
    dataset = dataset_label,
    channel = channel_label,
    n_channel_pair_rows = nrow(x),
    n_solver_called = sum(x$solver_called),
    n_empty_pairs = sum(x$empty_either),
    n_invalid_input_pairs = sum(x$invalid_input),
    n_constant_cost_pairs = sum(x$constant_cost %in% TRUE),
    n_finite_ab_and_ba = sum(x$finite_ab & x$finite_ba),
    n_nonfinite_either_direction = sum(x$any_nonfinite),
    n_plan_accepted_both = sum(x$ab_plan_accepted & x$ba_plan_accepted),
    n_finite_direction_gaps = length(gaps),
    mean_absolute_direction_gap = if (length(gaps)) mean(gaps) else NA_real_,
    median_absolute_direction_gap = if (length(gaps)) median(gaps) else NA_real_,
    p95_absolute_direction_gap = quantile_or_na(0.95),
    p99_absolute_direction_gap = quantile_or_na(0.99),
    max_absolute_direction_gap = if (length(gaps)) max(gaps) else NA_real_,
    n_gap_gt_1e_12 = sum(gaps > 1e-12),
    n_gap_gt_1e_10 = sum(gaps > 1e-10),
    n_gap_gt_1e_8 = sum(gaps > 1e-8),
    n_gap_gt_1e_6 = sum(gaps > 1e-6),
    n_gap_gt_1e_4 = sum(gaps > 1e-4),
    diagnostic_only = TRUE,
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
summary_index <- 0L
for (dataset in datasets) {
  dataset_rows <- summary_data[summary_data$dataset == dataset, , drop = FALSE]
  for (channel in c("fragment", "pooled_derived")) {
    summary_index <- summary_index + 1L
    summary_rows[[summary_index]] <- summarize_group(
      dataset_rows[dataset_rows$channel == channel, , drop = FALSE],
      dataset, channel
    )
  }
  summary_index <- summary_index + 1L
  summary_rows[[summary_index]] <- summarize_group(
    dataset_rows, dataset, "all_channels"
  )
}
summary_index <- summary_index + 1L
summary_rows[[summary_index]] <- summarize_group(
  summary_data, "all_datasets", "all_channels"
)
sinkhorn_summary <- do.call(rbind, summary_rows)
utils::write.csv(
  sinkhorn_summary,
  file.path(output_dir, "ot_sinkhorn100_directionality_summary.csv"),
  row.names = FALSE
)

if (length(top_parts)) {
  top_candidates <- do.call(rbind, top_parts)
  top_split <- split(
    top_candidates,
    interaction(top_candidates$dataset, top_candidates$channel, drop = TRUE)
  )
  top_extremes <- do.call(rbind, lapply(top_split, function(x) {
    head(x[order(x$absolute_direction_gap, decreasing = TRUE), , drop = FALSE],
         top_n)
  }))
  rownames(top_extremes) <- NULL
} else {
  top_extremes <- utils::read.csv(
    gzfile(all_pair_file), nrows = 0L, stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
utils::write.csv(
  top_extremes,
  file.path(output_dir, "ot_sinkhorn100_directionality_top_extremes.csv"),
  row.names = FALSE
)

if (length(nonfinite_parts)) {
  nonfinite_rows <- do.call(rbind, nonfinite_parts)
} else {
  nonfinite_rows <- top_extremes[0, , drop = FALSE]
}
utils::write.csv(
  nonfinite_rows,
  file.path(output_dir, "ot_sinkhorn100_nonfinite.csv"),
  row.names = FALSE
)

safe_exact <- function(a, b) {
  solver_error <- NA_character_
  value <- tryCatch(
    compute_one(
      a, b, method = "ppm_wasserstein", ppm = ppm,
      align_wasserstein = FALSE,
      wasserstein_transition_mult = transition_mult,
      ot_method = "exact", sinkhorn_epsilon = sinkhorn_epsilon,
      sinkhorn_niter = sinkhorn_niter
    ),
    error = function(e) {
      solver_error <<- conditionMessage(e)
      NA_real_
    }
  )
  list(value = value, error = solver_error)
}

exact_sample_parts <- list()
for (di in seq_along(datasets)) {
  dataset <- datasets[[di]]
  spectra <- loaded[[dataset]]$spectra
  main_matrix <- loaded[[dataset]]$matrix
  ids <- as.character(spectra$df_spec$id)
  n <- length(ids)
  all_pairs <- utils::combn(n, 2L)
  n_population <- ncol(all_pairs)
  n_sample <- min(exact_sample_size, n_population)
  set.seed(seed + di - 1L)
  selected_ranks <- sort(sample.int(n_population, n_sample, replace = FALSE))

  exact_one <- function(sample_position) {
    rank <- selected_ranks[[sample_position]]
    i <- all_pairs[1L, rank]
    j <- all_pairs[2L, rank]
    frag_prepared <- prepare_channel_pair(
      spectra$frag_list[[i]], spectra$frag_list[[j]]
    )
    loss_prepared <- prepare_channel_pair(
      spectra$loss_list[[i]], spectra$loss_list[[j]]
    )
    frag_ab <- safe_exact(spectra$frag_list[[i]], spectra$frag_list[[j]])
    frag_ba <- safe_exact(spectra$frag_list[[j]], spectra$frag_list[[i]])
    loss_ab <- safe_exact(spectra$loss_list[[i]], spectra$loss_list[[j]])
    loss_ba <- safe_exact(spectra$loss_list[[j]], spectra$loss_list[[i]])
    combined_ab <- sqrt(
      w_frag * frag_ab$value^2 + w_loss * loss_ab$value^2
    )
    combined_ba <- sqrt(
      w_frag * frag_ba$value^2 + w_loss * loss_ba$value^2
    )
    component_values <- c(
      frag_ab$value, frag_ba$value, loss_ab$value, loss_ba$value,
      combined_ab, combined_ba
    )
    all_finite <- all(is.finite(component_values))
    fragment_gap <- abs(frag_ab$value - frag_ba$value)
    loss_gap <- abs(loss_ab$value - loss_ba$value)
    combined_gap <- abs(combined_ab - combined_ba)
    main_ab <- main_matrix[i, j]
    main_ba <- main_matrix[j, i]
    main_diff_ab <- abs(combined_ab - main_ab)
    main_diff_ba <- abs(combined_ba - main_ba)
    exact_symmetry_pass <- all_finite &&
      fragment_gap <= exact_tolerance && loss_gap <= exact_tolerance &&
      combined_gap <= exact_tolerance
    main_match_pass <- all_finite && is.finite(main_ab) && is.finite(main_ba) &&
      main_diff_ab <= main_tolerance && main_diff_ba <= main_tolerance
    data.frame(
      dataset = dataset,
      sample_position = sample_position,
      pair_index = rank,
      spectrum_index_a = i,
      spectrum_index_b = j,
      spectrum_id_a = ids[[i]],
      spectrum_id_b = ids[[j]],
      fragment_input_status = frag_prepared$status,
      fragment_empty_either = frag_prepared$empty_a || frag_prepared$empty_b,
      fragment_constant_cost = frag_prepared$constant_cost,
      pooled_derived_input_status = loss_prepared$status,
      pooled_derived_empty_either = loss_prepared$empty_a || loss_prepared$empty_b,
      pooled_derived_constant_cost = loss_prepared$constant_cost,
      exact_fragment_ab = frag_ab$value,
      exact_fragment_ba = frag_ba$value,
      exact_fragment_absolute_direction_gap = fragment_gap,
      exact_pooled_derived_ab = loss_ab$value,
      exact_pooled_derived_ba = loss_ba$value,
      exact_pooled_derived_absolute_direction_gap = loss_gap,
      combined_exact_ab = combined_ab,
      combined_exact_ba = combined_ba,
      combined_exact_absolute_direction_gap = combined_gap,
      main_matrix_ab = main_ab,
      main_matrix_ba = main_ba,
      combined_exact_vs_main_ab_absolute_difference = main_diff_ab,
      combined_exact_vs_main_ba_absolute_difference = main_diff_ba,
      fragment_ab_error = frag_ab$error,
      fragment_ba_error = frag_ba$error,
      pooled_derived_ab_error = loss_ab$error,
      pooled_derived_ba_error = loss_ba$error,
      all_exact_values_finite = all_finite,
      exact_symmetry_pass = exact_symmetry_pass,
      main_matrix_match_pass = main_match_pass,
      sample_pass = exact_symmetry_pass && main_match_pass,
      exact_tolerance = exact_tolerance,
      main_tolerance = main_tolerance,
      stringsAsFactors = FALSE
    )
  }

  exact_results <- parallel_map(seq_len(n_sample), exact_one, n_cores_effective)
  exact_sample_parts[[di]] <- do.call(rbind, exact_results)
  message("[", dataset, "] exact AB/BA sample comparisons: ", n_sample)
}

exact_sample <- do.call(rbind, exact_sample_parts)
exact_sample_file <- file.path(
  output_dir, "ot_exact_primary_sample_comparison.csv"
)
utils::write.csv(exact_sample, exact_sample_file, row.names = FALSE)

max_finite_or_na <- function(x) {
  finite <- x[is.finite(x)]
  if (length(finite)) max(finite) else NA_real_
}

exact_sample_summary <- do.call(rbind, lapply(
  split(exact_sample, exact_sample$dataset),
  function(x) data.frame(
    dataset = x$dataset[[1L]],
    n_sampled_unordered_pairs = nrow(x),
    n_all_exact_values_finite = sum(x$all_exact_values_finite),
    max_fragment_exact_direction_gap =
      max_finite_or_na(x$exact_fragment_absolute_direction_gap),
    max_pooled_derived_exact_direction_gap =
      max_finite_or_na(x$exact_pooled_derived_absolute_direction_gap),
    max_combined_exact_direction_gap =
      max_finite_or_na(x$combined_exact_absolute_direction_gap),
    max_combined_exact_vs_main_difference = max_finite_or_na(
      c(x$combined_exact_vs_main_ab_absolute_difference,
        x$combined_exact_vs_main_ba_absolute_difference)
    ),
    n_exact_symmetry_pass = sum(x$exact_symmetry_pass),
    n_main_matrix_match_pass = sum(x$main_matrix_match_pass),
    n_sample_pass = sum(x$sample_pass),
    all_sample_pass = all(x$sample_pass),
    stringsAsFactors = FALSE
  )
))
utils::write.csv(
  exact_sample_summary,
  file.path(output_dir, "ot_exact_primary_sample_summary.csv"),
  row.names = FALSE
)

# Approximate-backend problems above are expected diagnostic observations and
# do not enter this terminal condition.  Only the exact publication contract
# can fail the run.
if (!all(exact_sample$sample_pass)) {
  stop(
    "Exact-primary sample comparison failed; see ", exact_sample_file,
    call. = FALSE
  )
}

message("Exact-primary health and sample comparison passed.")
message("Wrote raw Sinkhorn diagnostic (nonfinite/direction gaps do not fail): ",
        all_pair_file)
