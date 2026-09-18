#!/usr/bin/env Rscript

# Publication perturbation benchmark for ppmWass.
#
# This driver deliberately keeps two distinct experiments:
#
# 1. score_level
#    Perturb the already-preprocessed fragment spectrum and the already-built
#    pooled derived mass-difference representation directly. This reproduces
#    the interpretation of the historical robustness analysis: sensitivity of
#    the scoring representations, not a complete acquisition pipeline.
#
# 2. regenerated_representation
#    Perturb the raw fragment peak string, then rerun peak filtering/centroiding,
#    high-mass-reference estimation, anchored differences, pairwise differences,
#    and their pooled representation before scoring.
#
# Queries, condition seeds, and fragment-query seeds are paired across the two
# experiments. The original, unperturbed RECETOX spectra remain the library, so
# the endpoint is a query-only self-retention stress test rather than unknown-
# sample identification.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]])
}

has_arg <- function(name) {
  any(grepl(paste0("^--", name, "="), args))
}

flag_is_true <- function(name, default = FALSE) {
  value <- get_arg(name, if (default) "true" else "false")
  tolower(value) %in% c("1", "true", "yes", "y")
}

split_character <- function(name, default) {
  value <- trimws(strsplit(get_arg(name, default), ",", fixed = TRUE)[[1]])
  value[nzchar(value)]
}

split_numeric <- function(name, default) {
  value <- suppressWarnings(as.numeric(split_character(name, default)))
  if (any(!is.finite(value))) stop("--", name, " must be a comma-separated numeric vector.")
  value
}

as_single_integer <- function(name, default, lower = NULL) {
  value <- suppressWarnings(as.integer(get_arg(name, default)))
  if (length(value) != 1L || is.na(value)) stop("--", name, " must be an integer.")
  if (!is.null(lower) && value < lower) stop("--", name, " must be >= ", lower, ".")
  value
}

as_single_numeric <- function(name, default, lower = NULL,
                              lower_inclusive = TRUE) {
  value <- suppressWarnings(as.numeric(get_arg(name, default)))
  if (length(value) != 1L || !is.finite(value)) stop("--", name, " must be numeric.")
  if (!is.null(lower)) {
    ok <- if (lower_inclusive) value >= lower else value > lower
    if (!ok) stop("--", name, " is outside its allowed range.")
  }
  value
}

safe_normalize <- function(path, must_work = TRUE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

repo_dir <- safe_normalize(get_arg("repo", getwd()))
bundle_root <- safe_normalize(get_arg(
  "bundle-root", file.path(repo_dir, "..", "..")
))
recetox_msp <- safe_normalize(get_arg(
  "recetox-msp", file.path(bundle_root, "inputs", "RECETOX_merged.msp")
))
spectra_rds <- get_arg("spectra-rds", "")
if (nzchar(spectra_rds)) spectra_rds <- safe_normalize(spectra_rds)

output_dir <- get_arg(
  "output-dir", file.path(bundle_root, "results", "current", "perturbation")
)
manifest_dir <- get_arg(
  "manifest-dir", file.path(bundle_root, "manifest", "perturbation")
)
log_file <- get_arg("log-file", file.path(output_dir, "run_perturbation_publication.log"))

methods <- split_character(
  "methods",
  paste(c(
    "ppm_wasserstein", "composite", "entropy_weighted",
    "entropy_unweighted", "cosine", "weighted_cosine", "hellinger"
  ), collapse = ",")
)
expected_methods <- c(
  "ppm_wasserstein", "composite", "entropy_weighted",
  "entropy_unweighted", "cosine", "weighted_cosine", "hellinger"
)
if (!length(methods) || anyDuplicated(methods)) {
  stop("--methods must contain one or more unique methods.")
}

jitter_levels <- split_numeric("jitter-levels", "30,100")
drift_levels <- split_numeric("drift-levels", "30,100")
uniform_shift_levels <- split_numeric("uniform-shift-levels", "15,30,45,100")
reference_anchor_error_levels <- split_numeric(
  "reference-anchor-error-ppm", "0,5,10,20,50,100"
)
if (any(c(jitter_levels, drift_levels, uniform_shift_levels) <= 0)) {
  stop("All perturbation magnitudes must be > 0 ppm.")
}
if (any(reference_anchor_error_levels < 0) || anyDuplicated(reference_anchor_error_levels)) {
  stop("--reference-anchor-error-ppm must contain unique values >= 0.")
}

n_replicates <- as_single_integer("replicates", "3", lower = 1L)
base_seed <- as_single_integer("seed", "20260718", lower = 1L)
n_cores <- as_single_integer("n-cores", "8", lower = 1L)
max_spectra <- as_single_integer("max-spectra", "0", lower = 0L)
smoke_spectra <- as_single_integer("smoke-spectra", "6", lower = 2L)
smoke <- flag_is_true("smoke", FALSE)
overwrite <- flag_is_true("overwrite", FALSE)
save_distances <- flag_is_true("save-distances", FALSE)

tol_ppm <- as_single_numeric("tol-ppm", "15", lower = 0, lower_inclusive = FALSE)
transition_mult <- as_single_numeric(
  "transition-mult", "3", lower = 0, lower_inclusive = FALSE
)
sinkhorn_epsilon <- as_single_numeric(
  "sinkhorn-epsilon", "0.05", lower = 0, lower_inclusive = FALSE
)
sinkhorn_iterations <- as_single_integer("sinkhorn-iterations", "100", lower = 1L)
ot_method <- get_arg("ot-method", "exact")
use_parallel <- flag_is_true("parallel", TRUE)

if (smoke) {
  n_replicates <- 1L
  max_spectra <- if (max_spectra > 0L) min(max_spectra, smoke_spectra) else smoke_spectra
  use_parallel <- FALSE
}

score_dir <- file.path(output_dir, "score_level")
regenerated_dir <- file.path(output_dir, "regenerated_representation")
comparison_dir <- file.path(output_dir, "comparison")
reference_anchor_dir <- file.path(output_dir, "reference_anchor_error")

final_outputs <- c(
  file.path(score_dir, "per_query_tie_aware.csv"),
  file.path(score_dir, "summary_tie_aware.csv"),
  file.path(score_dir, "summary_across_replicates_tie_aware.csv"),
  file.path(regenerated_dir, "per_query_tie_aware.csv"),
  file.path(regenerated_dir, "summary_tie_aware.csv"),
  file.path(regenerated_dir, "summary_across_replicates_tie_aware.csv"),
  file.path(regenerated_dir, "representation_regeneration_qc.csv"),
  file.path(comparison_dir, "paired_per_query.csv"),
  file.path(comparison_dir, "paired_summary.csv"),
  file.path(comparison_dir, "paired_summary_across_replicates.csv"),
  file.path(reference_anchor_dir, "per_query_tie_aware.csv"),
  file.path(reference_anchor_dir, "summary_tie_aware.csv"),
  file.path(reference_anchor_dir, "pairwise_invariance_qc.csv"),
  file.path(manifest_dir, "run_manifest.rds")
)
existing_outputs <- final_outputs[file.exists(final_outputs)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite an existing perturbation run. Use a new --output-dir/",
    "--manifest-dir, or explicitly pass --overwrite=true. Existing file: ",
    existing_outputs[[1]], call. = FALSE
  )
}

for (path in c(
  output_dir, score_dir, regenerated_dir, comparison_dir,
  reference_anchor_dir, manifest_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

if (nzchar(log_file)) {
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log_connection <- file(log_file, open = "wt")
  sink(log_connection, split = TRUE)
  sink(log_connection, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_connection)
  }, add = TRUE)
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this script from the development tree.")
}
pkgload::load_all(repo_dir, quiet = TRUE)

statistics_helper <- file.path(
  repo_dir, "inst", "scripts", "lib", "publication_retrieval_statistics.R"
)
if (!file.exists(statistics_helper)) stop("Missing statistics helper: ", statistics_helper)
source(statistics_helper, local = TRUE)

parse_ei_internal <- getFromNamespace("parse_ei", "ppmWass")
process_spectrum_internal <- getFromNamespace("process_single_spectrum", "ppmWass")
build_derived_internal <- getFromNamespace("build_loss_peaks", "ppmWass")

RNGkind("L'Ecuyer-CMRG")

params <- eihrms_default_params()
params$min_mz <- 35
params$max_mz <- 650
params$centroid_ppm <- 10
params$noise_thr <- 0.01
params$topK <- 200L
params$use_sqrt <- FALSE
params$loss_top_peaks <- 40L
params$loss_max_peaks <- 250L
params$loss_min <- 5
params$loss_max <- 350
params$use_typical_loss <- FALSE
params$use_split_loss <- FALSE
params$use_mref_confidence <- FALSE
params$use_extended_losses <- FALSE
params$use_derivatization_losses <- FALSE
params$derivatization_mode <- "off"
params$tol_ppm <- tol_ppm
params$wasserstein_transition_mult <- transition_mult
params$ot_method <- ot_method
params$sinkhorn_epsilon <- sinkhorn_epsilon
params$sinkhorn_niter <- sinkhorn_iterations
params$mass_power <- 3
params$intensity_power <- 0.5
params$w_frag <- 0.70
params$w_loss <- 0.30
params$use_parallel <- use_parallel
params$n_cores <- n_cores
params$class_detection_ppm <- 15

for (method in methods) {
  p <- params
  p$distance_method <- method
  validate_params(p)
}

if (!smoke && !identical(methods, expected_methods)) {
  warning(
    "This is not the complete seven-method publication set. Methods: ",
    paste(methods, collapse = ", ")
  )
}

subset_spectra_local <- function(spectra, idx) {
  old_n <- nrow(spectra$df_spec)
  out <- spectra
  out$df_spec <- spectra$df_spec[idx, , drop = FALSE]
  component_names <- c(
    "frag_list", "loss_list", "derived_list", "loss_typ_list",
    "loss_anchor_list", "loss_pair_list", "derived_anchor_list",
    "derived_pair_list", "loss_anchor_typ_list", "loss_pair_typ_list"
  )
  for (name in component_names) {
    value <- spectra[[name]]
    if (!is.null(value) && length(value) == old_n) out[[name]] <- value[idx]
  }
  for (name in c("mref_confidence", "ri")) {
    value <- spectra[[name]]
    if (!is.null(value) && length(value) == old_n) out[[name]] <- value[idx]
  }
  out
}

filter_publication_spectra <- function(spectra) {
  df <- spectra$df_spec
  metadata_ok <- !is.na(df$RI) & !is.na(df$inchikey) & nchar(df$inchikey) >= 14L
  fragment_ok <- vapply(spectra$frag_list, function(fragment) {
    is.matrix(fragment) && nrow(fragment) >= 2L && ncol(fragment) >= 2L &&
      all(is.finite(fragment[, 1:2, drop = FALSE])) && sum(fragment[, 2]) > 0
  }, logical(1))
  keep <- which(metadata_ok & fragment_ok)
  message(
    "RECETOX filtering: input=", nrow(df), " retained=", length(keep),
    " removed=", nrow(df) - length(keep)
  )
  subset_spectra_local(spectra, keep)
}

load_recetox <- function() {
  spectra <- if (nzchar(spectra_rds)) {
    message("Loading prepared spectra: ", spectra_rds)
    readRDS(spectra_rds)
  } else {
    message("Building RECETOX spectra from: ", recetox_msp)
    build_spectra_from_msp(recetox_msp, params, require_ri = FALSE, progress = TRUE)
  }
  required <- c("df_spec", "frag_list", "loss_list")
  missing <- setdiff(required, names(spectra))
  if (length(missing)) stop("Prepared spectra are missing: ", paste(missing, collapse = ", "))
  if (!all(c("id", "ei", "RI", "inchikey") %in% names(spectra$df_spec))) {
    stop("df_spec must contain id, ei, RI, and inchikey for the perturbation benchmark.")
  }
  spectra <- filter_publication_spectra(spectra)
  if (max_spectra > 0L && nrow(spectra$df_spec) > max_spectra) {
    spectra <- subset_spectra_local(spectra, seq_len(max_spectra))
    message("Deterministic --max-spectra subset retained: ", nrow(spectra$df_spec))
  }
  ids <- as.character(spectra$df_spec$id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Filtered RECETOX query IDs must be unique, nonmissing strings.")
  }
  names(spectra$frag_list) <- ids
  names(spectra$loss_list) <- ids
  spectra$derived_list <- spectra$loss_list
  spectra
}

spectra <- load_recetox()
ids <- as.character(spectra$df_spec$id)
raw_ei <- stats::setNames(as.character(spectra$df_spec$ei), ids)
library_frag <- spectra$frag_list
library_derived <- spectra$loss_list
n_queries <- length(ids)
if (n_queries < 2L) stop("At least two valid spectra are required.")

make_conditions <- function() {
  base <- rbind(
    data.frame(
      perturbation_kind = "random_jitter", magnitude_ppm = jitter_levels,
      kind_code = 1L, stringsAsFactors = FALSE
    ),
    data.frame(
      perturbation_kind = "systematic_drift", magnitude_ppm = drift_levels,
      kind_code = 2L, stringsAsFactors = FALSE
    ),
    data.frame(
      perturbation_kind = "uniform_shift", magnitude_ppm = uniform_shift_levels,
      kind_code = 3L, stringsAsFactors = FALSE
    )
  )
  base <- base[order(base$kind_code, base$magnitude_ppm), , drop = FALSE]
  rows <- lapply(seq_len(nrow(base)), function(i) {
    data.frame(base[i, , drop = FALSE], replicate = seq_len(n_replicates))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$condition_seed <- base_seed + out$kind_code * 100000000 +
    as.integer(round(out$magnitude_ppm * 1000)) * 100L + out$replicate
  if (any(out$condition_seed > .Machine$integer.max)) {
    stop("Derived condition seed exceeds R's integer seed range; use a smaller --seed.")
  }
  magnitude_label <- gsub("\\.", "p", format(out$magnitude_ppm, trim = TRUE, scientific = FALSE))
  out$condition_id <- sprintf(
    "%s_%sppm_rep%02d", out$perturbation_kind, magnitude_label, out$replicate
  )
  if (smoke) out <- out[1L, , drop = FALSE]
  out
}

conditions <- make_conditions()

make_reference_anchor_conditions <- function() {
  out <- data.frame(
    perturbation_kind = "reference_anchor_error",
    magnitude_ppm = reference_anchor_error_levels,
    replicate = 1L,
    stringsAsFactors = FALSE
  )
  out$condition_seed <- base_seed + 900000000 +
    as.integer(round(out$magnitude_ppm * 1000)) * 100L
  if (any(out$condition_seed > .Machine$integer.max)) {
    stop("Reference-anchor condition seed exceeds R's integer seed range.")
  }
  magnitude_label <- gsub(
    "\\.", "p", format(out$magnitude_ppm, trim = TRUE, scientific = FALSE)
  )
  out$condition_id <- paste0("reference_anchor_error_", magnitude_label, "ppm")
  out
}

reference_anchor_conditions <- make_reference_anchor_conditions()

query_fragment_seeds <- function(condition_seed) {
  as.integer(condition_seed + seq_len(n_queries) * 1009L)
}

query_derived_seeds <- function(condition_seed) {
  as.integer(condition_seed + 500000000 + seq_len(n_queries) * 1009L)
}

query_reference_anchor_seeds <- function(condition_seed) {
  # The systematic reference-anchor errors do not draw random values. This
  # deterministic seed ledger keeps query/method pairing explicit and leaves a
  # stable seed available if a random-error sensitivity is added separately.
  as.integer(condition_seed + seq_len(n_queries) * 1009L)
}

empty_spectrum <- function() {
  matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("mz", "intensity")))
}

coerce_peak_matrix <- function(spectrum) {
  if (is.null(spectrum) || !length(spectrum) || nrow(spectrum) == 0L) {
    return(empty_spectrum())
  }
  out <- as.matrix(spectrum[, 1:2, drop = FALSE])
  storage.mode(out) <- "double"
  colnames(out) <- c("mz", "intensity")
  out
}

summarize_derived_exact <- function(mz_values, intensity_values) {
  if (!length(mz_values) || sum(intensity_values, na.rm = TRUE) <= 0) {
    return(empty_spectrum())
  }
  derived <- tibble::tibble(mz = mz_values, intensity = intensity_values) |>
    dplyr::mutate(mz_bin = round(mz, 3)) |>
    dplyr::group_by(mz_bin) |>
    dplyr::summarise(intensity = sum(intensity), .groups = "drop") |>
    dplyr::arrange(mz_bin)
  if (nrow(derived) > params$loss_max_peaks) {
    derived <- derived |>
      dplyr::arrange(dplyr::desc(intensity)) |>
      utils::head(params$loss_max_peaks) |>
      dplyr::arrange(mz_bin)
  }
  derived <- derived |>
    dplyr::mutate(intensity = intensity / sum(intensity))
  cbind(mz = derived$mz_bin, intensity = derived$intensity)
}

raw_pairwise_components <- function(fragment) {
  fragment <- coerce_peak_matrix(fragment)
  if (nrow(fragment) < 2L) {
    return(list(mz = numeric(), intensity = numeric()))
  }
  mz <- fragment[, 1]
  intensity <- fragment[, 2]
  order_by_intensity <- order(intensity, decreasing = TRUE)
  k <- min(params$loss_top_peaks, length(order_by_intensity))
  top_mz <- mz[order_by_intensity[seq_len(k)]]
  top_intensity <- intensity[order_by_intensity[seq_len(k)]]
  difference_matrix <- abs(outer(top_mz, top_mz, "-"))
  intensity_matrix <- outer(top_intensity, top_intensity, FUN = pmin)
  upper <- upper.tri(difference_matrix, diag = FALSE)
  pair_mz <- as.numeric(difference_matrix[upper])
  pair_intensity <- as.numeric(intensity_matrix[upper])
  keep <- is.finite(pair_mz) & pair_mz >= params$loss_min &
    pair_mz <= params$loss_max
  list(mz = pair_mz[keep], intensity = pair_intensity[keep])
}

max_matrix_component_difference <- function(current, reference, column) {
  if (identical(current, reference)) return(0)
  if (!identical(dim(current), dim(reference)) || !nrow(current)) return(Inf)
  max(abs(current[, column] - reference[, column]), na.rm = TRUE)
}

prepare_reference_anchor_baseline <- function() {
  anchored <- vector("list", n_queries)
  pairwise <- vector("list", n_queries)
  pooled <- vector("list", n_queries)
  raw_pairwise <- vector("list", n_queries)
  estimated_mref <- numeric(n_queries)
  mref_confidence <- numeric(n_queries)

  for (i in seq_len(n_queries)) {
    info <- build_derived_internal(
      library_frag[[i]], params, return_info = TRUE, deriv_type = "none"
    )
    raw_pair <- raw_pairwise_components(library_frag[[i]])
    pair_from_raw <- summarize_derived_exact(raw_pair$mz, raw_pair$intensity)
    if (!identical(pair_from_raw, info$lossB_peaks)) {
      stop(
        "Reference-anchor baseline could not exactly reproduce the pairwise ",
        "derived representation for query ", ids[[i]], call. = FALSE
      )
    }
    if (!identical(info$loss_peaks, library_derived[[i]])) {
      stop(
        "Reference-anchor baseline does not match the publication pooled ",
        "derived representation for query ", ids[[i]], call. = FALSE
      )
    }
    anchored[[i]] <- info$lossA_peaks
    pairwise[[i]] <- info$lossB_peaks
    pooled[[i]] <- info$loss_peaks
    raw_pairwise[[i]] <- raw_pair
    estimated_mref[[i]] <- info$Mref
    mref_confidence[[i]] <- info$mref_confidence
  }
  for (object_name in c("anchored", "pairwise", "pooled", "raw_pairwise")) {
    object <- get(object_name)
    names(object) <- ids
    assign(object_name, object)
  }
  names(estimated_mref) <- ids
  names(mref_confidence) <- ids
  list(
    anchored = anchored,
    pairwise = pairwise,
    pooled = pooled,
    raw_pairwise = raw_pairwise,
    estimated_mref = estimated_mref,
    mref_confidence = mref_confidence
  )
}

build_reference_anchor_queries <- function(condition, baseline) {
  query_seeds <- query_reference_anchor_seeds(condition$condition_seed)
  fragment <- library_frag
  anchored <- vector("list", n_queries)
  pairwise <- baseline$pairwise
  pooled <- vector("list", n_queries)
  qc <- vector("list", n_queries)

  for (i in seq_len(n_queries)) {
    fragment_i <- library_frag[[i]]
    original_mref <- baseline$estimated_mref[[i]]
    perturbed_mref <- original_mref * (1 + condition$magnitude_ppm * 1e-6)
    anchor_mz <- perturbed_mref - fragment_i[, 1]
    anchor_intensity <- fragment_i[, 2]
    if (isTRUE(params$use_mref_confidence)) {
      confidence_power <- if (is.null(params$mref_conf_power)) {
        1
      } else {
        params$mref_conf_power
      }
      anchor_intensity <- anchor_intensity *
        (baseline$mref_confidence[[i]]^confidence_power)
    }
    keep_anchor <- is.finite(anchor_mz) & anchor_mz >= params$loss_min &
      anchor_mz <= params$loss_max
    anchor_mz <- anchor_mz[keep_anchor]
    anchor_intensity <- anchor_intensity[keep_anchor]

    raw_pair <- baseline$raw_pairwise[[i]]
    anchored[[i]] <- summarize_derived_exact(anchor_mz, anchor_intensity)
    pooled[[i]] <- summarize_derived_exact(
      c(anchor_mz, raw_pair$mz),
      c(anchor_intensity, raw_pair$intensity)
    )

    fragment_identical <- identical(fragment_i, library_frag[[i]])
    pairwise_identical <- identical(pairwise[[i]], baseline$pairwise[[i]])
    pairwise_serialized_identical <- identical(
      serialize(pairwise[[i]], NULL, version = 3),
      serialize(baseline$pairwise[[i]], NULL, version = 3)
    )
    error_zero <- identical(as.numeric(condition$magnitude_ppm), 0)
    anchored_baseline_identical <- identical(anchored[[i]], baseline$anchored[[i]])
    pooled_baseline_identical <- identical(pooled[[i]], baseline$pooled[[i]])

    if (!fragment_identical || !pairwise_identical || !pairwise_serialized_identical) {
      stop(
        "Reference-anchor invariance failure for query ", ids[[i]],
        " at ", condition$magnitude_ppm, " ppm.", call. = FALSE
      )
    }
    if (error_zero && (!anchored_baseline_identical || !pooled_baseline_identical)) {
      stop(
        "The 0-ppm reference-anchor reconstruction is not bit-identical to ",
        "baseline for query ", ids[[i]], call. = FALSE
      )
    }

    qc[[i]] <- data.frame(
      condition_id = condition$condition_id,
      reference_anchor_error_ppm = condition$magnitude_ppm,
      condition_seed = condition$condition_seed,
      query_id = ids[[i]],
      reference_anchor_seed = query_seeds[[i]],
      estimated_mref = original_mref,
      perturbed_mref = perturbed_mref,
      reference_error_da = perturbed_mref - original_mref,
      realized_reference_error_ppm =
        (perturbed_mref - original_mref) / original_mref * 1e6,
      fragment_bit_identical = fragment_identical,
      pairwise_bit_identical = pairwise_identical,
      pairwise_serialized_bit_identical = pairwise_serialized_identical,
      pairwise_mz_max_abs_difference = max_matrix_component_difference(
        pairwise[[i]], baseline$pairwise[[i]], 1L
      ),
      pairwise_intensity_max_abs_difference = max_matrix_component_difference(
        pairwise[[i]], baseline$pairwise[[i]], 2L
      ),
      anchored_bit_identical_to_baseline = anchored_baseline_identical,
      pooled_bit_identical_to_baseline = pooled_baseline_identical,
      fragment_count = nrow(fragment_i),
      anchored_count = nrow(anchored[[i]]),
      pairwise_count = nrow(pairwise[[i]]),
      pooled_count = nrow(pooled[[i]]),
      stringsAsFactors = FALSE
    )
  }
  names(anchored) <- ids
  names(pooled) <- ids
  list(
    frag_list = fragment,
    derived_list = pooled,
    derived_anchor_list = anchored,
    derived_pair_list = pairwise,
    fragment_seeds = query_seeds,
    derived_seeds = query_seeds,
    qc = do.call(rbind, qc)
  )
}

perturb_coordinates <- function(spectrum, kind, magnitude_ppm, seed) {
  spectrum <- coerce_peak_matrix(spectrum)
  n <- nrow(spectrum)
  if (!n) {
    return(list(
      spectrum = spectrum,
      realized_ppm = numeric(),
      theta = NA_real_, drift_intercept_ppm = NA_real_,
      drift_slope_ppm_per_da = NA_real_, uniform_sign = NA_integer_
    ))
  }
  set.seed(seed)
  theta <- NA_real_
  drift_intercept <- NA_real_
  drift_slope <- NA_real_
  uniform_sign <- NA_integer_

  if (kind == "random_jitter") {
    shift_ppm <- stats::rnorm(n, mean = 0, sd = magnitude_ppm)
  } else if (kind == "systematic_drift") {
    theta <- stats::runif(1L, 0, 2 * pi)
    drift_intercept <- magnitude_ppm * cos(theta)
    drift_slope <- magnitude_ppm / 100 * sin(theta)
    shift_ppm <- drift_intercept + drift_slope * (spectrum[, 1] - 100)
  } else if (kind == "uniform_shift") {
    uniform_sign <- sample(c(-1L, 1L), 1L)
    shift_ppm <- rep(uniform_sign * magnitude_ppm, n)
  } else {
    stop("Unknown perturbation kind: ", kind)
  }

  out <- spectrum
  out[, 1] <- out[, 1] * (1 + shift_ppm * 1e-6)
  out <- out[order(out[, 1]), , drop = FALSE]
  list(
    spectrum = out,
    realized_ppm = shift_ppm,
    theta = theta,
    drift_intercept_ppm = drift_intercept,
    drift_slope_ppm_per_da = drift_slope,
    uniform_sign = uniform_sign
  )
}

draw_qc_row <- function(draw) {
  shift <- draw$realized_ppm
  finite <- is.finite(shift)
  shift <- shift[finite]
  data.frame(
    realized_mean_ppm = if (length(shift)) mean(shift) else NA_real_,
    realized_sd_ppm = if (length(shift) > 1L) stats::sd(shift) else 0,
    realized_rms_ppm = if (length(shift)) sqrt(mean(shift^2)) else NA_real_,
    realized_min_ppm = if (length(shift)) min(shift) else NA_real_,
    realized_max_ppm = if (length(shift)) max(shift) else NA_real_,
    drift_theta_rad = draw$theta,
    drift_intercept_ppm = draw$drift_intercept_ppm,
    drift_slope_ppm_per_da = draw$drift_slope_ppm_per_da,
    uniform_sign = draw$uniform_sign,
    stringsAsFactors = FALSE
  )
}

serialize_peak_string <- function(spectrum) {
  spectrum <- coerce_peak_matrix(spectrum)
  if (!nrow(spectrum)) return("")
  paste(
    sprintf("%.17g:%.17g", spectrum[, 1], spectrum[, 2]),
    collapse = " "
  )
}

build_score_level_queries <- function(condition) {
  fragment_seeds <- query_fragment_seeds(condition$condition_seed)
  derived_seeds <- query_derived_seeds(condition$condition_seed)
  fragment <- vector("list", n_queries)
  derived <- vector("list", n_queries)
  qc <- vector("list", n_queries * 2L)

  for (i in seq_len(n_queries)) {
    fragment_draw <- perturb_coordinates(
      library_frag[[i]], condition$perturbation_kind,
      condition$magnitude_ppm, fragment_seeds[[i]]
    )
    derived_draw <- perturb_coordinates(
      library_derived[[i]], condition$perturbation_kind,
      condition$magnitude_ppm, derived_seeds[[i]]
    )
    fragment[[i]] <- fragment_draw$spectrum
    derived[[i]] <- derived_draw$spectrum
    qc[[2L * i - 1L]] <- cbind(
      data.frame(
        query_id = ids[[i]], representation = "preprocessed_fragment",
        seed = fragment_seeds[[i]], n_peaks = nrow(fragment[[i]]),
        stringsAsFactors = FALSE
      ),
      draw_qc_row(fragment_draw)
    )
    qc[[2L * i]] <- cbind(
      data.frame(
        query_id = ids[[i]], representation = "precomputed_pooled_derived",
        seed = derived_seeds[[i]], n_peaks = nrow(derived[[i]]),
        stringsAsFactors = FALSE
      ),
      draw_qc_row(derived_draw)
    )
  }
  names(fragment) <- ids
  names(derived) <- ids
  list(
    frag_list = fragment,
    derived_list = derived,
    fragment_seeds = fragment_seeds,
    derived_seeds = derived_seeds,
    qc = do.call(rbind, qc)
  )
}

build_regenerated_queries <- function(condition) {
  fragment_seeds <- query_fragment_seeds(condition$condition_seed)
  fragment <- vector("list", n_queries)
  pooled <- vector("list", n_queries)
  anchored <- vector("list", n_queries)
  pairwise <- vector("list", n_queries)
  mref_confidence <- numeric(n_queries)
  qc <- vector("list", n_queries)

  for (i in seq_len(n_queries)) {
    raw_peaks <- parse_ei_internal(raw_ei[[i]])
    draw <- perturb_coordinates(
      raw_peaks, condition$perturbation_kind,
      condition$magnitude_ppm, fragment_seeds[[i]]
    )
    rebuilt <- process_spectrum_internal(
      serialize_peak_string(draw$spectrum), params, deriv_mode = "off"
    )
    fragment[[i]] <- rebuilt$frag
    pooled[[i]] <- rebuilt$loss
    anchored[[i]] <- rebuilt$loss_anchor
    pairwise[[i]] <- rebuilt$loss_pair
    mref_confidence[[i]] <- rebuilt$mref_confidence
    qc[[i]] <- cbind(
      data.frame(
        query_id = ids[[i]], fragment_seed = fragment_seeds[[i]],
        raw_peak_count = nrow(raw_peaks),
        perturbed_raw_peak_count = nrow(draw$spectrum),
        processed_fragment_count = nrow(rebuilt$frag),
        estimated_high_mass_reference = rebuilt$mref,
        mref_confidence = rebuilt$mref_confidence,
        anchored_difference_count = nrow(rebuilt$loss_anchor),
        pairwise_difference_count = nrow(rebuilt$loss_pair),
        pooled_derived_count = nrow(rebuilt$loss),
        stringsAsFactors = FALSE
      ),
      draw_qc_row(draw)
    )
  }
  for (value in c("fragment", "pooled", "anchored", "pairwise")) {
    object <- get(value)
    names(object) <- ids
    assign(value, object)
  }
  names(mref_confidence) <- ids
  list(
    frag_list = fragment,
    derived_list = pooled,
    derived_anchor_list = anchored,
    derived_pair_list = pairwise,
    mref_confidence = mref_confidence,
    fragment_seeds = fragment_seeds,
    derived_seeds = rep(NA_integer_, n_queries),
    qc = do.call(rbind, qc)
  )
}

per_query_self_tie_aware <- function(distance_matrix, condition, pipeline,
                                     method, fragment_seeds, derived_seeds,
                                     tie_seed) {
  query_ids <- rownames(distance_matrix)
  library_ids <- colnames(distance_matrix)
  if (!identical(query_ids, ids) || !identical(library_ids, ids)) {
    stop("Distance matrix IDs are not aligned to the paired query/library IDs.")
  }
  rows <- vector("list", length(query_ids))
  for (i in seq_along(query_ids)) {
    target <- match(query_ids[[i]], library_ids)
    relevant <- seq_along(library_ids) == target
    metrics <- first_relevant_tie_metrics(
      distance_matrix[i, ], relevant, ks = c(1L, 5L, 10L), tolerance = 0
    )
    set.seed(as.integer(tie_seed + i * 17L))
    random_rank <- random_first_relevant_rank(
      metrics$n_better, metrics$tie_size, metrics$relevant_in_tie
    )
    rows[[i]] <- data.frame(
      dataset = "RECETOX",
      pipeline = pipeline,
      method = method,
      condition_id = condition$condition_id,
      perturbation_kind = condition$perturbation_kind,
      magnitude_ppm = condition$magnitude_ppm,
      replicate = condition$replicate,
      condition_seed = condition$condition_seed,
      query_index = i,
      query_id = query_ids[[i]],
      target_library_id = library_ids[[target]],
      fragment_seed = fragment_seeds[[i]],
      derived_seed = derived_seeds[[i]],
      target_distance = distance_matrix[i, target],
      n_better = metrics$n_better,
      tie_size = metrics$tie_size,
      relevant_in_tie = metrics$relevant_in_tie,
      optimistic_rank = metrics$optimistic_rank,
      pessimistic_rank = metrics$pessimistic_rank,
      random_rank = random_rank,
      top1_optimistic = metrics$top1_optimistic,
      top1_fractional = metrics$top1_fractional,
      top1_pessimistic = metrics$top1_pessimistic,
      top1_random = as.numeric(random_rank <= 1L),
      top5_optimistic = metrics$top5_optimistic,
      top5_fractional = metrics$top5_fractional,
      top5_pessimistic = metrics$top5_pessimistic,
      top5_random = as.numeric(random_rank <= 5L),
      top10_optimistic = metrics$top10_optimistic,
      top10_fractional = metrics$top10_fractional,
      top10_pessimistic = metrics$top10_pessimistic,
      top10_random = as.numeric(random_rank <= 10L),
      rr_optimistic = 1 / metrics$optimistic_rank,
      rr_fractional = metrics$expected_rr,
      rr_pessimistic = 1 / metrics$pessimistic_rank,
      rr_random = 1 / random_rank,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_self_ties <- function(per_query) {
  mean_finite <- function(column) {
    x <- per_query[[column]]
    if (!any(is.finite(x))) return(NA_real_)
    mean(x[is.finite(x)])
  }
  identity_columns <- c(
    "dataset", "pipeline", "method", "condition_id", "perturbation_kind",
    "magnitude_ppm", "replicate", "condition_seed"
  )
  out <- per_query[1L, identity_columns, drop = FALSE]
  metric_columns <- c(
    "top1_optimistic", "top1_fractional", "top1_pessimistic", "top1_random",
    "top5_optimistic", "top5_fractional", "top5_pessimistic", "top5_random",
    "top10_optimistic", "top10_fractional", "top10_pessimistic", "top10_random",
    "rr_optimistic", "rr_fractional", "rr_pessimistic", "rr_random"
  )
  for (column in metric_columns) out[[column]] <- mean_finite(column)
  out$n_queries <- nrow(per_query)
  out$n_rank1_ties <- sum(per_query$n_better == 0 & per_query$tie_size > 1L, na.rm = TRUE)
  out$n_any_target_ties <- sum(per_query$tie_size > 1L, na.rm = TRUE)
  out$mean_target_distance <- mean_finite("target_distance")
  out
}

summarize_across_replicates <- function(per_query) {
  group_columns <- c(
    "dataset", "pipeline", "method", "perturbation_kind", "magnitude_ppm"
  )
  metric_columns <- c(
    "top1_optimistic", "top1_fractional", "top1_pessimistic", "top1_random",
    "top5_optimistic", "top5_fractional", "top5_pessimistic", "top5_random",
    "top10_optimistic", "top10_fractional", "top10_pessimistic", "top10_random",
    "rr_optimistic", "rr_fractional", "rr_pessimistic", "rr_random"
  )
  key <- interaction(per_query[group_columns], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(per_query)), key)
  rows <- lapply(groups, function(index) {
    data <- per_query[index, , drop = FALSE]
    out <- data[1L, group_columns, drop = FALSE]
    for (column in metric_columns) {
      value <- data[[column]]
      out[[column]] <- if (any(is.finite(value))) {
        mean(value[is.finite(value)])
      } else {
        NA_real_
      }
    }
    target_distance <- data$target_distance
    out$n_unique_queries <- length(unique(data$query_id))
    out$n_replicates <- length(unique(data$replicate))
    out$n_query_replicates <- nrow(data)
    out$n_rank1_ties <- sum(data$n_better == 0 & data$tie_size > 1L, na.rm = TRUE)
    out$n_any_target_ties <- sum(data$tie_size > 1L, na.rm = TRUE)
    out$mean_target_distance <- if (any(is.finite(target_distance))) {
      mean(target_distance[is.finite(target_distance)])
    } else {
      NA_real_
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

check_finite_distances <- function(distance_matrix, condition, pipeline, method) {
  bad <- which(!is.finite(distance_matrix), arr.ind = TRUE)
  if (!nrow(bad)) return(invisible(TRUE))
  diagnostic <- data.frame(
    condition_id = condition$condition_id,
    pipeline = pipeline,
    method = method,
    query_id = rownames(distance_matrix)[bad[, 1]],
    library_id = colnames(distance_matrix)[bad[, 2]],
    value = distance_matrix[bad],
    stringsAsFactors = FALSE
  )
  diagnostic_path <- file.path(
    output_dir,
    paste0("nonfinite_", condition$condition_id, "_", pipeline, "_", method, ".csv")
  )
  utils::write.csv(diagnostic, diagnostic_path, row.names = FALSE)
  stop(
    "Nonfinite distances remain after package fallback for ",
    condition$condition_id, "/", pipeline, "/", method,
    ". Diagnostic: ", diagnostic_path, call. = FALSE
  )
}

compute_one <- function(query, condition, pipeline, method, method_index) {
  p <- params
  p$distance_method <- method
  p <- validate_params(p)
  message(
    "[", condition$condition_id, "] ", pipeline, " | ", method,
    " | queries=", n_queries
  )
  distance_matrix <- compute_distance_matrix_search(
    query_frag_list = query$frag_list,
    query_loss_list = query$derived_list,
    lib_frag_list = library_frag,
    lib_loss_list = library_derived,
    params = p,
    progress = FALSE
  )
  check_finite_distances(distance_matrix, condition, pipeline, method)
  if (save_distances) {
    pipeline_dir <- switch(
      pipeline,
      score_level = score_dir,
      regenerated_representation = regenerated_dir,
      reference_anchor_error = reference_anchor_dir,
      stop("Unknown output pipeline: ", pipeline)
    )
    distance_dir <- file.path(
      pipeline_dir,
      "distance_matrices"
    )
    dir.create(distance_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      distance_matrix,
      file.path(distance_dir, paste0(condition$condition_id, "__", method, ".rds"))
    )
  }
  tie_seed <- as.integer(condition$condition_seed + method_index * 100000L)
  per_query <- per_query_self_tie_aware(
    distance_matrix, condition, pipeline, method,
    query$fragment_seeds, query$derived_seeds, tie_seed
  )
  list(per_query = per_query, summary = summarize_self_ties(per_query))
}

flatten_parameter <- function(value) {
  if (is.null(value)) return("<NULL>")
  if (is.function(value)) return("<function>")
  if (is.atomic(value)) return(paste(as.character(value), collapse = ";"))
  paste(capture.output(str(value, give.attr = FALSE)), collapse = " ")
}

git_output <- suppressWarnings(tryCatch(
  system2(
    "git", c("-C", repo_dir, "rev-parse", "HEAD"),
    stdout = TRUE, stderr = FALSE
  ),
  error = function(e) character()
))
git_commit <- if (length(git_output) && is.null(attr(git_output, "status"))) {
  git_output[[1]]
} else {
  NA_character_
}
git_status <- suppressWarnings(tryCatch(
  system2(
    "git", c("-C", repo_dir, "status", "--porcelain"),
    stdout = TRUE, stderr = FALSE
  ),
  error = function(e) character()
))
package_tree_dirty <- length(git_status) > 0L

definitions <- data.frame(
  perturbation_kind = c(
    "random_jitter", "systematic_drift", "uniform_shift",
    "reference_anchor_error"
  ),
  definition = c(
    paste0(
      "Independent per-peak Gaussian displacement: epsilon_i ~ N(0, magnitude_ppm); ",
      "mz_i' = mz_i * (1 + epsilon_i * 1e-6)."
    ),
    paste0(
      "One theta ~ Uniform(0, 2*pi) per query; a = magnitude_ppm*cos(theta), ",
      "b = magnitude_ppm/100*sin(theta), epsilon_i = a + b*(mz_i - 100), ",
      "mz_i' = mz_i*(1 + epsilon_i*1e-6). Magnitude is the legacy scale ",
      "parameter, not the maximum realized ppm shift."
    ),
    paste0(
      "One sign sampled uniformly from {-1,+1} per query; all peaks receive ",
      "epsilon_i = sign*magnitude_ppm and mz_i' = mz_i*(1 + epsilon_i*1e-6)."
    ),
    paste0(
      "Systematic positive reference error: Mref' = Mref*(1 + error_ppm*1e-6). ",
      "Fragment peaks are unchanged. Anchored coordinates are rebuilt as ",
      "Mref' - fragment_mz; the cached pairwise object is retained bit-for-bit; ",
      "pooled derived peaks are rebuilt from corrected unsummarized anchored ",
      "components plus the original unsummarized pairwise components using the ",
      "package's 0.001-Da binning, truncation, and normalization order."
    )
  ),
  score_level_application = c(
    rep(
      paste0(
        "Applied directly and independently to the preprocessed fragment and ",
        "precomputed pooled-derived coordinates; no downstream regeneration."
      ),
      3L
    ),
    "Not applicable: evaluated in the separate reference_anchor_error pipeline."
  ),
  regenerated_application = c(
    rep(
      paste0(
        "Applied to the raw fragment peak string, followed by preprocessing, ",
        "high-mass-reference re-estimation, anchored and pairwise difference ",
        "regeneration, pooling, and scoring."
      ),
      3L
    ),
    "Not applicable: fragment acquisition coordinates remain unchanged."
  ),
  reference_anchor_application = c(
    rep("Not applicable: reference Mref is not externally perturbed.", 3L),
    paste0(
      "Only the estimated query Mref and anchored derived component change; ",
      "fragment and pairwise derived components are required to be bit-identical ",
      "to baseline before scoring."
    )
  ),
  stringsAsFactors = FALSE
)

condition_manifest <- conditions
condition_manifest$dataset <- "RECETOX"
condition_manifest$fragment_seed_rule <- "condition_seed + query_index * 1009"
condition_manifest$score_level_derived_seed_rule <-
  "condition_seed + 500000000 + query_index * 1009"
condition_manifest$regenerated_derived_seed_rule <-
  "none: derived representations are regenerated from the fragment perturbation"

reference_anchor_condition_manifest <- reference_anchor_conditions
reference_anchor_condition_manifest$dataset <- "RECETOX"
reference_anchor_condition_manifest$reference_error_formula <-
  "Mref_perturbed = Mref_estimated * (1 + error_ppm * 1e-6)"
reference_anchor_condition_manifest$query_seed_rule <-
  "condition_seed + query_index * 1009 (ledger only; systematic error uses no RNG draw)"
reference_anchor_condition_manifest$fragment_policy <- "bit-identical to baseline"
reference_anchor_condition_manifest$pairwise_policy <- "bit-identical to baseline"
reference_anchor_condition_manifest$pooled_policy <- paste(
  "rebuild from corrected unsummarized anchored components plus original",
  "unsummarized pairwise components"
)

effective_parameters <- data.frame(
  parameter = names(params),
  value = vapply(params, flatten_parameter, character(1)),
  stringsAsFactors = FALSE
)

input_manifest <- data.frame(
  dataset = "RECETOX",
  msp_path = recetox_msp,
  msp_md5 = unname(tools::md5sum(recetox_msp)),
  prepared_spectra_rds = if (nzchar(spectra_rds)) spectra_rds else NA_character_,
  prepared_spectra_md5 = if (nzchar(spectra_rds)) {
    unname(tools::md5sum(spectra_rds))
  } else {
    NA_character_
  },
  n_filtered_queries = n_queries,
  stringsAsFactors = FALSE
)

run_metadata <- list(
  script = safe_normalize(file.path(repo_dir, "inst", "scripts", "run_perturbation_publication.R")),
  command_args = args,
  timestamp_started = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  package_git_commit = git_commit,
  package_tree_dirty = package_tree_dirty,
  dataset = "RECETOX",
  input_manifest = input_manifest,
  methods = methods,
  conditions = condition_manifest,
  reference_anchor_conditions = reference_anchor_condition_manifest,
  perturbation_definitions = definitions,
  effective_parameters = params,
  pipeline_definitions = list(
    score_level = paste(
      "Direct m/z-coordinate perturbation of preprocessed fragment and pooled",
      "derived representations. The fragment and derived branches use separate",
      "documented query seeds, matching the historical score-level design."
    ),
    regenerated_representation = paste(
      "Raw EI fragment coordinates are perturbed with the paired fragment seed;",
      "process_single_spectrum() then reruns preprocessing, high-mass-reference",
      "estimation, anchored/pairwise construction, and pooled-derived construction."
    ),
    reference_anchor_error = paste(
      "Fragment spectra remain unchanged. A systematic positive ppm error is",
      "applied only to each query's estimated high-mass reference; anchored",
      "differences and the pooled derived representation are rebuilt while the",
      "cached pairwise derived object must remain bit-identical."
    )
  ),
  retrieval_endpoint = paste(
    "Query-only self-retention against the original unperturbed library;",
    "fractional expected tie handling is primary, with optimistic, pessimistic,",
    "and fixed-seed randomized sensitivity columns."
  ),
  smoke = smoke,
  max_spectra = max_spectra,
  save_distances = save_distances
)

utils::write.csv(
  definitions, file.path(manifest_dir, "perturbation_definitions.csv"), row.names = FALSE
)
utils::write.csv(
  condition_manifest, file.path(manifest_dir, "perturbation_conditions_and_seeds.csv"),
  row.names = FALSE
)
utils::write.csv(
  reference_anchor_condition_manifest,
  file.path(manifest_dir, "reference_anchor_error_conditions_and_seeds.csv"),
  row.names = FALSE
)
utils::write.csv(
  effective_parameters, file.path(manifest_dir, "effective_parameters.csv"), row.names = FALSE
)
utils::write.csv(
  input_manifest, file.path(manifest_dir, "input_manifest.csv"), row.names = FALSE
)
writeLines(
  c(
    paste("Rscript", safe_normalize(file.path(repo_dir, "inst", "scripts", "run_perturbation_publication.R"))),
    paste(args, collapse = " ")
  ),
  file.path(manifest_dir, "command.txt"), useBytes = TRUE
)
writeLines(capture.output(sessionInfo()), file.path(manifest_dir, "sessionInfo.txt"))
saveRDS(run_metadata, file.path(manifest_dir, "run_manifest.rds"))

status_file <- file.path(manifest_dir, "run_status.txt")
run_completed <- FALSE
writeLines(
  paste("started", run_metadata$timestamp_started), status_file, useBytes = TRUE
)
on.exit({
  if (!run_completed) {
    writeLines(
      paste("incomplete", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      status_file, useBytes = TRUE
    )
  }
}, add = TRUE)

message("Repository: ", repo_dir)
message("RECETOX input: ", recetox_msp)
message("Queries: ", n_queries)
message("Methods: ", paste(methods, collapse = ", "))
message("Conditions: ", nrow(conditions))
message("Output: ", output_dir)
message("Smoke mode: ", smoke)

score_per_query <- list()
score_summary <- list()
score_qc <- list()
regenerated_per_query <- list()
regenerated_summary <- list()
regenerated_qc <- list()
result_index <- 0L

for (condition_index in seq_len(nrow(conditions))) {
  condition <- conditions[condition_index, , drop = FALSE]
  message("\n=== ", condition$condition_id, " ===")
  score_query <- build_score_level_queries(condition)
  regenerated_query <- build_regenerated_queries(condition)

  score_qc[[condition_index]] <- cbind(
    condition[, c(
      "condition_id", "perturbation_kind", "magnitude_ppm", "replicate",
      "condition_seed"
    ), drop = FALSE][rep(1L, nrow(score_query$qc)), , drop = FALSE],
    score_query$qc
  )
  regenerated_qc[[condition_index]] <- cbind(
    condition[, c(
      "condition_id", "perturbation_kind", "magnitude_ppm", "replicate",
      "condition_seed"
    ), drop = FALSE][rep(1L, nrow(regenerated_query$qc)), , drop = FALSE],
    regenerated_query$qc
  )

  for (method_index in seq_along(methods)) {
    method <- methods[[method_index]]
    result_index <- result_index + 1L
    score_result <- compute_one(
      score_query, condition, "score_level", method, method_index
    )
    regenerated_result <- compute_one(
      regenerated_query, condition, "regenerated_representation", method, method_index
    )
    score_per_query[[result_index]] <- score_result$per_query
    score_summary[[result_index]] <- score_result$summary
    regenerated_per_query[[result_index]] <- regenerated_result$per_query
    regenerated_summary[[result_index]] <- regenerated_result$summary
  }
}

message("\n=== Reference-anchor error analysis ===")
reference_anchor_baseline <- prepare_reference_anchor_baseline()
reference_per_query <- list()
reference_summary <- list()
reference_invariance_qc <- list()
reference_result_index <- 0L

for (condition_index in seq_len(nrow(reference_anchor_conditions))) {
  condition <- reference_anchor_conditions[condition_index, , drop = FALSE]
  reference_query <- build_reference_anchor_queries(
    condition, reference_anchor_baseline
  )
  reference_invariance_qc[[condition_index]] <- reference_query$qc

  for (method_index in seq_along(methods)) {
    method <- methods[[method_index]]
    reference_result_index <- reference_result_index + 1L
    result <- compute_one(
      reference_query, condition, "reference_anchor_error", method, method_index
    )
    query_match <- match(result$per_query$query_id, reference_query$qc$query_id)
    result$per_query$reference_anchor_error_ppm <- condition$magnitude_ppm
    result$per_query$reference_anchor_seed <-
      reference_query$qc$reference_anchor_seed[query_match]
    result$per_query$estimated_mref <-
      reference_query$qc$estimated_mref[query_match]
    result$per_query$perturbed_mref <-
      reference_query$qc$perturbed_mref[query_match]
    result$per_query$reference_error_da <-
      reference_query$qc$reference_error_da[query_match]
    reference_per_query[[reference_result_index]] <- result$per_query
    reference_summary[[reference_result_index]] <- result$summary
  }
}

score_per_query <- do.call(rbind, score_per_query)
score_summary <- do.call(rbind, score_summary)
score_qc <- do.call(rbind, score_qc)
regenerated_per_query <- do.call(rbind, regenerated_per_query)
regenerated_summary <- do.call(rbind, regenerated_summary)
regenerated_qc <- do.call(rbind, regenerated_qc)
reference_per_query <- do.call(rbind, reference_per_query)
reference_summary <- do.call(rbind, reference_summary)
reference_invariance_qc <- do.call(rbind, reference_invariance_qc)
rownames(score_per_query) <- NULL
rownames(score_summary) <- NULL
rownames(score_qc) <- NULL
rownames(regenerated_per_query) <- NULL
rownames(regenerated_summary) <- NULL
rownames(regenerated_qc) <- NULL
rownames(reference_per_query) <- NULL
rownames(reference_summary) <- NULL
rownames(reference_invariance_qc) <- NULL
score_across_replicates <- summarize_across_replicates(score_per_query)
regenerated_across_replicates <- summarize_across_replicates(regenerated_per_query)

utils::write.csv(
  score_per_query, file.path(score_dir, "per_query_tie_aware.csv"), row.names = FALSE
)
utils::write.csv(
  score_summary, file.path(score_dir, "summary_tie_aware.csv"), row.names = FALSE
)
utils::write.csv(
  score_across_replicates,
  file.path(score_dir, "summary_across_replicates_tie_aware.csv"),
  row.names = FALSE
)
utils::write.csv(
  score_qc, file.path(score_dir, "direct_perturbation_qc.csv"), row.names = FALSE
)
saveRDS(
  list(
    per_query = score_per_query,
    summary = score_summary,
    summary_across_replicates = score_across_replicates,
    qc = score_qc
  ),
  file.path(score_dir, "score_level_results.rds")
)

utils::write.csv(
  regenerated_per_query,
  file.path(regenerated_dir, "per_query_tie_aware.csv"), row.names = FALSE
)
utils::write.csv(
  regenerated_summary,
  file.path(regenerated_dir, "summary_tie_aware.csv"), row.names = FALSE
)
utils::write.csv(
  regenerated_across_replicates,
  file.path(regenerated_dir, "summary_across_replicates_tie_aware.csv"),
  row.names = FALSE
)
utils::write.csv(
  regenerated_qc,
  file.path(regenerated_dir, "representation_regeneration_qc.csv"), row.names = FALSE
)
saveRDS(
  list(
    per_query = regenerated_per_query,
    summary = regenerated_summary,
    summary_across_replicates = regenerated_across_replicates,
    representation_qc = regenerated_qc
  ),
  file.path(regenerated_dir, "regenerated_representation_results.rds")
)

utils::write.csv(
  reference_per_query,
  file.path(reference_anchor_dir, "per_query_tie_aware.csv"), row.names = FALSE
)
utils::write.csv(
  reference_summary,
  file.path(reference_anchor_dir, "summary_tie_aware.csv"), row.names = FALSE
)
utils::write.csv(
  reference_invariance_qc,
  file.path(reference_anchor_dir, "pairwise_invariance_qc.csv"), row.names = FALSE
)
saveRDS(
  list(
    per_query = reference_per_query,
    summary = reference_summary,
    pairwise_invariance_qc = reference_invariance_qc
  ),
  file.path(reference_anchor_dir, "reference_anchor_error_results.rds")
)

comparison_keys <- c(
  "dataset", "method", "condition_id", "perturbation_kind", "magnitude_ppm",
  "replicate", "condition_seed", "query_index", "query_id", "target_library_id",
  "fragment_seed"
)
comparison_metrics <- c(
  "target_distance", "n_better", "tie_size", "optimistic_rank", "pessimistic_rank",
  "top1_optimistic", "top1_fractional", "top1_pessimistic", "top1_random",
  "top5_optimistic", "top5_fractional", "top5_pessimistic", "top5_random",
  "top10_optimistic", "top10_fractional", "top10_pessimistic", "top10_random",
  "rr_optimistic", "rr_fractional", "rr_pessimistic", "rr_random"
)
score_compare <- score_per_query[, c(comparison_keys, comparison_metrics), drop = FALSE]
regenerated_compare <- regenerated_per_query[, c(comparison_keys, comparison_metrics), drop = FALSE]
names(score_compare)[match(comparison_metrics, names(score_compare))] <-
  paste0(comparison_metrics, "_score_level")
names(regenerated_compare)[match(comparison_metrics, names(regenerated_compare))] <-
  paste0(comparison_metrics, "_regenerated")
paired_per_query <- merge(
  score_compare, regenerated_compare, by = comparison_keys, all = FALSE, sort = FALSE
)
for (metric in comparison_metrics) {
  paired_per_query[[paste0("delta_regenerated_minus_score_", metric)]] <-
    paired_per_query[[paste0(metric, "_regenerated")]] -
    paired_per_query[[paste0(metric, "_score_level")]]
}

summary_keys <- c(
  "dataset", "method", "condition_id", "perturbation_kind", "magnitude_ppm",
  "replicate", "condition_seed"
)
summary_metrics <- c(
  "top1_optimistic", "top1_fractional", "top1_pessimistic", "top1_random",
  "top5_optimistic", "top5_fractional", "top5_pessimistic", "top5_random",
  "top10_optimistic", "top10_fractional", "top10_pessimistic", "top10_random",
  "rr_optimistic", "rr_fractional", "rr_pessimistic", "rr_random",
  "n_rank1_ties", "n_any_target_ties", "mean_target_distance"
)
score_summary_compare <- score_summary[, c(summary_keys, summary_metrics), drop = FALSE]
regenerated_summary_compare <-
  regenerated_summary[, c(summary_keys, summary_metrics), drop = FALSE]
names(score_summary_compare)[match(summary_metrics, names(score_summary_compare))] <-
  paste0(summary_metrics, "_score_level")
names(regenerated_summary_compare)[match(summary_metrics, names(regenerated_summary_compare))] <-
  paste0(summary_metrics, "_regenerated")
paired_summary <- merge(
  score_summary_compare, regenerated_summary_compare,
  by = summary_keys, all = FALSE, sort = FALSE
)
for (metric in summary_metrics) {
  paired_summary[[paste0("delta_regenerated_minus_score_", metric)]] <-
    paired_summary[[paste0(metric, "_regenerated")]] -
    paired_summary[[paste0(metric, "_score_level")]]
}

aggregate_keys <- c("dataset", "method", "perturbation_kind", "magnitude_ppm")
aggregate_metrics <- c(
  "top1_optimistic", "top1_fractional", "top1_pessimistic", "top1_random",
  "top5_optimistic", "top5_fractional", "top5_pessimistic", "top5_random",
  "top10_optimistic", "top10_fractional", "top10_pessimistic", "top10_random",
  "rr_optimistic", "rr_fractional", "rr_pessimistic", "rr_random",
  "n_rank1_ties", "n_any_target_ties", "mean_target_distance"
)
score_aggregate_compare <-
  score_across_replicates[, c(aggregate_keys, aggregate_metrics), drop = FALSE]
regenerated_aggregate_compare <-
  regenerated_across_replicates[, c(aggregate_keys, aggregate_metrics), drop = FALSE]
names(score_aggregate_compare)[match(aggregate_metrics, names(score_aggregate_compare))] <-
  paste0(aggregate_metrics, "_score_level")
names(regenerated_aggregate_compare)[
  match(aggregate_metrics, names(regenerated_aggregate_compare))
] <- paste0(aggregate_metrics, "_regenerated")
paired_across_replicates <- merge(
  score_aggregate_compare, regenerated_aggregate_compare,
  by = aggregate_keys, all = FALSE, sort = FALSE
)
for (metric in aggregate_metrics) {
  paired_across_replicates[[paste0("delta_regenerated_minus_score_", metric)]] <-
    paired_across_replicates[[paste0(metric, "_regenerated")]] -
    paired_across_replicates[[paste0(metric, "_score_level")]]
}

utils::write.csv(
  paired_per_query, file.path(comparison_dir, "paired_per_query.csv"), row.names = FALSE
)
utils::write.csv(
  paired_summary, file.path(comparison_dir, "paired_summary.csv"), row.names = FALSE
)
utils::write.csv(
  paired_across_replicates,
  file.path(comparison_dir, "paired_summary_across_replicates.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    per_query = paired_per_query,
    summary = paired_summary,
    summary_across_replicates = paired_across_replicates
  ),
  file.path(comparison_dir, "paired_pipeline_comparison.rds")
)

run_metadata$timestamp_completed <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
run_metadata$output_files <- c(
  score_per_query = file.path(score_dir, "per_query_tie_aware.csv"),
  score_summary = file.path(score_dir, "summary_tie_aware.csv"),
  score_summary_across_replicates = file.path(
    score_dir, "summary_across_replicates_tie_aware.csv"
  ),
  regenerated_per_query = file.path(regenerated_dir, "per_query_tie_aware.csv"),
  regenerated_summary = file.path(regenerated_dir, "summary_tie_aware.csv"),
  regenerated_summary_across_replicates = file.path(
    regenerated_dir, "summary_across_replicates_tie_aware.csv"
  ),
  regeneration_qc = file.path(regenerated_dir, "representation_regeneration_qc.csv"),
  paired_per_query = file.path(comparison_dir, "paired_per_query.csv"),
  paired_summary = file.path(comparison_dir, "paired_summary.csv"),
  paired_summary_across_replicates = file.path(
    comparison_dir, "paired_summary_across_replicates.csv"
  ),
  reference_anchor_per_query = file.path(
    reference_anchor_dir, "per_query_tie_aware.csv"
  ),
  reference_anchor_summary = file.path(
    reference_anchor_dir, "summary_tie_aware.csv"
  ),
  reference_anchor_pairwise_invariance_qc = file.path(
    reference_anchor_dir, "pairwise_invariance_qc.csv"
  )
)
saveRDS(run_metadata, file.path(manifest_dir, "run_manifest.rds"))
run_completed <- TRUE
writeLines(
  paste("completed", run_metadata$timestamp_completed), status_file, useBytes = TRUE
)

message("Perturbation benchmark complete.")
message("Score-level output: ", score_dir)
message("Regenerated-representation output: ", regenerated_dir)
message("Paired comparison: ", comparison_dir)
message("Reference-anchor error output: ", reference_anchor_dir)
message("Manifest: ", manifest_dir)
