#!/usr/bin/env Rscript

# Publication-grade supplementary stress analyses for ppmWass.
#
# All query perturbations are applied to raw EI fragment peaks. The package
# preprocessing pipeline is then rerun to regenerate fragments, anchored and
# pairwise mass differences, and the pooled derived representation. In
# particular, this script never shifts or perturbs a cached pooled loss_list.
#
# Default tasks reproduce the legacy v4 grids where those existed and add:
#   - Sinkhorn epsilon sensitivity, and
#   - an absolute saturation-width sweep at fixed ppmWass base cost.
#
# Examples:
#   Rscript inst/scripts/run_supplementary_stress_publication.R --smoke
#   Rscript inst/scripts/run_supplementary_stress_publication.R \
#     --tasks=factorial,background,siloxane --replicates=3
#   Rscript inst/scripts/run_supplementary_stress_publication.R \
#     --resume --output-dir=/path/to/results --manifest-dir=/path/to/manifest

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  substring(hit[[1]], nchar(prefix) + 1L)
}

flag_value <- function(name, default = FALSE) {
  if (any(args == paste0("--", name))) return(TRUE)
  raw <- get_arg(name, if (default) "true" else "false")
  tolower(raw) %in% c("1", "true", "yes", "y", "on")
}

split_character <- function(name, default) {
  value <- trimws(strsplit(get_arg(name, default), ",", fixed = TRUE)[[1]])
  value[nzchar(value)]
}

split_numeric <- function(name, default) {
  value <- suppressWarnings(as.numeric(split_character(name, default)))
  if (!length(value) || any(!is.finite(value))) {
    stop("--", name, " must be a comma-separated finite numeric vector.")
  }
  value
}

single_integer <- function(name, default, lower = NULL) {
  value <- suppressWarnings(as.integer(get_arg(name, default)))
  if (length(value) != 1L || is.na(value)) {
    stop("--", name, " must be one integer.")
  }
  if (!is.null(lower) && value < lower) {
    stop("--", name, " must be >= ", lower, ".")
  }
  value
}

single_numeric <- function(name, default, lower = NULL,
                           lower_inclusive = TRUE) {
  value <- suppressWarnings(as.numeric(get_arg(name, default)))
  if (length(value) != 1L || !is.finite(value)) {
    stop("--", name, " must be one finite number.")
  }
  if (!is.null(lower)) {
    ok <- if (lower_inclusive) value >= lower else value > lower
    if (!ok) stop("--", name, " is outside its allowed range.")
  }
  value
}

normalize_existing <- function(path) {
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

normalize_future <- function(path) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(getwd(), path)
  missing_components <- character()
  ancestor <- path
  while (!file.exists(ancestor)) {
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) break
    missing_components <- c(basename(ancestor), missing_components)
    ancestor <- parent
  }
  ancestor <- normalizePath(ancestor, winslash = "/", mustWork = TRUE)
  if (!length(missing_components)) return(ancestor)
  do.call(file.path, c(list(ancestor), as.list(missing_components)))
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalize_existing(sub("^--file=", "", script_arg[[1]]))
} else {
  NA_character_
}

repo_dir <- normalize_existing(get_arg("repo", getwd()))
bundle_root <- normalize_existing(get_arg(
  "bundle-root", file.path(repo_dir, "..", "..")
))
spectra_rds <- get_arg("spectra-rds", "")
if (nzchar(spectra_rds)) spectra_rds <- normalize_existing(spectra_rds)
use_recetox_msp <- !nzchar(spectra_rds)
recetox_msp_candidate <- get_arg(
  "recetox-msp", file.path(bundle_root, "inputs", "RECETOX_merged.msp")
)
recetox_msp_exists <- file.exists(recetox_msp_candidate)
if (use_recetox_msp && !recetox_msp_exists) {
  stop("RECETOX MSP does not exist and --spectra-rds was not supplied: ",
       recetox_msp_candidate)
}
recetox_msp <- if (recetox_msp_exists) {
  normalize_existing(recetox_msp_candidate)
} else {
  normalize_future(recetox_msp_candidate)
}

output_dir <- normalize_future(get_arg(
  "output-dir",
  file.path(bundle_root, "results", "current", "supplementary_stress")
))
manifest_dir <- normalize_future(get_arg(
  "manifest-dir",
  file.path(bundle_root, "manifest", "supplementary_stress")
))
log_file <- normalize_future(get_arg(
  "log-file", file.path(output_dir, "run_supplementary_stress_publication.log")
))

all_tasks <- c(
  "factorial", "background", "siloxane", "breakdown",
  "tolerance", "sinkhorn", "transition_width"
)
tasks <- split_character("tasks", paste(all_tasks, collapse = ","))
task_aliases <- c(
  tier1 = "factorial", bg = "background", class = "breakdown",
  tol = "tolerance", transition = "transition_width"
)
tasks <- unname(ifelse(tasks %in% names(task_aliases), task_aliases[tasks], tasks))
if (!length(tasks) || anyDuplicated(tasks) || any(!tasks %in% all_tasks)) {
  stop("--tasks must contain unique values from: ", paste(all_tasks, collapse = ", "))
}

expected_methods <- c(
  "ppm_wasserstein", "entropy_weighted", "entropy_unweighted",
  "cosine", "weighted_cosine", "composite", "hellinger"
)
methods <- split_character("methods", paste(expected_methods, collapse = ","))
if (!length(methods) || anyDuplicated(methods) ||
    any(!methods %in% expected_methods)) {
  stop("--methods must contain unique supported publication methods.")
}

smoke <- flag_value("smoke", FALSE)
resume <- flag_value("resume", FALSE)
overwrite <- flag_value("overwrite", FALSE)
if (resume && overwrite) stop("--resume and --overwrite are mutually exclusive.")

n_replicates <- single_integer("replicates", "3", lower = 1L)
base_seed <- single_integer("seed", "42", lower = 1L)
n_cores <- single_integer("n-cores", "4", lower = 1L)
max_spectra <- single_integer("max-spectra", "0", lower = 0L)
smoke_spectra <- single_integer("smoke-spectra", "6", lower = 2L)
use_parallel <- flag_value("parallel", TRUE)
write_query_level <- flag_value("write-query-level", TRUE)
require_clean <- flag_value("require-clean", !smoke)

hard_match_tolerance_ppm <- single_numeric(
  "hard-match-tolerance-ppm", get_arg("tol-ppm", "15"),
  lower = 0, lower_inclusive = FALSE
)
ppmWass_base_cost_ppm <- single_numeric(
  "ppmWass-base-cost-ppm", as.character(hard_match_tolerance_ppm),
  lower = 0, lower_inclusive = FALSE
)
transition_multiplier <- single_numeric(
  "transition-mult", "3", lower = 0, lower_inclusive = FALSE
)
sinkhorn_epsilon <- single_numeric(
  "sinkhorn-epsilon", "0.05", lower = 0, lower_inclusive = FALSE
)
sinkhorn_iterations <- single_integer("sinkhorn-iterations", "50", lower = 1L)
ot_method <- get_arg("ot-method", "exact")
if (!ot_method %in% c("sinkhorn", "greenkhorn", "exact")) {
  stop("--ot-method must be sinkhorn, greenkhorn, or exact.")
}

mz_levels <- split_numeric("mz-levels", "0,10,30,100,300,500")
intensity_levels <- split_numeric("intensity-levels", "0,0.1,0.3")
dropout_levels <- split_numeric("dropout-levels", "0,0.1,0.3")
background_levels <- split_numeric("background-levels", "0,0.01,0.02,0.05,0.1")
background_peak_counts <- as.integer(split_numeric(
  "background-peak-counts", "100,500"
))
siloxane_levels <- split_numeric("siloxane-levels", "0,0.05,0.2,0.5")
tolerance_levels <- split_numeric("tolerance-levels", "5,10,15,20,30")
tolerance_mz_levels <- split_numeric("tolerance-mz-levels", "0,10,30,100")
sinkhorn_niter_levels <- as.integer(split_numeric(
  "sinkhorn-niter-levels", "30,50,100,200"
))
sinkhorn_epsilon_levels <- split_numeric(
  "sinkhorn-epsilon-levels", "0.02,0.05,0.1"
)
sinkhorn_transition_levels <- split_numeric(
  "sinkhorn-transition-levels", "2,3,5"
)
sinkhorn_mz_levels <- split_numeric("sinkhorn-mz-levels", "30,100")
absolute_width_levels <- split_numeric("absolute-width-levels", "30,45,75")
calibration_shift_levels <- split_numeric(
  "calibration-shift-levels", "0,10,20,30,45,60,100,150"
)
failure_mz_ppm <- single_numeric("failure-mz-ppm", "100", lower = 0)
failure_intensity_sd <- single_numeric("failure-intensity-sd", "0.1", lower = 0)
failure_examples_per_method <- single_integer(
  "failure-examples-per-method", "200", lower = 1L
)
include_exact_reference <- flag_value("include-exact-reference", TRUE)

check_unique_nonnegative <- function(x, name, upper = Inf, positive = FALSE) {
  if (anyDuplicated(x)) stop("--", name, " must not contain duplicates.")
  lower_ok <- if (positive) x > 0 else x >= 0
  if (any(!lower_ok) || any(x > upper)) {
    stop("--", name, " contains values outside its allowed range.")
  }
}

check_unique_nonnegative(mz_levels, "mz-levels")
check_unique_nonnegative(intensity_levels, "intensity-levels")
check_unique_nonnegative(dropout_levels, "dropout-levels", upper = 1)
check_unique_nonnegative(background_levels, "background-levels", upper = 1)
check_unique_nonnegative(background_peak_counts, "background-peak-counts", positive = TRUE)
check_unique_nonnegative(siloxane_levels, "siloxane-levels", upper = 1)
check_unique_nonnegative(tolerance_levels, "tolerance-levels", positive = TRUE)
check_unique_nonnegative(tolerance_mz_levels, "tolerance-mz-levels")
check_unique_nonnegative(sinkhorn_niter_levels, "sinkhorn-niter-levels", positive = TRUE)
check_unique_nonnegative(sinkhorn_epsilon_levels, "sinkhorn-epsilon-levels", positive = TRUE)
check_unique_nonnegative(
  sinkhorn_transition_levels, "sinkhorn-transition-levels", positive = TRUE
)
check_unique_nonnegative(sinkhorn_mz_levels, "sinkhorn-mz-levels")
check_unique_nonnegative(absolute_width_levels, "absolute-width-levels", positive = TRUE)
check_unique_nonnegative(calibration_shift_levels, "calibration-shift-levels")

if (smoke) {
  n_replicates <- 1L
  max_spectra <- if (max_spectra > 0L) {
    min(max_spectra, smoke_spectra)
  } else {
    smoke_spectra
  }
  use_parallel <- FALSE
  require_clean <- FALSE
}

git_lines <- function(arguments) {
  suppressWarnings(tryCatch(
    system2("git", c("-C", repo_dir, arguments), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  ))
}

git_commit_lines <- git_lines(c("rev-parse", "HEAD"))
git_commit <- if (length(git_commit_lines)) git_commit_lines[[1]] else NA_character_
git_status <- git_lines(c("status", "--porcelain"))
package_tree_dirty <- length(git_status) > 0L
if (require_clean && package_tree_dirty) {
  stop(
    "Publication mode requires a clean package tree. Commit the audited source ",
    "or pass --require-clean=false for a non-final diagnostic run."
  )
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to execute this driver from the development tree.")
}
pkgload::load_all(repo_dir, quiet = TRUE)

parse_ei_internal <- getFromNamespace("parse_ei", "ppmWass")
process_spectrum_internal <- getFromNamespace("process_single_spectrum", "ppmWass")
requested_approx_combined_internal <- getFromNamespace(
  "combined_distance_requested_approx", "ppmWass"
)
na_safe_group_key_internal <- getFromNamespace("na_safe_group_key", "ppmWass")
summarize_publication_metric_internal <- getFromNamespace(
  "summarize_publication_metric", "ppmWass"
)

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
params$tol_ppm <- hard_match_tolerance_ppm
params$wasserstein_transition_mult <- transition_multiplier
params$ot_method <- ot_method
params$sinkhorn_epsilon <- sinkhorn_epsilon
params$sinkhorn_niter <- sinkhorn_iterations
params$mass_power <- 3
params$intensity_power <- 0.5
params$w_frag <- 0.70
params$w_loss <- 0.30
params$backend <- "pair_loop"
params$use_parallel <- use_parallel
params$n_cores <- n_cores
params$class_detection_ppm <- 15

for (method in methods) {
  p <- params
  p$distance_method <- method
  validate_params(p)
}

format_id_number <- function(x) {
  out <- format(x, trim = TRUE, scientific = FALSE, digits = 12)
  out <- sub("(\\.[0-9]*?)0+$", "\\1", out)
  out <- sub("\\.$", "", out)
  out <- gsub("-", "m", out, fixed = TRUE)
  gsub(".", "p", out, fixed = TRUE)
}

add_seed_metadata <- function(df, task, task_code, paired_key) {
  df$task <- task
  df$condition_index <- seq_len(nrow(df))
  df$paired_seed_group <- paired_key
  seed_group_index <- match(paired_key, unique(paired_key))
  seed_value <- base_seed + task_code * 10000000 +
    seed_group_index * 1000 + df$replicate
  if (any(seed_value > .Machine$integer.max)) {
    stop("Derived condition seed exceeds R's integer seed range.")
  }
  df$condition_seed <- as.integer(seed_value)
  df
}

make_factorial_conditions <- function() {
  grid <- expand.grid(
    mz_noise_ppm = mz_levels,
    intensity_noise_sd = intensity_levels,
    dropout_prob = dropout_levels,
    replicate = seq_len(n_replicates),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  key <- paste(
    "mz", grid$mz_noise_ppm, "int", grid$intensity_noise_sd,
    "drop", grid$dropout_prob, sep = "_"
  )
  grid <- add_seed_metadata(grid, "factorial", 1L, key)
  grid$condition_id <- paste0(
    "mz", format_id_number(grid$mz_noise_ppm),
    "_int", format_id_number(grid$intensity_noise_sd),
    "_drop", format_id_number(grid$dropout_prob),
    "_rep", sprintf("%02d", grid$replicate)
  )
  if (smoke) {
    hit <- which(
      grid$mz_noise_ppm == 30 & grid$intensity_noise_sd == 0.1 &
        grid$dropout_prob == 0.1
    )
    grid <- grid[if (length(hit)) hit[[1]] else 1L, , drop = FALSE]
  }
  grid
}

make_background_conditions <- function() {
  grid <- expand.grid(
    bg_level_frac = background_levels,
    n_bg_peaks = background_peak_counts,
    replicate = seq_len(n_replicates),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid <- grid[
    grid$bg_level_frac != 0 |
      grid$n_bg_peaks == background_peak_counts[[1]],
    , drop = FALSE
  ]
  key <- paste("bg", grid$bg_level_frac, "n", grid$n_bg_peaks, sep = "_")
  grid <- add_seed_metadata(grid, "background", 2L, key)
  grid$condition_id <- paste0(
    "bg", format_id_number(grid$bg_level_frac),
    "_n", grid$n_bg_peaks,
    "_rep", sprintf("%02d", grid$replicate)
  )
  if (smoke) {
    hit <- which(grid$bg_level_frac == 0.02 & grid$n_bg_peaks == 100)
    grid <- grid[if (length(hit)) hit[[1]] else 1L, , drop = FALSE]
  }
  grid
}

make_siloxane_conditions <- function() {
  grid <- expand.grid(
    level_frac = siloxane_levels,
    replicate = seq_len(n_replicates),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  key <- paste("siloxane", grid$level_frac, sep = "_")
  grid <- add_seed_metadata(grid, "siloxane", 3L, key)
  grid$condition_id <- paste0(
    "level", format_id_number(grid$level_frac),
    "_rep", sprintf("%02d", grid$replicate)
  )
  if (smoke) {
    hit <- which(grid$level_frac == 0.2)
    grid <- grid[if (length(hit)) hit[[1]] else 1L, , drop = FALSE]
  }
  grid
}

make_breakdown_conditions <- function() {
  grid <- data.frame(
    mz_noise_ppm = c(0, 30, 100, 300),
    intensity_noise_sd = c(0, 0.1, 0.1, 0.1),
    dropout_prob = 0,
    stringsAsFactors = FALSE
  )
  grid <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    cbind(grid[rep(i, n_replicates), , drop = FALSE],
          replicate = seq_len(n_replicates))
  }))
  rownames(grid) <- NULL
  key <- paste(
    "mz", grid$mz_noise_ppm, "int", grid$intensity_noise_sd,
    "drop", grid$dropout_prob, sep = "_"
  )
  grid <- add_seed_metadata(grid, "breakdown", 4L, key)
  grid$condition_id <- paste0(
    "mz", format_id_number(grid$mz_noise_ppm),
    "_int", format_id_number(grid$intensity_noise_sd),
    "_drop", format_id_number(grid$dropout_prob),
    "_rep", sprintf("%02d", grid$replicate)
  )
  if (smoke) {
    hit <- which(
      grid$mz_noise_ppm == 100 & grid$intensity_noise_sd == 0.1
    )
    grid <- grid[if (length(hit)) hit[[1]] else 1L, , drop = FALSE]
  }
  grid
}

make_tolerance_conditions <- function() {
  rows <- list()
  index <- 0L
  for (mz in tolerance_mz_levels) {
    for (replicate in seq_len(n_replicates)) {
      for (tol in tolerance_levels) {
        index <- index + 1L
        rows[[index]] <- data.frame(
          tol_ppm = tol,
          mz_noise_ppm = mz,
          intensity_noise_sd = 0.1,
          dropout_prob = 0,
          replicate = replicate,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  grid <- do.call(rbind, rows)
  key <- paste(
    "mz", grid$mz_noise_ppm, "int", grid$intensity_noise_sd,
    "drop", grid$dropout_prob, sep = "_"
  )
  grid <- add_seed_metadata(grid, "tolerance", 5L, key)
  grid$condition_id <- paste0(
    "tol", format_id_number(grid$tol_ppm),
    "_mz", format_id_number(grid$mz_noise_ppm),
    "_rep", sprintf("%02d", grid$replicate)
  )
  if (smoke) {
    hit <- which(grid$tol_ppm == 15 & grid$mz_noise_ppm == 30)
    grid <- grid[if (length(hit)) hit[[1]] else 1L, , drop = FALSE]
  }
  grid
}

make_sinkhorn_conditions <- function() {
  rows <- list()
  index <- 0L
  for (mz in sinkhorn_mz_levels) {
    for (replicate in seq_len(n_replicates)) {
      for (niter in sinkhorn_niter_levels) {
        for (epsilon in sinkhorn_epsilon_levels) {
          for (mult in sinkhorn_transition_levels) {
            index <- index + 1L
            rows[[index]] <- data.frame(
              ot_method_condition = "sinkhorn",
              sinkhorn_niter = as.integer(niter),
              sinkhorn_epsilon = epsilon,
              transition_mult = mult,
              mz_noise_ppm = mz,
              intensity_noise_sd = 0.1,
              dropout_prob = 0,
              replicate = replicate,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  if (include_exact_reference) {
    index <- index + 1L
    rows[[index]] <- data.frame(
      ot_method_condition = "exact",
      sinkhorn_niter = sinkhorn_iterations,
      sinkhorn_epsilon = sinkhorn_epsilon,
      transition_mult = 3,
      mz_noise_ppm = 30,
      intensity_noise_sd = 0.1,
      dropout_prob = 0,
      replicate = 1L,
      stringsAsFactors = FALSE
    )
  }
  grid <- do.call(rbind, rows)
  key <- paste(
    "mz", grid$mz_noise_ppm, "int", grid$intensity_noise_sd,
    "drop", grid$dropout_prob, sep = "_"
  )
  grid <- add_seed_metadata(grid, "sinkhorn", 6L, key)
  sinkhorn_id <- paste0(
    "n", grid$sinkhorn_niter,
    "_eps", format_id_number(grid$sinkhorn_epsilon),
    "_tm", format_id_number(grid$transition_mult),
    "_mz", format_id_number(grid$mz_noise_ppm),
    "_rep", sprintf("%02d", grid$replicate)
  )
  grid$condition_id <- ifelse(
    grid$ot_method_condition == "exact",
    paste0("exact_reference_mz30_rep", sprintf("%02d", grid$replicate)),
    sinkhorn_id
  )
  if (smoke) {
    hit <- which(
      grid$ot_method_condition == "sinkhorn" &
      grid$sinkhorn_niter == 50 & grid$sinkhorn_epsilon == 0.05 &
        grid$transition_mult == 3 & grid$mz_noise_ppm == 30
    )
    requested_index <- if (length(hit)) hit[[1]] else
      which(grid$ot_method_condition == "sinkhorn")[[1L]]
    exact_indices <- which(grid$ot_method_condition == "exact")
    grid <- grid[unique(c(requested_index, exact_indices)), , drop = FALSE]
  }
  grid
}

make_transition_width_conditions <- function() {
  rows <- list()
  index <- 0L
  for (shift in calibration_shift_levels) {
    for (replicate in seq_len(n_replicates)) {
      for (width in absolute_width_levels) {
        index <- index + 1L
        rows[[index]] <- data.frame(
          saturation_width_ppm_requested = width,
          calibration_shift_ppm = shift,
          normalized_shift = shift / width,
          replicate = replicate,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  grid <- do.call(rbind, rows)
  key <- paste("shift", grid$calibration_shift_ppm, sep = "_")
  grid <- add_seed_metadata(grid, "transition_width", 7L, key)
  grid$condition_id <- paste0(
    "width", format_id_number(grid$saturation_width_ppm_requested),
    "_shift", format_id_number(grid$calibration_shift_ppm),
    "_rep", sprintf("%02d", grid$replicate)
  )
  if (smoke) {
    hit <- which(
      grid$saturation_width_ppm_requested == 45 &
        grid$calibration_shift_ppm == 45
    )
    grid <- grid[if (length(hit)) hit[[1]] else 1L, , drop = FALSE]
  }
  grid
}

condition_grids <- list(
  factorial = make_factorial_conditions(),
  background = make_background_conditions(),
  siloxane = make_siloxane_conditions(),
  breakdown = make_breakdown_conditions(),
  tolerance = make_tolerance_conditions(),
  sinkhorn = make_sinkhorn_conditions(),
  transition_width = make_transition_width_conditions()
)

task_specs <- list(
  factorial = list(
    methods = methods,
    analysis_cols = c("mz_noise_ppm", "intensity_noise_sd", "dropout_prob"),
    builder = "factorial"
  ),
  background = list(
    methods = methods,
    analysis_cols = c("bg_level_frac", "n_bg_peaks"),
    builder = "background"
  ),
  siloxane = list(
    methods = methods,
    analysis_cols = "level_frac",
    builder = "siloxane"
  ),
  breakdown = list(
    methods = methods,
    analysis_cols = c("mz_noise_ppm", "intensity_noise_sd", "dropout_prob"),
    builder = "factorial"
  ),
  tolerance = list(
    methods = methods,
    analysis_cols = c(
      "tol_ppm", "mz_noise_ppm", "intensity_noise_sd", "dropout_prob"
    ),
    builder = "factorial"
  ),
  sinkhorn = list(
    methods = "ppm_wasserstein",
    analysis_cols = c(
      "ot_method_condition", "sinkhorn_niter", "sinkhorn_epsilon",
      "transition_mult",
      "mz_noise_ppm", "intensity_noise_sd", "dropout_prob"
    ),
    builder = "factorial"
  ),
  transition_width = list(
    methods = "ppm_wasserstein",
    analysis_cols = c(
      "saturation_width_ppm_requested", "calibration_shift_ppm",
      "normalized_shift"
    ),
    builder = "calibration"
  )
)

if (any(vapply(condition_grids, function(x) anyDuplicated(x$condition_id),
               integer(1)) > 0L)) {
  stop("Internal error: duplicate condition IDs.")
}

flatten_parameter <- function(value) {
  if (is.null(value)) return("<NULL>")
  if (is.function(value)) return("<function>")
  if (is.atomic(value)) return(paste(as.character(value), collapse = ";"))
  paste(capture.output(str(value, give.attr = FALSE)), collapse = " ")
}

object_fingerprint <- function(object) {
  path <- tempfile("supplementary-stress-fingerprint-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 3)
  unname(tools::md5sum(path))
}

input_manifest <- data.frame(
  dataset = c(
    if (use_recetox_msp) "RECETOX_MSP" else NULL,
    if (nzchar(spectra_rds)) "PREPARED_SPECTRA_RDS" else NULL
  ),
  path = c(
    if (use_recetox_msp) recetox_msp else NULL,
    if (nzchar(spectra_rds)) spectra_rds else NULL
  ),
  md5 = c(
    if (use_recetox_msp) unname(tools::md5sum(recetox_msp)) else NULL,
    if (nzchar(spectra_rds)) unname(tools::md5sum(spectra_rds)) else NULL
  ),
  stringsAsFactors = FALSE
)

package_version_or_na <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}

runtime_versions <- data.frame(
  component = c(
    "R", "R_platform", "operating_system", "ppmWass", "pkgload",
    "approxOT", "transport", "msentropy"
  ),
  version = c(
    R.version.string,
    R.version$platform,
    paste(Sys.info()[c("sysname", "release", "machine")], collapse = ";"),
    package_version_or_na("ppmWass"),
    package_version_or_na("pkgload"),
    package_version_or_na("approxOT"),
    package_version_or_na("transport"),
    package_version_or_na("msentropy")
  ),
  stringsAsFactors = FALSE
)

configuration <- list(
  schema_version = 3L,
  script_path = script_path,
  script_md5 = if (!is.na(script_path)) unname(tools::md5sum(script_path)) else NA_character_,
  repo_dir = repo_dir,
  bundle_root = bundle_root,
  output_dir = output_dir,
  manifest_dir = manifest_dir,
  input_manifest = input_manifest,
  package_git_commit = git_commit,
  package_git_status = git_status,
  runtime_versions = runtime_versions,
  tasks = tasks,
  methods = methods,
  smoke = smoke,
  n_replicates = n_replicates,
  base_seed = base_seed,
  n_cores = n_cores,
  max_spectra = max_spectra,
  use_parallel = use_parallel,
  write_query_level = write_query_level,
  effective_parameters = params,
  condition_grids = condition_grids[tasks],
  failure_mz_ppm = failure_mz_ppm,
  failure_intensity_sd = failure_intensity_sd,
  failure_examples_per_method = failure_examples_per_method,
  include_exact_reference = include_exact_reference,
  named_sinkhorn_execution_policy =
    "raw_requested_only_no_retry_no_exact_fallback",
  named_sinkhorn_orientation = "query_to_library",
  named_sinkhorn_nonfinite_policy = paste(
    "for constant ground cost, use its coupling-independent analytic value",
    "while retaining the one raw requested solver attempt in the audit;",
    "otherwise retain invalid output as NA, report exhaustive counts and",
    "examples, and never retry, fall back, impute, or compute retrieval",
    "metrics for an incomplete row"
  )
)
configuration_fingerprint <- object_fingerprint(configuration)

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Atomic rename failed for: ", path)
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Atomic rename failed for: ", path)
  invisible(path)
}

atomic_write_lines <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(object, temporary, useBytes = TRUE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Atomic rename failed for: ", path)
  invisible(path)
}

directory_has_files <- function(path) {
  dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0L
}

guard_removal_target <- function(path) {
  resolved <- normalize_future(path)
  forbidden <- unique(vapply(
    c("/", path.expand("~"), repo_dir, bundle_root),
    normalize_future, character(1)
  ))
  if (resolved %in% forbidden || nchar(resolved) < 12L) {
    stop("Refusing destructive overwrite of unsafe path: ", resolved)
  }
  resolved
}

configuration_path <- file.path(manifest_dir, "run_configuration.rds")
if (overwrite) {
  for (path in unique(c(output_dir, manifest_dir))) {
    if (dir.exists(path)) unlink(guard_removal_target(path), recursive = TRUE)
  }
} else if (resume) {
  if (!file.exists(configuration_path)) {
    stop("--resume requested but run_configuration.rds is absent: ", configuration_path)
  }
  previous <- readRDS(configuration_path)
  if (!identical(previous$fingerprint, configuration_fingerprint)) {
    stop(
      "Resume configuration mismatch. Use the original arguments/source state ",
      "or start a new output directory."
    )
  }
} else if (directory_has_files(output_dir) || directory_has_files(manifest_dir)) {
  stop(
    "Refusing to overwrite an existing supplementary run. Use a new output ",
    "directory, --resume, or explicit --overwrite."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
if (!resume) {
  atomic_save_rds(
    list(fingerprint = configuration_fingerprint, configuration = configuration),
    configuration_path
  )
}

log_connection <- file(log_file, open = if (resume) "at" else "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
log_closed <- FALSE
close_run_log <- function() {
  if (log_closed) return(invisible(NULL))
  sink(type = "message")
  sink()
  close(log_connection)
  log_closed <<- TRUE
  invisible(NULL)
}
on.exit(close_run_log(), add = TRUE)

timestamp_started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
command_text <- paste(
  shQuote(c("Rscript", script_path, args), type = "sh"),
  collapse = " "
)
atomic_write_lines(command_text, file.path(manifest_dir, "command.txt"))
atomic_write_lines(git_commit, file.path(manifest_dir, "git_commit.txt"))
atomic_write_lines(
  if (length(git_status)) git_status else "<clean>",
  file.path(manifest_dir, "git_status.txt")
)
atomic_write_csv(input_manifest, file.path(manifest_dir, "input_checksums.csv"))
atomic_write_csv(runtime_versions, file.path(manifest_dir, "runtime_versions.csv"))
atomic_write_csv(
  data.frame(
    parameter = names(params),
    value = vapply(params, flatten_parameter, character(1)),
    stringsAsFactors = FALSE
  ),
  file.path(manifest_dir, "effective_parameters.csv")
)
for (task in tasks) {
  atomic_write_csv(
    condition_grids[[task]],
    file.path(manifest_dir, paste0("grid_", task, ".csv"))
  )
}

analysis_definitions <- data.frame(
  analysis = all_tasks,
  definition = c(
    "Legacy Tier-1 random per-peak m/z jitter x multiplicative intensity noise x dropout.",
    "Legacy continuum background: uniform m/z locations and capped exponential intensities.",
    "Legacy fixed eight-peak siloxane/GC contaminant series.",
    "Representative class, failure, and nearest-score correctness AUROC analysis.",
    "All methods across hard-match/base-cost settings with paired query perturbations.",
    paste(
      "Raw requested-only ppmWass Sinkhorn iterations x epsilon x transition",
      "multiplier with paired queries; one approxOT call per nonempty channel,",
      "no retry or exact fallback, explicit nonfinite/residual gates, plus an",
      "exact-OT reference."
    ),
    paste(
      "ppmWass absolute saturation width at fixed base cost; calibration shift is",
      "applied to raw peaks and the derived representation is regenerated."
    )
  ),
  stringsAsFactors = FALSE
)
atomic_write_csv(
  analysis_definitions[analysis_definitions$analysis %in% tasks, , drop = FALSE],
  file.path(manifest_dir, "analysis_definitions.csv")
)

message("Supplementary stress run fingerprint: ", configuration_fingerprint)
message("Tasks: ", paste(tasks, collapse = ", "))
message("Methods: ", paste(methods, collapse = ", "))
message("Replicates: ", n_replicates, " | base Sinkhorn iterations: ", sinkhorn_iterations)

subset_spectra_local <- function(spectra, idx) {
  old_n <- nrow(spectra$df_spec)
  out <- spectra
  out$df_spec <- spectra$df_spec[idx, , drop = FALSE]
  components <- c(
    "frag_list", "loss_list", "derived_list", "loss_typ_list",
    "loss_anchor_list", "loss_pair_list", "derived_anchor_list",
    "derived_pair_list", "loss_anchor_typ_list", "loss_pair_typ_list"
  )
  for (name in components) {
    value <- spectra[[name]]
    if (!is.null(value) && length(value) == old_n) out[[name]] <- value[idx]
  }
  for (name in c("mref_confidence", "ri")) {
    value <- spectra[[name]]
    if (!is.null(value) && length(value) == old_n) out[[name]] <- value[idx]
  }
  out
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
  required_metadata <- c("id", "ei", "RI", "inchikey")
  missing_metadata <- setdiff(required_metadata, names(spectra$df_spec))
  if (length(missing_metadata)) {
    stop("df_spec is missing: ", paste(missing_metadata, collapse = ", "))
  }
  df <- spectra$df_spec
  metadata_ok <- !is.na(df$RI) & !is.na(df$inchikey) & nchar(df$inchikey) >= 14L
  fragment_ok <- vapply(spectra$frag_list, function(fragment) {
    is.matrix(fragment) && ncol(fragment) >= 2L && nrow(fragment) >= 2L &&
      all(is.finite(fragment[, 1:2, drop = FALSE])) &&
      sum(fragment[, 2], na.rm = TRUE) > 0
  }, logical(1))
  keep <- which(metadata_ok & fragment_ok)
  message(
    "RECETOX filter: input=", nrow(df), " retained=", length(keep),
    " removed=", nrow(df) - length(keep)
  )
  spectra <- subset_spectra_local(spectra, keep)
  if (max_spectra > 0L && nrow(spectra$df_spec) > max_spectra) {
    spectra <- subset_spectra_local(spectra, seq_len(max_spectra))
    message("Deterministic max-spectra subset: ", nrow(spectra$df_spec))
  }
  ids <- as.character(spectra$df_spec$id)
  if (length(ids) < 2L || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Filtered RECETOX IDs must be at least two unique nonmissing strings.")
  }
  names(spectra$frag_list) <- ids
  names(spectra$loss_list) <- ids
  spectra$derived_list <- spectra$loss_list

  if (nzchar(spectra_rds)) {
    compatibility <- do.call(rbind, lapply(seq_along(ids), function(i) {
      rebuilt <- process_spectrum_internal(
        as.character(spectra$df_spec$ei[[i]]), params, deriv_mode = "off"
      )
      data.frame(
        query_id = ids[[i]],
        fragment_bit_identical = identical(rebuilt$frag, spectra$frag_list[[i]]),
        derived_bit_identical = identical(rebuilt$loss, spectra$loss_list[[i]]),
        stringsAsFactors = FALSE
      )
    }))
    atomic_write_csv(
      compatibility,
      file.path(manifest_dir, "prepared_spectra_compatibility.csv")
    )
    if (!all(compatibility$fragment_bit_identical) ||
        !all(compatibility$derived_bit_identical)) {
      stop(
        "--spectra-rds is incompatible with the current preprocessing code ",
        "or effective parameters; see prepared_spectra_compatibility.csv."
      )
    }
  }

  if ("breakdown" %in% tasks) {
    if (!"compound_class" %in% names(spectra$df_spec)) {
      stop("The breakdown task requires df_spec$compound_class.")
    }
    observed_class <- as.character(spectra$df_spec$compound_class)
    informative <- !is.na(observed_class) & nzchar(observed_class) &
      !tolower(observed_class) %in% c("unknown", "unclassified")
    if (!any(informative)) {
      stop("The breakdown task requires at least one informative compound class.")
    }
  }
  spectra
}

spectra <- load_recetox()
ids <- as.character(spectra$df_spec$id)
raw_ei <- stats::setNames(as.character(spectra$df_spec$ei), ids)
library_frag <- spectra$frag_list
library_derived <- spectra$loss_list
n_queries <- length(ids)
compound_class <- if ("compound_class" %in% names(spectra$df_spec)) {
  as.character(spectra$df_spec$compound_class)
} else {
  rep("unknown", n_queries)
}
compound_class[is.na(compound_class) | !nzchar(compound_class)] <- "unknown"
names(compound_class) <- ids

empty_spectrum <- function() {
  matrix(numeric(), ncol = 2L, dimnames = list(NULL, c("mz", "intensity")))
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

serialize_peak_string <- function(spectrum) {
  spectrum <- coerce_peak_matrix(spectrum)
  keep <- is.finite(spectrum[, 1]) & is.finite(spectrum[, 2]) &
    spectrum[, 1] > 0 & spectrum[, 2] > 0
  spectrum <- spectrum[keep, , drop = FALSE]
  if (!nrow(spectrum)) return("")
  paste(sprintf("%.17g:%.17g", spectrum[, 1], spectrum[, 2]), collapse = " ")
}

SILOXANE_MZ_DEFAULT <- c(
  73.0473, 147.0654, 207.0327, 221.0848,
  281.0518, 295.1075, 355.0708, 429.0884
)

transform_raw <- function(raw, condition, builder_kind, query_seed) {
  raw <- coerce_peak_matrix(raw)
  out <- raw
  realized_shift <- numeric()
  injected_count <- 0L
  calibration_sign <- NA_integer_

  if (builder_kind == "factorial") {
    set.seed(query_seed)
    if (condition$mz_noise_ppm > 0 && nrow(out)) {
      realized_shift <- stats::rnorm(nrow(out), 0, condition$mz_noise_ppm)
      out[, 1] <- out[, 1] * (1 + realized_shift * 1e-6)
    }
    if (condition$intensity_noise_sd > 0 && nrow(out)) {
      out[, 2] <- out[, 2] *
        (1 + stats::rnorm(nrow(out), 0, condition$intensity_noise_sd))
      out[, 2] <- pmax(out[, 2], 0)
    }
    if (condition$dropout_prob > 0 && nrow(out)) {
      keep <- stats::runif(nrow(out)) > condition$dropout_prob
      out <- out[keep, , drop = FALSE]
      if (length(realized_shift)) realized_shift <- realized_shift[keep]
    }
  } else if (builder_kind == "background") {
    set.seed(query_seed)
    if (condition$bg_level_frac > 0 && condition$n_bg_peaks > 0 && nrow(out)) {
      max_background <- condition$bg_level_frac * max(out[, 2], na.rm = TRUE)
      bg_mz <- stats::runif(condition$n_bg_peaks, params$min_mz, params$max_mz)
      bg_int <- pmin(
        stats::rexp(condition$n_bg_peaks, rate = 3 / max_background),
        max_background
      )
      out <- rbind(out, cbind(mz = bg_mz, intensity = bg_int))
      injected_count <- condition$n_bg_peaks
    }
  } else if (builder_kind == "siloxane") {
    set.seed(query_seed)
    if (condition$level_frac > 0 && nrow(out)) {
      in_range <- SILOXANE_MZ_DEFAULT >= min(out[, 1]) - 1 &
        SILOXANE_MZ_DEFAULT <= max(out[, 1]) + 1
      injection_mz <- SILOXANE_MZ_DEFAULT[in_range]
      if (length(injection_mz)) {
        injection_intensity <- rep(
          condition$level_frac * max(out[, 2], na.rm = TRUE),
          length(injection_mz)
        )
        out <- rbind(
          out, cbind(mz = injection_mz, intensity = injection_intensity)
        )
        injected_count <- length(injection_mz)
      }
    }
  } else if (builder_kind == "calibration") {
    set.seed(query_seed)
    if (condition$calibration_shift_ppm > 0 && nrow(out)) {
      calibration_sign <- sample(c(-1L, 1L), 1L)
      realized_shift <- rep(
        calibration_sign * condition$calibration_shift_ppm, nrow(out)
      )
      out[, 1] <- out[, 1] * (1 + realized_shift * 1e-6)
    }
  } else {
    stop("Unknown query builder: ", builder_kind)
  }

  keep <- is.finite(out[, 1]) & is.finite(out[, 2]) &
    out[, 1] > 0 & out[, 2] > 0
  out <- out[keep, , drop = FALSE]
  if (nrow(out)) out <- out[order(out[, 1]), , drop = FALSE]
  finite_shift <- realized_shift[is.finite(realized_shift)]
  list(
    spectrum = out,
    shift_mean_ppm = if (length(finite_shift)) mean(finite_shift) else NA_real_,
    shift_sd_ppm = if (length(finite_shift) > 1L) stats::sd(finite_shift) else
      if (length(finite_shift)) 0 else NA_real_,
    shift_min_ppm = if (length(finite_shift)) min(finite_shift) else NA_real_,
    shift_max_ppm = if (length(finite_shift)) max(finite_shift) else NA_real_,
    injected_peak_count = injected_count,
    calibration_sign = calibration_sign
  )
}

build_regenerated_query <- function(condition, builder_kind) {
  fragment <- vector("list", n_queries)
  derived <- vector("list", n_queries)
  qc <- vector("list", n_queries)
  query_seeds <- as.integer(condition$condition_seed + seq_len(n_queries) * 1009L)
  if (any(query_seeds > .Machine$integer.max)) {
    stop("Query seed exceeds R's integer seed range.")
  }

  for (i in seq_len(n_queries)) {
    raw <- parse_ei_internal(raw_ei[[i]])
    transformed <- transform_raw(raw, condition, builder_kind, query_seeds[[i]])
    rebuilt <- process_spectrum_internal(
      serialize_peak_string(transformed$spectrum),
      params,
      deriv_mode = "off"
    )
    fragment[[i]] <- coerce_peak_matrix(rebuilt$frag)
    derived[[i]] <- coerce_peak_matrix(rebuilt$loss)
    if (any(!is.finite(fragment[[i]])) || any(!is.finite(derived[[i]]))) {
      stop("Regenerated representation is nonfinite for query: ", ids[[i]])
    }
    qc[[i]] <- data.frame(
      query_index = i,
      query_id = ids[[i]],
      query_seed = query_seeds[[i]],
      raw_peak_count = nrow(raw),
      transformed_raw_peak_count = nrow(transformed$spectrum),
      processed_fragment_count = nrow(fragment[[i]]),
      anchored_difference_count = nrow(rebuilt$loss_anchor),
      pairwise_difference_count = nrow(rebuilt$loss_pair),
      pooled_derived_count = nrow(derived[[i]]),
      estimated_high_mass_reference = rebuilt$mref,
      mref_confidence = rebuilt$mref_confidence,
      realized_shift_mean_ppm = transformed$shift_mean_ppm,
      realized_shift_sd_ppm = transformed$shift_sd_ppm,
      realized_shift_min_ppm = transformed$shift_min_ppm,
      realized_shift_max_ppm = transformed$shift_max_ppm,
      injected_peak_count = transformed$injected_peak_count,
      calibration_sign = transformed$calibration_sign,
      fragment_bit_identical_to_library = identical(
        fragment[[i]], library_frag[[i]]
      ),
      derived_bit_identical_to_library = identical(
        derived[[i]], library_derived[[i]]
      ),
      stringsAsFactors = FALSE
    )
  }
  names(fragment) <- ids
  names(derived) <- ids
  clean_condition <- switch(
    builder_kind,
    factorial = condition$mz_noise_ppm == 0 &&
      condition$intensity_noise_sd == 0 && condition$dropout_prob == 0,
    background = condition$bg_level_frac == 0,
    siloxane = condition$level_frac == 0,
    calibration = condition$calibration_shift_ppm == 0,
    FALSE
  )
  qc_table <- do.call(rbind, qc)
  if (isTRUE(clean_condition) &&
      (!all(qc_table$fragment_bit_identical_to_library) ||
       !all(qc_table$derived_bit_identical_to_library))) {
    stop(
      "Clean raw-peak regeneration does not reproduce the loaded library ",
      "representations exactly; prepared spectra or preprocessing settings ",
      "are incompatible."
    )
  }
  list(
    frag_list = fragment,
    derived_list = derived,
    query_seeds = query_seeds,
    qc = qc_table
  )
}

settings_for <- function(task, condition, method) {
  p <- params
  p$distance_method <- method
  if (task == "tolerance") {
    p$tol_ppm <- condition$tol_ppm
  } else if (task == "sinkhorn") {
    p$tol_ppm <- ppmWass_base_cost_ppm
    p$ot_method <- condition$ot_method_condition
    p$sinkhorn_niter <- as.integer(condition$sinkhorn_niter)
    p$sinkhorn_epsilon <- condition$sinkhorn_epsilon
    p$wasserstein_transition_mult <- condition$transition_mult
  } else if (task == "transition_width") {
    p$tol_ppm <- ppmWass_base_cost_ppm
    p$ot_method <- ot_method
    p$sinkhorn_niter <- sinkhorn_iterations
    p$sinkhorn_epsilon <- sinkhorn_epsilon
    p$wasserstein_transition_mult <-
      condition$saturation_width_ppm_requested / ppmWass_base_cost_ppm
  } else {
    p$tol_ppm <- if (method == "ppm_wasserstein") {
      ppmWass_base_cost_ppm
    } else {
      hard_match_tolerance_ppm
    }
  }
  p <- validate_params(p)
  is_ppm_wasserstein <- identical(method, "ppm_wasserstein")
  raw_requested_only <- identical(task, "sinkhorn") &&
    is_ppm_wasserstein &&
    identical(as.character(condition$ot_method_condition), "sinkhorn")
  list(
    params = p,
    raw_requested_only = raw_requested_only,
    row = data.frame(
      hard_match_tolerance_ppm = if (method == "ppm_wasserstein") {
        hard_match_tolerance_ppm
      } else {
        p$tol_ppm
      },
      ppmWass_base_cost_ppm = if (method == "ppm_wasserstein") {
        p$tol_ppm
      } else {
        ppmWass_base_cost_ppm
      },
      transition_multiplier = p$wasserstein_transition_mult,
      saturation_width_ppm = if (method == "ppm_wasserstein") {
        p$tol_ppm * p$wasserstein_transition_mult
      } else {
        ppmWass_base_cost_ppm * transition_multiplier
      },
      sinkhorn_epsilon_effective = if (!is_ppm_wasserstein ||
        identical(p$ot_method, "exact")) NA_real_ else p$sinkhorn_epsilon,
      sinkhorn_iterations_effective = if (!is_ppm_wasserstein ||
        identical(p$ot_method, "exact")) NA_integer_ else p$sinkhorn_niter,
      solver_backend = if (is_ppm_wasserstein) p$ot_method else
        "not_applicable",
      ot_estimand = if (!is_ppm_wasserstein) {
        "not_applicable"
      } else if (identical(p$ot_method, "exact")) {
        "unregularized_exact_transport_cost"
      } else if (raw_requested_only) {
        "entropic_finite_iteration_raw_requested_only"
      } else {
        "entropic_approximate_transport_with_validated_fallback_policy"
      },
      solver_execution_policy = if (!is_ppm_wasserstein) {
        "not_applicable"
      } else if (raw_requested_only) {
        "raw_requested_only_no_retry_no_exact_fallback"
      } else if (identical(p$ot_method, "exact")) {
        "exact_primary_no_approximate_fallback"
      } else {
        "validated_retry_then_exact_fallback"
      },
      solver_orientation = if (is_ppm_wasserstein) "query_to_library" else
        "not_applicable",
      selected_path_contract = if (!is_ppm_wasserstein) {
        "not_applicable"
      } else if (raw_requested_only) {
        "approx_requested_only_or_explicit_analytic_input_contract"
      } else if (identical(p$ot_method, "exact")) {
        "exact_primary_or_explicit_analytic_input_contract"
      } else {
        "validated_backend_path"
      },
      distance_matrix_backend = p$backend,
      stringsAsFactors = FALSE
    )
  )
}

write_nonfinite_diagnostic <- function(distance_matrix, query, task,
                                       condition, method) {
  bad <- which(!is.finite(distance_matrix), arr.ind = TRUE)
  diagnostic <- data.frame(
    task = task,
    condition_id = condition$condition_id,
    method = method,
    query_id = rownames(distance_matrix)[bad[, 1]],
    library_id = colnames(distance_matrix)[bad[, 2]],
    value = distance_matrix[bad],
    query_fragment_count = vapply(
      query$frag_list[bad[, 1]], nrow, integer(1)
    ),
    query_derived_count = vapply(
      query$derived_list[bad[, 1]], nrow, integer(1)
    ),
    library_fragment_count = vapply(
      library_frag[bad[, 2]], nrow, integer(1)
    ),
    library_derived_count = vapply(
      library_derived[bad[, 2]], nrow, integer(1)
    ),
    stringsAsFactors = FALSE
  )
  path <- file.path(
    output_dir, "nonfinite_diagnostics",
    paste0(task, "__", condition$condition_id, "__", method, ".csv")
  )
  atomic_write_csv(diagnostic, path)
  path
}

validate_distance_matrix <- function(distance_matrix, query, task,
                                     condition, method,
                                     allow_nonfinite = FALSE) {
  if (!is.matrix(distance_matrix) ||
      !identical(dim(distance_matrix), c(n_queries, n_queries))) {
    stop("Unexpected distance-matrix dimensions for ", task, "/", condition$condition_id)
  }
  if (!identical(rownames(distance_matrix), ids) ||
      !identical(colnames(distance_matrix), ids)) {
    stop("Distance-matrix IDs are not aligned for ", task, "/", condition$condition_id)
  }
  if (any(!is.finite(distance_matrix))) {
    if (isTRUE(allow_nonfinite)) return(invisible(FALSE))
    diagnostic <- write_nonfinite_diagnostic(
      distance_matrix, query, task, condition, method
    )
    stop(
      "Strict nonfinite stop for ", task, "/", condition$condition_id, "/",
      method, ". Diagnostic: ", diagnostic
    )
  }
  invisible(TRUE)
}

requested_only_count_names <- c(
  "n_pair_distances", "n_channel_records", "n_solver_eligible_channel_records",
  "n_solver_attempts", "n_requested_setting_attempts",
  "n_unexpected_setting_attempts", "n_retry_attempts", "n_exact_attempts",
  "n_records_with_multiple_attempts", "n_record_contract_mismatches",
  "n_fallback_records",
  "n_selected_requested_only", "n_selected_requested_only_invalid",
  "n_selected_analytic_constant_raw_diagnostic",
  "n_analytic_input_contracts", "n_empty_channel_analytic_contracts",
  "n_unexpected_selected_paths",
  "n_nonfinite_channel_distances", "n_nonfinite_combined_distances",
  "n_raw_invalid_channel_records", "n_raw_plan_validation_failures",
  "n_raw_plan_nonfinite_attempts", "n_raw_plan_nonfinite_mass",
  "n_marginal_residual_exceeds_tolerance"
)

empty_requested_only_counts <- function() {
  # A full production scan can accumulate more than 2^31-1 plan cells across
  # all query-library/channel attempts. Keep counters as doubles to avoid
  # integer overflow while retaining exact integer representation at this scale.
  stats::setNames(numeric(length(requested_only_count_names)),
                  requested_only_count_names)
}

finite_max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

compute_requested_only_matrix <- function(query, effective, condition) {
  requested_niter <- as.integer(effective$params$sinkhorn_niter)
  requested_epsilon <- as.numeric(effective$params$sinkhorn_epsilon)
  requested_method <- as.character(effective$params$ot_method)

  compute_row <- function(i) {
    distances <- rep(NA_real_, n_queries)
    counts <- empty_requested_only_counts()
    selected_paths <- character()
    validation_errors <- character()
    nonfinite_examples <- list()
    raw_invalid_examples <- list()
    residual_extremes <- NULL
    nonfinite_index <- 0L
    raw_invalid_index <- 0L
    current_j <- NA_integer_

    collect_record <- function(record) {
      counts[["n_channel_records"]] <<-
        counts[["n_channel_records"]] + 1L
      path <- as.character(record$selected_path)
      selected_paths <<- unique(c(selected_paths, path))
      record_contract_match <- isTRUE(record$requested_only) &&
        identical(as.character(record$requested_method), requested_method) &&
        identical(as.integer(record$requested_niter), requested_niter) &&
        isTRUE(all.equal(
          as.numeric(record$requested_epsilon), requested_epsilon,
          tolerance = 0
        )) &&
        identical(as.character(record$context$task), "sinkhorn") &&
        identical(
          as.character(record$context$condition_id),
          as.character(condition$condition_id)
        ) &&
        identical(as.character(record$context$orientation), "query_to_library") &&
        identical(as.integer(record$context$query_index), as.integer(i)) &&
        identical(as.character(record$context$query_id), ids[[i]]) &&
        identical(
          as.integer(record$context$library_index), as.integer(current_j)
        ) &&
        identical(as.character(record$context$library_id), ids[[current_j]]) &&
        as.character(record$context$channel) %in% c("fragment", "derived")
      if (!record_contract_match) {
        counts[["n_record_contract_mismatches"]] <<-
          counts[["n_record_contract_mismatches"]] + 1L
      }
      if (identical(path, "approx_requested_only")) {
        counts[["n_selected_requested_only"]] <<-
          counts[["n_selected_requested_only"]] + 1L
        counts[["n_solver_eligible_channel_records"]] <<-
          counts[["n_solver_eligible_channel_records"]] + 1L
      } else if (identical(path, "approx_requested_only_invalid")) {
        counts[["n_selected_requested_only_invalid"]] <<-
          counts[["n_selected_requested_only_invalid"]] + 1L
        counts[["n_solver_eligible_channel_records"]] <<-
          counts[["n_solver_eligible_channel_records"]] + 1L
      } else if (identical(
        path, "analytic_constant_cost_with_raw_requested_diagnostic"
      )) {
        counts[["n_selected_analytic_constant_raw_diagnostic"]] <<-
          counts[["n_selected_analytic_constant_raw_diagnostic"]] + 1L
        counts[["n_solver_eligible_channel_records"]] <<-
          counts[["n_solver_eligible_channel_records"]] + 1L
        counts[["n_analytic_input_contracts"]] <<-
          counts[["n_analytic_input_contracts"]] + 1L
      } else if (path %in% c(
        "empty_spectrum", "empty_aligned_spectrum", "empty_positive_marginal"
      )) {
        counts[["n_analytic_input_contracts"]] <<-
          counts[["n_analytic_input_contracts"]] + 1L
        counts[["n_empty_channel_analytic_contracts"]] <<-
          counts[["n_empty_channel_analytic_contracts"]] + 1L
      } else {
        counts[["n_unexpected_selected_paths"]] <<-
          counts[["n_unexpected_selected_paths"]] + 1L
      }
      if (isTRUE(record$fallback_used)) {
        counts[["n_fallback_records"]] <<-
          counts[["n_fallback_records"]] + 1L
      }
      if (!is.finite(record$selected_total)) {
        counts[["n_nonfinite_channel_distances"]] <<-
          counts[["n_nonfinite_channel_distances"]] + 1L
      }
      if (record$status %in% c(
        "raw_invalid", "analytic_constant_cost_raw_diagnostic_invalid"
      )) {
        counts[["n_raw_invalid_channel_records"]] <<-
          counts[["n_raw_invalid_channel_records"]] + 1L
      }

      attempts <- record$attempts
      if (is.null(attempts)) attempts <- data.frame()
      counts[["n_solver_attempts"]] <<-
        counts[["n_solver_attempts"]] + nrow(attempts)
      if (nrow(attempts) > 1L) {
        counts[["n_records_with_multiple_attempts"]] <<-
          counts[["n_records_with_multiple_attempts"]] + 1L
      }
      if (!nrow(attempts)) return(invisible(NULL))

      setting_match <- attempts$step == "approx_requested_only" &
        attempts$backend == paste0("approxOT_", requested_method) &
        attempts$niter == requested_niter &
        is.finite(attempts$epsilon) &
        abs(attempts$epsilon - requested_epsilon) <=
          .Machine$double.eps * max(1, abs(requested_epsilon))
      counts[["n_requested_setting_attempts"]] <<-
        counts[["n_requested_setting_attempts"]] + sum(setting_match)
      counts[["n_unexpected_setting_attempts"]] <<-
        counts[["n_unexpected_setting_attempts"]] + sum(!setting_match)
      counts[["n_retry_attempts"]] <<- counts[["n_retry_attempts"]] +
        sum(grepl("retry", attempts$step, fixed = TRUE))
      counts[["n_exact_attempts"]] <<- counts[["n_exact_attempts"]] +
        sum(grepl("transport_exact", attempts$backend, fixed = TRUE) |
            grepl("exact", attempts$step, fixed = TRUE))
      counts[["n_raw_plan_validation_failures"]] <<-
        counts[["n_raw_plan_validation_failures"]] +
        sum(!attempts$accepted | is.na(attempts$accepted))
      nonfinite_plan <-
        (!is.na(attempts$nonfinite_mass) & attempts$nonfinite_mass > 0L) |
        (!is.na(attempts$validation_error) &
         attempts$validation_error %in% c(
           "nonfinite_plan_mass", "nonfinite_total_cost"
         ))
      counts[["n_raw_plan_nonfinite_attempts"]] <<-
        counts[["n_raw_plan_nonfinite_attempts"]] + sum(nonfinite_plan)
      counts[["n_raw_plan_nonfinite_mass"]] <<-
        counts[["n_raw_plan_nonfinite_mass"]] +
        sum(attempts$nonfinite_mass, na.rm = TRUE)
      raw_invalid_attempt <- !attempts$accepted | is.na(attempts$accepted) |
        (!is.na(attempts$solver_error) & nzchar(attempts$solver_error)) |
        nonfinite_plan
      if (any(raw_invalid_attempt) && raw_invalid_index < 5L) {
        attempt_index <- which(raw_invalid_attempt)[[1L]]
        raw_invalid_index <<- raw_invalid_index + 1L
        raw_invalid_examples[[raw_invalid_index]] <<- data.frame(
          task = "sinkhorn", condition_id = condition$condition_id,
          query_id = ids[[i]], library_id = ids[[current_j]],
          channel = as.character(record$context$channel),
          orientation = "query_to_library",
          selected_path = path, selected_distance = record$selected_total,
          analytic_identity_used = identical(
            path, "analytic_constant_cost_with_raw_requested_diagnostic"
          ),
          record_status = record$status,
          raw_step = attempts$step[[attempt_index]],
          raw_backend = attempts$backend[[attempt_index]],
          sinkhorn_epsilon_requested = requested_epsilon,
          sinkhorn_iterations_requested = requested_niter,
          sinkhorn_epsilon_executed = attempts$epsilon[[attempt_index]],
          sinkhorn_iterations_executed = attempts$niter[[attempt_index]],
          raw_total_cost = attempts$total_cost[[attempt_index]],
          raw_accepted_at_marginal_tolerance =
            attempts$accepted[[attempt_index]],
          raw_validation_error = attempts$validation_error[[attempt_index]],
          raw_solver_error = attempts$solver_error[[attempt_index]],
          raw_nonfinite_plan_mass = attempts$nonfinite_mass[[attempt_index]],
          raw_row_residual_linf = attempts$row_residual_linf[[attempt_index]],
          raw_col_residual_linf = attempts$col_residual_linf[[attempt_index]],
          stringsAsFactors = FALSE
        )
      }
      residual_failure <- !is.na(attempts$validation_error) &
        attempts$validation_error == "marginal_residual_exceeds_tolerance"
      counts[["n_marginal_residual_exceeds_tolerance"]] <<-
        counts[["n_marginal_residual_exceeds_tolerance"]] +
        sum(residual_failure)
      observed_errors <- attempts$validation_error[
        !is.na(attempts$validation_error) & nzchar(attempts$validation_error)
      ]
      validation_errors <<- unique(c(validation_errors, observed_errors))

      candidate <- data.frame(
        task = "sinkhorn", condition_id = condition$condition_id,
        query_id = ids[[i]], library_id = ids[[current_j]],
        channel = as.character(record$context$channel),
        orientation = "query_to_library",
        selected_path = path,
        solver_backend = attempts$backend[[1]],
        sinkhorn_epsilon_requested = requested_epsilon,
        sinkhorn_iterations_requested = requested_niter,
        sinkhorn_epsilon_executed = attempts$epsilon[[1]],
        sinkhorn_iterations_executed = attempts$niter[[1]],
        total_cost = attempts$total_cost[[1]],
        accepted_at_marginal_tolerance = attempts$accepted[[1]],
        validation_error = attempts$validation_error[[1]],
        row_residual_linf = attempts$row_residual_linf[[1]],
        col_residual_linf = attempts$col_residual_linf[[1]],
        max_marginal_residual_linf = finite_max_or_na(c(
          attempts$row_residual_linf[[1]], attempts$col_residual_linf[[1]]
        )),
        nonfinite_plan_mass = attempts$nonfinite_mass[[1]],
        stringsAsFactors = FALSE
      )
      residual_extremes <<- if (is.null(residual_extremes)) {
        candidate
      } else {
        rbind(residual_extremes, candidate)
      }
      residual_extremes <<- head(
        residual_extremes[
          order(residual_extremes$max_marginal_residual_linf,
                decreasing = TRUE, na.last = TRUE), , drop = FALSE
        ],
        5L
      )
      invisible(NULL)
    }

    for (j in seq_len(n_queries)) {
      current_j <- j
      distance <- requested_approx_combined_internal(
        query$frag_list[[i]], library_frag[[j]],
        query$derived_list[[i]], library_derived[[j]],
        params = effective$params,
        diagnostics = collect_record,
        context = list(
          task = "sinkhorn", condition_id = condition$condition_id,
          orientation = "query_to_library", query_index = i,
          query_id = ids[[i]], library_index = j, library_id = ids[[j]]
        )
      )
      distances[[j]] <- distance
      counts[["n_pair_distances"]] <- counts[["n_pair_distances"]] + 1L

      if (!is.finite(distance)) {
        counts[["n_nonfinite_combined_distances"]] <-
          counts[["n_nonfinite_combined_distances"]] + 1L
        if (nonfinite_index < 5L) {
          nonfinite_index <- nonfinite_index + 1L
          nonfinite_examples[[nonfinite_index]] <- data.frame(
            task = "sinkhorn", condition_id = condition$condition_id,
            orientation = "query_to_library",
            query_id = ids[[i]], library_id = ids[[j]],
            value = distance,
            query_fragment_count = nrow(query$frag_list[[i]]),
            query_derived_count = nrow(query$derived_list[[i]]),
            library_fragment_count = nrow(library_frag[[j]]),
            library_derived_count = nrow(library_derived[[j]]),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    residual_frame <- if (is.null(residual_extremes)) data.frame() else
      residual_extremes
    list(
      i = i, distances = distances, counts = counts,
      selected_paths = unique(selected_paths),
      validation_errors = unique(validation_errors),
      nonfinite_examples = if (length(nonfinite_examples)) {
        do.call(rbind, nonfinite_examples)
      } else data.frame(),
      raw_invalid_examples = if (length(raw_invalid_examples)) {
        do.call(rbind, raw_invalid_examples)
      } else data.frame(),
      residual_extremes = residual_frame
    )
  }

  indices <- seq_len(n_queries)
  if (isTRUE(effective$params$use_parallel) &&
      effective$params$n_cores > 1L && .Platform$OS.type != "windows") {
    rows <- parallel::mclapply(
      indices, compute_row,
      mc.cores = min(effective$params$n_cores, length(indices))
    )
  } else {
    rows <- lapply(indices, compute_row)
  }
  distance_matrix <- matrix(
    NA_real_, n_queries, n_queries, dimnames = list(ids, ids)
  )
  for (row in rows) distance_matrix[row$i, ] <- row$distances

  counts <- Reduce(`+`, lapply(rows, `[[`, "counts"))
  selected_paths <- sort(unique(unlist(lapply(rows, `[[`, "selected_paths"))))
  validation_errors <- sort(unique(unlist(lapply(rows, `[[`, "validation_errors"))))
  residual_frames <- lapply(rows, `[[`, "residual_extremes")
  residual_frames <- residual_frames[vapply(residual_frames, nrow, integer(1)) > 0L]
  residual_extremes <- if (length(residual_frames)) {
    frame <- do.call(rbind, residual_frames)
    frame <- frame[
      order(frame$max_marginal_residual_linf, decreasing = TRUE,
            na.last = TRUE), , drop = FALSE
    ]
    rownames(frame) <- NULL
    head(frame, 100L)
  } else data.frame()
  nonfinite_frames <- lapply(rows, `[[`, "nonfinite_examples")
  nonfinite_frames <- nonfinite_frames[vapply(nonfinite_frames, nrow, integer(1)) > 0L]
  nonfinite_examples <- if (length(nonfinite_frames)) {
    frame <- do.call(rbind, nonfinite_frames)
    rownames(frame) <- NULL
    head(frame, 100L)
  } else data.frame()
  raw_invalid_frames <- lapply(rows, `[[`, "raw_invalid_examples")
  raw_invalid_frames <- raw_invalid_frames[
    vapply(raw_invalid_frames, nrow, integer(1)) > 0L
  ]
  raw_invalid_examples <- if (length(raw_invalid_frames)) {
    frame <- do.call(rbind, raw_invalid_frames)
    rownames(frame) <- NULL
    head(frame, 100L)
  } else data.frame()

  expected_pair_distances <- n_queries * n_queries
  expected_channel_records <- 2L * expected_pair_distances
  selected_path_category_records <-
    counts[["n_selected_requested_only"]] +
    counts[["n_selected_requested_only_invalid"]] +
    counts[["n_selected_analytic_constant_raw_diagnostic"]] +
    counts[["n_empty_channel_analytic_contracts"]] +
    counts[["n_unexpected_selected_paths"]]
  provenance_gate_pass <-
    counts[["n_pair_distances"]] == expected_pair_distances &&
    counts[["n_channel_records"]] == expected_channel_records &&
    selected_path_category_records == counts[["n_channel_records"]] &&
    counts[["n_solver_attempts"]] ==
      counts[["n_solver_eligible_channel_records"]] &&
    counts[["n_requested_setting_attempts"]] ==
      counts[["n_solver_attempts"]] &&
    counts[["n_unexpected_setting_attempts"]] == 0L &&
    counts[["n_retry_attempts"]] == 0L &&
    counts[["n_exact_attempts"]] == 0L &&
    counts[["n_records_with_multiple_attempts"]] == 0L &&
    counts[["n_record_contract_mismatches"]] == 0L &&
    counts[["n_fallback_records"]] == 0L &&
    counts[["n_unexpected_selected_paths"]] == 0L
  selected_output_finite_gate_pass <-
    counts[["n_nonfinite_channel_distances"]] == 0L &&
    counts[["n_nonfinite_combined_distances"]] == 0L
  raw_diagnostic_plan_valid_gate_pass <-
    counts[["n_raw_plan_validation_failures"]] == 0L &&
    counts[["n_raw_invalid_channel_records"]] == 0L
  marginal_residual_gate_pass <-
    counts[["n_marginal_residual_exceeds_tolerance"]] == 0L

  residual_values <- if (nrow(residual_extremes)) {
    residual_extremes$max_marginal_residual_linf
  } else numeric()
  audit <- cbind(
    condition[1L, , drop = FALSE],
    data.frame(
      method = "ppm_wasserstein",
      ot_estimand = "entropic_finite_iteration_raw_requested_only",
      solver_execution_policy =
        "raw_requested_only_no_retry_no_exact_fallback",
      solver_orientation = "query_to_library",
      solver_api = "approxOT::transport_plan_given_C",
      selected_paths_observed = paste(selected_paths, collapse = ";"),
      validation_errors_observed = if (length(validation_errors)) {
        paste(validation_errors, collapse = ";")
      } else "<none>",
      sinkhorn_epsilon_requested = requested_epsilon,
      sinkhorn_iterations_requested = requested_niter,
      n_pair_distances = counts[["n_pair_distances"]],
      n_pair_distances_expected = expected_pair_distances,
      n_channel_records = counts[["n_channel_records"]],
      n_channel_records_expected = expected_channel_records,
      n_selected_path_category_records = selected_path_category_records,
      n_solver_eligible_channel_records =
        counts[["n_solver_eligible_channel_records"]],
      n_solver_attempts = counts[["n_solver_attempts"]],
      n_requested_setting_attempts =
        counts[["n_requested_setting_attempts"]],
      n_unexpected_setting_attempts =
        counts[["n_unexpected_setting_attempts"]],
      n_retry_attempts = counts[["n_retry_attempts"]],
      n_exact_attempts = counts[["n_exact_attempts"]],
      n_records_with_multiple_attempts =
        counts[["n_records_with_multiple_attempts"]],
      n_record_contract_mismatches =
        counts[["n_record_contract_mismatches"]],
      n_fallback_records = counts[["n_fallback_records"]],
      n_selected_requested_only = counts[["n_selected_requested_only"]],
      n_selected_requested_only_invalid =
        counts[["n_selected_requested_only_invalid"]],
      n_selected_analytic_constant_raw_diagnostic =
        counts[["n_selected_analytic_constant_raw_diagnostic"]],
      n_analytic_input_contracts = counts[["n_analytic_input_contracts"]],
      n_empty_channel_analytic_contracts =
        counts[["n_empty_channel_analytic_contracts"]],
      n_unexpected_selected_paths =
        counts[["n_unexpected_selected_paths"]],
      n_nonfinite_channel_distances =
        counts[["n_nonfinite_channel_distances"]],
      n_selected_nonfinite_channel_distances =
        counts[["n_nonfinite_channel_distances"]],
      n_nonfinite_combined_distances =
        counts[["n_nonfinite_combined_distances"]],
      n_selected_nonfinite_combined_distances =
        counts[["n_nonfinite_combined_distances"]],
      n_raw_invalid_channel_records =
        counts[["n_raw_invalid_channel_records"]],
      n_raw_plan_validation_failures =
        counts[["n_raw_plan_validation_failures"]],
      n_raw_plan_nonfinite_attempts =
        counts[["n_raw_plan_nonfinite_attempts"]],
      n_raw_diagnostic_nonfinite_attempts =
        counts[["n_raw_plan_nonfinite_attempts"]],
      n_raw_plan_nonfinite_mass = counts[["n_raw_plan_nonfinite_mass"]],
      raw_plan_validation_failure_rate = if (counts[["n_solver_attempts"]]) {
        counts[["n_raw_plan_validation_failures"]] /
          counts[["n_solver_attempts"]]
      } else NA_real_,
      n_marginal_residual_exceeds_tolerance =
        counts[["n_marginal_residual_exceeds_tolerance"]],
      max_marginal_residual_linf = finite_max_or_na(residual_values),
      marginal_residual_tolerance =
        getFromNamespace("PPMWASS_OT_MARGINAL_TOLERANCE", "ppmWass"),
      requested_only_provenance_gate_pass = provenance_gate_pass,
      selected_output_finite_gate_pass = selected_output_finite_gate_pass,
      finite_output_gate_pass = selected_output_finite_gate_pass,
      raw_diagnostic_plan_valid_gate_pass =
        raw_diagnostic_plan_valid_gate_pass,
      raw_solver_health_gate_pass = raw_diagnostic_plan_valid_gate_pass,
      marginal_residual_gate_pass = marginal_residual_gate_pass,
      retrieval_metric_gate_pass =
        provenance_gate_pass && selected_output_finite_gate_pass,
      strict_plan_gate_pass = provenance_gate_pass &&
        selected_output_finite_gate_pass &&
        raw_diagnostic_plan_valid_gate_pass && marginal_residual_gate_pass,
      stringsAsFactors = FALSE
    )
  )
  rownames(audit) <- NULL
  list(
    distance_matrix = distance_matrix, solver_audit = audit,
    nonfinite_examples = nonfinite_examples,
    raw_invalid_examples = raw_invalid_examples,
    residual_extremes = residual_extremes
  )
}

single_target_metrics <- function(distance_matrix, condition, method,
                                  setting_row, query_seeds) {
  top_k <- c(1L, 3L, 5L, 10L)
  rows <- vector("list", n_queries)
  condition_block <- condition[rep(1L, n_queries), , drop = FALSE]
  setting_block <- setting_row[rep(1L, n_queries), , drop = FALSE]

  for (i in seq_len(n_queries)) {
    d <- distance_matrix[i, ]
    target <- i
    target_distance <- d[[target]]
    retrieval_metric_eligible <- all(is.finite(d))
    n_nonfinite_library_distances <- sum(!is.finite(d))
    if (retrieval_metric_eligible) {
      ranked <- order(d, seq_along(d))
      old_rank <- match(target, ranked)
      n_better <- sum(d < target_distance)
      tied <- d == target_distance
      tie_size <- sum(tied)
      optimistic_rank <- n_better + 1L
      pessimistic_rank <- n_better + tie_size
      fractional_topk <- vapply(top_k, function(k) {
        max(0, min(1, (k - n_better) / tie_size))
      }, numeric(1))
      names(fractional_topk) <- paste0("top", top_k, "_fractional")
      nn_library_id <- ids[[ranked[[1]]]]
      nn_compound_class <- compound_class[[ranked[[1]]]]
      nn_distance <- d[[ranked[[1]]]]
      rr_fractional_value <- mean(1 / (n_better + seq_len(tie_size)))
    } else {
      old_rank <- n_better <- tie_size <- optimistic_rank <-
        pessimistic_rank <- NA_integer_
      fractional_topk <- stats::setNames(
        rep(NA_real_, length(top_k)), paste0("top", top_k, "_fractional")
      )
      nn_library_id <- NA_character_
      nn_compound_class <- NA_character_
      nn_distance <- NA_real_
      rr_fractional_value <- NA_real_
    }

    metric <- data.frame(
      query_index = i,
      query_id = ids[[i]],
      query_compound_class = compound_class[[i]],
      target_library_id = ids[[target]],
      query_seed = query_seeds[[i]],
      target_distance = target_distance,
      retrieval_metric_eligible = retrieval_metric_eligible,
      n_nonfinite_library_distances = n_nonfinite_library_distances,
      nn_library_id = nn_library_id,
      nn_compound_class = nn_compound_class,
      nn_distance = nn_distance,
      rank_old_order = old_rank,
      n_better = n_better,
      tie_size = tie_size,
      optimistic_rank = optimistic_rank,
      pessimistic_rank = pessimistic_rank,
      top_1 = as.numeric(old_rank <= 1L),
      top_3 = as.numeric(old_rank <= 3L),
      top_5 = as.numeric(old_rank <= 5L),
      top_10 = as.numeric(old_rank <= 10L),
      top1_optimistic = as.numeric(optimistic_rank <= 1L),
      top1_fractional = fractional_topk[["top1_fractional"]],
      top1_pessimistic = as.numeric(pessimistic_rank <= 1L),
      top3_optimistic = as.numeric(optimistic_rank <= 3L),
      top3_fractional = fractional_topk[["top3_fractional"]],
      top3_pessimistic = as.numeric(pessimistic_rank <= 3L),
      top5_optimistic = as.numeric(optimistic_rank <= 5L),
      top5_fractional = fractional_topk[["top5_fractional"]],
      top5_pessimistic = as.numeric(pessimistic_rank <= 5L),
      top10_optimistic = as.numeric(optimistic_rank <= 10L),
      top10_fractional = fractional_topk[["top10_fractional"]],
      top10_pessimistic = as.numeric(pessimistic_rank <= 10L),
      rr_old_order = 1 / old_rank,
      rr_optimistic = 1 / optimistic_rank,
      rr_fractional = rr_fractional_value,
      rr_pessimistic = 1 / pessimistic_rank,
      stringsAsFactors = FALSE
    )
    rows[[i]] <- cbind(
      condition_block[i, , drop = FALSE],
      data.frame(method = method, stringsAsFactors = FALSE),
      setting_block[i, , drop = FALSE],
      metric
    )
  }
  do.call(rbind, rows)
}

nearest_score_auroc <- function(per_query) {
  eligible <- !is.na(per_query$retrieval_metric_eligible) &
    per_query$retrieval_metric_eligible
  correct <- !is.na(per_query$top_1) & per_query$top_1 == 1
  d_correct <- per_query$nn_distance[correct & is.finite(per_query$nn_distance)]
  d_incorrect <- per_query$nn_distance[
    eligible & !correct & is.finite(per_query$nn_distance)
  ]
  if (!length(d_correct) || !length(d_incorrect)) return(NA_real_)
  comparison <- outer(d_correct, d_incorrect, "<")
  tied <- outer(d_correct, d_incorrect, "==")
  mean(comparison + 0.5 * tied)
}

summarize_one <- function(per_query) {
  identity_cols <- c(
    names(condition_grids[[per_query$task[[1]]]]),
    "method", "hard_match_tolerance_ppm", "ppmWass_base_cost_ppm",
    "transition_multiplier", "saturation_width_ppm",
    "sinkhorn_epsilon_effective", "sinkhorn_iterations_effective",
    "solver_backend", "ot_estimand", "solver_execution_policy",
    "solver_orientation", "selected_path_contract", "distance_matrix_backend"
  )
  identity_cols <- unique(identity_cols[identity_cols %in% names(per_query)])
  out <- per_query[1L, identity_cols, drop = FALSE]
  metric_cols <- c(
    "top_1", "top_3", "top_5", "top_10",
    "top1_optimistic", "top1_fractional", "top1_pessimistic",
    "top3_optimistic", "top3_fractional", "top3_pessimistic",
    "top5_optimistic", "top5_fractional", "top5_pessimistic",
    "top10_optimistic", "top10_fractional", "top10_pessimistic",
    "rr_old_order", "rr_optimistic", "rr_fractional", "rr_pessimistic"
  )
  for (column in metric_cols) {
    value <- per_query[[column]]
    out[[column]] <- if (all(is.finite(value))) mean(value) else NA_real_
  }
  out$score_auroc <- nearest_score_auroc(per_query)
  out$n_queries <- nrow(per_query)
  out$n_metric_eligible_queries <- sum(per_query$retrieval_metric_eligible)
  out$n_metric_ineligible_queries <- sum(!per_query$retrieval_metric_eligible)
  out$n_nonfinite_library_distances <-
    sum(per_query$n_nonfinite_library_distances)
  out$n_correct <- sum(per_query$top_1 == 1, na.rm = TRUE)
  out$n_incorrect <- sum(per_query$top_1 == 0, na.rm = TRUE)
  out$n_rank1_target_ties <- sum(
    per_query$n_better == 0 & per_query$tie_size > 1, na.rm = TRUE
  )
  out$n_any_target_ties <- sum(per_query$tie_size > 1, na.rm = TRUE)
  out
}

checkpoint_path_for <- function(task, condition_id, method) {
  file.path(
    output_dir, "checkpoints", task,
    paste0(condition_id, "__", method, ".rds")
  )
}

read_valid_checkpoint <- function(path, task, condition_id, method) {
  if (!file.exists(path)) return(NULL)
  object <- readRDS(path)
  valid <- identical(object$schema_version, 3L) &&
    identical(object$configuration_fingerprint, configuration_fingerprint) &&
    identical(object$task, task) &&
    identical(object$condition_id, condition_id) &&
    identical(object$method, method)
  if (!valid) stop("Incompatible checkpoint: ", path)
  object
}

run_task <- function(task) {
  spec <- task_specs[[task]]
  grid <- condition_grids[[task]]
  message("=== ", task, ": ", nrow(grid), " conditions ===")
  cached_seed <- NA_integer_
  cached_query <- NULL

  for (row_index in seq_len(nrow(grid))) {
    condition <- grid[row_index, , drop = FALSE]
    pending <- vapply(spec$methods, function(method) {
      is.null(read_valid_checkpoint(
        checkpoint_path_for(task, condition$condition_id, method),
        task, condition$condition_id, method
      ))
    }, logical(1))
    pending_methods <- spec$methods[pending]
    if (!length(pending_methods)) {
      message("[", task, "] resume skip: ", condition$condition_id)
      next
    }

    if (is.null(cached_query) ||
        !identical(cached_seed, as.integer(condition$condition_seed))) {
      cached_query <- build_regenerated_query(condition, spec$builder)
      cached_seed <- as.integer(condition$condition_seed)
    }

    for (method in pending_methods) {
      effective <- settings_for(task, condition, method)
      message(
        "[", task, " ", row_index, "/", nrow(grid), "] ",
        condition$condition_id, " | ", method
      )
      started <- proc.time()[["elapsed"]]
      requested_result <- NULL
      if (isTRUE(effective$raw_requested_only)) {
        requested_result <- compute_requested_only_matrix(
          cached_query, effective, condition
        )
        distance_matrix <- requested_result$distance_matrix
      } else {
        distance_matrix <- compute_distance_matrix_search(
          query_frag_list = cached_query$frag_list,
          query_loss_list = cached_query$derived_list,
          lib_frag_list = library_frag,
          lib_loss_list = library_derived,
          params = effective$params,
          progress = FALSE
        )
      }
      elapsed <- proc.time()[["elapsed"]] - started
      validate_distance_matrix(
        distance_matrix, cached_query, task, condition, method,
        allow_nonfinite = isTRUE(effective$raw_requested_only)
      )
      per_query <- single_target_metrics(
        distance_matrix, condition, method, effective$row,
        cached_query$query_seeds
      )
      summary <- summarize_one(per_query)
      summary$elapsed_seconds <- elapsed
      audit_gate_cols <- c(
        "requested_only_provenance_gate_pass",
        "selected_output_finite_gate_pass", "finite_output_gate_pass",
        "raw_diagnostic_plan_valid_gate_pass", "raw_solver_health_gate_pass",
        "marginal_residual_gate_pass",
        "retrieval_metric_gate_pass",
        "strict_plan_gate_pass", "n_nonfinite_channel_distances",
        "n_nonfinite_combined_distances",
        "n_raw_plan_validation_failures", "n_raw_plan_nonfinite_attempts",
        "raw_plan_validation_failure_rate",
        "n_marginal_residual_exceeds_tolerance",
        "max_marginal_residual_linf"
      )
      for (column in audit_gate_cols) summary[[column]] <- NA
      if (!is.null(requested_result)) {
        for (column in audit_gate_cols) {
          summary[[column]] <- requested_result$solver_audit[[column]]
        }
      }
      checkpoint <- list(
        schema_version = 3L,
        configuration_fingerprint = configuration_fingerprint,
        task = task,
        condition_id = condition$condition_id,
        method = method,
        per_query = per_query,
        summary = summary,
        requested_only_solver_audit = if (!is.null(requested_result)) {
          requested_result$solver_audit
        } else NULL,
        requested_only_nonfinite_examples = if (!is.null(requested_result)) {
          requested_result$nonfinite_examples
        } else NULL,
        requested_only_raw_invalid_examples = if (!is.null(requested_result)) {
          requested_result$raw_invalid_examples
        } else NULL,
        requested_only_residual_extremes = if (!is.null(requested_result)) {
          requested_result$residual_extremes
        } else NULL,
        regeneration_qc = cbind(
          condition[rep(1L, nrow(cached_query$qc)), , drop = FALSE],
          cached_query$qc
        )
      )
      atomic_save_rds(
        checkpoint, checkpoint_path_for(task, condition$condition_id, method)
      )
    }
  }
  invisible(TRUE)
}

metric_columns <- c(
  "top_1", "top_3", "top_5", "top_10",
  "top1_optimistic", "top1_fractional", "top1_pessimistic",
  "top3_optimistic", "top3_fractional", "top3_pessimistic",
  "top5_optimistic", "top5_fractional", "top5_pessimistic",
  "top10_optimistic", "top10_fractional", "top10_pessimistic",
  "rr_old_order", "rr_optimistic", "rr_fractional", "rr_pessimistic",
  "score_auroc", "elapsed_seconds"
)

group_summary <- function(data, group_cols, metrics = metric_columns) {
  group_cols <- unique(group_cols[group_cols %in% names(data)])
  metrics <- metrics[metrics %in% names(data)]
  key <- na_safe_group_key_internal(data, group_cols)
  groups <- split(seq_len(nrow(data)), key)
  rows <- lapply(groups, function(index) {
    block <- data[index, , drop = FALSE]
    out <- block[1L, group_cols, drop = FALSE]
    for (metric in metrics) {
      value <- block[[metric]]
      require_complete_sinkhorn_metric <-
        identical(as.character(block$task[[1L]]), "sinkhorn") &&
        grepl("^(top_|top[0-9]+_|rr_)", metric)
      metric_summary <- summarize_publication_metric_internal(
        value, require_complete = require_complete_sinkhorn_metric
      )
      out[[paste0("n_finite_", metric, "_replicates")]] <-
        metric_summary$n_finite
      out[[paste0("mean_", metric)]] <- metric_summary$mean
      out[[paste0("sd_", metric)]] <- metric_summary$sd
    }
    out$n_replicates <- length(unique(block$replicate))
    out
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

class_summary <- function(per_query, spec, include_replicate = TRUE) {
  group_cols <- c(
    "task", "method", spec$analysis_cols,
    if (include_replicate) "replicate" else NULL,
    "query_compound_class"
  )
  key <- interaction(per_query[group_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(per_query)), key)
  rows <- lapply(groups, function(index) {
    block <- per_query[index, , drop = FALSE]
    out <- block[1L, group_cols, drop = FALSE]
    for (column in c(
      "top_1", "top_3", "top_5", "top_10",
      "top1_fractional", "top3_fractional",
      "top5_fractional", "top10_fractional"
    )) {
      out[[column]] <- mean(block[[column]])
    }
    out$n_query_replicates <- nrow(block)
    out$n_unique_queries <- length(unique(block$query_id))
    out
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

failure_outputs <- function(per_query, task_dir, spec) {
  is_failure_condition <- per_query$mz_noise_ppm == failure_mz_ppm &
    per_query$intensity_noise_sd == failure_intensity_sd
  failures <- per_query[is_failure_condition & per_query$top_1 == 0, , drop = FALSE]
  if (!nrow(failures)) {
    atomic_write_csv(failures, file.path(task_dir, "failure_cases_all.csv"))
    atomic_write_csv(failures, file.path(task_dir, "failure_cases_examples.csv"))
    atomic_write_csv(
      data.frame(
        method = character(), query_compound_class = character(),
        nn_compound_class = character(), n = integer(),
        same_compound_class = logical(), stringsAsFactors = FALSE
      ),
      file.path(task_dir, "failure_confusion.csv")
    )
    atomic_write_csv(
      data.frame(
        method = character(), same_compound_class = logical(), n = integer(),
        total_failures = integer(), fraction = numeric(),
        stringsAsFactors = FALSE
      ),
      file.path(task_dir, "failure_same_cross_fraction.csv")
    )
    return(invisible(NULL))
  }
  failures$same_compound_class <-
    failures$query_compound_class == failures$nn_compound_class
  atomic_write_csv(failures, file.path(task_dir, "failure_cases_all.csv"))

  by_method <- split(failures, failures$method)
  examples <- do.call(rbind, lapply(by_method, function(block) {
    block <- block[
      order(block$rank_old_order, decreasing = TRUE, block$query_id),
      , drop = FALSE
    ]
    head(block, failure_examples_per_method)
  }))
  rownames(examples) <- NULL
  atomic_write_csv(examples, file.path(task_dir, "failure_cases_examples.csv"))

  confusion_cols <- c(
    "method", "query_compound_class", "nn_compound_class"
  )
  key <- interaction(failures[confusion_cols], drop = TRUE, lex.order = TRUE)
  confusion <- do.call(rbind, lapply(split(seq_len(nrow(failures)), key), function(index) {
    out <- failures[index[[1]], confusion_cols, drop = FALSE]
    out$n <- length(index)
    out
  }))
  rownames(confusion) <- NULL
  confusion$same_compound_class <-
    confusion$query_compound_class == confusion$nn_compound_class
  atomic_write_csv(confusion, file.path(task_dir, "failure_confusion.csv"))

  key_same <- interaction(
    failures$method, failures$same_compound_class,
    drop = TRUE, lex.order = TRUE
  )
  same_cross <- do.call(rbind, lapply(split(seq_len(nrow(failures)), key_same), function(index) {
    block <- failures[index, , drop = FALSE]
    data.frame(
      method = block$method[[1]],
      same_compound_class = block$same_compound_class[[1]],
      n = nrow(block),
      stringsAsFactors = FALSE
    )
  }))
  totals <- tapply(same_cross$n, same_cross$method, sum)
  same_cross$total_failures <- unname(totals[same_cross$method])
  same_cross$fraction <- same_cross$n / same_cross$total_failures
  atomic_write_csv(
    same_cross, file.path(task_dir, "failure_same_cross_fraction.csv")
  )
  invisible(NULL)
}

write_task_outputs <- function(task) {
  spec <- task_specs[[task]]
  grid <- condition_grids[[task]]
  checkpoints <- list()
  index <- 0L
  for (row_index in seq_len(nrow(grid))) {
    condition <- grid[row_index, , drop = FALSE]
    for (method in spec$methods) {
      path <- checkpoint_path_for(task, condition$condition_id, method)
      object <- read_valid_checkpoint(
        path, task, condition$condition_id, method
      )
      if (is.null(object)) stop("Missing completed checkpoint: ", path)
      index <- index + 1L
      checkpoints[[index]] <- object
    }
  }

  task_dir <- file.path(output_dir, task)
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  summaries <- do.call(rbind, lapply(checkpoints, function(x) x[["summary"]]))
  rownames(summaries) <- NULL
  atomic_write_csv(summaries, file.path(task_dir, "summary_by_replicate.csv"))

  setting_cols <- c(
    "hard_match_tolerance_ppm", "ppmWass_base_cost_ppm",
    "transition_multiplier", "saturation_width_ppm",
    "sinkhorn_epsilon_effective", "sinkhorn_iterations_effective",
    "solver_backend", "ot_estimand", "solver_execution_policy",
    "solver_orientation", "selected_path_contract", "distance_matrix_backend"
  )
  aggregate_cols <- c("task", "method", spec$analysis_cols, setting_cols)
  across <- group_summary(summaries, aggregate_cols)
  atomic_write_csv(
    across, file.path(task_dir, "summary_across_replicates.csv")
  )

  regeneration_qc <- unique(do.call(
    rbind, lapply(checkpoints, function(x) x[["regeneration_qc"]])
  ))
  rownames(regeneration_qc) <- NULL
  atomic_write_csv(
    regeneration_qc, file.path(task_dir, "representation_regeneration_qc.csv")
  )

  per_query <- do.call(rbind, lapply(checkpoints, function(x) x[["per_query"]]))
  rownames(per_query) <- NULL
  if (write_query_level) {
    atomic_write_csv(per_query, file.path(task_dir, "per_query.csv"))
  }
  requested_audits <- Filter(
    Negate(is.null),
    lapply(checkpoints, `[[`, "requested_only_solver_audit")
  )
  requested_nonfinite <- Filter(
    function(x) !is.null(x) && nrow(x) > 0L,
    lapply(checkpoints, `[[`, "requested_only_nonfinite_examples")
  )
  requested_raw_invalid <- Filter(
    function(x) !is.null(x) && nrow(x) > 0L,
    lapply(checkpoints, `[[`, "requested_only_raw_invalid_examples")
  )
  requested_residuals <- Filter(
    function(x) !is.null(x) && nrow(x) > 0L,
    lapply(checkpoints, `[[`, "requested_only_residual_extremes")
  )
  requested_audit_frame <- if (length(requested_audits)) {
    do.call(rbind, requested_audits)
  } else NULL
  requested_nonfinite_frame <- if (length(requested_nonfinite)) {
    do.call(rbind, requested_nonfinite)
  } else NULL
  requested_raw_invalid_frame <- if (length(requested_raw_invalid)) {
    do.call(rbind, requested_raw_invalid)
  } else NULL
  requested_residual_frame <- if (length(requested_residuals)) {
    do.call(rbind, requested_residuals)
  } else NULL
  atomic_save_rds(
    list(
      per_query = per_query,
      summary_by_replicate = summaries,
      summary_across_replicates = across,
      representation_regeneration_qc = regeneration_qc,
      requested_only_solver_audit = requested_audit_frame,
      requested_only_nonfinite_examples = requested_nonfinite_frame,
      requested_only_raw_invalid_examples = requested_raw_invalid_frame,
      requested_only_residual_extremes = requested_residual_frame
    ),
    file.path(task_dir, "results.rds")
  )

  auroc_cols <- c(
    "task", "method", spec$analysis_cols, "replicate",
    "score_auroc", "n_correct", "n_incorrect"
  )
  atomic_write_csv(
    summaries[, unique(auroc_cols[auroc_cols %in% names(summaries)]), drop = FALSE],
    file.path(task_dir, "nearest_score_correctness_auroc.csv")
  )

  if (task == "breakdown") {
    atomic_write_csv(
      class_summary(per_query, spec, include_replicate = TRUE),
      file.path(task_dir, "class_summary_by_replicate.csv")
    )
    atomic_write_csv(
      class_summary(per_query, spec, include_replicate = FALSE),
      file.path(task_dir, "class_summary_pooled_replicates.csv")
    )
    failure_outputs(per_query, task_dir, spec)
  }

  if (task == "tolerance") {
    optimal_key <- interaction(
      across[c("method", "mz_noise_ppm")],
      drop = TRUE, lex.order = TRUE
    )
    optimal <- do.call(rbind, lapply(split(seq_len(nrow(across)), optimal_key), function(ii) {
      block <- across[ii, , drop = FALSE]
      block <- block[
        order(-block$mean_top_1, block$tol_ppm),
        , drop = FALSE
      ]
      block[1L, , drop = FALSE]
    }))
    rownames(optimal) <- NULL
    atomic_write_csv(optimal, file.path(task_dir, "optimal_tolerance.csv"))
  }

  if (task == "sinkhorn") {
    if (is.null(requested_audit_frame) || !nrow(requested_audit_frame)) {
      stop("Sinkhorn task is missing requested-only solver audit rows.")
    }
    rownames(requested_audit_frame) <- NULL
    atomic_write_csv(
      requested_audit_frame,
      file.path(task_dir, "requested_only_solver_audit_by_condition.csv")
    )

    if (is.null(requested_nonfinite_frame)) {
      requested_nonfinite_frame <- data.frame(
        task = character(), condition_id = character(),
        orientation = character(), query_id = character(),
        library_id = character(), value = numeric(),
        query_fragment_count = integer(), query_derived_count = integer(),
        library_fragment_count = integer(), library_derived_count = integer(),
        stringsAsFactors = FALSE
      )
    }
    rownames(requested_nonfinite_frame) <- NULL
    atomic_write_csv(
      requested_nonfinite_frame,
      file.path(task_dir, "requested_only_nonfinite_examples.csv")
    )

    if (is.null(requested_raw_invalid_frame)) {
      requested_raw_invalid_frame <- data.frame(
        task = character(), condition_id = character(),
        query_id = character(), library_id = character(),
        channel = character(), orientation = character(),
        selected_path = character(), selected_distance = numeric(),
        analytic_identity_used = logical(), record_status = character(),
        raw_step = character(), raw_backend = character(),
        sinkhorn_epsilon_requested = numeric(),
        sinkhorn_iterations_requested = integer(),
        sinkhorn_epsilon_executed = numeric(),
        sinkhorn_iterations_executed = integer(), raw_total_cost = numeric(),
        raw_accepted_at_marginal_tolerance = logical(),
        raw_validation_error = character(), raw_solver_error = character(),
        raw_nonfinite_plan_mass = integer(), raw_row_residual_linf = numeric(),
        raw_col_residual_linf = numeric(), stringsAsFactors = FALSE
      )
    }
    rownames(requested_raw_invalid_frame) <- NULL
    atomic_write_csv(
      requested_raw_invalid_frame,
      file.path(task_dir, "requested_only_raw_invalid_examples.csv")
    )

    if (is.null(requested_residual_frame)) {
      requested_residual_frame <- data.frame(
        task = character(), condition_id = character(),
        query_id = character(), library_id = character(),
        channel = character(), orientation = character(),
        selected_path = character(), solver_backend = character(),
        sinkhorn_epsilon_requested = numeric(),
        sinkhorn_iterations_requested = integer(),
        sinkhorn_epsilon_executed = numeric(),
        sinkhorn_iterations_executed = integer(),
        total_cost = numeric(), accepted_at_marginal_tolerance = logical(),
        validation_error = character(), row_residual_linf = numeric(),
        col_residual_linf = numeric(), max_marginal_residual_linf = numeric(),
        nonfinite_plan_mass = integer(), stringsAsFactors = FALSE
      )
    }
    if (nrow(requested_residual_frame)) {
      requested_residual_frame <- requested_residual_frame[
        order(requested_residual_frame$max_marginal_residual_linf,
              decreasing = TRUE, na.last = TRUE), , drop = FALSE
      ]
    }
    rownames(requested_residual_frame) <- NULL
    atomic_write_csv(
      requested_residual_frame,
      file.path(task_dir, "requested_only_plan_residual_extremes.csv")
    )

    requested_solver_manifest <- data.frame(
      schema_version = 3L,
      task = "sinkhorn",
      ot_estimand = "entropic_finite_iteration_raw_requested_only",
      execution_policy = "raw_requested_only_no_retry_no_exact_fallback",
      orientation = "query_to_library",
      solver_api = "approxOT::transport_plan_given_C",
      n_conditions = nrow(requested_audit_frame),
      n_pair_distances = sum(requested_audit_frame$n_pair_distances),
      n_pair_distances_expected =
        sum(requested_audit_frame$n_pair_distances_expected),
      n_channel_records = sum(requested_audit_frame$n_channel_records),
      n_channel_records_expected =
        sum(requested_audit_frame$n_channel_records_expected),
      n_selected_path_category_records =
        sum(requested_audit_frame$n_selected_path_category_records),
      n_solver_attempts = sum(requested_audit_frame$n_solver_attempts),
      n_requested_setting_attempts =
        sum(requested_audit_frame$n_requested_setting_attempts),
      n_retry_attempts = sum(requested_audit_frame$n_retry_attempts),
      n_exact_attempts = sum(requested_audit_frame$n_exact_attempts),
      n_fallback_records = sum(requested_audit_frame$n_fallback_records),
      n_record_contract_mismatches =
        sum(requested_audit_frame$n_record_contract_mismatches),
      n_selected_analytic_constant_raw_diagnostic = sum(
        requested_audit_frame$n_selected_analytic_constant_raw_diagnostic
      ),
      n_empty_channel_analytic_contracts =
        sum(requested_audit_frame$n_empty_channel_analytic_contracts),
      n_nonfinite_channel_distances =
        sum(requested_audit_frame$n_nonfinite_channel_distances),
      n_selected_nonfinite_channel_distances =
        sum(requested_audit_frame$n_selected_nonfinite_channel_distances),
      n_nonfinite_combined_distances =
        sum(requested_audit_frame$n_nonfinite_combined_distances),
      n_selected_nonfinite_combined_distances =
        sum(requested_audit_frame$n_selected_nonfinite_combined_distances),
      n_raw_invalid_channel_records =
        sum(requested_audit_frame$n_raw_invalid_channel_records),
      n_raw_plan_validation_failures =
        sum(requested_audit_frame$n_raw_plan_validation_failures),
      n_raw_plan_nonfinite_attempts =
        sum(requested_audit_frame$n_raw_plan_nonfinite_attempts),
      n_raw_diagnostic_nonfinite_attempts =
        sum(requested_audit_frame$n_raw_diagnostic_nonfinite_attempts),
      n_raw_plan_nonfinite_mass =
        sum(requested_audit_frame$n_raw_plan_nonfinite_mass),
      raw_plan_validation_failure_rate =
        sum(requested_audit_frame$n_raw_plan_validation_failures) /
        sum(requested_audit_frame$n_solver_attempts),
      n_marginal_residual_exceeds_tolerance = sum(
        requested_audit_frame$n_marginal_residual_exceeds_tolerance
      ),
      max_marginal_residual_linf = finite_max_or_na(
        requested_audit_frame$max_marginal_residual_linf
      ),
      all_requested_only_provenance_gate_pass =
        all(requested_audit_frame$requested_only_provenance_gate_pass),
      all_selected_output_finite_gate_pass =
        all(requested_audit_frame$selected_output_finite_gate_pass),
      all_finite_output_gate_pass =
        all(requested_audit_frame$finite_output_gate_pass),
      all_raw_diagnostic_plan_valid_gate_pass =
        all(requested_audit_frame$raw_diagnostic_plan_valid_gate_pass),
      all_raw_solver_health_gate_pass =
        all(requested_audit_frame$raw_solver_health_gate_pass),
      all_marginal_residual_gate_pass =
        all(requested_audit_frame$marginal_residual_gate_pass),
      n_retrieval_metric_eligible_conditions =
        sum(requested_audit_frame$retrieval_metric_gate_pass),
      n_retrieval_metric_ineligible_conditions =
        sum(!requested_audit_frame$retrieval_metric_gate_pass),
      all_retrieval_metric_gate_pass =
        all(requested_audit_frame$retrieval_metric_gate_pass),
      stringsAsFactors = FALSE
    )
    atomic_write_csv(
      requested_solver_manifest,
      file.path(task_dir, "requested_only_solver_manifest.csv")
    )
    if (!isTRUE(requested_solver_manifest$all_requested_only_provenance_gate_pass)) {
      stop(
        "Requested-only Sinkhorn provenance gate failed; no retry/exact-fallback ",
        "sensitivity result may be published."
      )
    }
    if (!isTRUE(requested_solver_manifest$all_retrieval_metric_gate_pass)) {
      stop(
        "At least one requested-only Sinkhorn condition lacks a complete finite ",
        "selected distance matrix; retrieval sensitivity must not be published."
      )
    }
  }
  invisible(TRUE)
}

for (task in tasks) {
  run_task(task)
  write_task_outputs(task)
}

session_text <- capture.output(sessionInfo())
atomic_write_lines(session_text, file.path(manifest_dir, "sessionInfo.txt"))

message("All task outputs are complete; closing the run log before checksumming.")
close_run_log()
output_files <- list.files(
  output_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE
)
output_files <- output_files[file.info(output_files)$isdir %in% FALSE]
if (file.exists(log_file)) output_files <- unique(c(output_files, log_file))
output_checksums <- data.frame(
  path = normalizePath(output_files, winslash = "/", mustWork = TRUE),
  md5 = unname(tools::md5sum(output_files)),
  stringsAsFactors = FALSE
)
atomic_write_csv(output_checksums, file.path(manifest_dir, "output_checksums.csv"))

task_counts <- do.call(rbind, lapply(tasks, function(task) {
  data.frame(
    task = task,
    n_conditions = nrow(condition_grids[[task]]),
    n_methods = length(task_specs[[task]]$methods),
    n_completed_matrices =
      nrow(condition_grids[[task]]) * length(task_specs[[task]]$methods),
    stringsAsFactors = FALSE
  )
}))

run_manifest <- list(
  schema_version = 3L,
  configuration_fingerprint = configuration_fingerprint,
  timestamp_started = timestamp_started,
  timestamp_completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script_path = script_path,
  command = command_text,
  package_git_commit = git_commit,
  package_tree_dirty = package_tree_dirty,
  package_git_status = git_status,
  configuration = configuration,
  n_filtered_queries = n_queries,
  input_manifest = input_manifest,
  task_counts = task_counts,
  output_checksums = output_checksums,
  session_info = sessionInfo()
)
atomic_save_rds(run_manifest, file.path(manifest_dir, "run_manifest.rds"))
atomic_write_csv(task_counts, file.path(manifest_dir, "task_counts.csv"))
atomic_write_lines(
  c(
    paste("configuration_fingerprint:", configuration_fingerprint),
    paste("completed:", run_manifest$timestamp_completed),
    paste("n_filtered_queries:", n_queries),
    paste("tasks:", paste(tasks, collapse = ",")),
    paste("smoke:", smoke),
    "status: complete"
  ),
  file.path(manifest_dir, "COMPLETED.txt")
)

message("Supplementary stress analyses completed.")
message("Results: ", output_dir)
message("Manifest: ", manifest_dir)
