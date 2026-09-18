#!/usr/bin/env Rscript

# Publication-grade explicit-TIC mixture self-retention benchmark.
#
# The fragment-only and default combined (fragment / pooled derived
# mass-difference = 70 / 30) modes use exactly the same query spectra,
# interferent sets, seeds, and mixed fragment spectra. In combined mode the
# secondary representation is regenerated from the mixed fragment spectrum;
# the clean library representation is never copied into a mixed query.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

flag_is_true <- function(name, default = FALSE) {
  value <- tolower(get_arg(name, if (default) "true" else "false"))
  value %in% c("1", "true", "yes", "y")
}

split_arg <- function(name, default) {
  value <- trimws(strsplit(get_arg(name, default), ",", fixed = TRUE)[[1L]])
  value[nzchar(value)]
}

script_path <- function() {
  full <- commandArgs(trailingOnly = FALSE)
  hit <- sub("^--file=", "", grep("^--file=", full, value = TRUE))
  if (!length(hit)) return(NA_character_)
  normalizePath(hit[[1L]], mustWork = TRUE)
}

script_file <- script_path()
repo_default <- if (is.na(script_file)) getwd() else dirname(dirname(dirname(script_file)))
repo_dir <- normalizePath(get_arg("repo", repo_default), mustWork = TRUE)
bundle_default <- if (basename(dirname(repo_dir)) == "code") {
  dirname(dirname(repo_dir))
} else {
  dirname(repo_dir)
}
bundle_dir <- normalizePath(get_arg("bundle-dir", bundle_default), mustWork = TRUE)
input_dir <- normalizePath(get_arg("input-dir", file.path(bundle_dir, "inputs")), mustWork = TRUE)

canonical_dataset <- function(x) {
  key <- gsub("[^a-z0-9]", "", tolower(x))
  if (key == "recetox") return("RECETOX")
  if (key %in% c("hrei", "hreimsdb")) return("HREI-MSDB")
  stop("Unknown dataset: ", x, ". Use RECETOX or HREI-MSDB.")
}

datasets <- unique(vapply(
  split_arg("datasets", "RECETOX,HREI-MSDB"), canonical_dataset, character(1L)
))
methods <- unique(split_arg(
  "methods",
  "ppm_wasserstein,composite,entropy_weighted,entropy_unweighted,cosine,weighted_cosine,hellinger"
))
valid_methods <- c(
  "ppm_wasserstein", "composite", "entropy_weighted",
  "entropy_unweighted", "cosine", "weighted_cosine", "hellinger"
)
bad_methods <- setdiff(methods, valid_methods)
if (length(bad_methods)) stop("Unsupported method(s): ", paste(bad_methods, collapse = ", "))

modes <- unique(split_arg("modes", "fragment_only,combined_default"))
valid_modes <- c("fragment_only", "combined_default")
bad_modes <- setdiff(modes, valid_modes)
if (length(bad_modes)) stop("Unsupported mode(s): ", paste(bad_modes, collapse = ", "))

signal_weights <- as.numeric(split_arg("signal-weights", "1,0.5,0.2,0.1,0.05"))
if (!length(signal_weights) || any(!is.finite(signal_weights)) ||
    any(signal_weights < 0 | signal_weights > 1)) {
  stop("signal-weights must be finite values between 0 and 1.")
}
signal_weights <- unique(signal_weights)

replicates <- as.integer(get_arg("replicates", "3"))
query_n <- as.integer(get_arg("query-n", "60"))
noise_components <- as.integer(get_arg("noise-components", "3"))
seed <- as.integer(get_arg("seed", "20260718"))
ri_window <- as.numeric(get_arg("interferent-ri-window", "50"))
tol_ppm <- as.numeric(get_arg("tol-ppm", "15"))
merge_ppm <- as.numeric(get_arg("mixture-merge-ppm", "10"))
transition_mult <- as.numeric(get_arg("transition-mult", "3"))
ot_method <- tolower(get_arg("ot-method", "exact"))
sinkhorn_epsilon <- as.numeric(get_arg("sinkhorn-epsilon", "0.05"))
sinkhorn_niter <- as.integer(get_arg("sinkhorn-niter", "100"))
n_cores <- as.integer(get_arg("n-cores", "8"))
use_parallel <- flag_is_true("parallel", TRUE)
show_progress <- flag_is_true("progress", FALSE)
save_distances <- flag_is_true("save-distances", FALSE)
resume <- flag_is_true("resume", FALSE)
tie_tolerance <- as.numeric(get_arg("tie-tolerance", "0"))

integer_args <- c(
  replicates = replicates, query_n = query_n,
  noise_components = noise_components, sinkhorn_niter = sinkhorn_niter,
  n_cores = n_cores
)
if (any(!is.finite(integer_args)) || any(integer_args < 1L)) {
  stop("replicates, query-n, noise-components, sinkhorn-niter, and n-cores must be >= 1.")
}
if (!is.finite(seed) || seed < 1L) stop("seed must be a positive integer.")
if (!is.finite(ri_window) || ri_window < 0) stop("interferent-ri-window must be non-negative.")
if (!is.finite(tol_ppm) || tol_ppm <= 0) stop("tol-ppm must be positive.")
if (!is.finite(merge_ppm) || merge_ppm < 0) stop("mixture-merge-ppm must be non-negative.")
if (!is.finite(transition_mult) || transition_mult <= 0) stop("transition-mult must be positive.")
if (!ot_method %in% c("exact", "sinkhorn", "greenkhorn")) {
  stop("ot-method must be exact, sinkhorn, or greenkhorn.")
}
if (!is.finite(sinkhorn_epsilon) || sinkhorn_epsilon <= 0) stop("sinkhorn-epsilon must be positive.")
if (!is.finite(tie_tolerance) || tie_tolerance < 0) stop("tie-tolerance must be non-negative.")

recetox_msp <- normalizePath(
  get_arg("recetox-msp", file.path(input_dir, "RECETOX_merged.msp")),
  mustWork = TRUE
)
hrei_msp <- normalizePath(
  get_arg("hrei-msp", file.path(input_dir, "HREI-MSDB.msp")),
  mustWork = TRUE
)
recetox_spectra_rds <- get_arg("recetox-spectra-rds", "")
hrei_spectra_rds <- get_arg("hrei-spectra-rds", "")
if (nzchar(recetox_spectra_rds)) {
  recetox_spectra_rds <- normalizePath(recetox_spectra_rds, mustWork = TRUE)
}
if (nzchar(hrei_spectra_rds)) {
  hrei_spectra_rds <- normalizePath(hrei_spectra_rds, mustWork = TRUE)
}

output_dir <- normalizePath(
  get_arg("output-dir", file.path(bundle_dir, "results", "current", "mixture_explicit_tic")),
  mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
audit_dir <- file.path(output_dir, "audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
for (mode in modes) dir.create(file.path(output_dir, mode), recursive = TRUE, showWarnings = FALSE)
if (save_distances) dir.create(file.path(output_dir, "distance_matrices"), recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this script from the development tree.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required for the publication manifest.")
}
pkgload::load_all(repo_dir, quiet = TRUE)
if (identical(ot_method, "exact") &&
    !requireNamespace("transport", quietly = TRUE)) {
  stop("ot-method=exact requires transport; publication runs fail closed.")
}

stats_helper <- file.path(repo_dir, "inst", "scripts", "lib",
                          "publication_retrieval_statistics.R")
if (!file.exists(stats_helper)) stop("Missing statistics helper: ", stats_helper)
source(stats_helper, local = TRUE)

RNGkind("L'Ecuyer-CMRG")

empty_spectrum <- function() {
  matrix(numeric(0), ncol = 2L,
         dimnames = list(NULL, c("mz", "intensity")))
}

empty_spectrum_list <- function(ids) {
  out <- rep(list(empty_spectrum()), length(ids))
  names(out) <- ids
  out
}

prefix14 <- function(x) {
  x <- as.character(x)
  out <- substr(x, 1L, 14L)
  out[is.na(x) | !nzchar(x) | nchar(x) < 14L] <- NA_character_
  out
}

hash_seed <- function(base_seed, ...) {
  modulus <- 2147483646
  key <- paste(c(base_seed, ...), collapse = "|")
  acc <- as.double(base_seed %% modulus)
  for (value in utf8ToInt(enc2utf8(key))) {
    acc <- (acc * 131 + value) %% modulus
  }
  as.integer(acc + 1)
}

sample_fixed <- function(values, size, sample_seed, replace = FALSE) {
  if (!length(values) || size < 1L) return(values[integer()])
  set.seed(sample_seed)
  values[sample.int(length(values), size = size, replace = replace)]
}

md5_text <- function(x) {
  path <- tempfile("ppmWass-manifest-")
  on.exit(unlink(path), add = TRUE)
  writeLines(enc2utf8(x), path, useBytes = TRUE)
  unname(tools::md5sum(path))
}

git_capture <- function(extra_args) {
  out <- suppressWarnings(tryCatch(
    system2("git", c("-C", shQuote(repo_dir), extra_args), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  ))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) return(character())
  out
}

base_params <- function() {
  params <- eihrms_default_params()
  params$min_mz <- 35
  params$max_mz <- 650
  params$centroid_ppm <- 10
  params$noise_thr <- 0.01
  params$topK <- 200L
  params$class_detection_ppm <- 15
  params$tol_ppm <- tol_ppm
  params$wasserstein_transition_mult <- transition_mult
  params$ot_method <- ot_method
  params$sinkhorn_epsilon <- sinkhorn_epsilon
  params$sinkhorn_niter <- sinkhorn_niter
  params$use_typical_loss <- FALSE
  params$use_split_loss <- FALSE
  params$use_mref_confidence <- FALSE
  params$w_frag <- 0.70
  params$w_loss <- 0.30
  params$backend <- "pair_loop"
  params$use_parallel <- use_parallel
  params$n_cores <- n_cores
  params
}

params0 <- base_params()

mode_weights <- list(
  fragment_only = c(fragment = 1.00, pooled_derived_mass_difference = 0.00),
  combined_default = c(fragment = 0.70, pooled_derived_mass_difference = 0.30)
)

config <- list(
  schema_version = "ppmWass-explicit-tic-mixture-v2-exact-ot",
  datasets = datasets,
  methods = methods,
  modes = modes,
  signal_weights = signal_weights,
  replicates = replicates,
  query_n = query_n,
  interferent_components = noise_components,
  base_seed = seed,
  interferent_selection = list(
    exclude_query_itself = TRUE,
    require_different_14_character_inchikey_prefix_when_possible = TRUE,
    prefer_absolute_ri_difference_at_most = ri_window,
    sample_without_replacement_when_possible = TRUE
  ),
  mixture = list(
    signal_weight_term = "assigned pre-merge TIC mixing weight",
    component_normalization = "signal and each interferent are independently TIC-normalized",
    interferent_pool = "equal-TIC average of selected interferent spectra",
    peak_merge_ppm = merge_ppm,
    signal_peak_dropout_probability = 0,
    source_attribution = paste(
      "Source labels are retained diagnostically through the same anchor-based",
      "coalescence bins; merged-bin signal/noise contributions are the sums of",
      "their labelled pre-merge intensities."
    )
  ),
  representations = list(
    fragment_only = paste(
      "The explicit-TIC mixed fragment spectrum is scored with fragment weight 1",
      "and pooled-derived weight 0."
    ),
    combined_default = paste(
      "The pooled derived mass-difference representation is regenerated from",
      "the explicit-TIC mixed fragment spectrum, then scored 70:30 against the",
      "corresponding clean library fragment and derived representations."
    )
  ),
  tie_handling = list(
    definition = if (tie_tolerance == 0) {
      "full-precision exact equality"
    } else {
      paste0("absolute distance tolerance <= ", format(tie_tolerance, scientific = TRUE))
    },
    tolerance = tie_tolerance,
    primary_policy = "fractional_expected_under_uniform_random_order_within_tie_block",
    sensitivity_policies = c("optimistic", "pessimistic", "fixed-seed_randomized"),
    cutoffs = c(1L, 5L, 10L),
    candidate_recall_definition = "correct unperturbed library entry retained at or above cutoff K"
  ),
  nonfinite_policy = "abort; do not replace or silently discard nonfinite distances",
  ot_estimand = if (identical(ot_method, "exact"))
    "unregularized_exact_transport_cost" else
    "finite_iteration_regularized_plan_transport_cost",
  effective_base_parameters = params0,
  mode_weights = lapply(mode_weights, as.list),
  execution = list(
    use_parallel = use_parallel,
    n_cores = n_cores,
    show_progress = show_progress,
    save_distance_matrices = save_distances
  )
)

input_sources <- c(
  RECETOX = if (nzchar(recetox_spectra_rds)) recetox_spectra_rds else recetox_msp,
  `HREI-MSDB` = if (nzchar(hrei_spectra_rds)) hrei_spectra_rds else hrei_msp
)
input_sources <- input_sources[names(input_sources) %in% datasets]
input_checksums <- data.frame(
  dataset = names(input_sources),
  path = unname(input_sources),
  source_type = ifelse(grepl("[.]rds$", input_sources, ignore.case = TRUE),
                       "prepared_spectra_rds", "raw_msp"),
  md5 = unname(tools::md5sum(input_sources)),
  bytes = unname(file.info(input_sources)$size),
  stringsAsFactors = FALSE
)
# Input content and runner source are part of the resume signature; a run can
# never be resumed silently after either changes.
config$input_files <- input_checksums
config$script_md5 <- if (!is.na(script_file)) {
  unname(tools::md5sum(script_file))
} else {
  NA_character_
}
config_json <- jsonlite::toJSON(config, auto_unbox = TRUE, null = "null", digits = NA)
config_signature <- md5_text(config_json)

manifest_path <- file.path(output_dir, "manifest.json")
if (file.exists(manifest_path)) {
  if (!resume) {
    stop("Output manifest already exists. Choose a new --output-dir or pass --resume=true: ",
         manifest_path)
  }
  previous_manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  if (is.null(previous_manifest$config_signature) ||
      !identical(as.character(previous_manifest$config_signature), config_signature)) {
    stop("Cannot resume: the existing manifest has a different configuration signature.")
  }
}

utils::write.csv(input_checksums, file.path(output_dir, "input_checksums.csv"), row.names = FALSE)

git_commit <- git_capture(c("rev-parse", "HEAD"))
git_status <- git_capture(c("status", "--porcelain"))
manifest_base <- list(
  schema_version = config$schema_version,
  config_signature = config_signature,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  repository = repo_dir,
  bundle = bundle_dir,
  output_directory = output_dir,
  script = list(
    path = script_file,
    md5 = if (!is.na(script_file)) unname(tools::md5sum(script_file)) else NA_character_
  ),
  git = list(
    available = length(git_commit) > 0L,
    commit = if (length(git_commit)) git_commit[[1L]] else NA_character_,
    dirty = if (length(git_commit)) length(git_status) > 0L else NA,
    status_porcelain = git_status
  ),
  input_files = input_checksums,
  command = paste(commandArgs(), collapse = " "),
  config = config,
  seed_derivation = paste(
    "Every stochastic unit receives hash_seed(base_seed, dataset, stage, replicate,",
    "query position, condition). The polynomial character hash uses multiplier 131",
    "modulo 2,147,483,646 and returns a positive R seed."
  ),
  outputs = list(
    design = "mixture_query_interferent_design.csv",
    diagnostics = "audit/mixture_weight_diagnostics.csv",
    per_query = "mixture_per_query_metrics.csv",
    per_replicate = "mixture_per_replicate_metrics.csv",
    summary = "mixture_summary.csv",
    timings = "condition_timings.csv",
    mode_subdirectories = modes
  )
)

write_manifest <- function(status) {
  manifest <- manifest_base
  manifest$status <- status
  manifest$updated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  jsonlite::write_json(
    manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE,
    null = "null", digits = NA, na = "null"
  )
}

write_manifest("running")
run_complete <- FALSE
on.exit({
  if (!run_complete && exists("manifest_base", inherits = FALSE)) {
    try(write_manifest("incomplete"), silent = TRUE)
  }
}, add = TRUE)

log_path <- file.path(output_dir, "run.log")
log_connection <- file(log_path, open = if (resume) "at" else "wt")
sink(log_connection, split = TRUE)
sink(log_connection, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_connection)
}, add = TRUE)

message("Repository: ", repo_dir)
message("Output: ", output_dir)
message("Datasets: ", paste(datasets, collapse = ", "))
message("Methods: ", paste(methods, collapse = ", "))
message("Modes: ", paste(modes, collapse = ", "))
message("Assigned pre-merge TIC signal weights: ", paste(signal_weights, collapse = ", "))
message("Replicates: ", replicates, "; queries per dataset: ", query_n)

subset_spectra_local <- function(spectra_obj, idx) {
  out <- list(
    df_spec = spectra_obj$df_spec[idx, , drop = FALSE],
    frag_list = spectra_obj$frag_list[idx],
    loss_list = if (!is.null(spectra_obj$loss_list)) spectra_obj$loss_list[idx] else NULL,
    ri = if (!is.null(spectra_obj$ri)) spectra_obj$ri[idx] else spectra_obj$df_spec$RI[idx]
  )
  out
}

normalize_hrei_msp_names <- function(input, output) {
  lines <- readLines(input, warn = FALSE)
  starts <- grep("^Name:", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  for (i in seq_along(starts)) {
    block <- lines[starts[[i]]:ends[[i]]]
    db_line <- grep("^DB#:", block, value = TRUE)
    db_id <- if (length(db_line)) {
      trimws(sub("^DB#:[[:space:]]*", "", db_line[[1L]]))
    } else {
      as.character(i)
    }
    old_name <- trimws(sub("^Name:[[:space:]]*", "", lines[starts[[i]]]))
    lines[starts[[i]]] <- paste0("Name: ", old_name, " __HREI_DB", db_id)
  }
  writeLines(lines, output, useBytes = TRUE)
  output
}

prepare_dataset <- function(spectra_obj, params, dataset) {
  required <- c("df_spec", "frag_list")
  missing <- required[!required %in% names(spectra_obj)]
  if (length(missing)) stop(dataset, " spectra object lacks: ", paste(missing, collapse = ", "))
  df <- spectra_obj$df_spec
  if (!all(c("id", "RI", "inchikey") %in% names(df))) {
    stop(dataset, " df_spec must contain id, RI, and inchikey columns.")
  }
  valid_frag <- vapply(spectra_obj$frag_list, function(x) {
    is.matrix(x) && ncol(x) >= 2L && nrow(x) >= 2L &&
      all(is.finite(x[, 1L:2L])) && sum(x[, 2L]) > 0
  }, logical(1L))
  valid <- valid_frag & is.finite(suppressWarnings(as.numeric(df$RI))) &
    !is.na(df$inchikey) & nchar(as.character(df$inchikey)) >= 14L
  if (!any(valid)) stop(dataset, " has no valid spectra after RI/InChIKey/peak filtering.")
  out <- subset_spectra_local(spectra_obj, which(valid))
  ids <- as.character(out$df_spec$id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(dataset, " library IDs must be non-missing and unique for self-retention.")
  }
  names(out$frag_list) <- ids
  if (is.null(out$loss_list) || length(out$loss_list) != length(out$frag_list)) {
    build_loss <- getFromNamespace("build_loss_peaks", "ppmWass")
    out$loss_list <- lapply(out$frag_list, build_loss, params = params)
  }
  names(out$loss_list) <- ids
  out$ri <- suppressWarnings(as.numeric(out$df_spec$RI))
  names(out$ri) <- ids
  out
}

load_dataset <- function(dataset, params) {
  prepared <- if (dataset == "RECETOX") recetox_spectra_rds else hrei_spectra_rds
  if (nzchar(prepared)) {
    message("[", dataset, "] loading prepared spectra RDS (explicit override): ", prepared)
    return(prepare_dataset(readRDS(prepared), params, dataset))
  }
  if (dataset == "RECETOX") {
    message("[RECETOX] parsing and preprocessing ", recetox_msp)
    raw <- build_spectra_from_msp(recetox_msp, params, require_ri = FALSE, progress = show_progress)
  } else {
    prepared_dir <- file.path(output_dir, "prepared_inputs")
    dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)
    unique_msp <- file.path(prepared_dir, "HREI-MSDB_unique_names.msp")
    normalize_hrei_msp_names(hrei_msp, unique_msp)
    message("[HREI-MSDB] parsing and preprocessing ", unique_msp)
    raw <- build_spectra_from_msp(unique_msp, params, require_ri = FALSE, progress = show_progress)
  }
  prepare_dataset(raw, params, dataset)
}

select_design <- function(spectra, dataset, dataset_index) {
  n_library <- length(spectra$frag_list)
  n_select <- min(query_n, n_library)
  query_seed <- hash_seed(seed, dataset, "query_selection")
  query_indices <- sort(sample_fixed(seq_len(n_library), n_select, query_seed))
  ids <- names(spectra$frag_list)
  prefixes <- prefix14(spectra$df_spec$inchikey)
  ri <- spectra$ri
  plan <- vector("list", replicates)
  plan_rows <- list()
  row_number <- 0L

  for (replicate in seq_len(replicates)) {
    plan[[replicate]] <- vector("list", length(query_indices))
    for (query_position in seq_along(query_indices)) {
      query_index <- query_indices[[query_position]]
      interferent_seed <- hash_seed(
        seed, dataset, "interferent_selection", replicate, query_position
      )
      eligible_self_excluded <- setdiff(seq_len(n_library), query_index)
      eligible <- eligible_self_excluded
      prefix_fallback <- FALSE
      if (!is.na(prefixes[[query_index]])) {
        different <- eligible[prefixes[eligible] != prefixes[[query_index]] |
                                is.na(prefixes[eligible])]
        if (length(different) >= min(noise_components, length(eligible_self_excluded))) {
          eligible <- different
        } else {
          prefix_fallback <- TRUE
        }
      }
      if (!length(eligible)) stop("No eligible interferent for query ", ids[[query_index]])

      set.seed(interferent_seed)
      near <- integer()
      if (is.finite(ri[[query_index]])) {
        near <- eligible[is.finite(ri[eligible]) &
                           abs(ri[eligible] - ri[[query_index]]) <= ri_window]
      }
      chosen <- integer()
      if (length(near)) {
        take_near <- min(noise_components, length(near))
        chosen <- near[sample.int(length(near), take_near, replace = FALSE)]
      }
      remaining_n <- noise_components - length(chosen)
      replacement <- FALSE
      if (remaining_n > 0L) {
        remaining <- setdiff(eligible, chosen)
        if (length(remaining) >= remaining_n) {
          chosen <- c(chosen, remaining[sample.int(length(remaining), remaining_n,
                                                   replace = FALSE)])
        } else {
          replacement <- TRUE
          pool <- if (length(remaining)) remaining else eligible
          chosen <- c(chosen, pool[sample.int(length(pool), remaining_n, replace = TRUE)])
        }
      }

      mix_seed <- hash_seed(seed, dataset, "mixture", replicate, query_position)
      plan[[replicate]][[query_position]] <- list(
        query_position = query_position,
        query_index = query_index,
        interferent_indices = chosen,
        query_selection_seed = query_seed,
        interferent_seed = interferent_seed,
        mixture_seed = mix_seed
      )
      row_number <- row_number + 1L
      plan_rows[[row_number]] <- data.frame(
        dataset = dataset,
        replicate = replicate,
        query_position = query_position,
        query_library_index = query_index,
        query_id = ids[[query_index]],
        query_inchikey = as.character(spectra$df_spec$inchikey[[query_index]]),
        query_inchikey_prefix = prefixes[[query_index]],
        query_ri = ri[[query_index]],
        interferent_library_indices = paste(chosen, collapse = ";"),
        interferent_ids = paste(ids[chosen], collapse = ";"),
        interferent_inchikey_prefixes = paste(prefixes[chosen], collapse = ";"),
        interferent_ris = paste(format(ri[chosen], digits = 16), collapse = ";"),
        n_eligible_interferents = length(eligible),
        n_ri_window_interferents = length(near),
        n_chosen_within_ri_window = sum(chosen %in% near),
        prefix_exclusion_fallback = prefix_fallback,
        sampling_with_replacement = replacement,
        query_selection_seed = query_seed,
        interferent_selection_seed = interferent_seed,
        mixture_seed = mix_seed,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    query_indices = query_indices,
    plan = plan,
    rows = do.call(rbind, plan_rows)
  )
}

source_attribution <- function(signal, noise, signal_weight, merge_ppm) {
  normalize_tic <- getFromNamespace("normalize_tic_spectrum", "ppmWass")
  signal <- normalize_tic(signal)
  noise <- normalize_tic(noise)
  signal_rows <- if (nrow(signal) && signal_weight > 0) {
    data.frame(mz = signal[, 1L], intensity = signal[, 2L] * signal_weight,
               source = "signal", stringsAsFactors = FALSE)
  } else data.frame(mz = numeric(), intensity = numeric(), source = character())
  noise_rows <- if (nrow(noise) && signal_weight < 1) {
    data.frame(mz = noise[, 1L], intensity = noise[, 2L] * (1 - signal_weight),
               source = "interferent", stringsAsFactors = FALSE)
  } else data.frame(mz = numeric(), intensity = numeric(), source = character())
  x <- rbind(signal_rows, noise_rows)
  if (!nrow(x)) {
    return(c(
      source_labelled_signal_tic_postmerge = NA_real_,
      source_labelled_interferent_tic_postmerge = NA_real_,
      n_source_coalesced_bins = 0,
      n_bins_with_both_sources = 0
    ))
  }
  x <- x[order(x$mz), , drop = FALSE]
  if (merge_ppm <= 0) {
    total <- sum(x$intensity)
    return(c(
      source_labelled_signal_tic_postmerge =
        sum(x$intensity[x$source == "signal"]) / total,
      source_labelled_interferent_tic_postmerge =
        sum(x$intensity[x$source == "interferent"]) / total,
      n_source_coalesced_bins = nrow(x),
      n_bins_with_both_sources = 0
    ))
  }
  signal_bin <- numeric(nrow(x))
  noise_bin <- numeric(nrow(x))
  k <- 0L
  i <- 1L
  while (i <= nrow(x)) {
    anchor <- x$mz[[i]]
    tolerance <- anchor * merge_ppm * 1e-6
    j <- i
    while (j <= nrow(x) && abs(x$mz[[j]] - anchor) <= tolerance) j <- j + 1L
    idx <- i:(j - 1L)
    k <- k + 1L
    signal_bin[[k]] <- sum(x$intensity[idx][x$source[idx] == "signal"])
    noise_bin[[k]] <- sum(x$intensity[idx][x$source[idx] == "interferent"])
    i <- j
  }
  signal_bin <- signal_bin[seq_len(k)]
  noise_bin <- noise_bin[seq_len(k)]
  total <- sum(signal_bin) + sum(noise_bin)
  c(
    source_labelled_signal_tic_postmerge = sum(signal_bin) / total,
    source_labelled_interferent_tic_postmerge = sum(noise_bin) / total,
    n_source_coalesced_bins = k,
    n_bins_with_both_sources = sum(signal_bin > 0 & noise_bin > 0)
  )
}

build_mixed_queries <- function(spectra, design, dataset, replicate, signal_weight) {
  combine_noise <- getFromNamespace("combine_noise_spectra", "ppmWass")
  mix_impl <- getFromNamespace("mix_signal_noise_impl", "ppmWass")
  build_loss <- getFromNamespace("build_loss_peaks", "ppmWass")
  ids <- names(spectra$frag_list)
  query_ids <- ids[design$query_indices]
  fragment <- vector("list", length(query_ids))
  derived <- vector("list", length(query_ids))
  diagnostics <- vector("list", length(query_ids))
  names(fragment) <- names(derived) <- query_ids

  for (query_position in seq_along(query_ids)) {
    item <- design$plan[[replicate]][[query_position]]
    signal <- spectra$frag_list[[item$query_index]]
    noise <- combine_noise(spectra$frag_list[item$interferent_indices], merge_ppm)
    mixed <- mix_impl(
      signal = signal,
      noise = noise,
      signal_fraction = signal_weight,
      dropout_prob = 0,
      merge_ppm = merge_ppm,
      seed = item$mixture_seed
    )
    if (!nrow(mixed$spectrum) || abs(sum(mixed$spectrum[, 2L]) - 1) > 1e-10) {
      stop("Invalid TIC-normalized mixture for ", dataset, "/", query_ids[[query_position]])
    }
    fragment[[query_position]] <- mixed$spectrum
    derived[[query_position]] <- build_loss(mixed$spectrum, params0)
    attribution <- source_attribution(signal, noise, signal_weight, merge_ppm)
    d <- as.list(mixed$diagnostics)
    requested_weight <- d$requested_premerge_signal_tic_weight
    if (is.null(requested_weight)) requested_weight <- signal_weight
    assigned_weight <- d$assigned_premerge_signal_tic_weight
    if (is.null(assigned_weight)) assigned_weight <- d$realized_signal_fraction
    package_source_weight <- d$source_labelled_postmerge_signal_tic_fraction
    if (is.null(package_source_weight)) package_source_weight <- NA_real_
    diagnostics[[query_position]] <- data.frame(
      dataset = dataset,
      replicate = replicate,
      query_position = query_position,
      query_id = query_ids[[query_position]],
      interferent_ids = paste(ids[item$interferent_indices], collapse = ";"),
      query_selection_seed = item$query_selection_seed,
      interferent_selection_seed = item$interferent_seed,
      mixture_seed = item$mixture_seed,
      requested_premerge_signal_tic_weight = requested_weight,
      assigned_premerge_signal_tic_weight = assigned_weight,
      assigned_premerge_interferent_tic_weight = 1 - assigned_weight,
      package_source_labelled_signal_tic_postmerge = package_source_weight,
      source_labelled_signal_tic_postmerge = attribution[["source_labelled_signal_tic_postmerge"]],
      source_labelled_interferent_tic_postmerge = attribution[["source_labelled_interferent_tic_postmerge"]],
      n_bins_with_both_sources = attribution[["n_bins_with_both_sources"]],
      n_signal_peaks_before = d$n_signal_peaks_before,
      n_signal_peaks_after = d$n_signal_peaks_after,
      n_interferent_pool_peaks = d$n_noise_peaks,
      n_mixed_fragment_peaks = d$n_mixture_peaks,
      n_regenerated_pooled_derived_peaks = nrow(derived[[query_position]]),
      modes_using_identical_mixed_fragment = paste(modes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  list(fragment = fragment, derived = derived,
       diagnostics = do.call(rbind, diagnostics))
}

per_query_tie_metrics <- function(distance_matrix, spectra, design, dataset,
                                  mode, method, signal_weight, replicate) {
  if (any(!is.finite(distance_matrix))) {
    stop(
      "Publication mixture run aborted: ", dataset, "/", mode, "/", method,
      "/signal=", signal_weight, "/replicate=", replicate,
      " produced ", sum(!is.finite(distance_matrix)), " nonfinite distance cells."
    )
  }
  ids <- names(spectra$frag_list)
  rows <- vector("list", nrow(distance_matrix))
  n_library <- ncol(distance_matrix)
  for (query_position in seq_len(nrow(distance_matrix))) {
    item <- design$plan[[replicate]][[query_position]]
    relevant <- rep(FALSE, n_library)
    relevant[[item$query_index]] <- TRUE
    metrics <- first_relevant_tie_metrics(
      distance_matrix[query_position, ], relevant,
      ks = c(1L, 5L, 10L), tolerance = tie_tolerance
    )
    tie_random_seed <- hash_seed(
      seed, dataset, "tie_randomization", mode, method,
      format(signal_weight, digits = 17), replicate, query_position
    )
    set.seed(tie_random_seed)
    random_rank <- random_first_relevant_rank(
      metrics$n_better, metrics$tie_size, metrics$relevant_in_tie
    )
    plan_item <- design$plan[[replicate]][[query_position]]
    rows[[query_position]] <- data.frame(
      dataset = dataset,
      mode = mode,
      method = method,
      signal_weight = signal_weight,
      interferent_weight = 1 - signal_weight,
      replicate = replicate,
      query_position = query_position,
      query_library_index = item$query_index,
      query_id = ids[[item$query_index]],
      query_inchikey = as.character(spectra$df_spec$inchikey[[item$query_index]]),
      query_inchikey_prefix = prefix14(spectra$df_spec$inchikey[[item$query_index]]),
      interferent_library_indices = paste(plan_item$interferent_indices, collapse = ";"),
      interferent_ids = paste(ids[plan_item$interferent_indices], collapse = ";"),
      query_selection_seed = item$query_selection_seed,
      interferent_selection_seed = item$interferent_seed,
      mixture_seed = item$mixture_seed,
      tie_random_seed = tie_random_seed,
      n_library = n_library,
      correct_distance = metrics$best_distance,
      n_strictly_better = metrics$n_better,
      tie_size = metrics$tie_size,
      relevant_in_tie = metrics$relevant_in_tie,
      rank_optimistic = metrics$optimistic_rank,
      rank_pessimistic = metrics$pessimistic_rank,
      rank_random = random_rank,
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
      mrr_optimistic = 1 / metrics$optimistic_rank,
      mrr_fractional = metrics$expected_rr,
      mrr_pessimistic = 1 / metrics$pessimistic_rank,
      mrr_random = 1 / random_rank,
      tie_crosses_top1 = metrics$optimistic_rank <= 1L & metrics$pessimistic_rank > 1L,
      tie_crosses_top5 = metrics$optimistic_rank <= 5L & metrics$pessimistic_rank > 5L,
      tie_crosses_top10 = metrics$optimistic_rank <= 10L & metrics$pessimistic_rank > 10L,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

metric_columns <- list(
  top1 = c(
    optimistic = "top1_optimistic", fractional_expected = "top1_fractional",
    pessimistic = "top1_pessimistic", fixed_seed_randomized = "top1_random"
  ),
  top5 = c(
    optimistic = "top5_optimistic", fractional_expected = "top5_fractional",
    pessimistic = "top5_pessimistic", fixed_seed_randomized = "top5_random"
  ),
  top10 = c(
    optimistic = "top10_optimistic", fractional_expected = "top10_fractional",
    pessimistic = "top10_pessimistic", fixed_seed_randomized = "top10_random"
  ),
  mrr = c(
    optimistic = "mrr_optimistic", fractional_expected = "mrr_fractional",
    pessimistic = "mrr_pessimistic", fixed_seed_randomized = "mrr_random"
  )
)

summarize_replicates <- function(per_query) {
  group_cols <- c("dataset", "mode", "method", "signal_weight",
                  "interferent_weight", "replicate")
  group_key <- interaction(per_query[group_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(per_query, group_key)
  rows <- list()
  row_number <- 0L
  for (group in groups) {
    n_library <- unique(group$n_library)
    if (length(n_library) != 1L) stop("Internal error: inconsistent library sizes.")
    for (metric in names(metric_columns)) {
      cutoff <- switch(metric, top1 = 1L, top5 = 5L, top10 = 10L, mrr = NA_integer_)
      random_baseline <- if (metric == "mrr") {
        sum(1 / seq_len(n_library)) / n_library
      } else {
        min(cutoff, n_library) / n_library
      }
      for (policy in names(metric_columns[[metric]])) {
        value <- group[[metric_columns[[metric]][[policy]]]]
        row_number <- row_number + 1L
        rows[[row_number]] <- data.frame(
          dataset = group$dataset[[1L]],
          mode = group$mode[[1L]],
          method = group$method[[1L]],
          signal_weight = group$signal_weight[[1L]],
          interferent_weight = group$interferent_weight[[1L]],
          replicate = group$replicate[[1L]],
          metric = metric,
          candidate_recall_cutoff = cutoff,
          tie_policy = policy,
          primary = policy == "fractional_expected",
          estimate = mean(value, na.rm = TRUE),
          query_sd = if (sum(is.finite(value)) > 1L) stats::sd(value, na.rm = TRUE) else 0,
          n_queries = sum(is.finite(value)),
          n_library = n_library,
          random_ranking_baseline = random_baseline,
          difference_above_random = mean(value, na.rm = TRUE) - random_baseline,
          fold_over_random = mean(value, na.rm = TRUE) / random_baseline,
          n_any_tie = sum(group$tie_size > 1L, na.rm = TRUE),
          n_ties_crossing_top1 = sum(group$tie_crosses_top1, na.rm = TRUE),
          n_ties_crossing_top5 = sum(group$tie_crosses_top5, na.rm = TRUE),
          n_ties_crossing_top10 = sum(group$tie_crosses_top10, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

summarize_across_replicates <- function(per_replicate) {
  group_cols <- c("dataset", "mode", "method", "signal_weight", "interferent_weight",
                  "metric", "candidate_recall_cutoff", "tie_policy", "primary")
  grouping_data <- per_replicate[group_cols]
  grouping_data$candidate_recall_cutoff[
    is.na(grouping_data$candidate_recall_cutoff)
  ] <- -1L
  group_key <- interaction(grouping_data, drop = TRUE, lex.order = TRUE)
  groups <- split(per_replicate, group_key)
  rows <- lapply(groups, function(group) {
    value <- group$estimate
    data.frame(
      dataset = group$dataset[[1L]],
      mode = group$mode[[1L]],
      method = group$method[[1L]],
      signal_weight = group$signal_weight[[1L]],
      interferent_weight = group$interferent_weight[[1L]],
      metric = group$metric[[1L]],
      candidate_recall_cutoff = group$candidate_recall_cutoff[[1L]],
      tie_policy = group$tie_policy[[1L]],
      primary = group$primary[[1L]],
      estimate = mean(value, na.rm = TRUE),
      replicate_sd = if (sum(is.finite(value)) > 1L) stats::sd(value, na.rm = TRUE) else 0,
      replicate_se = if (sum(is.finite(value)) > 1L) {
        stats::sd(value, na.rm = TRUE) / sqrt(sum(is.finite(value)))
      } else 0,
      replicate_min = min(value, na.rm = TRUE),
      replicate_max = max(value, na.rm = TRUE),
      n_replicates = sum(is.finite(value)),
      n_queries_per_replicate = unique(group$n_queries)[[1L]],
      n_library = unique(group$n_library)[[1L]],
      random_ranking_baseline = unique(group$random_ranking_baseline)[[1L]],
      difference_above_random = mean(value, na.rm = TRUE) -
        unique(group$random_ranking_baseline)[[1L]],
      fold_over_random = mean(value, na.rm = TRUE) /
        unique(group$random_ranking_baseline)[[1L]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

append_rows <- function(existing, new_rows) {
  if (!nrow(existing)) new_rows else rbind(existing, new_rows)
}

checkpoint_per_query <- file.path(output_dir, ".checkpoint_per_query.csv")
checkpoint_timings <- file.path(output_dir, ".checkpoint_timings.csv")
per_query_all <- if (resume && file.exists(checkpoint_per_query)) {
  utils::read.csv(checkpoint_per_query, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
timing_all <- if (resume && file.exists(checkpoint_timings)) {
  utils::read.csv(checkpoint_timings, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()

condition_complete <- function(dataset, mode, method, signal_weight, replicate, n_query) {
  if (!nrow(per_query_all)) return(FALSE)
  hit <- per_query_all$dataset == dataset & per_query_all$mode == mode &
    per_query_all$method == method & per_query_all$replicate == replicate &
    abs(per_query_all$signal_weight - signal_weight) <= 1e-12
  sum(hit, na.rm = TRUE) == n_query
}

design_rows_all <- list()
diagnostics_rows_all <- list()
dataset_qc_rows <- list()

for (dataset_index in seq_along(datasets)) {
  dataset <- datasets[[dataset_index]]
  message("\n=== ", dataset, " ===")
  spectra <- load_dataset(dataset, params0)
  design <- select_design(spectra, dataset, dataset_index)
  design_rows_all[[dataset]] <- design$rows
  dataset_qc_rows[[dataset]] <- data.frame(
    dataset = dataset,
    n_library = length(spectra$frag_list),
    n_queries = length(design$query_indices),
    n_replicates = replicates,
    n_interferents_per_query = noise_components,
    stringsAsFactors = FALSE
  )
  diagnostic_index <- 0L
  dataset_diagnostics <- list()

  for (replicate in seq_len(replicates)) {
    for (signal_weight in signal_weights) {
      mixed <- build_mixed_queries(
        spectra, design, dataset, replicate, signal_weight
      )
      diagnostic_index <- diagnostic_index + 1L
      dataset_diagnostics[[diagnostic_index]] <- mixed$diagnostics

      for (mode in modes) {
        weights <- mode_weights[[mode]]
        query_loss <- if (mode == "combined_default") {
          mixed$derived
        } else {
          empty_spectrum_list(names(mixed$fragment))
        }
        library_loss <- if (mode == "combined_default") {
          spectra$loss_list
        } else {
          empty_spectrum_list(names(spectra$frag_list))
        }

        for (method in methods) {
          if (condition_complete(
              dataset, mode, method, signal_weight, replicate,
              length(design$query_indices))) {
            message("[resume] skipping completed ", dataset, "/", mode, "/", method,
                    "/signal=", signal_weight, "/replicate=", replicate)
            next
          }
          message("[score] ", dataset, " | ", mode, " | ", method,
                  " | assigned signal TIC=", signal_weight,
                  " | replicate=", replicate)
          params <- params0
          params$distance_method <- method
          params$w_frag <- unname(weights[["fragment"]])
          params$w_loss <- unname(weights[["pooled_derived_mass_difference"]])
          params <- getFromNamespace("validate_params", "ppmWass")(params)
          started <- proc.time()[["elapsed"]]
          distance_matrix <- compute_distance_matrix_search(
            query_frag_list = mixed$fragment,
            query_loss_list = query_loss,
            lib_frag_list = spectra$frag_list,
            lib_loss_list = library_loss,
            params = params,
            progress = show_progress
          )
          elapsed <- proc.time()[["elapsed"]] - started
          per_query <- per_query_tie_metrics(
            distance_matrix, spectra, design, dataset, mode, method,
            signal_weight, replicate
          )
          per_query_all <- append_rows(per_query_all, per_query)
          timing_all <- append_rows(timing_all, data.frame(
            dataset = dataset,
            mode = mode,
            method = method,
            signal_weight = signal_weight,
            replicate = replicate,
            n_query = nrow(distance_matrix),
            n_library = ncol(distance_matrix),
            elapsed_seconds = elapsed,
            stringsAsFactors = FALSE
          ))
          utils::write.csv(per_query_all, checkpoint_per_query, row.names = FALSE)
          utils::write.csv(timing_all, checkpoint_timings, row.names = FALSE)
          if (save_distances) {
            weight_label <- gsub("[.]", "p", format(signal_weight, trim = TRUE))
            saveRDS(
              distance_matrix,
              file.path(
                output_dir, "distance_matrices",
                paste(dataset, mode, method, paste0("signal_", weight_label),
                      paste0("rep_", replicate), sep = "__") |> paste0(".rds")
              )
            )
          }
        }
      }
    }
  }
  diagnostics_rows_all[[dataset]] <- do.call(rbind, dataset_diagnostics)
}

design_all <- do.call(rbind, design_rows_all)
diagnostics_all <- do.call(rbind, diagnostics_rows_all)
dataset_qc <- do.call(rbind, dataset_qc_rows)

sort_per_query <- order(
  match(per_query_all$dataset, datasets), match(per_query_all$mode, modes),
  match(per_query_all$method, methods), -per_query_all$signal_weight,
  per_query_all$replicate, per_query_all$query_position
)
per_query_all <- per_query_all[sort_per_query, , drop = FALSE]
per_replicate <- summarize_replicates(per_query_all)
summary_all <- summarize_across_replicates(per_replicate)

utils::write.csv(design_all, file.path(output_dir, "mixture_query_interferent_design.csv"), row.names = FALSE)
utils::write.csv(diagnostics_all, file.path(audit_dir, "mixture_weight_diagnostics.csv"), row.names = FALSE)
utils::write.csv(dataset_qc, file.path(output_dir, "dataset_qc.csv"), row.names = FALSE)
utils::write.csv(per_query_all, file.path(output_dir, "mixture_per_query_metrics.csv"), row.names = FALSE)
utils::write.csv(per_replicate, file.path(output_dir, "mixture_per_replicate_metrics.csv"), row.names = FALSE)
utils::write.csv(summary_all, file.path(output_dir, "mixture_summary.csv"), row.names = FALSE)
utils::write.csv(timing_all, file.path(output_dir, "condition_timings.csv"), row.names = FALSE)

effective_rows <- do.call(rbind, lapply(modes, function(mode) {
  do.call(rbind, lapply(methods, function(method) {
    data.frame(
      mode = mode,
      method = method,
      fragment_weight = mode_weights[[mode]][["fragment"]],
      pooled_derived_mass_difference_weight =
        mode_weights[[mode]][["pooled_derived_mass_difference"]],
      min_mz = params0$min_mz,
      max_mz = params0$max_mz,
      centroid_ppm = params0$centroid_ppm,
      noise_threshold_relative_to_base_peak = params0$noise_thr,
      top_k_peaks = params0$topK,
      hard_match_tolerance_ppm = params0$tol_ppm,
      ppmWass_base_cost_ppm = params0$tol_ppm,
      transition_multiplier = params0$wasserstein_transition_mult,
      saturation_width_ppm = params0$tol_ppm * params0$wasserstein_transition_mult,
      solver_backend = params0$ot_method,
      ot_estimand = if (identical(params0$ot_method, "exact"))
        "unregularized_exact_transport_cost" else
        "finite_iteration_regularized_plan_transport_cost",
      sinkhorn_epsilon = if (identical(params0$ot_method, "exact"))
        NA_real_ else params0$sinkhorn_epsilon,
      sinkhorn_iterations = if (identical(params0$ot_method, "exact"))
        NA_integer_ else params0$sinkhorn_niter,
      stringsAsFactors = FALSE
    )
  }))
}))
utils::write.csv(effective_rows, file.path(output_dir, "effective_parameters.csv"), row.names = FALSE)

for (mode in modes) {
  mode_dir <- file.path(output_dir, mode)
  utils::write.csv(
    per_query_all[per_query_all$mode == mode, , drop = FALSE],
    file.path(mode_dir, "per_query_metrics.csv"), row.names = FALSE
  )
  utils::write.csv(
    per_replicate[per_replicate$mode == mode, , drop = FALSE],
    file.path(mode_dir, "per_replicate_metrics.csv"), row.names = FALSE
  )
  utils::write.csv(
    summary_all[summary_all$mode == mode, , drop = FALSE],
    file.path(mode_dir, "summary.csv"), row.names = FALSE
  )
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"), useBytes = TRUE)
write_manifest("complete")
run_complete <- TRUE
message("Explicit-TIC mixture benchmark complete: ", output_dir)
