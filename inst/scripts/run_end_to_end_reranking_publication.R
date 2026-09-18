#!/usr/bin/env Rscript

# Publication-grade end-to-end Top-K candidate reranking benchmark.
#
# The full first-pass query-by-library search, Top-K extraction, candidate-only
# ppm-Wasserstein scoring, and result formatting are executed through
# search_topk_rerank(). Input loading and preprocessing are timed separately.
# First-pass baseline and final Top-1/MRR use the fractional expectation under
# uniform ordering within an exact full-precision distance-tie block.

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
first_pass_methods <- unique(split_arg(
  "first-pass-methods", "entropy_weighted,composite,hellinger,cosine"
))
valid_first_pass <- c("entropy_weighted", "composite", "hellinger", "cosine")
bad_methods <- setdiff(first_pass_methods, valid_first_pass)
if (length(bad_methods)) {
  stop("Unsupported first-pass method(s): ", paste(bad_methods, collapse = ", "))
}

k_values <- unique(as.integer(split_arg("k-values", "5,10,20,50")))
if (!length(k_values) || any(!is.finite(k_values)) || any(k_values < 1L)) {
  stop("k-values must be positive integers.")
}
k_values <- sort(k_values)

execution_modes <- unique(split_arg("execution-modes", "serial,parallel_8core"))
valid_execution_modes <- c("serial", "parallel_8core")
bad_execution <- setdiff(execution_modes, valid_execution_modes)
if (length(bad_execution)) {
  stop("Unsupported execution mode(s): ", paste(bad_execution, collapse = ", "))
}

parallel_cores <- as.integer(get_arg("parallel-cores", "8"))
query_limit <- as.integer(get_arg("query-limit", "0"))
seed <- as.integer(get_arg("seed", "20260718"))
tol_ppm <- as.numeric(get_arg("tol-ppm", "15"))
transition_mult <- as.numeric(get_arg("transition-mult", "3"))
ot_method <- tolower(get_arg("ot-method", "exact"))
sinkhorn_epsilon <- as.numeric(get_arg("sinkhorn-epsilon", "0.05"))
sinkhorn_niter <- as.integer(get_arg("sinkhorn-niter", "100"))
tie_tolerance <- as.numeric(get_arg("tie-tolerance", "0"))
show_progress <- flag_is_true("progress", FALSE)
resume <- flag_is_true("resume", FALSE)

if (!is.finite(parallel_cores) || parallel_cores < 1L) {
  stop("parallel-cores must be >= 1.")
}
if (!is.finite(query_limit) || query_limit < 0L) stop("query-limit must be >= 0.")
if (!is.finite(seed) || seed < 1L) stop("seed must be a positive integer.")
if (!is.finite(tol_ppm) || tol_ppm <= 0) stop("tol-ppm must be positive.")
if (!is.finite(transition_mult) || transition_mult <= 0) stop("transition-mult must be positive.")
if (!ot_method %in% c("exact", "sinkhorn", "greenkhorn")) {
  stop("ot-method must be exact, sinkhorn, or greenkhorn.")
}
if (!is.finite(sinkhorn_epsilon) || sinkhorn_epsilon <= 0) {
  stop("sinkhorn-epsilon must be positive.")
}
if (!is.finite(sinkhorn_niter) || sinkhorn_niter < 1L) {
  stop("sinkhorn-niter must be >= 1.")
}
if (!is.finite(tie_tolerance) || tie_tolerance < 0) {
  stop("tie-tolerance must be non-negative.")
}

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
  get_arg(
    "output-dir",
    file.path(bundle_dir, "results", "current", "end_to_end_reranking")
  ),
  mustWork = FALSE
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

stats_helper <- file.path(
  repo_dir, "inst", "scripts", "lib", "publication_retrieval_statistics.R"
)
if (!file.exists(stats_helper)) stop("Missing statistics helper: ", stats_helper)
source(stats_helper, local = TRUE)
RNGkind("L'Ecuyer-CMRG")

hash_seed <- function(base_seed, ...) {
  modulus <- 2147483646
  key <- paste(c(base_seed, ...), collapse = "|")
  acc <- as.double(base_seed %% modulus)
  for (value in utf8ToInt(enc2utf8(key))) {
    acc <- (acc * 131 + value) %% modulus
  }
  as.integer(acc + 1)
}

md5_text <- function(x) {
  path <- tempfile("ppmWass-rerank-manifest-")
  on.exit(unlink(path), add = TRUE)
  writeLines(enc2utf8(x), path, useBytes = TRUE)
  unname(tools::md5sum(path))
}

git_capture <- function(extra_args) {
  out <- suppressWarnings(tryCatch(
    system2("git", c("-C", shQuote(repo_dir), extra_args),
            stdout = TRUE, stderr = TRUE),
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
  params$use_parallel <- FALSE
  params$n_cores <- 1L
  params
}

params0 <- base_params()
resolved_parallel_cores <- getFromNamespace("resolve_parallel_cores", "ppmWass")(
  parallel_cores
)

input_sources <- c(
  RECETOX = if (nzchar(recetox_spectra_rds)) recetox_spectra_rds else recetox_msp,
  `HREI-MSDB` = if (nzchar(hrei_spectra_rds)) hrei_spectra_rds else hrei_msp
)
input_sources <- input_sources[names(input_sources) %in% datasets]
input_checksums <- data.frame(
  dataset = names(input_sources),
  path = unname(input_sources),
  source_type = ifelse(
    grepl("[.]rds$", input_sources, ignore.case = TRUE),
    "prepared_spectra_rds", "raw_msp"
  ),
  md5 = unname(tools::md5sum(input_sources)),
  bytes = unname(file.info(input_sources)$size),
  stringsAsFactors = FALSE
)

config <- list(
  schema_version = "ppmWass-end-to-end-reranking-v3-exact-ot",
  datasets = datasets,
  query_definition = paste(
    "All spectra whose full InChIKey occurs more than once; the identical",
    "query ID is excluded from the library and other full-InChIKey entries",
    "are relevant. query_limit=0 uses all eligible queries."
  ),
  query_limit = query_limit,
  first_pass_methods = first_pass_methods,
  candidate_k = k_values,
  rerank_method = "ppm_wasserstein",
  execution_modes = execution_modes,
  parallel_cores_requested = parallel_cores,
  parallel_cores_effective_at_start = resolved_parallel_cores,
  combined_representation_weights = list(
    fragment = 0.70,
    pooled_derived_mass_difference = 0.30
  ),
  tie_handling = list(
    definition = if (tie_tolerance == 0) {
      "full-precision exact equality"
    } else {
      paste0("absolute distance tolerance <= ", format(tie_tolerance, scientific = TRUE))
    },
    tolerance = tie_tolerance,
    primary = "fractional_expected_under_uniform_random_order_within_tie_block",
    sensitivity = c("optimistic", "pessimistic", "fixed-seed_randomized")
  ),
  candidate_selection_ties = paste(
    "search_topk_rerank() uses deterministic library-index order when a",
    "first-pass tie crosses the Top-K boundary; candidate retention is the",
    "query-level indicator that at least one relevant entry is retained."
  ),
  candidate_retention_definition = paste(
    "Binary query-level hit: at least one relevant full-InChIKey library entry",
    "is present in the realized Top-K candidate set; this is not item-level recall."
  ),
  first_pass_baseline = list(
    search_scope = paste(
      "Full query-by-library first-pass distance matrix after excluding the",
      "identical query ID; the baseline is independent of candidate K."
    ),
    relevance = "full InChIKey equality among the self-excluded library",
    metrics = paste(
      "Top-1 and MRR from the first relevant full-InChIKey entry; fractional",
      "expectation within the configured full-precision distance-tie block,",
      "with optimistic and pessimistic sensitivity values."
    )
  ),
  ot_estimand = if (identical(ot_method, "exact"))
    "unregularized_exact_transport_cost" else
    "finite_iteration_regularized_plan_transport_cost",
  candidate_fraction = "effective_K / (N - 1), where the identical query ID is excluded",
  timing = list(
    input_loading = "MSP parsing (and HREI name normalization) or readRDS",
    preprocessing = paste(
      "MSP conversion, fragment/pooled-derived construction, and eligibility",
      "filtering; for an explicit RDS override, validation/filtering only"
    ),
    query_setup = "duplicate-full-InChIKey query selection and list subsetting",
    first_pass = "search_topk_rerank() timing$first_pass_sec",
    topk_extraction = "search_topk_rerank() timing$topk_extraction_sec",
    ppm_candidate_scoring = "search_topk_rerank() timing$rerank_sec",
    result_formatting = "search_topk_rerank() timing$result_formatting_sec",
    search_total = "search_topk_rerank() timing$total_sec",
    cold_pipeline_total = paste(
      "sum of separately measured shared setup stages (input loading,",
      "preprocessing, and query setup) plus condition-specific search total;",
      "this is not a directly timed single cold execution"
    )
  ),
  memory = list(
    measurement = "gc(reset=TRUE) immediately before search; gc() max used immediately after",
    reported_value = "sum of Ncells and Vcells max-used megabytes in the parent R process",
    scope = paste(
      "R-managed parent-process heap only; excludes native allocations, OS",
      "RSS, and fork/PSOCK worker-process heaps"
    )
  ),
  nonfinite_policy = paste(
    "Stop if any first-pass distance is nonfinite after self-exclusion or if",
    "any selected ppm-Wasserstein candidate fails to return a finite score."
  ),
  effective_base_parameters = params0,
  base_seed = seed,
  input_files = input_checksums,
  script_md5 = if (!is.na(script_file)) {
    unname(tools::md5sum(script_file))
  } else {
    NA_character_
  }
)
config_json <- jsonlite::toJSON(config, auto_unbox = TRUE, null = "null", digits = NA)
config_signature <- md5_text(config_json)

manifest_path <- file.path(output_dir, "manifest.json")
if (file.exists(manifest_path)) {
  if (!resume) {
    stop(
      "Output manifest already exists. Choose a new --output-dir or pass --resume=true: ",
      manifest_path
    )
  }
  previous_manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  if (is.null(previous_manifest$config_signature) ||
      !identical(as.character(previous_manifest$config_signature), config_signature)) {
    stop("Cannot resume: the existing manifest has a different configuration signature.")
  }
}

utils::write.csv(
  input_checksums, file.path(output_dir, "input_checksums.csv"), row.names = FALSE
)
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
    md5 = config$script_md5
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
  outputs = list(
    query_design = "reranking_query_design.csv",
    dataset_stages = "dataset_stage_timings.csv",
    per_query = "reranking_per_query.csv",
    condition_summary = "reranking_condition_summary.csv",
    runtime_stages = "runtime_stage_timings.csv",
    serial_parallel_equivalence = "serial_parallel_equivalence.csv",
    effective_parameters = "effective_parameters.csv"
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
message("First-pass methods: ", paste(first_pass_methods, collapse = ", "))
message("K: ", paste(k_values, collapse = ", "))
message("Execution modes: ", paste(execution_modes, collapse = ", "))
message("Rerank method: ppm_wasserstein; combined fragment/derived weights: 70/30")

normalize_hrei_msp_names <- function(input, output) {
  lines <- readLines(input, warn = FALSE)
  starts <- grep("^Name:", lines)
  if (!length(starts)) stop("No Name records found in HREI MSP: ", input)
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

subset_spectra_local <- function(spectra_obj, idx) {
  list(
    df_spec = spectra_obj$df_spec[idx, , drop = FALSE],
    frag_list = spectra_obj$frag_list[idx],
    loss_list = spectra_obj$loss_list[idx],
    ri = if (!is.null(spectra_obj$ri)) {
      spectra_obj$ri[idx]
    } else {
      spectra_obj$df_spec$RI[idx]
    }
  )
}

prepare_dataset <- function(spectra_obj, dataset) {
  if (!all(c("df_spec", "frag_list", "loss_list") %in% names(spectra_obj))) {
    stop(dataset, " spectra object must contain df_spec, frag_list, and loss_list.")
  }
  df <- spectra_obj$df_spec
  if (!all(c("id", "RI", "inchikey") %in% names(df))) {
    stop(dataset, " df_spec must contain id, RI, and inchikey.")
  }
  if (length(spectra_obj$frag_list) != nrow(df) ||
      length(spectra_obj$loss_list) != nrow(df)) {
    stop(dataset, " spectra list lengths do not match df_spec.")
  }
  valid_frag <- vapply(spectra_obj$frag_list, function(x) {
    is.matrix(x) && ncol(x) >= 2L && nrow(x) >= 2L &&
      all(is.finite(x[, 1L:2L])) && sum(x[, 2L]) > 0
  }, logical(1L))
  valid <- valid_frag & is.finite(suppressWarnings(as.numeric(df$RI))) &
    !is.na(df$inchikey) & nchar(as.character(df$inchikey)) >= 14L
  if (!any(valid)) stop(dataset, " has no valid spectra after eligibility filtering.")
  out <- subset_spectra_local(spectra_obj, which(valid))
  ids <- as.character(out$df_spec$id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(dataset, " library IDs must be non-missing and unique.")
  }
  names(out$frag_list) <- names(out$loss_list) <- ids
  out$ri <- suppressWarnings(as.numeric(out$df_spec$RI))
  names(out$ri) <- ids
  out
}

preprocess_msp_library <- function(msp_library, params, dataset) {
  df_spec <- convert_msp_to_internal(msp_library, require_ri = FALSE)
  process_one <- getFromNamespace("process_single_spectrum", "ppmWass")
  fragment <- vector("list", nrow(df_spec))
  derived <- vector("list", nrow(df_spec))
  for (i in seq_len(nrow(df_spec))) {
    processed <- process_one(df_spec$ei[[i]], params)
    fragment[[i]] <- processed$frag
    derived[[i]] <- processed$loss
    if (show_progress && (i %% 100L == 0L || i == nrow(df_spec))) {
      message("[", dataset, "] preprocessing ", i, " / ", nrow(df_spec))
    }
  }
  names(fragment) <- names(derived) <- as.character(df_spec$id)
  prepare_dataset(
    list(
      df_spec = df_spec,
      frag_list = fragment,
      loss_list = derived,
      ri = df_spec$RI
    ),
    dataset
  )
}

load_dataset_timed <- function(dataset, params) {
  rds_path <- if (dataset == "RECETOX") recetox_spectra_rds else hrei_spectra_rds
  input_start <- proc.time()[["elapsed"]]
  if (nzchar(rds_path)) {
    message("[", dataset, "] reading prepared spectra RDS override: ", rds_path)
    loaded <- readRDS(rds_path)
    input_seconds <- proc.time()[["elapsed"]] - input_start
    preprocess_start <- proc.time()[["elapsed"]]
    spectra <- prepare_dataset(loaded, dataset)
    preprocessing_seconds <- proc.time()[["elapsed"]] - preprocess_start
    preprocessing_kind <- "prepared RDS validation and eligibility filtering"
  } else {
    msp_path <- if (dataset == "RECETOX") recetox_msp else hrei_msp
    if (dataset == "HREI-MSDB") {
      prepared_dir <- file.path(output_dir, "prepared_inputs")
      dir.create(prepared_dir, recursive = TRUE, showWarnings = FALSE)
      normalized <- file.path(prepared_dir, "HREI-MSDB_unique_names.msp")
      normalize_hrei_msp_names(msp_path, normalized)
      msp_path <- normalized
    }
    message("[", dataset, "] parsing MSP: ", msp_path)
    loaded <- read_msp(msp_path, progress = if (show_progress) 1000L else 0L)
    input_seconds <- proc.time()[["elapsed"]] - input_start
    preprocess_start <- proc.time()[["elapsed"]]
    spectra <- preprocess_msp_library(loaded, params, dataset)
    preprocessing_seconds <- proc.time()[["elapsed"]] - preprocess_start
    preprocessing_kind <- "MSP conversion, fragment/derived construction, and eligibility filtering"
  }
  list(
    spectra = spectra,
    input_loading_seconds = input_seconds,
    preprocessing_seconds = preprocessing_seconds,
    preprocessing_kind = preprocessing_kind
  )
}

build_query_set <- function(spectra, dataset) {
  start <- proc.time()[["elapsed"]]
  inchikey <- as.character(spectra$df_spec$inchikey)
  names(inchikey) <- names(spectra$frag_list)
  counts <- table(inchikey[!is.na(inchikey) & nzchar(inchikey)])
  eligible_keys <- names(counts[counts > 1L])
  query_index <- which(inchikey %in% eligible_keys)
  if (query_limit > 0L) query_index <- utils::head(query_index, query_limit)
  if (!length(query_index)) stop(dataset, " has no eligible duplicate-full-InChIKey queries.")
  query_ids <- names(spectra$frag_list)[query_index]
  query_inchikey <- inchikey[query_index]
  n_relevant <- vapply(query_inchikey, function(key) sum(inchikey == key) - 1L, integer(1L))
  if (any(n_relevant < 1L)) stop("Internal query construction error in ", dataset)
  rows <- data.frame(
    dataset = dataset,
    query_position = seq_along(query_index),
    query_library_index = query_index,
    query_id = query_ids,
    full_inchikey = unname(query_inchikey),
    inchikey_prefix = prefix_inchikey(unname(query_inchikey)),
    n_relevant_library_entries_excluding_self = n_relevant,
    stringsAsFactors = FALSE
  )
  list(
    query_index = query_index,
    query_ids = query_ids,
    query_inchikey = unname(query_inchikey),
    fragment = spectra$frag_list[query_index],
    derived = spectra$loss_list[query_index],
    rows = rows,
    setup_seconds = proc.time()[["elapsed"]] - start
  )
}

heap_peak_after_search <- function(gc_result) {
  mb_columns <- which(colnames(gc_result) == "(Mb)")
  max_mb_column <- utils::tail(mb_columns, 1L)
  c(
    r_heap_peak_ncells = unname(gc_result["Ncells", "max used"]),
    r_heap_peak_vcells = unname(gc_result["Vcells", "max used"]),
    r_heap_peak_ncells_mb = unname(gc_result["Ncells", max_mb_column]),
    r_heap_peak_vcells_mb = unname(gc_result["Vcells", max_mb_column]),
    r_heap_peak_total_mb = sum(gc_result[, max_mb_column])
  )
}

evaluate_condition <- function(search_result, spectra, queries, dataset,
                               first_pass_method, requested_k, effective_k,
                               execution_mode, effective_cores) {
  library_ids <- names(spectra$frag_list)
  library_inchikey <- as.character(spectra$df_spec$inchikey)
  names(library_inchikey) <- library_ids
  result_rows <- search_result$results
  first_pass_matrix <- search_result$first_pass_distance_matrix
  expected_first_pass_dim <- c(length(queries$query_ids), length(library_ids))
  if (!is.matrix(first_pass_matrix) || !is.numeric(first_pass_matrix) ||
      !identical(dim(first_pass_matrix), expected_first_pass_dim)) {
    stop(
      "Missing or malformed full first-pass distance matrix for ", dataset, "/",
      first_pass_method, "/K=", requested_k, "/", execution_mode, "."
    )
  }
  if (!identical(rownames(first_pass_matrix), queries$query_ids) ||
      !identical(colnames(first_pass_matrix), library_ids)) {
    stop(
      "First-pass distance matrix IDs are not aligned for ", dataset, "/",
      first_pass_method, "/K=", requested_k, "/", execution_mode, "."
    )
  }
  rows <- vector("list", length(queries$query_ids))

  for (query_position in seq_along(queries$query_ids)) {
    query_id <- queries$query_ids[[query_position]]
    query_key <- queries$query_inchikey[[query_position]]
    if (sum(library_ids == query_id) != 1L) {
      stop("Query ID must occur exactly once in the library: ", query_id)
    }
    first_pass_eligible <- library_ids != query_id
    first_pass_distances <- first_pass_matrix[query_id, first_pass_eligible]
    first_pass_relevant <-
      library_inchikey[first_pass_eligible] == query_key
    first_pass_relevant[is.na(first_pass_relevant)] <- FALSE
    if (any(!is.finite(first_pass_distances))) {
      stop(
        "Nonfinite full-library first-pass distance: ", dataset, "/",
        first_pass_method, "/K=", requested_k, "/", execution_mode, "/",
        query_id, "."
      )
    }
    if (!any(first_pass_relevant)) {
      stop("No self-excluded full-InChIKey relevant entry for query ", query_id)
    }
    first_pass_metrics <- first_relevant_tie_metrics(
      first_pass_distances, first_pass_relevant, ks = 1L,
      tolerance = tie_tolerance
    )
    first_pass_top1_optimistic <- first_pass_metrics$top1_optimistic
    first_pass_top1_fractional <- first_pass_metrics$top1_fractional
    first_pass_top1_pessimistic <- first_pass_metrics$top1_pessimistic
    first_pass_mrr_optimistic <- 1 / first_pass_metrics$optimistic_rank
    first_pass_mrr_fractional <- first_pass_metrics$expected_rr
    first_pass_mrr_pessimistic <- 1 / first_pass_metrics$pessimistic_rank

    current <- result_rows[result_rows$query_id == query_id, , drop = FALSE]
    if (nrow(current) != effective_k) {
      stop(
        "Nonfinite/short rerank result: ", dataset, "/", first_pass_method,
        "/K=", requested_k, "/", execution_mode, "/", query_id,
        " returned ", nrow(current), " of ", effective_k, " candidates."
      )
    }
    candidate_index <- match(current$lib_id, library_ids)
    if (anyNA(candidate_index)) stop("Unknown library ID returned for query ", query_id)
    relevant <- library_inchikey[candidate_index] == query_key & current$lib_id != query_id
    relevant[is.na(relevant)] <- FALSE
    candidate_retention <- as.numeric(any(relevant))
    tie_random_seed <- hash_seed(
      seed, dataset, first_pass_method, requested_k, query_position, "final_tie"
    )

    if (candidate_retention > 0) {
      metrics <- first_relevant_tie_metrics(
        current$rerank_distance, relevant, ks = 1L, tolerance = tie_tolerance
      )
      set.seed(tie_random_seed)
      random_rank <- random_first_relevant_rank(
        metrics$n_better, metrics$tie_size, metrics$relevant_in_tie
      )
      top1_optimistic <- metrics$top1_optimistic
      top1_fractional <- metrics$top1_fractional
      top1_pessimistic <- metrics$top1_pessimistic
      top1_random <- as.numeric(random_rank == 1L)
      mrr_optimistic <- 1 / metrics$optimistic_rank
      mrr_fractional <- metrics$expected_rr
      mrr_pessimistic <- 1 / metrics$pessimistic_rank
      mrr_random <- 1 / random_rank
      best_distance <- metrics$best_distance
      n_better <- metrics$n_better
      tie_size <- metrics$tie_size
      relevant_in_tie <- metrics$relevant_in_tie
      rank_optimistic <- metrics$optimistic_rank
      rank_pessimistic <- metrics$pessimistic_rank
    } else {
      random_rank <- NA_real_
      top1_optimistic <- top1_fractional <- top1_pessimistic <- top1_random <- 0
      mrr_optimistic <- mrr_fractional <- mrr_pessimistic <- mrr_random <- 0
      best_distance <- NA_real_
      n_better <- tie_size <- relevant_in_tie <- NA_integer_
      rank_optimistic <- rank_pessimistic <- NA_real_
    }

    first_pass_order <- order(current$first_pass_rank)
    rerank_order <- order(current$rerank_rank)
    rows[[query_position]] <- data.frame(
      dataset = dataset,
      first_pass_method = first_pass_method,
      rerank_method = "ppm_wasserstein",
      requested_k = requested_k,
      effective_k = effective_k,
      execution_mode = execution_mode,
      effective_cores = effective_cores,
      query_position = query_position,
      query_id = query_id,
      full_inchikey = query_key,
      n_library = length(library_ids),
      eligible_library_size_excluding_self = length(library_ids) - 1L,
      n_relevant_library_entries_excluding_self = sum(
        library_inchikey == query_key & library_ids != query_id,
        na.rm = TRUE
      ),
      n_relevant_candidates = sum(relevant),
      candidate_retention = candidate_retention,
      candidate_recall = candidate_retention,
      candidate_ids_first_pass_order = paste(current$lib_id[first_pass_order], collapse = ";"),
      candidate_ids_rerank_order = paste(current$lib_id[rerank_order], collapse = ";"),
      first_pass_best_relevant_distance = first_pass_metrics$best_distance,
      first_pass_n_strictly_better = first_pass_metrics$n_better,
      first_pass_tie_size = first_pass_metrics$tie_size,
      first_pass_relevant_in_tie = first_pass_metrics$relevant_in_tie,
      first_pass_rank_optimistic = first_pass_metrics$optimistic_rank,
      first_pass_rank_pessimistic = first_pass_metrics$pessimistic_rank,
      first_pass_top1_optimistic = first_pass_top1_optimistic,
      first_pass_top1_fractional = first_pass_top1_fractional,
      first_pass_top1_pessimistic = first_pass_top1_pessimistic,
      first_pass_mrr_optimistic = first_pass_mrr_optimistic,
      first_pass_mrr_fractional = first_pass_mrr_fractional,
      first_pass_mrr_pessimistic = first_pass_mrr_pessimistic,
      first_pass_tie_crosses_top1 = isTRUE(
        first_pass_metrics$optimistic_rank <= 1L &&
          first_pass_metrics$pessimistic_rank > 1L
      ),
      final_best_relevant_distance = best_distance,
      final_n_strictly_better = n_better,
      final_tie_size = tie_size,
      final_relevant_in_tie = relevant_in_tie,
      final_rank_optimistic = rank_optimistic,
      final_rank_pessimistic = rank_pessimistic,
      final_rank_random = random_rank,
      final_top1_optimistic = top1_optimistic,
      final_top1_fractional = top1_fractional,
      final_top1_pessimistic = top1_pessimistic,
      final_top1_random = top1_random,
      final_mrr_optimistic = mrr_optimistic,
      final_mrr_fractional = mrr_fractional,
      final_mrr_pessimistic = mrr_pessimistic,
      final_mrr_random = mrr_random,
      final_tie_crosses_top1 = isTRUE(
        is.finite(rank_optimistic) && is.finite(rank_pessimistic) &&
          rank_optimistic <= 1L && rank_pessimistic > 1L
      ),
      tie_random_seed = tie_random_seed,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_conditions <- function(per_query, timings) {
  group_columns <- c(
    "dataset", "first_pass_method", "rerank_method", "requested_k",
    "effective_k", "execution_mode", "effective_cores", "n_library",
    "eligible_library_size_excluding_self"
  )
  key <- interaction(per_query[group_columns], drop = TRUE, lex.order = TRUE)
  groups <- split(per_query, key)
  rows <- lapply(groups, function(group) {
    se <- function(x) {
      if (sum(is.finite(x)) > 1L) stats::sd(x, na.rm = TRUE) /
        sqrt(sum(is.finite(x))) else 0
    }
    data.frame(
      dataset = group$dataset[[1L]],
      first_pass_method = group$first_pass_method[[1L]],
      rerank_method = group$rerank_method[[1L]],
      requested_k = group$requested_k[[1L]],
      effective_k = group$effective_k[[1L]],
      execution_mode = group$execution_mode[[1L]],
      effective_cores = group$effective_cores[[1L]],
      n_queries = nrow(group),
      n_library = group$n_library[[1L]],
      eligible_library_size_excluding_self =
        group$eligible_library_size_excluding_self[[1L]],
      candidate_retention = mean(group$candidate_retention),
      candidate_retention_se = se(group$candidate_retention),
      candidate_recall = mean(group$candidate_retention),
      candidate_recall_se = se(group$candidate_retention),
      first_pass_top1_fractional = mean(group$first_pass_top1_fractional),
      first_pass_top1_fractional_se = se(group$first_pass_top1_fractional),
      first_pass_mrr_fractional = mean(group$first_pass_mrr_fractional),
      first_pass_mrr_fractional_se = se(group$first_pass_mrr_fractional),
      first_pass_top1_optimistic = mean(group$first_pass_top1_optimistic),
      first_pass_top1_pessimistic = mean(group$first_pass_top1_pessimistic),
      first_pass_mrr_optimistic = mean(group$first_pass_mrr_optimistic),
      first_pass_mrr_pessimistic = mean(group$first_pass_mrr_pessimistic),
      n_queries_with_first_pass_tie = sum(
        group$first_pass_tie_size > 1L, na.rm = TRUE
      ),
      n_first_pass_ties_crossing_top1 = sum(
        group$first_pass_tie_crosses_top1, na.rm = TRUE
      ),
      final_top1_fractional = mean(group$final_top1_fractional),
      final_top1_fractional_se = se(group$final_top1_fractional),
      final_mrr_fractional = mean(group$final_mrr_fractional),
      final_mrr_fractional_se = se(group$final_mrr_fractional),
      final_top1_optimistic = mean(group$final_top1_optimistic),
      final_top1_pessimistic = mean(group$final_top1_pessimistic),
      final_mrr_optimistic = mean(group$final_mrr_optimistic),
      final_mrr_pessimistic = mean(group$final_mrr_pessimistic),
      n_queries_with_final_tie = sum(group$final_tie_size > 1L, na.rm = TRUE),
      n_final_ties_crossing_top1 = sum(group$final_tie_crosses_top1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  merge(
    summary, timings,
    by = c(
      "dataset", "first_pass_method", "rerank_method", "requested_k",
      "effective_k", "execution_mode", "effective_cores", "n_queries",
      "n_library", "eligible_library_size_excluding_self"
    ),
    all.x = TRUE, sort = FALSE
  )
}

serial_parallel_audit <- function(per_query) {
  if (!all(c("serial", "parallel_8core") %in% unique(per_query$execution_mode))) {
    return(data.frame())
  }
  key_columns <- c("dataset", "first_pass_method", "requested_k", "query_id")
  value_columns <- c(
    "candidate_ids_first_pass_order", "candidate_ids_rerank_order",
    "candidate_retention",
    "first_pass_top1_optimistic", "first_pass_top1_fractional",
    "first_pass_top1_pessimistic", "first_pass_mrr_optimistic",
    "first_pass_mrr_fractional", "first_pass_mrr_pessimistic",
    "final_top1_fractional", "final_mrr_fractional",
    "tie_random_seed", "final_rank_random", "final_top1_random",
    "final_mrr_random"
  )
  serial <- per_query[per_query$execution_mode == "serial",
                      c(key_columns, value_columns), drop = FALSE]
  parallel <- per_query[per_query$execution_mode == "parallel_8core",
                        c(key_columns, value_columns), drop = FALSE]
  names(serial)[match(value_columns, names(serial))] <- paste0(value_columns, "_serial")
  names(parallel)[match(value_columns, names(parallel))] <- paste0(value_columns, "_parallel")
  joined <- merge(serial, parallel, by = key_columns, all = TRUE)
  joined$candidate_set_equal <-
    joined$candidate_ids_first_pass_order_serial ==
      joined$candidate_ids_first_pass_order_parallel
  joined$rerank_order_equal <-
    joined$candidate_ids_rerank_order_serial ==
      joined$candidate_ids_rerank_order_parallel
  joined$candidate_retention_abs_difference <- abs(
    joined$candidate_retention_serial - joined$candidate_retention_parallel
  )
  joined$candidate_recall_abs_difference <-
    joined$candidate_retention_abs_difference
  joined$first_pass_top1_optimistic_abs_difference <- abs(
    joined$first_pass_top1_optimistic_serial -
      joined$first_pass_top1_optimistic_parallel
  )
  joined$first_pass_top1_fractional_abs_difference <- abs(
    joined$first_pass_top1_fractional_serial -
      joined$first_pass_top1_fractional_parallel
  )
  joined$first_pass_top1_pessimistic_abs_difference <- abs(
    joined$first_pass_top1_pessimistic_serial -
      joined$first_pass_top1_pessimistic_parallel
  )
  joined$first_pass_mrr_optimistic_abs_difference <- abs(
    joined$first_pass_mrr_optimistic_serial -
      joined$first_pass_mrr_optimistic_parallel
  )
  joined$first_pass_mrr_fractional_abs_difference <- abs(
    joined$first_pass_mrr_fractional_serial -
      joined$first_pass_mrr_fractional_parallel
  )
  joined$first_pass_mrr_pessimistic_abs_difference <- abs(
    joined$first_pass_mrr_pessimistic_serial -
      joined$first_pass_mrr_pessimistic_parallel
  )
  joined$final_top1_abs_difference <- abs(
    joined$final_top1_fractional_serial - joined$final_top1_fractional_parallel
  )
  joined$final_mrr_abs_difference <- abs(
    joined$final_mrr_fractional_serial - joined$final_mrr_fractional_parallel
  )
  joined$tie_random_seed_equal <-
    joined$tie_random_seed_serial == joined$tie_random_seed_parallel
  joined$final_rank_random_equal <-
    (joined$final_rank_random_serial == joined$final_rank_random_parallel) |
    (is.na(joined$final_rank_random_serial) &
       is.na(joined$final_rank_random_parallel))
  joined$final_top1_random_abs_difference <- abs(
    joined$final_top1_random_serial - joined$final_top1_random_parallel
  )
  joined$final_mrr_random_abs_difference <- abs(
    joined$final_mrr_random_serial - joined$final_mrr_random_parallel
  )
  joined
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

condition_complete <- function(dataset, first_pass_method, requested_k,
                               execution_mode, n_query) {
  if (!nrow(per_query_all)) return(FALSE)
  hit <- per_query_all$dataset == dataset &
    per_query_all$first_pass_method == first_pass_method &
    per_query_all$requested_k == requested_k &
    per_query_all$execution_mode == execution_mode
  count <- sum(hit, na.rm = TRUE)
  if (count == 0L) return(FALSE)
  if (count != n_query) {
    stop(
      "Checkpoint has an incomplete or duplicated condition: ", dataset, "/",
      first_pass_method, "/K=", requested_k, "/", execution_mode,
      " (", count, " rows; expected ", n_query, ")."
    )
  }
  TRUE
}

query_design_rows <- list()
dataset_stage_rows <- list()

for (dataset_index in seq_along(datasets)) {
  dataset <- datasets[[dataset_index]]
  message("\n=== ", dataset, " ===")
  loaded <- load_dataset_timed(dataset, params0)
  spectra <- loaded$spectra
  queries <- build_query_set(spectra, dataset)
  query_design_rows[[dataset]] <- queries$rows
  dataset_stage_rows[[dataset]] <- data.frame(
    dataset = dataset,
    source_type = input_checksums$source_type[input_checksums$dataset == dataset],
    n_library = length(spectra$frag_list),
    n_eligible_queries = length(queries$query_ids),
    input_loading_seconds = loaded$input_loading_seconds,
    preprocessing_seconds = loaded$preprocessing_seconds,
    query_setup_seconds = queries$setup_seconds,
    preprocessing_kind = loaded$preprocessing_kind,
    stringsAsFactors = FALSE
  )

  for (first_pass_method in first_pass_methods) {
    for (requested_k in k_values) {
      effective_k <- min(requested_k, length(spectra$frag_list) - 1L)
      if (effective_k < 1L) stop(dataset, " library is too small after self-exclusion.")
      for (execution_mode in execution_modes) {
        if (condition_complete(
            dataset, first_pass_method, requested_k, execution_mode,
            length(queries$query_ids))) {
          message("[resume] skipping ", dataset, "/", first_pass_method,
                  "/K=", requested_k, "/", execution_mode)
          next
        }
        use_parallel <- execution_mode == "parallel_8core"
        effective_cores <- if (use_parallel) resolved_parallel_cores else 1L
        params <- params0
        params$use_parallel <- use_parallel
        params$n_cores <- effective_cores
        params <- getFromNamespace("validate_params", "ppmWass")(params)

        message(
          "[rerank] ", dataset, " | first=", first_pass_method,
          " | K=", requested_k,
          if (effective_k != requested_k) paste0(" (effective ", effective_k, ")") else "",
          " | ", execution_mode, " | cores=", effective_cores
        )
        invisible(gc(reset = TRUE))
        search_result <- search_topk_rerank(
          query_frag_list = queries$fragment,
          query_loss_list = queries$derived,
          lib_frag_list = spectra$frag_list,
          lib_loss_list = spectra$loss_list,
          params = params,
          first_pass_method = first_pass_method,
          rerank_method = "ppm_wasserstein",
          first_pass_top_k = effective_k,
          final_top_k = effective_k,
          query_ids = queries$query_ids,
          lib_ids = names(spectra$frag_list),
          exclude_self = TRUE,
          return_first_pass = TRUE,
          use_parallel = use_parallel,
          n_cores = effective_cores,
          progress = show_progress
        )
        heap <- heap_peak_after_search(gc())

        expected_finite <- length(spectra$frag_list) - 1L
        if (any(search_result$candidate_summary$n_finite_first_pass != expected_finite)) {
          stop(
            "Nonfinite first-pass distances detected for ", dataset, "/",
            first_pass_method, "/K=", requested_k, "/", execution_mode, "."
          )
        }
        if (any(search_result$candidate_summary$n_candidates != effective_k) ||
            any(search_result$candidate_summary$n_returned != effective_k)) {
          stop(
            "Nonfinite or missing ppm-Wasserstein candidate scores detected for ",
            dataset, "/", first_pass_method, "/K=", requested_k, "/",
            execution_mode, "."
          )
        }
        required_timing <- c(
          "first_pass_sec", "topk_extraction_sec", "rerank_sec",
          "result_formatting_sec", "total_sec"
        )
        if (!all(required_timing %in% names(search_result$timing))) {
          stop(
            "search_topk_rerank() lacks required publication timing fields: ",
            paste(setdiff(required_timing, names(search_result$timing)), collapse = ", ")
          )
        }

        per_query <- evaluate_condition(
          search_result, spectra, queries, dataset, first_pass_method,
          requested_k, effective_k, execution_mode, effective_cores
        )
        per_query_all <- append_rows(per_query_all, per_query)

        ppm_evaluations <- sum(search_result$candidate_summary$n_candidates)
        eligible_cells <- length(queries$query_ids) * expected_finite
        search_stage_sum <- search_result$timing$first_pass_sec +
          search_result$timing$topk_extraction_sec +
          search_result$timing$rerank_sec +
          search_result$timing$result_formatting_sec
        timing_row <- data.frame(
          dataset = dataset,
          first_pass_method = first_pass_method,
          rerank_method = "ppm_wasserstein",
          requested_k = requested_k,
          effective_k = effective_k,
          execution_mode = execution_mode,
          effective_cores = effective_cores,
          n_queries = length(queries$query_ids),
          n_library = length(spectra$frag_list),
          eligible_library_size_excluding_self = expected_finite,
          input_loading_seconds = loaded$input_loading_seconds,
          preprocessing_seconds = loaded$preprocessing_seconds,
          query_setup_seconds = queries$setup_seconds,
          first_pass_seconds = search_result$timing$first_pass_sec,
          topk_extraction_seconds = search_result$timing$topk_extraction_sec,
          ppm_candidate_scoring_seconds = search_result$timing$rerank_sec,
          result_formatting_seconds = search_result$timing$result_formatting_sec,
          search_total_seconds = search_result$timing$total_sec,
          search_unattributed_overhead_seconds =
            search_result$timing$total_sec - search_stage_sum,
          cold_pipeline_total_seconds = loaded$input_loading_seconds +
            loaded$preprocessing_seconds + queries$setup_seconds +
            search_result$timing$total_sec,
          cold_pipeline_total_definition = paste(
            "sum of separately measured shared setup stages plus",
            "condition-specific search; not a directly timed single cold run"
          ),
          ppm_candidate_evaluations = ppm_evaluations,
          eligible_query_library_cells = eligible_cells,
          ppm_candidate_fraction = ppm_evaluations / eligible_cells,
          k_over_n_minus_one = effective_k / expected_finite,
          r_heap_peak_ncells = heap[["r_heap_peak_ncells"]],
          r_heap_peak_vcells = heap[["r_heap_peak_vcells"]],
          r_heap_peak_ncells_mb = heap[["r_heap_peak_ncells_mb"]],
          r_heap_peak_vcells_mb = heap[["r_heap_peak_vcells_mb"]],
          r_heap_peak_total_mb = heap[["r_heap_peak_total_mb"]],
          r_heap_peak_definition = paste(
            "gc(reset=TRUE) max used; parent R Ncells+Vcells MB; excludes",
            "native/RSS and worker heaps"
          ),
          stringsAsFactors = FALSE
        )
        if (abs(timing_row$ppm_candidate_fraction - timing_row$k_over_n_minus_one) > 1e-15) {
          stop("Internal candidate-fraction accounting error.")
        }
        timing_all <- append_rows(timing_all, timing_row)
        utils::write.csv(per_query_all, checkpoint_per_query, row.names = FALSE)
        utils::write.csv(timing_all, checkpoint_timings, row.names = FALSE)
      }
    }
  }
}

query_design <- do.call(rbind, query_design_rows)
dataset_stages <- do.call(rbind, dataset_stage_rows)
condition_summary <- summarize_conditions(per_query_all, timing_all)
equivalence <- serial_parallel_audit(per_query_all)

if (nrow(equivalence)) {
  first_pass_difference_columns <- c(
    "first_pass_top1_optimistic_abs_difference",
    "first_pass_top1_fractional_abs_difference",
    "first_pass_top1_pessimistic_abs_difference",
    "first_pass_mrr_optimistic_abs_difference",
    "first_pass_mrr_fractional_abs_difference",
    "first_pass_mrr_pessimistic_abs_difference"
  )
  first_pass_equivalence_failed <- Reduce(
    `|`,
    lapply(first_pass_difference_columns, function(column) {
      is.na(equivalence[[column]]) | equivalence[[column]] > 1e-12
    })
  )
  failed <- is.na(equivalence$candidate_set_equal) |
    is.na(equivalence$rerank_order_equal) |
    !equivalence$candidate_set_equal | !equivalence$rerank_order_equal |
    equivalence$candidate_retention_abs_difference > 0 |
    first_pass_equivalence_failed |
    equivalence$final_top1_abs_difference > 1e-12 |
    equivalence$final_mrr_abs_difference > 1e-12 |
    is.na(equivalence$tie_random_seed_equal) |
    !equivalence$tie_random_seed_equal |
    is.na(equivalence$final_rank_random_equal) |
    !equivalence$final_rank_random_equal |
    equivalence$final_top1_random_abs_difference > 1e-12 |
    equivalence$final_mrr_random_abs_difference > 1e-12
  if (any(failed)) {
    utils::write.csv(
      equivalence, file.path(output_dir, "serial_parallel_equivalence.csv"),
      row.names = FALSE
    )
    stop("Serial/parallel reranking equivalence failed for ", sum(failed), " query-condition rows.")
  }
}

sort_query <- order(
  match(per_query_all$dataset, datasets),
  match(per_query_all$first_pass_method, first_pass_methods),
  per_query_all$requested_k,
  match(per_query_all$execution_mode, execution_modes),
  per_query_all$query_position
)
per_query_all <- per_query_all[sort_query, , drop = FALSE]
sort_timing <- order(
  match(timing_all$dataset, datasets),
  match(timing_all$first_pass_method, first_pass_methods),
  timing_all$requested_k,
  match(timing_all$execution_mode, execution_modes)
)
timing_all <- timing_all[sort_timing, , drop = FALSE]

effective_parameters <- do.call(rbind, lapply(first_pass_methods, function(method) {
  do.call(rbind, lapply(execution_modes, function(execution_mode) {
    data.frame(
      first_pass_method = method,
      rerank_method = "ppm_wasserstein",
      execution_mode = execution_mode,
      effective_cores = if (execution_mode == "serial") 1L else resolved_parallel_cores,
      fragment_weight = 0.70,
      pooled_derived_mass_difference_weight = 0.30,
      min_mz = params0$min_mz,
      max_mz = params0$max_mz,
      centroid_ppm = params0$centroid_ppm,
      noise_threshold_relative_to_base_peak = params0$noise_thr,
      top_k_preprocessed_peaks = params0$topK,
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
      backend = params0$backend,
      stringsAsFactors = FALSE
    )
  }))
}))

utils::write.csv(
  query_design, file.path(output_dir, "reranking_query_design.csv"), row.names = FALSE
)
utils::write.csv(
  dataset_stages, file.path(output_dir, "dataset_stage_timings.csv"), row.names = FALSE
)
utils::write.csv(
  per_query_all, file.path(output_dir, "reranking_per_query.csv"), row.names = FALSE
)
utils::write.csv(
  condition_summary, file.path(output_dir, "reranking_condition_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  timing_all, file.path(output_dir, "runtime_stage_timings.csv"), row.names = FALSE
)
utils::write.csv(
  equivalence, file.path(output_dir, "serial_parallel_equivalence.csv"), row.names = FALSE
)
utils::write.csv(
  effective_parameters, file.path(output_dir, "effective_parameters.csv"), row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"), useBytes = TRUE
)

write_manifest("complete")
run_complete <- TRUE
message("End-to-end reranking benchmark complete: ", output_dir)
