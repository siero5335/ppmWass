#!/usr/bin/env Rscript

# Publication-grade secondary-representation ablation.
#
# The seven configurations are evaluated with identical preprocessing and all
# seven publication methods. Split derived representations are enabled while
# spectra are built, then each ablation configuration passes exactly one
# selected derived representation through the legacy loss slot. This avoids
# invoking the package's internal split-channel weighting during the ablation.
#
# Primary inference uses the exact full-precision tie definition, fractional
# random-order expectations, and clustered bootstrap confidence intervals from
# inst/scripts/lib/publication_retrieval_statistics.R.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]])
}

flag_is_true <- function(name, default = FALSE) {
  value <- get_arg(name, if (default) "true" else "false")
  tolower(value) %in% c("1", "true", "yes", "y")
}

split_arg <- function(name, default) {
  value <- get_arg(name, default)
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
}

`%||%` <- function(x, y) if (is.null(x)) y else x

canonical_methods <- c(
  "ppm_wasserstein",
  "composite",
  "entropy_weighted",
  "entropy_unweighted",
  "cosine",
  "weighted_cosine",
  "hellinger"
)

configuration_definitions <- data.frame(
  configuration = c(
    "fragment_only",
    "anchored_only",
    "pairwise_only",
    "pooled_only",
    "fragment_plus_anchored",
    "fragment_plus_pairwise",
    "fragment_plus_pooled"
  ),
  label = c(
    "Fragment only",
    "Anchored differences only",
    "Pairwise differences only",
    "Pooled derived differences only",
    "Fragment + anchored differences",
    "Fragment + pairwise differences",
    "Fragment + pooled derived differences"
  ),
  fragment_source = c(
    "fragment", "empty", "empty", "empty",
    "fragment", "fragment", "fragment"
  ),
  derived_source = c(
    "empty", "anchored", "pairwise", "pooled",
    "anchored", "pairwise", "pooled"
  ),
  w_frag = c(1, 0, 0, 0, 0.70, 0.70, 0.70),
  w_derived = c(0, 1, 1, 1, 0.30, 0.30, 0.30),
  stringsAsFactors = FALSE
)

repo_dir <- normalizePath(get_arg("repo", getwd()), mustWork = TRUE)
bundle_root <- normalizePath(
  get_arg("bundle-root", file.path(repo_dir, "..", "..")),
  mustWork = TRUE
)
input_dir <- normalizePath(
  get_arg("input-dir", file.path(bundle_root, "inputs")),
  mustWork = TRUE
)
recetox_msp <- normalizePath(
  get_arg("recetox-msp", file.path(input_dir, "RECETOX_merged.msp")),
  mustWork = TRUE
)
hrei_msp <- normalizePath(
  get_arg("hrei-msp", file.path(input_dir, "HREI-MSDB.msp")),
  mustWork = TRUE
)
output_dir <- get_arg(
  "output-dir",
  file.path(bundle_root, "results", "current", "secondary_ablation")
)
log_file <- get_arg("log-file", file.path(output_dir, "run.log"))

datasets <- split_arg("datasets", "RECETOX,HREI-MSDB")
methods <- split_arg("methods", paste(canonical_methods, collapse = ","))
configurations <- split_arg(
  "configurations",
  paste(configuration_definitions$configuration, collapse = ",")
)

query_boot <- as.integer(get_arg("query-boot", "1000"))
cluster_boot <- as.integer(get_arg("cluster-boot", "20000"))
seed <- as.integer(get_arg("seed", "20260718"))
tol_ppm <- as.numeric(get_arg("tol-ppm", "15"))
n_cores <- as.integer(get_arg("n-cores", "8"))
use_parallel <- flag_is_true("parallel", TRUE)
write_distance_csv <- flag_is_true("write-distance-csv", FALSE)
resume <- flag_is_true("resume", FALSE)
overwrite <- flag_is_true("overwrite", FALSE)
smoke <- flag_is_true("smoke", FALSE)
smoke_clusters <- as.integer(get_arg("smoke-clusters", "2"))
smoke_per_cluster <- as.integer(get_arg("smoke-per-cluster", "2"))
recetox_spectra_rds <- get_arg("recetox-spectra-rds", "")
hrei_spectra_rds <- get_arg("hrei-spectra-rds", "")

if (!all(datasets %in% c("RECETOX", "HREI-MSDB"))) {
  stop("Unknown dataset(s): ", paste(setdiff(datasets, c("RECETOX", "HREI-MSDB")), collapse = ", "))
}
if (!length(methods) || !all(methods %in% canonical_methods)) {
  stop("Unknown method(s): ", paste(setdiff(methods, canonical_methods), collapse = ", "))
}
if (!length(configurations) ||
    !all(configurations %in% configuration_definitions$configuration)) {
  stop(
    "Unknown configuration(s): ",
    paste(setdiff(configurations, configuration_definitions$configuration), collapse = ", ")
  )
}
if (!"fragment_only" %in% configurations) {
  stop("fragment_only must be included because paired differences use it as the reference.")
}
if (!is.finite(query_boot) || query_boot < 1L ||
    !is.finite(cluster_boot) || cluster_boot < 1L) {
  stop("query-boot and cluster-boot must be positive integers.")
}
if (!is.finite(seed) || !is.finite(tol_ppm) || tol_ppm <= 0 ||
    !is.finite(n_cores) || n_cores < 1L) {
  stop("Invalid seed, tol-ppm, or n-cores.")
}
if (smoke && (!is.finite(smoke_clusters) || smoke_clusters < 1L ||
              !is.finite(smoke_per_cluster) || smoke_per_cluster < 2L)) {
  stop("Smoke mode requires smoke-clusters >= 1 and smoke-per-cluster >= 2.")
}
if (resume && overwrite) stop("Choose at most one of --resume=true and --overwrite=true.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (nzchar(log_file)) {
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  log_connection <- file(log_file, open = if (file.exists(log_file)) "at" else "wt")
  sink(log_connection, split = TRUE)
  sink(log_connection, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_connection)
  }, add = TRUE)
  cat("\n=== secondary ablation run started ", format(Sys.time(), tz = "UTC"), " UTC ===\n", sep = "")
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to run this script from the development tree.")
}
pkgload::load_all(repo_dir, quiet = TRUE)

stats_helper <- file.path(
  repo_dir, "inst", "scripts", "lib", "publication_retrieval_statistics.R"
)
if (!file.exists(stats_helper)) stop("Missing statistics helper: ", stats_helper)
source(stats_helper, local = TRUE)
RNGkind("L'Ecuyer-CMRG")

message("Repository: ", repo_dir)
message("Bundle root: ", bundle_root)
message("Output dir: ", normalizePath(output_dir, mustWork = FALSE))
message("Datasets: ", paste(datasets, collapse = ", "))
message("Methods: ", paste(methods, collapse = ", "))
message("Configurations: ", paste(configurations, collapse = ", "))
message("Exact fractional tie policy; cluster bootstrap R: ", cluster_boot)
message("Smoke mode: ", smoke)

base_params <- function(build_secondary = FALSE) {
  params <- eihrms_default_params()
  params$min_mz <- 35
  params$max_mz <- 650
  params$noise_thr <- 0.01
  params$topK <- 200L
  params$class_detection_ppm <- 15
  params$tol_ppm <- tol_ppm
  params$ot_method <- "exact"
  params$sinkhorn_epsilon <- 0.05
  params$sinkhorn_niter <- 100L
  params$wasserstein_transition_mult <- 3
  params$use_typical_loss <- FALSE
  params$use_mref_confidence <- FALSE
  # TRUE only while building spectra so anchored and pairwise representations
  # are retained separately. Scoring configurations set this back to FALSE.
  params$use_split_loss <- isTRUE(build_secondary)
  params$use_parallel <- use_parallel
  params$n_cores <- n_cores
  params$backend <- "pair_loop"
  params
}

aligned_fields <- c(
  "frag_list", "loss_list", "derived_list", "loss_typ_list",
  "loss_anchor_list", "loss_pair_list", "derived_anchor_list",
  "derived_pair_list", "loss_anchor_typ_list", "loss_pair_typ_list"
)

subset_spectra_local <- function(spectra, idx) {
  n_before <- nrow(spectra$df_spec)
  out <- spectra
  out$df_spec <- spectra$df_spec[idx, , drop = FALSE]
  for (field in aligned_fields) {
    value <- spectra[[field]]
    if (!is.null(value) && length(value) == n_before) out[[field]] <- value[idx]
  }
  for (field in c("ri", "mref_confidence")) {
    value <- spectra[[field]]
    if (!is.null(value) && length(value) == n_before) out[[field]] <- value[idx]
  }
  ids <- out$df_spec$id
  for (field in aligned_fields) {
    if (!is.null(out[[field]]) && length(out[[field]]) == length(ids)) {
      names(out[[field]]) <- ids
    }
  }
  if (!is.null(out$ri) && length(out$ri) == length(ids)) names(out$ri) <- ids
  if (!is.null(out$mref_confidence) && length(out$mref_confidence) == length(ids)) {
    names(out$mref_confidence) <- ids
  }
  out
}

filter_valid_spectra_local <- function(spectra) {
  valid <- vapply(spectra$frag_list, function(x) {
    is.matrix(x) && ncol(x) >= 2L && nrow(x) >= 2L &&
      all(is.finite(x[, 1:2, drop = FALSE])) && sum(x[, 2]) > 0
  }, logical(1))
  if (!all(valid)) {
    message("Removed ", sum(!valid), " degenerate fragment spectra.")
    spectra <- subset_spectra_local(spectra, which(valid))
  }
  spectra
}

normalize_hrei_msp_names <- function(input, output) {
  lines <- readLines(input, warn = FALSE)
  starts <- grep("^Name:", lines)
  ends <- c(starts[-1] - 1L, length(lines))
  for (i in seq_along(starts)) {
    block <- lines[starts[[i]]:ends[[i]]]
    db_line <- grep("^DB#:", block, value = TRUE)
    db_id <- if (length(db_line)) {
      trimws(sub("^DB#:[[:space:]]*", "", db_line[[1]]))
    } else {
      as.character(i)
    }
    original_name <- trimws(sub("^Name:[[:space:]]*", "", lines[starts[[i]]]))
    lines[starts[[i]]] <- paste0("Name: ", original_name, " __HREI_DB", db_id)
  }
  writeLines(lines, output, useBytes = TRUE)
  output
}

validate_secondary_representations <- function(spectra, dataset) {
  required <- c("frag_list", "loss_list", "loss_anchor_list", "loss_pair_list")
  missing <- required[vapply(required, function(x) is.null(spectra[[x]]), logical(1))]
  if (length(missing)) {
    stop(
      dataset, " spectra do not contain required secondary representations: ",
      paste(missing, collapse = ", "),
      ". Rebuild with params$use_split_loss = TRUE."
    )
  }
  n <- nrow(spectra$df_spec)
  bad_length <- required[vapply(required, function(x) length(spectra[[x]]) != n, logical(1))]
  if (length(bad_length)) {
    stop(dataset, " has misaligned representation lists: ", paste(bad_length, collapse = ", "))
  }
  ids <- spectra$df_spec$id
  if (anyDuplicated(ids)) stop(dataset, " has duplicate spectrum IDs after preprocessing.")
  for (field in required) names(spectra[[field]]) <- ids
  spectra
}

load_dataset <- function(dataset, dataset_dir) {
  prepared_path <- if (dataset == "RECETOX") recetox_spectra_rds else hrei_spectra_rds
  if (nzchar(prepared_path)) {
    message("Loading prepared split representations: ", prepared_path)
    spectra <- readRDS(normalizePath(prepared_path, mustWork = TRUE))
  } else {
    params <- base_params(build_secondary = TRUE)
    input <- if (dataset == "RECETOX") {
      recetox_msp
    } else {
      unique_msp <- file.path(dataset_dir, "HREI-MSDB_unique_names_secondary_ablation.msp")
      if (!file.exists(unique_msp) || overwrite) {
        normalize_hrei_msp_names(hrei_msp, unique_msp)
      }
      unique_msp
    }
    spectra_all <- build_spectra_from_msp(
      input, params, require_ri = FALSE, progress = !smoke
    )
    df <- spectra_all$df_spec
    idx <- which(
      !is.na(df$RI) & !is.na(df$inchikey) & df$inchikey != "" &
        nchar(df$inchikey) >= 14L
    )
    spectra <- subset_spectra_local(spectra_all, idx)
    spectra$raw_msp_metadata <- spectra_all$msp_metadata
    spectra$msp_parser_diagnostics <- attr(
      spectra_all$msp_metadata, "msp_parser_diagnostics"
    )
  }
  spectra <- filter_valid_spectra_local(spectra)
  validate_secondary_representations(spectra, dataset)
}

select_smoke_subset <- function(spectra, dataset) {
  inchikey <- as.character(spectra$df_spec$inchikey)
  representation_complete <- vapply(seq_along(inchikey), function(i) {
    all(vapply(
      c("loss_list", "loss_anchor_list", "loss_pair_list"),
      function(field) {
        value <- spectra[[field]][[i]]
        is.matrix(value) && nrow(value) > 0L
      },
      logical(1)
    ))
  }, logical(1))

  complete_counts <- table(inchikey[representation_complete])
  groups <- names(complete_counts[complete_counts >= 2L])
  if (!length(groups)) {
    all_counts <- table(inchikey)
    groups <- names(all_counts[all_counts >= 2L])
    representation_complete[] <- TRUE
  }
  if (!length(groups)) stop(dataset, " has no duplicate full-InChIKey group for smoke testing.")

  groups <- head(groups, smoke_clusters)
  idx <- unlist(lapply(groups, function(group) {
    candidates <- which(inchikey == group & representation_complete)
    head(candidates, smoke_per_cluster)
  }), use.names = FALSE)
  idx <- sort(unique(idx))
  if (length(idx) < 2L) stop("Smoke subset selection produced fewer than two spectra.")
  message(
    "Smoke subset: ", length(idx), " spectra from ", length(unique(inchikey[idx])),
    " full-InChIKey cluster(s)."
  )
  subset_spectra_local(spectra, idx)
}

empty_spectrum <- function() {
  matrix(
    numeric(0), nrow = 0L, ncol = 2L,
    dimnames = list(NULL, c("mz", "intensity"))
  )
}

empty_spectrum_list <- function(ids) {
  out <- rep(list(empty_spectrum()), length(ids))
  names(out) <- ids
  out
}

configuration_inputs <- function(spectra, definition) {
  ids <- spectra$df_spec$id
  empty <- empty_spectrum_list(ids)
  fragments <- if (definition$fragment_source[[1]] == "fragment") {
    spectra$frag_list
  } else {
    empty
  }
  derived <- switch(
    definition$derived_source[[1]],
    empty = empty,
    anchored = spectra$loss_anchor_list,
    pairwise = spectra$loss_pair_list,
    pooled = spectra$loss_list,
    stop("Unhandled derived source: ", definition$derived_source[[1]])
  )
  names(fragments) <- ids
  names(derived) <- ids
  list(fragments = fragments, derived = derived)
}

validate_distance_matrix <- function(dm, ids, context) {
  if (!is.matrix(dm) || !identical(dim(dm), c(length(ids), length(ids)))) {
    stop(context, " returned a distance matrix with unexpected dimensions.")
  }
  if (!identical(rownames(dm), ids) || !identical(colnames(dm), ids)) {
    stop(context, " returned a distance matrix with misaligned IDs.")
  }
  if (any(!is.finite(dm))) {
    stop(
      context, " produced ", sum(!is.finite(dm)),
      " nonfinite cells after deterministic pair-level fallback."
    )
  }
  dm
}

write_distance_csv_local <- function(dm, path) {
  out <- data.frame(id = rownames(dm), dm, check.names = FALSE)
  utils::write.csv(out, path, row.names = FALSE)
}

checkpoint_distance <- function(path, compute, ids, context) {
  if (file.exists(path) && resume) {
    message("[resume] ", path)
    return(validate_distance_matrix(readRDS(path), ids, context))
  }
  if (file.exists(path) && !overwrite) {
    stop("Checkpoint exists; use --resume=true or --overwrite=true: ", path)
  }
  dm <- validate_distance_matrix(compute(), ids, context)
  saveRDS(dm, path, compress = TRUE)
  dm
}

representation_qc <- function(spectra, dataset) {
  fields <- c(
    fragment = "frag_list", pooled = "loss_list",
    anchored = "loss_anchor_list", pairwise = "loss_pair_list"
  )
  do.call(rbind, lapply(names(fields), function(representation) {
    values <- spectra[[fields[[representation]]]]
    peaks <- vapply(values, function(x) if (is.matrix(x)) nrow(x) else 0L, integer(1))
    data.frame(
      dataset = dataset,
      representation = representation,
      n_spectra = length(values),
      n_empty = sum(peaks == 0L),
      min_peaks = min(peaks),
      median_peaks = stats::median(peaks),
      max_peaks = max(peaks),
      stringsAsFactors = FALSE
    )
  }))
}

add_configuration_columns <- function(data, definition) {
  data$configuration <- definition$configuration[[1]]
  data$configuration_label <- definition$label[[1]]
  data$fragment_source <- definition$fragment_source[[1]]
  data$derived_source <- definition$derived_source[[1]]
  data$w_frag <- definition$w_frag[[1]]
  data$w_derived <- definition$w_derived[[1]]
  data
}

run_dataset <- function(dataset, dataset_index) {
  dataset_dir <- file.path(output_dir, gsub("[^A-Za-z0-9]+", "_", dataset))
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n=== Dataset: ", dataset, " ===")

  spectra <- load_dataset(dataset, dataset_dir)
  if (smoke) spectra <- select_smoke_subset(spectra, dataset)
  spectra <- validate_secondary_representations(spectra, dataset)
  spectra_path <- file.path(dataset_dir, "spectra_secondary_ablation.rds")
  if (!file.exists(spectra_path) || overwrite) {
    saveRDS(spectra, spectra_path, compress = TRUE)
  } else if (!resume) {
    stop("Prepared spectra checkpoint exists; use --resume=true or --overwrite=true: ", spectra_path)
  }

  ids <- spectra$df_spec$id
  inchikeys <- as.character(spectra$df_spec$inchikey)
  names(inchikeys) <- ids
  config_table <- configuration_definitions[
    match(configurations, configuration_definitions$configuration), , drop = FALSE
  ]

  per_configuration <- setNames(vector("list", nrow(config_table)), config_table$configuration)
  query_rows <- list()
  summary_rows <- list()
  tie_sensitivity_rows <- list()
  tie_diagnostic_rows <- list()
  matrix_qc_rows <- list()
  row_index <- 0L

  for (ci in seq_len(nrow(config_table))) {
    definition <- config_table[ci, , drop = FALSE]
    configuration <- definition$configuration[[1]]
    configuration_dir <- file.path(dataset_dir, configuration)
    dir.create(configuration_dir, recursive = TRUE, showWarnings = FALSE)
    inputs <- configuration_inputs(spectra, definition)
    per_configuration[[configuration]] <- list()

    message("\n[", dataset, "] configuration: ", configuration)
    for (method in methods) {
      method_index <- match(method, canonical_methods)
      context <- paste(dataset, configuration, method, sep = "/")
      message("  method: ", method)
      params <- base_params(build_secondary = FALSE)
      params$distance_method <- method
      params$w_frag <- definition$w_frag[[1]]
      params$w_loss <- definition$w_derived[[1]]
      params$use_split_loss <- FALSE
      params$use_typical_loss <- FALSE
      params$use_mref_confidence <- FALSE
      params <- validate_params(params)

      matrix_path <- file.path(configuration_dir, paste0("dist_", method, ".rds"))
      dm <- checkpoint_distance(
        matrix_path,
        compute = function() {
          compute_distance_matrix(
            inputs$fragments,
            inputs$derived,
            params,
            progress = FALSE
          )
        },
        ids = ids,
        context = context
      )
      if (write_distance_csv) {
        csv_path <- file.path(configuration_dir, paste0("dist_", method, ".csv"))
        if (!file.exists(csv_path) || overwrite || resume) {
          write_distance_csv_local(dm, csv_path)
        }
      }

      stat_seed <- seed + dataset_index * 1000000L + ci * 10000L + method_index * 100L
      per_query <- per_query_metrics_tie_aware(
        dm,
        inchikeys,
        tie_tolerance = 0,
        random_seed = stat_seed
      )
      per_configuration[[configuration]][[method]] <- per_query

      query_output <- add_configuration_columns(per_query, definition)
      query_output$dataset <- dataset
      query_output$method <- method
      row_index <- row_index + 1L
      query_rows[[row_index]] <- query_output
      utils::write.csv(
        query_output,
        file.path(configuration_dir, paste0("per_query_", method, ".csv")),
        row.names = FALSE
      )

      summary <- summarize_tie_aware_metrics(
        per_query,
        dataset,
        method,
        query_boot_R = query_boot,
        cluster_boot_R = cluster_boot,
        seed = stat_seed
      )
      summary$tie_definition <- "full_precision_exact_equality"
      summary$cluster_unit <- ifelse(
        summary$metric %in% c("p_at_1", "map", "p_at_5", "p_at_10"),
        "inchikey_prefix_14", "full_inchikey"
      )
      summary$estimand <- "query_weighted_mean_conditional_on_fixed_library"
      summary$query_boot_R <- query_boot
      summary$cluster_boot_R <- cluster_boot
      summary$seed <- stat_seed
      summary_rows[[row_index]] <- add_configuration_columns(summary, definition)

      tie_sensitivity <- summarize_tie_sensitivity(per_query, dataset, method)
      tie_sensitivity$tie_definition <- "full_precision_exact_equality"
      tie_sensitivity_rows[[row_index]] <- add_configuration_columns(
        tie_sensitivity, definition
      )

      diagnostic <- tie_diagnostics(per_query, dataset, method)
      diagnostic$tie_definition <- "full_precision_exact_equality"
      tie_diagnostic_rows[[row_index]] <- add_configuration_columns(
        diagnostic, definition
      )

      matrix_qc_rows[[row_index]] <- data.frame(
        dataset = dataset,
        configuration = configuration,
        method = method,
        n_rows = nrow(dm),
        n_cols = ncol(dm),
        n_nonfinite = sum(!is.finite(dm)),
        min_distance = min(dm),
        max_distance = max(dm),
        directional_full_matrix = identical(method, "composite"),
        stringsAsFactors = FALSE
      )
      rm(dm)
      invisible(gc(FALSE))
    }
  }

  metric_specs <- data.frame(
    metric = c("top1", "mrr", "p_at_1"),
    value_col = c("top1_fractional", "rr_fractional", "p_at_1_fractional"),
    cluster_col = c("inchikey", "inchikey", "inchikey_prefix"),
    stringsAsFactors = FALSE
  )
  paired_rows <- list()
  paired_index <- 0L
  for (ci in seq_len(nrow(config_table))) {
    definition <- config_table[ci, , drop = FALSE]
    configuration <- definition$configuration[[1]]
    for (method in methods) {
      current <- per_configuration[[configuration]][[method]]
      fragment <- per_configuration[["fragment_only"]][[method]]
      method_index <- match(method, canonical_methods)
      for (mi in seq_len(nrow(metric_specs))) {
        paired_index <- paired_index + 1L
        paired_seed <- seed + 5000000L + dataset_index * 1000000L +
          ci * 10000L + method_index * 100L + mi
        difference <- paired_cluster_bootstrap_difference(
          current,
          fragment,
          value_col = metric_specs$value_col[[mi]],
          cluster_col = metric_specs$cluster_col[[mi]],
          R = cluster_boot,
          seed = paired_seed
        )
        paired_rows[[paired_index]] <- data.frame(
          dataset = dataset,
          method = method,
          configuration = configuration,
          configuration_label = definition$label[[1]],
          reference_configuration = "fragment_only",
          metric = metric_specs$metric[[mi]],
          tie_policy = "fractional_expected",
          tie_definition = "full_precision_exact_equality",
          difference_direction = "configuration_minus_fragment_only",
          cluster_unit = if (
            metric_specs$metric[[mi]] == "p_at_1"
          ) "inchikey_prefix_14" else "full_inchikey",
          estimand = "query_weighted_mean_conditional_on_fixed_library",
          difference = difference[["estimate"]],
          ci_low = difference[["ci_low"]],
          ci_high = difference[["ci_high"]],
          n_queries = as.integer(difference[["n_queries"]]),
          n_clusters = as.integer(difference[["n_clusters"]]),
          cluster_boot_R = cluster_boot,
          seed = paired_seed,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    dataset = dataset,
    spectra_path = spectra_path,
    query_level = do.call(rbind, query_rows),
    summary = do.call(rbind, summary_rows),
    paired_vs_fragment = do.call(rbind, paired_rows),
    tie_sensitivity = do.call(rbind, tie_sensitivity_rows),
    tie_diagnostics = do.call(rbind, tie_diagnostic_rows),
    representation_qc = representation_qc(spectra, dataset),
    matrix_qc = do.call(rbind, matrix_qc_rows)
  )
}

selected_configuration_table <- configuration_definitions[
  match(configurations, configuration_definitions$configuration), , drop = FALSE
]
utils::write.csv(
  selected_configuration_table,
  file.path(output_dir, "secondary_ablation_configuration_definitions.csv"),
  row.names = FALSE
)

results <- lapply(seq_along(datasets), function(i) run_dataset(datasets[[i]], i))
names(results) <- datasets

combined <- function(field) do.call(rbind, lapply(results, `[[`, field))
query_level <- combined("query_level")
summary <- combined("summary")
paired_vs_fragment <- combined("paired_vs_fragment")
tie_sensitivity <- combined("tie_sensitivity")
tie_diagnostics_all <- combined("tie_diagnostics")
representation_qc_all <- combined("representation_qc")
matrix_qc_all <- combined("matrix_qc")

utils::write.csv(
  query_level,
  file.path(output_dir, "secondary_ablation_query_level.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary,
  file.path(output_dir, "secondary_ablation_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_vs_fragment,
  file.path(output_dir, "secondary_ablation_paired_vs_fragment.csv"),
  row.names = FALSE
)
utils::write.csv(
  tie_sensitivity,
  file.path(output_dir, "secondary_ablation_tie_sensitivity.csv"),
  row.names = FALSE
)
utils::write.csv(
  tie_diagnostics_all,
  file.path(output_dir, "secondary_ablation_tie_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  representation_qc_all,
  file.path(output_dir, "secondary_ablation_representation_qc.csv"),
  row.names = FALSE
)
utils::write.csv(
  matrix_qc_all,
  file.path(output_dir, "secondary_ablation_distance_matrix_qc.csv"),
  row.names = FALSE
)

run_parameters <- data.frame(
  parameter = c(
    "timestamp_utc", "package_commit", "package_tree_dirty", "command",
    "repo_dir", "bundle_root", "recetox_msp", "hrei_msp",
    "datasets", "methods", "configurations", "query_boot_R",
    "cluster_boot_R", "seed", "hard_match_tolerance_ppm",
    "ppmWass_base_cost_ppm", "transition_multiplier", "saturation_width_ppm",
    "ot_estimand", "sinkhorn_epsilon", "sinkhorn_iterations", "solver_backend", "parallel",
    "n_cores", "smoke", "resume", "overwrite", "write_distance_csv"
  ),
  value = c(
    format(Sys.time(), tz = "UTC"),
    tryCatch(trimws(system2("git", c("-C", repo_dir, "rev-parse", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
             error = function(e) NA_character_),
    length(tryCatch(system2("git", c("-C", repo_dir, "status", "--porcelain"),
                            stdout = TRUE, stderr = FALSE),
                    error = function(e) character())) > 0L,
    paste(commandArgs(), collapse = " "),
    repo_dir, bundle_root, recetox_msp, hrei_msp,
    paste(datasets, collapse = ","), paste(methods, collapse = ","),
    paste(configurations, collapse = ","), query_boot, cluster_boot, seed,
    tol_ppm, tol_ppm, 3, tol_ppm * 3,
    "unregularized_exact_transport_cost", NA, NA, "transport_exact_unregularized",
    use_parallel, n_cores, smoke, resume, overwrite, write_distance_csv
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  run_parameters,
  file.path(output_dir, "secondary_ablation_run_parameters.csv"),
  row.names = FALSE
)

input_checksums <- data.frame(
  file = c(recetox_msp, hrei_msp, stats_helper),
  md5 = unname(tools::md5sum(c(recetox_msp, hrei_msp, stats_helper))),
  stringsAsFactors = FALSE
)
utils::write.csv(
  input_checksums,
  file.path(output_dir, "secondary_ablation_input_checksums.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    query_level = query_level,
    summary = summary,
    paired_vs_fragment = paired_vs_fragment,
    tie_sensitivity = tie_sensitivity,
    tie_diagnostics = tie_diagnostics_all,
    representation_qc = representation_qc_all,
    matrix_qc = matrix_qc_all,
    configuration_definitions = selected_configuration_table,
    run_parameters = run_parameters
  ),
  file.path(output_dir, "secondary_ablation_results.rds"),
  compress = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "secondary_ablation_sessionInfo.txt"),
  useBytes = TRUE
)

message("Secondary ablation complete.")
message("Summary: ", file.path(output_dir, "secondary_ablation_summary.csv"))
message(
  "Paired differences: ",
  file.path(output_dir, "secondary_ablation_paired_vs_fragment.csv")
)
