#!/usr/bin/env Rscript

# Audit legacy June retrieval results against a publication-grade rerun.
#
# The primary comparison joins legacy published query-level outcomes to metrics
# recalculated from each new full-precision distance matrix.  Deltas are always
# evaluated on query-ID-matched rows.  For Composite, a separate counterfactual
# reconstructs the legacy upper-triangle mirroring from the new directional
# matrix, so parser changes are held fixed while the directionality fix is
# isolated.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  substring(hit[[1]], nchar(prefix) + 1L)
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

canonical_methods <- c(
  "ppm_wasserstein",
  "composite",
  "entropy_weighted",
  "entropy_unweighted",
  "cosine",
  "weighted_cosine",
  "hellinger"
)
canonical_datasets <- c("RECETOX", "HREI-MSDB")

repo_dir <- normalizePath(get_arg("repo", getwd()), mustWork = TRUE)
bundle_root <- normalizePath(
  get_arg("bundle-root", file.path(repo_dir, "..", "..")),
  mustWork = TRUE
)
legacy_dir <- normalizePath(
  get_arg(
    "legacy-dir",
    file.path(
      bundle_root, "results", "legacy_20260616", "paper_bundle",
      "results", "main_tol15"
    )
  ),
  mustWork = TRUE
)
new_dir <- normalizePath(
  get_arg("new-dir", file.path(bundle_root, "results", "current", "main_clean")),
  mustWork = TRUE
)
output_dir <- normalizePath(
  get_arg(
    "output-dir",
    file.path(bundle_root, "manifest", "audit", "retrieval_old_new")
  ),
  mustWork = FALSE
)

datasets <- split_arg("datasets", paste(canonical_datasets, collapse = ","))
methods <- split_arg("methods", paste(canonical_methods, collapse = ","))
seed <- as.integer(get_arg("seed", "20260718"))
change_tolerance <- as.numeric(get_arg("change-tolerance", "1e-12"))
allow_partial <- flag_is_true("allow-partial", FALSE)
overwrite <- flag_is_true("overwrite", FALSE)

unknown_datasets <- setdiff(datasets, canonical_datasets)
unknown_methods <- setdiff(methods, canonical_methods)
if (length(unknown_datasets)) {
  stop("Unknown dataset(s): ", paste(unknown_datasets, collapse = ", "))
}
if (!length(methods) || length(unknown_methods)) {
  stop("Unknown method(s): ", paste(unknown_methods, collapse = ", "))
}
if (!is.finite(seed)) stop("seed must be a finite integer.")
if (!is.finite(change_tolerance) || change_tolerance < 0) {
  stop("change-tolerance must be a finite nonnegative number.")
}

stats_helper <- file.path(
  repo_dir, "inst", "scripts", "lib", "publication_retrieval_statistics.R"
)
if (!file.exists(stats_helper)) stop("Missing statistics helper: ", stats_helper)
source(stats_helper, local = TRUE)

output_files <- c(
  "old_new_retrieval_query_level.csv",
  "old_new_retrieval_aggregate.csv",
  "old_new_retrieval_unmatched_counts.csv",
  "composite_directionality_impact.csv",
  "composite_directionality_impact_query_level.csv",
  "composite_directionality_impact_aggregate.csv",
  "old_new_retrieval_audit_input_files.csv",
  "old_new_retrieval_audit_run_parameters.csv",
  "old_new_retrieval_audit_sessionInfo.txt",
  "old_new_retrieval_audit_results.rds"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
existing_outputs <- file.path(output_dir, output_files)
if (any(file.exists(existing_outputs)) && !overwrite) {
  stop(
    "Audit output already exists; use --overwrite=true or a fresh --output-dir: ",
    paste(basename(existing_outputs[file.exists(existing_outputs)]), collapse = ", ")
  )
}

message("Repository: ", repo_dir)
message("Legacy results: ", legacy_dir)
message("New results: ", new_dir)
message("Output: ", output_dir)
message("Datasets: ", paste(datasets, collapse = ", "))
message("Methods: ", paste(methods, collapse = ", "))

dataset_aliases <- function(dataset) {
  if (identical(dataset, "RECETOX")) {
    c("RECETOX")
  } else if (identical(dataset, "HREI-MSDB")) {
    c("HREI_MSDB", "HREI-MSDB")
  } else {
    stop("No directory aliases defined for dataset: ", dataset)
  }
}

find_dataset_dir <- function(root, dataset) {
  candidates <- file.path(root, dataset_aliases(dataset))
  hit <- candidates[dir.exists(candidates)]
  if (!length(hit)) {
    stop(
      "Missing ", dataset, " directory under ", root,
      "; tried: ", paste(candidates, collapse = ", ")
    )
  }
  normalizePath(hit[[1]], mustWork = TRUE)
}

first_existing <- function(paths, description) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop("Missing ", description, "; tried: ", paste(paths, collapse = ", "))
  }
  normalizePath(hit[[1]], mustWork = TRUE)
}

find_new_matrix <- function(dataset_dir, method) {
  first_existing(
    file.path(
      dataset_dir,
      c(
        paste0("dist_", method, "_tol15.rds"),
        paste0("dist_", method, ".rds"),
        paste0("dist_", method, "_tol15.csv"),
        paste0("dist_", method, ".csv")
      )
    ),
    paste0("new distance matrix for ", method)
  )
}

find_new_spectra <- function(dataset_dir) {
  first_existing(
    file.path(
      dataset_dir,
      c(
        "spectra_filtered_tol15.rds",
        "spectra_secondary_ablation.rds",
        "spectra_filtered.rds"
      )
    ),
    "new filtered spectra metadata"
  )
}

read_distance_matrix <- function(path) {
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    matrix <- readRDS(path)
    matrix <- as.matrix(matrix)
  } else {
    raw <- utils::read.csv(
      path, check.names = FALSE, stringsAsFactors = FALSE,
      na.strings = c("NA", "NaN", "Inf", "-Inf")
    )
    if (ncol(raw) < 2L) stop("Distance CSV has fewer than two columns: ", path)
    id_col <- if ("id" %in% names(raw)) "id" else names(raw)[[1]]
    ids <- as.character(raw[[id_col]])
    raw[[id_col]] <- NULL
    matrix <- suppressWarnings(as.matrix(data.frame(
      lapply(raw, as.numeric), check.names = FALSE
    )))
    rownames(matrix) <- ids
    colnames(matrix) <- names(raw)
  }

  if (!is.numeric(matrix)) storage.mode(matrix) <- "double"
  if (nrow(matrix) != ncol(matrix)) {
    stop("Distance matrix is not square: ", path, " [", paste(dim(matrix), collapse = " x "), "]")
  }
  if (is.null(rownames(matrix)) || is.null(colnames(matrix))) {
    stop("Distance matrix lacks row or column IDs: ", path)
  }
  if (anyNA(rownames(matrix)) || anyNA(colnames(matrix)) ||
      any(!nzchar(rownames(matrix))) || any(!nzchar(colnames(matrix))) ||
      anyDuplicated(rownames(matrix)) || anyDuplicated(colnames(matrix))) {
    stop("Distance matrix IDs are missing or duplicated: ", path)
  }
  if (!setequal(rownames(matrix), colnames(matrix))) {
    stop("Distance matrix row/column ID sets differ: ", path)
  }
  matrix <- matrix[rownames(matrix), rownames(matrix), drop = FALSE]
  if (any(!is.finite(matrix))) {
    stop("New distance matrix contains nonfinite cells: ", path)
  }
  matrix
}

read_spectra_metadata <- function(path) {
  spectra <- readRDS(path)
  if (is.null(spectra$df_spec) || !is.data.frame(spectra$df_spec)) {
    stop("Spectra RDS lacks a df_spec data frame: ", path)
  }
  required <- c("id", "inchikey")
  missing <- setdiff(required, names(spectra$df_spec))
  if (length(missing)) {
    stop("Spectra df_spec lacks column(s) ", paste(missing, collapse = ", "), ": ", path)
  }
  metadata <- data.frame(
    query_id = as.character(spectra$df_spec$id),
    inchikey = as.character(spectra$df_spec$inchikey),
    stringsAsFactors = FALSE
  )
  if (anyNA(metadata$query_id) || any(!nzchar(metadata$query_id)) ||
      anyDuplicated(metadata$query_id)) {
    stop("Spectra query IDs are missing or duplicated: ", path)
  }
  metadata$inchikey_prefix <- prefix_inchikey(metadata$inchikey)
  metadata
}

validate_matrix_metadata <- function(matrix, metadata, dataset, method) {
  missing_metadata <- setdiff(rownames(matrix), metadata$query_id)
  extra_metadata <- setdiff(metadata$query_id, rownames(matrix))
  if (length(missing_metadata) || length(extra_metadata)) {
    stop(
      dataset, "/", method, " matrix and spectra IDs differ; matrix-only=",
      length(missing_metadata), ", spectra-only=", length(extra_metadata)
    )
  }
  metadata[match(rownames(matrix), metadata$query_id), , drop = FALSE]
}

legacy_query_candidates <- function(dataset_dir, method) {
  file.path(
    dataset_dir,
    c(
      paste0("per_query_", method, "_tol15.csv"),
      paste0("per_query_", method, ".csv")
    )
  )
}

read_legacy_query <- function(dataset_dir, method, dataset) {
  candidates <- legacy_query_candidates(dataset_dir, method)
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) {
    path <- normalizePath(hit[[1]], mustWork = TRUE)
    legacy <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    path <- first_existing(
      file.path(
        dataset_dir,
        c("retrieval_per_query_metrics_tol15.csv", "retrieval_per_query_metrics.csv")
      ),
      paste0("legacy per-query results for ", method)
    )
    legacy <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"method" %in% names(legacy)) {
      stop("Combined legacy per-query file lacks method: ", path)
    }
    legacy <- legacy[as.character(legacy$method) == method, , drop = FALSE]
    if ("dataset" %in% names(legacy)) {
      accepted <- unique(c(dataset, dataset_aliases(dataset)))
      legacy <- legacy[as.character(legacy$dataset) %in% accepted, , drop = FALSE]
    }
  }

  required <- c("query_id", "top1", "rr", "p_at_1")
  missing <- setdiff(required, names(legacy))
  if (length(missing)) {
    stop("Legacy per-query file lacks column(s) ", paste(missing, collapse = ", "), ": ", path)
  }
  if (!nrow(legacy)) stop("Legacy per-query file has no rows for ", dataset, "/", method)
  legacy$query_id <- as.character(legacy$query_id)
  if (anyNA(legacy$query_id) || any(!nzchar(legacy$query_id)) ||
      anyDuplicated(legacy$query_id)) {
    stop("Legacy query IDs are missing or duplicated for ", dataset, "/", method)
  }

  inchikey <- if ("inchikey" %in% names(legacy)) {
    as.character(legacy$inchikey)
  } else {
    rep(NA_character_, nrow(legacy))
  }
  prefix <- if ("inchikey_prefix" %in% names(legacy)) {
    as.character(legacy$inchikey_prefix)
  } else {
    prefix_inchikey(inchikey)
  }
  standardized <- data.frame(
    query_id = legacy$query_id,
    old_inchikey = inchikey,
    old_inchikey_prefix = prefix,
    old_published_top1 = as.numeric(legacy$top1),
    old_published_mrr = as.numeric(legacy$rr),
    old_published_p_at_1 = as.numeric(legacy$p_at_1),
    stringsAsFactors = FALSE
  )
  list(data = standardized, path = path)
}

calculate_new_query <- function(matrix, metadata, random_seed) {
  aligned <- metadata[match(rownames(matrix), metadata$query_id), , drop = FALSE]
  inchikeys <- aligned$inchikey
  names(inchikeys) <- aligned$query_id
  metrics <- per_query_metrics_tie_aware(
    matrix, inchikeys, tie_tolerance = 0, random_seed = random_seed
  )
  if (!nrow(metrics)) stop("No eligible retrieval queries in new distance matrix.")
  if (anyDuplicated(metrics$query_id)) stop("New tie-aware query metrics contain duplicate IDs.")
  data.frame(
    query_id = as.character(metrics$query_id),
    new_inchikey = as.character(metrics$inchikey),
    new_inchikey_prefix = as.character(metrics$inchikey_prefix),
    new_optimistic_top1 = as.numeric(metrics$top1_optimistic),
    new_fractional_top1 = as.numeric(metrics$top1_fractional),
    new_optimistic_mrr = as.numeric(metrics$rr_optimistic),
    new_fractional_mrr = as.numeric(metrics$rr_fractional),
    new_optimistic_p_at_1 = as.numeric(metrics$p_at_1_optimistic),
    new_fractional_p_at_1 = as.numeric(metrics$p_at_1_fractional),
    stringsAsFactors = FALSE
  )
}

finite_delta <- function(new, old) {
  out <- rep(NA_real_, length(new))
  keep <- is.finite(new) & is.finite(old)
  out[keep] <- new[keep] - old[keep]
  out
}

same_identity <- function(old, new) {
  out <- rep(NA, length(old))
  keep <- !is.na(old) & nzchar(old) & !is.na(new) & nzchar(new)
  out[keep] <- old[keep] == new[keep]
  out
}

join_old_new_query <- function(old, new, dataset, method) {
  id_order <- unique(c(old$query_id, new$query_id))
  joined <- merge(old, new, by = "query_id", all = TRUE, sort = FALSE)
  joined <- joined[match(id_order, joined$query_id), , drop = FALSE]

  joined$old_query_present <- joined$query_id %in% old$query_id
  joined$new_query_present <- joined$query_id %in% new$query_id
  joined$query_match_status <- ifelse(
    joined$old_query_present & joined$new_query_present, "matched",
    ifelse(joined$old_query_present, "old_only", "new_only")
  )
  joined$inchikey_agreement <- same_identity(joined$old_inchikey, joined$new_inchikey)
  joined$inchikey_prefix_agreement <- same_identity(
    joined$old_inchikey_prefix, joined$new_inchikey_prefix
  )

  for (metric in c("top1", "mrr", "p_at_1")) {
    old_col <- paste0("old_published_", metric)
    optimistic_col <- paste0("new_optimistic_", metric)
    fractional_col <- paste0("new_fractional_", metric)
    joined[[paste0("old_eligible_", metric)]] <- is.finite(joined[[old_col]])
    joined[[paste0("new_eligible_", metric)]] <-
      is.finite(joined[[optimistic_col]]) & is.finite(joined[[fractional_col]])
    joined[[paste0("matched_eligible_", metric)]] <-
      joined[[paste0("old_eligible_", metric)]] &
      joined[[paste0("new_eligible_", metric)]]
    joined[[paste0("delta_new_optimistic_minus_old_published_", metric)]] <-
      finite_delta(joined[[optimistic_col]], joined[[old_col]])
    joined[[paste0("delta_new_fractional_minus_old_published_", metric)]] <-
      finite_delta(joined[[fractional_col]], joined[[old_col]])
  }

  joined$dataset <- dataset
  joined$method <- method
  joined$old_published_policy_top1_mrr <- "minimum-rank optimistic"
  joined$old_published_policy_p_at_1 <-
    "first sorted candidate; exact-tie order dependent"
  joined$new_tie_definition <- "full-precision exact equality"
  joined$new_primary_tie_policy <- "fractional random-order expectation"

  first <- c(
    "dataset", "method", "query_id", "query_match_status",
    "old_query_present", "new_query_present"
  )
  joined[, c(first, setdiff(names(joined), first)), drop = FALSE]
}

mean_finite <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

metric_specifications <- data.frame(
  metric = c("top1", "mrr", "p_at_1"),
  old_col = c(
    "old_published_top1", "old_published_mrr", "old_published_p_at_1"
  ),
  optimistic_col = c(
    "new_optimistic_top1", "new_optimistic_mrr", "new_optimistic_p_at_1"
  ),
  fractional_col = c(
    "new_fractional_top1", "new_fractional_mrr", "new_fractional_p_at_1"
  ),
  legacy_policy = c(
    "minimum-rank optimistic", "minimum-rank optimistic",
    "first sorted candidate; exact-tie order dependent"
  ),
  stringsAsFactors = FALSE
)

summarize_old_new <- function(joined, dataset, method) {
  rows <- lapply(seq_len(nrow(metric_specifications)), function(i) {
    spec <- metric_specifications[i, ]
    old <- joined[[spec$old_col]]
    optimistic <- joined[[spec$optimistic_col]]
    fractional <- joined[[spec$fractional_col]]
    matched <- is.finite(old) & is.finite(optimistic) & is.finite(fractional)
    optimistic_delta <- optimistic[matched] - old[matched]
    fractional_delta <- fractional[matched] - old[matched]
    data.frame(
      dataset = dataset,
      method = method,
      metric = spec$metric,
      legacy_policy = spec$legacy_policy,
      new_tie_definition = "full-precision exact equality",
      new_primary_tie_policy = "fractional random-order expectation",
      delta_direction = "new minus old published",
      change_tolerance = change_tolerance,
      n_old_all = sum(is.finite(old)),
      old_all_estimate = mean_finite(old),
      n_new_optimistic_all = sum(is.finite(optimistic)),
      new_optimistic_all_estimate = mean_finite(optimistic),
      n_new_fractional_all = sum(is.finite(fractional)),
      new_fractional_all_estimate = mean_finite(fractional),
      n_matched = sum(matched),
      old_matched_estimate = mean_finite(old[matched]),
      new_optimistic_matched_estimate = mean_finite(optimistic[matched]),
      new_fractional_matched_estimate = mean_finite(fractional[matched]),
      delta_optimistic_minus_old = mean_finite(optimistic_delta),
      delta_fractional_minus_old = mean_finite(fractional_delta),
      n_changed_optimistic_above_tolerance =
        sum(abs(optimistic_delta) > change_tolerance),
      n_changed_fractional_above_tolerance =
        sum(abs(fractional_delta) > change_tolerance),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_unmatched <- function(joined, dataset, method) {
  rows <- lapply(seq_len(nrow(metric_specifications)), function(i) {
    spec <- metric_specifications[i, ]
    old <- is.finite(joined[[spec$old_col]])
    new <- is.finite(joined[[spec$optimistic_col]]) &
      is.finite(joined[[spec$fractional_col]])
    identity_disagreement <- joined$old_query_present & joined$new_query_present &
      !is.na(joined$inchikey_agreement) & !joined$inchikey_agreement
    data.frame(
      dataset = dataset,
      method = method,
      metric = spec$metric,
      n_union_query_ids = nrow(joined),
      n_old_eligible = sum(old),
      n_new_eligible = sum(new),
      n_matched = sum(old & new),
      n_old_only = sum(old & !new),
      n_new_only = sum(!old & new),
      n_neither = sum(!old & !new),
      n_query_id_matched = sum(joined$old_query_present & joined$new_query_present),
      n_query_id_old_only = sum(joined$old_query_present & !joined$new_query_present),
      n_query_id_new_only = sum(!joined$old_query_present & joined$new_query_present),
      n_inchikey_disagreements_on_matched_ids = sum(identity_disagreement),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

reconstruct_legacy_mirroring <- function(directional) {
  mirrored <- directional
  lower <- lower.tri(mirrored)
  mirrored[lower] <- t(directional)[lower]
  diag(mirrored) <- 0
  if (!all(mirrored == t(mirrored))) {
    stop("Internal error: reconstructed triangular-mirroring matrix is not symmetric.")
  }
  mirrored
}

composite_metric_specifications <- data.frame(
  metric = c("top1", "mrr", "p_at_1"),
  optimistic_col = c("top1_optimistic", "rr_optimistic", "p_at_1_optimistic"),
  fractional_col = c("top1_fractional", "rr_fractional", "p_at_1_fractional"),
  stringsAsFactors = FALSE
)

composite_directionality_impact <- function(directional, metadata, dataset, random_seed) {
  mirrored <- reconstruct_legacy_mirroring(directional)
  aligned <- metadata[match(rownames(directional), metadata$query_id), , drop = FALSE]
  inchikeys <- aligned$inchikey
  names(inchikeys) <- aligned$query_id
  mirrored_query <- per_query_metrics_tie_aware(
    mirrored, inchikeys, tie_tolerance = 0, random_seed = random_seed
  )
  directional_query <- per_query_metrics_tie_aware(
    directional, inchikeys, tie_tolerance = 0, random_seed = random_seed
  )
  if (!setequal(mirrored_query$query_id, directional_query$query_id)) {
    stop(dataset, "/composite eligibility changed between matrices; this should be impossible.")
  }

  upper <- upper.tri(directional)
  absolute_asymmetry <- abs(directional[upper] - t(directional)[upper])
  max_asymmetry <- if (length(absolute_asymmetry)) max(absolute_asymmetry) else 0
  mean_asymmetry <- if (length(absolute_asymmetry)) mean(absolute_asymmetry) else 0
  n_asymmetric <- sum(absolute_asymmetry > 0)
  n_asymmetric_1e12 <- sum(absolute_asymmetry > 1e-12)

  query_rows <- list()
  aggregate_rows <- list()
  for (i in seq_len(nrow(composite_metric_specifications))) {
    spec <- composite_metric_specifications[i, ]
    rank_prefix <- if (identical(spec$metric, "p_at_1")) "prefix" else "full"
    add_rank_columns <- function(x) {
      n_better <- x[[paste0(rank_prefix, "_n_better")]]
      tie_size <- x[[paste0(rank_prefix, "_tie_size")]]
      relevant_in_tie <- x[[paste0(rank_prefix, "_relevant_in_tie")]]
      x$optimistic_rank_audit <- n_better + 1
      x$pessimistic_rank_audit <-
        n_better + (tie_size - relevant_in_tie) + 1
      x
    }
    mirrored_ranked <- add_rank_columns(mirrored_query)
    directional_ranked <- add_rank_columns(directional_query)
    left <- mirrored_query[, c(
      "query_id", "inchikey", "inchikey_prefix",
      spec$optimistic_col, spec$fractional_col
    ), drop = FALSE]
    left$mirrored_optimistic_rank <- mirrored_ranked$optimistic_rank_audit
    left$mirrored_pessimistic_rank <- mirrored_ranked$pessimistic_rank_audit
    names(left)[4:5] <- c("mirrored_optimistic", "mirrored_fractional")
    right <- directional_query[, c(
      "query_id", spec$optimistic_col, spec$fractional_col
    ), drop = FALSE]
    names(right)[2:3] <- c("directional_optimistic", "directional_fractional")
    right$directional_optimistic_rank <- directional_ranked$optimistic_rank_audit
    right$directional_pessimistic_rank <- directional_ranked$pessimistic_rank_audit
    query <- merge(left, right, by = "query_id", all = FALSE, sort = FALSE)
    query <- query[match(mirrored_query$query_id, query$query_id), , drop = FALSE]
    query$delta_optimistic <- finite_delta(
      query$directional_optimistic, query$mirrored_optimistic
    )
    query$delta_fractional <- finite_delta(
      query$directional_fractional, query$mirrored_fractional
    )

    query_rows[[i]] <- data.frame(
      record_level = "query",
      dataset = dataset,
      method = "composite",
      metric = spec$metric,
      query_id = query$query_id,
      inchikey = query$inchikey,
      inchikey_prefix = query$inchikey_prefix,
      comparison_scope =
        "corrected parser/same directional matrix; mirrored counterfactual versus directional",
      difference_direction = "directional minus mirrored",
      directional_orientation = "row=query; column=library",
      mirrored_counterfactual = "upper triangle copied to lower triangle",
      tie_definition = "full-precision exact equality",
      change_tolerance = change_tolerance,
      n_queries = NA_integer_,
      n_changed_optimistic_above_tolerance = NA_integer_,
      n_changed_fractional_above_tolerance = NA_integer_,
      mirrored_optimistic = query$mirrored_optimistic,
      directional_optimistic = query$directional_optimistic,
      delta_optimistic = query$delta_optimistic,
      mirrored_fractional = query$mirrored_fractional,
      directional_fractional = query$directional_fractional,
      delta_fractional = query$delta_fractional,
      mirrored_optimistic_rank = query$mirrored_optimistic_rank,
      mirrored_pessimistic_rank = query$mirrored_pessimistic_rank,
      directional_optimistic_rank = query$directional_optimistic_rank,
      directional_pessimistic_rank = query$directional_pessimistic_rank,
      n_asymmetric_unordered_pairs_exact = NA_integer_,
      n_asymmetric_unordered_pairs_above_1e_12 = NA_integer_,
      max_absolute_directional_asymmetry = NA_real_,
      mean_absolute_directional_asymmetry = NA_real_,
      stringsAsFactors = FALSE
    )

    keep_optimistic <- is.finite(query$mirrored_optimistic) &
      is.finite(query$directional_optimistic)
    keep_fractional <- is.finite(query$mirrored_fractional) &
      is.finite(query$directional_fractional)
    aggregate_rows[[i]] <- data.frame(
      record_level = "aggregate",
      dataset = dataset,
      method = "composite",
      metric = spec$metric,
      query_id = NA_character_,
      inchikey = NA_character_,
      inchikey_prefix = NA_character_,
      comparison_scope =
        "corrected parser/same directional matrix; mirrored counterfactual versus directional",
      difference_direction = "directional minus mirrored",
      directional_orientation = "row=query; column=library",
      mirrored_counterfactual = "upper triangle copied to lower triangle",
      tie_definition = "full-precision exact equality",
      change_tolerance = change_tolerance,
      n_queries = sum(keep_fractional),
      n_changed_optimistic_above_tolerance =
        sum(abs(query$delta_optimistic[keep_optimistic]) > change_tolerance),
      n_changed_fractional_above_tolerance =
        sum(abs(query$delta_fractional[keep_fractional]) > change_tolerance),
      mirrored_optimistic = mean_finite(query$mirrored_optimistic[keep_optimistic]),
      directional_optimistic = mean_finite(query$directional_optimistic[keep_optimistic]),
      delta_optimistic = mean_finite(query$delta_optimistic[keep_optimistic]),
      mirrored_fractional = mean_finite(query$mirrored_fractional[keep_fractional]),
      directional_fractional = mean_finite(query$directional_fractional[keep_fractional]),
      delta_fractional = mean_finite(query$delta_fractional[keep_fractional]),
      mirrored_optimistic_rank = mean_finite(query$mirrored_optimistic_rank),
      mirrored_pessimistic_rank = mean_finite(query$mirrored_pessimistic_rank),
      directional_optimistic_rank = mean_finite(query$directional_optimistic_rank),
      directional_pessimistic_rank = mean_finite(query$directional_pessimistic_rank),
      n_asymmetric_unordered_pairs_exact = n_asymmetric,
      n_asymmetric_unordered_pairs_above_1e_12 = n_asymmetric_1e12,
      max_absolute_directional_asymmetry = max_asymmetry,
      mean_absolute_directional_asymmetry = mean_asymmetry,
      stringsAsFactors = FALSE
    )
  }

  list(
    query = do.call(rbind, query_rows),
    aggregate = do.call(rbind, aggregate_rows)
  )
}

input_records <- list()
record_input <- function(role, dataset, method, path) {
  input_records[[length(input_records) + 1L]] <<- data.frame(
    role = role,
    dataset = dataset,
    method = method,
    path = normalizePath(path, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

query_rows <- list()
aggregate_rows <- list()
unmatched_rows <- list()
composite_query_rows <- list()
composite_aggregate_rows <- list()
row_index <- 0L

for (di in seq_along(datasets)) {
  dataset <- datasets[[di]]
  message("\n=== ", dataset, " ===")
  legacy_dataset_dir <- tryCatch(
    find_dataset_dir(legacy_dir, dataset),
    error = function(e) {
      if (!allow_partial) stop(e)
      message("Skipping ", dataset, ": ", conditionMessage(e))
      NULL
    }
  )
  new_dataset_dir <- tryCatch(
    find_dataset_dir(new_dir, dataset),
    error = function(e) {
      if (!allow_partial) stop(e)
      message("Skipping ", dataset, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(legacy_dataset_dir) || is.null(new_dataset_dir)) next

  spectra_path <- tryCatch(
    find_new_spectra(new_dataset_dir),
    error = function(e) {
      if (!allow_partial) stop(e)
      message("Skipping ", dataset, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(spectra_path)) next
  metadata <- read_spectra_metadata(spectra_path)
  record_input("new_spectra_metadata", dataset, NA_character_, spectra_path)

  for (mi in seq_along(methods)) {
    method <- methods[[mi]]
    message("[", dataset, "] ", method)
    components <- tryCatch({
      legacy <- read_legacy_query(legacy_dataset_dir, method, dataset)
      matrix_path <- find_new_matrix(new_dataset_dir, method)
      matrix <- read_distance_matrix(matrix_path)
      aligned_metadata <- validate_matrix_metadata(matrix, metadata, dataset, method)
      list(
        legacy = legacy,
        matrix_path = matrix_path,
        matrix = matrix,
        metadata = aligned_metadata
      )
    }, error = function(e) {
      if (!allow_partial) stop(e)
      message("Skipping ", dataset, "/", method, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(components)) next

    row_index <- row_index + 1L
    record_input("legacy_published_per_query", dataset, method, components$legacy$path)
    record_input("new_full_precision_distance_matrix", dataset, method, components$matrix_path)
    new_query <- calculate_new_query(
      components$matrix,
      components$metadata,
      random_seed = seed + di * 100000L + mi * 1000L
    )
    joined <- join_old_new_query(
      components$legacy$data, new_query, dataset, method
    )
    query_rows[[row_index]] <- joined
    aggregate_rows[[row_index]] <- summarize_old_new(joined, dataset, method)
    unmatched_rows[[row_index]] <- summarize_unmatched(joined, dataset, method)

    if (identical(method, "composite")) {
      impact <- composite_directionality_impact(
        components$matrix,
        components$metadata,
        dataset,
        random_seed = seed + di * 100000L + 90000L
      )
      composite_query_rows[[dataset]] <- impact$query
      composite_aggregate_rows[[dataset]] <- impact$aggregate
    }
  }
}

if (!length(query_rows)) stop("No dataset/method comparisons were completed.")
if (!allow_partial) {
  expected <- length(datasets) * length(methods)
  if (length(query_rows) != expected) {
    stop("Expected ", expected, " dataset/method comparisons but completed ", length(query_rows), ".")
  }
  if ("composite" %in% methods && length(composite_query_rows) != length(datasets)) {
    stop("Composite directionality output is incomplete.")
  }
}

query_output <- do.call(rbind, query_rows)
aggregate_output <- do.call(rbind, aggregate_rows)
unmatched_output <- do.call(rbind, unmatched_rows)

empty_composite <- data.frame(
  record_level = character(), dataset = character(), method = character(),
  metric = character(), query_id = character(), inchikey = character(),
  inchikey_prefix = character(), comparison_scope = character(),
  difference_direction = character(), directional_orientation = character(),
  mirrored_counterfactual = character(), tie_definition = character(),
  change_tolerance = numeric(), n_queries = integer(),
  n_changed_optimistic_above_tolerance = integer(),
  n_changed_fractional_above_tolerance = integer(), mirrored_optimistic = numeric(),
  directional_optimistic = numeric(), delta_optimistic = numeric(),
  mirrored_fractional = numeric(), directional_fractional = numeric(),
  delta_fractional = numeric(), mirrored_optimistic_rank = numeric(),
  mirrored_pessimistic_rank = numeric(), directional_optimistic_rank = numeric(),
  directional_pessimistic_rank = numeric(),
  n_asymmetric_unordered_pairs_exact = integer(),
  n_asymmetric_unordered_pairs_above_1e_12 = integer(),
  max_absolute_directional_asymmetry = numeric(),
  mean_absolute_directional_asymmetry = numeric(),
  stringsAsFactors = FALSE
)
composite_query_output <- if (length(composite_query_rows)) {
  do.call(rbind, composite_query_rows)
} else {
  empty_composite
}
composite_aggregate_output <- if (length(composite_aggregate_rows)) {
  do.call(rbind, composite_aggregate_rows)
} else {
  empty_composite
}
composite_output <- rbind(composite_query_output, composite_aggregate_output)

input_manifest <- do.call(rbind, input_records)
input_manifest$bytes <- as.numeric(file.info(input_manifest$path)$size)
input_manifest$md5 <- unname(tools::md5sum(input_manifest$path))

run_parameters <- data.frame(
  parameter = c(
    "timestamp_utc", "package_commit", "package_tree_dirty", "command",
    "repo_dir", "bundle_root", "legacy_dir", "new_dir",
    "output_dir", "datasets", "methods", "seed", "change_tolerance", "allow_partial",
    "legacy_top1_mrr_policy", "legacy_p_at_1_policy", "new_tie_definition",
    "new_primary_tie_policy", "aggregate_delta_population",
    "composite_counterfactual"
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
    repo_dir, bundle_root, legacy_dir, new_dir,
    output_dir, paste(datasets, collapse = ","), paste(methods, collapse = ","),
    seed, change_tolerance, allow_partial, "minimum-rank optimistic",
    "first sorted candidate; exact-tie order dependent",
    "full-precision exact equality", "fractional random-order expectation",
    "query-ID and metric-eligibility matched rows only",
    "same corrected-parser directional matrix; upper triangle copied to lower"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  query_output, file.path(output_dir, "old_new_retrieval_query_level.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  aggregate_output, file.path(output_dir, "old_new_retrieval_aggregate.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  unmatched_output, file.path(output_dir, "old_new_retrieval_unmatched_counts.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  composite_output, file.path(output_dir, "composite_directionality_impact.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  composite_query_output,
  file.path(output_dir, "composite_directionality_impact_query_level.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  composite_aggregate_output,
  file.path(output_dir, "composite_directionality_impact_aggregate.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  input_manifest,
  file.path(output_dir, "old_new_retrieval_audit_input_files.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  run_parameters,
  file.path(output_dir, "old_new_retrieval_audit_run_parameters.csv"),
  row.names = FALSE, na = ""
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "old_new_retrieval_audit_sessionInfo.txt"),
  useBytes = TRUE
)
saveRDS(
  list(
    query_level = query_output,
    aggregate = aggregate_output,
    unmatched = unmatched_output,
    composite_directionality = composite_output,
    inputs = input_manifest,
    parameters = run_parameters
  ),
  file.path(output_dir, "old_new_retrieval_audit_results.rds"),
  compress = TRUE
)

message("\nAudit complete.")
message("Dataset/method comparisons: ", length(query_rows))
message("Query-level rows: ", nrow(query_output))
message("Aggregate rows: ", nrow(aggregate_output))
message("Composite impact rows: ", nrow(composite_output))
message("Output: ", output_dir)
