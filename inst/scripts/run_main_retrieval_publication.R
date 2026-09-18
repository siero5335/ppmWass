#!/usr/bin/env Rscript

# Publication-grade clean replicate-retrieval rerun. This driver uses the
# full-precision MSP parser, explicit query-to-library orientation for the
# directional composite score, exact/fractional tie handling, and clustered
# bootstrap confidence intervals.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]])
}

flag_is_true <- function(name, default = FALSE) {
  val <- get_arg(name, if (default) "true" else "false")
  tolower(val) %in% c("1", "true", "yes", "y")
}

split_arg <- function(name, default) {
  trimws(strsplit(get_arg(name, default), ",", fixed = TRUE)[[1]])
}

`%||%` <- function(x, y) if (is.null(x)) y else x

repo_dir <- normalizePath(get_arg("repo", getwd()), mustWork = TRUE)
publication_bundle_root <- normalizePath(
  file.path(repo_dir, "..", ".."), mustWork = TRUE
)
provenance_path <- function(path) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  normalized <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(publication_bundle_root, .Platform$file.sep)
  if (startsWith(normalized, prefix)) {
    substring(normalized, nchar(prefix) + 1L)
  } else normalized
}
gc_dir <- get_arg("gc-dir", "/Users/a_eguchi/Desktop/gc")
n_boot <- as.integer(get_arg("boot", "1000"))
cluster_boot <- as.integer(get_arg("cluster-boot", "20000"))
seed <- as.integer(get_arg("seed", "20260430"))
tol_ppm <- as.numeric(get_arg("tol-ppm", "15"))
n_cores <- as.integer(get_arg("n-cores", "8"))
ot_method <- tolower(get_arg("ot-method", "exact"))
if (!ot_method %in% c("exact", "sinkhorn", "greenkhorn")) {
  stop("--ot-method must be exact, sinkhorn, or greenkhorn.")
}
spectra_rds <- get_arg("spectra-rds", "")
dataset_label_override <- get_arg("dataset-label", "")
reuse_distance_dir <- get_arg("reuse-distance-dir", "")
default_file_tag <- paste0(
  "tol", gsub("[^0-9A-Za-z]+", "p", format(tol_ppm, scientific = FALSE, trim = TRUE))
)
file_tag <- get_arg("file-tag", default_file_tag)
if (!grepl("^[A-Za-z0-9._-]+$", file_tag)) {
  stop("--file-tag may contain only letters, digits, dot, underscore, and hyphen.")
}
tagged_name <- function(stem, extension = ".csv") {
  paste0(stem, "_", file_tag, extension)
}
output_dir <- get_arg(
  "output-dir", file.path(gc_dir, paste0("benchmark_output_main_", file_tag))
)
output_dir <- normalizePath(output_dir, mustWork = FALSE)
if (nzchar(reuse_distance_dir)) {
  reuse_distance_dir <- normalizePath(reuse_distance_dir, mustWork = TRUE)
  if (identical(reuse_distance_dir, output_dir)) {
    stop("--reuse-distance-dir must differ from --output-dir.")
  }
}
log_file <- get_arg("log-file", file.path(output_dir, tagged_name("run", ".log")))

if (nzchar(spectra_rds)) {
  spectra_rds <- normalizePath(spectra_rds, mustWork = TRUE)
  datasets <- if (!nzchar(dataset_label_override)) "PREPARED" else dataset_label_override
} else {
  datasets <- split_arg("datasets", "RECETOX,HREI-MSDB")
  unknown_datasets <- setdiff(datasets, c("RECETOX", "HREI-MSDB"))
  if (length(unknown_datasets)) {
    stop("Unknown --datasets value(s): ", paste(unknown_datasets, collapse = ", "))
  }
}

resolve_msp_input <- function(path, required) {
  if (isTRUE(required)) return(normalizePath(path, mustWork = TRUE))
  # An unselected MSP is not an input merely because a machine happens to have
  # the default file. Keep prepared-mode provenance independent of host state.
  NA_character_
}
recetox_msp <- resolve_msp_input(
  get_arg("recetox-msp", file.path(gc_dir, "RECETOX_merged.msp")),
  !nzchar(spectra_rds) && "RECETOX" %in% datasets
)
hrei_msp <- resolve_msp_input(
  get_arg("hrei-msp", file.path(gc_dir, "HREI-MSDB.msp")),
  !nzchar(spectra_rds) && "HREI-MSDB" %in% datasets
)

methods_all <- split_arg(
  "methods",
  "ppm_wasserstein,composite,entropy_weighted,entropy_unweighted,cosine,weighted_cosine,hellinger"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (nzchar(log_file)) {
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
if (identical(ot_method, "exact") &&
    !requireNamespace("transport", quietly = TRUE)) {
  stop("--ot-method=exact requires transport; publication runs fail closed.")
}

stats_helper <- file.path(repo_dir, "inst", "scripts", "lib",
                          "publication_retrieval_statistics.R")
if (!file.exists(stats_helper)) stop("Missing statistics helper: ", stats_helper)
source(stats_helper, local = TRUE)
RNGkind("L'Ecuyer-CMRG")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required for figure generation.")
}

message("Repository: ", repo_dir)
message("Output dir: ", output_dir)
message("tol_ppm: ", tol_ppm)
message("Methods: ", paste(methods_all, collapse = ", "))
message("Query bootstrap R: ", n_boot)
message("Cluster bootstrap R: ", cluster_boot)
message(
  "Distance matrices: ",
  if (nzchar(reuse_distance_dir)) paste0("verified reuse from ", reuse_distance_dir) else
    "computed in this run"
)

git_commit <- tryCatch(
  suppressWarnings(trimws(system2(
    "git", c("-C", repo_dir, "rev-parse", "HEAD"),
    stdout = TRUE, stderr = FALSE
  ))),
  error = function(e) NA_character_
)
if (length(git_commit) != 1L || is.na(git_commit) || !nzchar(git_commit)) {
  git_commit <- NA_character_
}
git_status <- tryCatch(
  system2(
    "git", c("-C", repo_dir, "status", "--porcelain"),
    stdout = TRUE, stderr = FALSE
  ),
  error = function(e) character()
)
package_tree_dirty <- length(git_status) > 0L
distance_generation_commit <- git_commit
distance_source_parameter_file <- NA_character_
if (nzchar(reuse_distance_dir)) {
  distance_source_parameter_file <- file.path(
    reuse_distance_dir, tagged_name("main_run_parameters", ".csv")
  )
  if (!file.exists(distance_source_parameter_file)) {
    stop("Missing source run parameters for reused matrices: ",
         distance_source_parameter_file)
  }
  source_parameters <- utils::read.csv(
    distance_source_parameter_file, stringsAsFactors = FALSE,
    check.names = FALSE
  )
  source_commit_values <- unique(as.character(
    source_parameters$value[source_parameters$parameter == "package_commit"]
  ))
  source_dirty_values <- unique(tolower(as.character(
    source_parameters$value[source_parameters$parameter == "package_tree_dirty"]
  )))
  source_solver_values <- unique(as.character(
    source_parameters$value[source_parameters$parameter == "ot_method"]
  ))
  if (length(source_commit_values) != 1L || is.na(source_commit_values) ||
      !nzchar(source_commit_values) || length(source_dirty_values) != 1L ||
      !source_dirty_values %in% c("false", "0") ||
      length(source_solver_values) != 1L ||
      !identical(source_solver_values, ot_method)) {
    stop(
      "Reused distance matrices lack a unique clean source commit or matching OT backend."
    )
  }
  distance_generation_commit <- source_commit_values
}
exact_command <- paste(
  c(file.path(R.home("bin"), "Rscript"), commandArgs()), collapse = " "
)
has_approxOT <- requireNamespace("approxOT", quietly = TRUE)
has_transport <- requireNamespace("transport", quietly = TRUE)
effective_ot_backend <- if (identical(ot_method, "exact")) {
  "transport exact unregularized OT"
} else {
  paste0("approxOT finite-iteration ", ot_method,
         " with exact fallback on invalid plans")
}
nonfinite_fallback_sequence <- if (identical(ot_method, "exact")) {
  paste(
    "Use the analytic identity for constant ground-cost matrices; otherwise",
    "use transport exact unregularized OT. Abort on solver error, malformed",
    "plan, marginal validation failure, or any nonfinite selected value.",
    "No approximate or Hellinger substitution is permitted."
  )
} else {
  paste(
    "Use the analytic constant-cost shortcut; run the requested approximate",
    "backend at 100, 500, and 2000 iterations; use exact OT only if every",
    "approximate plan fails validation; otherwise stop."
  )
}
run_parameters <- data.frame(
  parameter = c(
    "timestamp_utc", "package_commit", "package_tree_dirty", "command",
    "repo_dir", "file_tag", "recetox_msp",
    "hrei_msp", "prepared_spectra_rds", "datasets", "methods",
    "query_bootstrap_R", "cluster_bootstrap_R", "seed",
    "hard_match_tolerance_ppm", "ppmWass_base_cost_ppm",
    "transition_multiplier", "saturation_width_ppm",
    "ot_method", "ot_estimand", "sinkhorn_epsilon", "sinkhorn_iterations",
    "ot_marginal_residual_tolerance", "approxOT_available",
    "transport_available", "solver_backend",
    "nonfinite_fallback", "nonfinite_matrix_policy", "parallel", "n_cores",
    "tie_primary", "tie_sensitivity", "tie_affected_definition",
    "cluster_ci", "cluster_estimand", "distance_matrix_mode",
    "distance_generation_commit", "statistics_recompute_commit",
    "distance_reuse_source", "distance_reuse_validation"
  ),
  value = c(
    format(Sys.time(), tz = "UTC"), git_commit, package_tree_dirty,
    exact_command, repo_dir, file_tag, recetox_msp,
    hrei_msp, if (nzchar(spectra_rds)) spectra_rds else NA_character_,
    paste(datasets, collapse = ","),
    paste(methods_all, collapse = ","), n_boot, cluster_boot, seed,
    tol_ppm, tol_ppm, 3, tol_ppm * 3, ot_method,
    if (identical(ot_method, "exact")) "unregularized_exact_transport_cost" else
      "finite_iteration_regularized_plan_transport_cost",
    if (identical(ot_method, "exact")) NA else 0.05,
    if (identical(ot_method, "exact")) NA else 100, 1e-8,
    has_approxOT, has_transport, effective_ot_backend,
    nonfinite_fallback_sequence,
    "abort on any nonfinite matrix cell; no matrix-cell replacement",
    flag_is_true("parallel", TRUE), n_cores,
    "full-precision exact equality; fractional expected value",
    "distance rounded to 12 decimal digits",
    paste(
      "optimistic hit differs from pessimistic hit at Top-1, Top-5, Top-10,",
      "or prefix P@1; tie_size > 1 alone is not classified as affected"
    ),
    "95% percentile, quantile type 8, R=20000 by default",
    paste(
      "query-weighted mean conditional on fixed candidate library;",
      "full-InChIKey clusters for Top-1/MRR and 14-character InChIKey-prefix",
      "clusters for P@1"
    ),
    if (nzchar(reuse_distance_dir)) "verified_reuse" else "computed_current_run",
    distance_generation_commit, git_commit,
    if (nzchar(reuse_distance_dir)) reuse_distance_dir else NA_character_,
    if (nzchar(reuse_distance_dir)) paste(
      "source and destination RDS/CSV MD5 identical; RDS square numeric,",
      "spectrum IDs exact, all cells finite, diagonal zero"
    ) else "computed under statistics_recompute_commit"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  run_parameters, file.path(output_dir, tagged_name("main_run_parameters", ".csv")),
  row.names = FALSE
)
checksum_files <- c(recetox_msp, hrei_msp, stats_helper)
if (nzchar(spectra_rds)) checksum_files <- c(checksum_files, spectra_rds)
if (nzchar(reuse_distance_dir)) {
  reused_distance_inputs <- unlist(lapply(datasets, function(dataset) {
    dataset_subdir <- gsub("[^A-Za-z0-9]+", "_", dataset)
    unlist(lapply(methods_all, function(method) {
      file.path(
        reuse_distance_dir, dataset_subdir,
        c(
          tagged_name(paste0("dist_", method), ".rds"),
          tagged_name(paste0("dist_", method), ".csv")
        )
      )
    }), use.names = FALSE)
  }), use.names = FALSE)
  missing_reused_inputs <- reused_distance_inputs[!file.exists(reused_distance_inputs)]
  if (length(missing_reused_inputs)) {
    stop("Missing reused distance input(s):\n  ",
         paste(missing_reused_inputs, collapse = "\n  "))
  }
  checksum_files <- c(checksum_files, distance_source_parameter_file,
                      reused_distance_inputs)
}
checksum_files <- unique(checksum_files[
  !is.na(checksum_files) & nzchar(checksum_files) & file.exists(checksum_files)
])
utils::write.csv(
  data.frame(
    file = checksum_files,
    md5 = unname(tools::md5sum(checksum_files)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, tagged_name("main_input_checksums", ".csv")),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, tagged_name("main_sessionInfo", ".txt")),
  useBytes = TRUE
)
writeLines(
  exact_command,
  file.path(output_dir, tagged_name("main_command", ".txt")), useBytes = TRUE
)

method_labels <- c(
  ppm_wasserstein = "ppm-Wasserstein",
  composite = "Composite",
  entropy_weighted = "Entropy (weighted)",
  entropy_unweighted = "Entropy (unweighted)",
  cosine = "Cosine",
  weighted_cosine = "Weighted cosine",
  hellinger = "Hellinger"
)

method_colors <- c(
  ppm_wasserstein = "#D62728",
  composite = "#9467BD",
  entropy_weighted = "#1F77B4",
  entropy_unweighted = "#17BECF",
  cosine = "#7F7F7F",
  weighted_cosine = "#8C564B",
  hellinger = "#2CA02C"
)

base_params <- function() {
  params <- eihrms_default_params()
  params$min_mz <- 35
  params$max_mz <- 650
  params$noise_thr <- 0.01
  params$topK <- 200L
  params$class_detection_ppm <- 15
  params$tol_ppm <- tol_ppm
  params$ot_method <- ot_method
  params$sinkhorn_epsilon <- 0.05
  params$sinkhorn_niter <- 100L
  params$wasserstein_transition_mult <- 3
  params$use_parallel <- flag_is_true("parallel", TRUE)
  params$n_cores <- n_cores
  params
}

subset_spectra_local <- function(spectra_obj, idx) {
  out <- list(
    df_spec = spectra_obj$df_spec[idx, , drop = FALSE],
    frag_list = spectra_obj$frag_list[idx],
    loss_list = spectra_obj$loss_list[idx],
    ri = spectra_obj$ri[idx]
  )
  for (nm in c("loss_typ_list", "mref_confidence", "loss_anchor_list", "loss_pair_list",
               "loss_anchor_typ_list", "loss_pair_typ_list")) {
    if (!is.null(spectra_obj[[nm]])) out[[nm]] <- spectra_obj[[nm]][idx]
  }
  out
}

filter_valid_spectra_local <- function(spectra_obj) {
  valid <- vapply(spectra_obj$frag_list, function(f) {
    if (is.null(f) || !is.matrix(f) || nrow(f) < 2) return(FALSE)
    if (ncol(f) >= 2 && all(f[, 2] == 0)) return(FALSE)
    TRUE
  }, logical(1))
  if (all(valid)) return(spectra_obj)
  message("filter_valid_spectra_local: removed ", sum(!valid), " degenerate entries")
  subset_spectra_local(spectra_obj, which(valid))
}

normalize_hrei_msp_names <- function(input, output) {
  lines <- readLines(input, warn = FALSE)
  starts <- grep("^Name:", lines)
  ends <- c(starts[-1] - 1L, length(lines))

  for (i in seq_along(starts)) {
    block <- lines[starts[i]:ends[i]]
    db_line <- grep("^DB#:", block, value = TRUE)
    db_id <- if (length(db_line)) {
      trimws(sub("^DB#:[[:space:]]*", "", db_line[[1]]))
    } else {
      as.character(i)
    }
    name <- trimws(sub("^Name:[[:space:]]*", "", lines[starts[i]]))
    lines[starts[i]] <- paste0("Name: ", name, " __HREI_DB", db_id)
  }

  writeLines(lines, output, useBytes = TRUE)
  output
}

load_dataset <- function(dataset, params) {
  if (dataset == "RECETOX") {
    sp_all <- build_spectra_from_msp(recetox_msp, params, require_ri = FALSE, progress = TRUE)
  } else if (dataset == "HREI-MSDB") {
    unique_msp <- file.path(output_dir, tagged_name("HREI-MSDB_unique_names", ".msp"))
    normalize_hrei_msp_names(hrei_msp, unique_msp)
    sp_all <- build_spectra_from_msp(unique_msp, params, require_ri = FALSE, progress = TRUE)
  } else {
    stop("Unknown dataset: ", dataset)
  }

  df <- sp_all$df_spec
  idx <- which(!is.na(df$RI) & !is.na(df$inchikey) & nchar(df$inchikey) >= 14)
  sp <- filter_valid_spectra_local(subset_spectra_local(sp_all, idx))
  sp$raw_msp_metadata <- sp_all$msp_metadata
  sp$msp_parser_diagnostics <- attr(sp_all$msp_metadata, "msp_parser_diagnostics")
  names(sp$frag_list) <- sp$df_spec$id
  names(sp$loss_list) <- sp$df_spec$id
  sp
}

sanitize_dist_matrix_local <- function(dm, dataset, method) {
  bad <- !is.finite(dm)
  if (!any(bad)) {
    return(list(mat = dm, n_nonfinite = 0L, replacement = NA_real_))
  }
  stop(
    "Publication run aborted: ", dataset, "/", method, " produced ",
    sum(bad), " nonfinite matrix cells after deterministic pair fallback.",
    call. = FALSE
  )
}

write_dist_csv <- function(dm, path) {
  out <- data.frame(id = rownames(dm), dm, check.names = FALSE)
  utils::write.csv(out, path, row.names = FALSE)
}

prefix14 <- function(x, n = 14) {
  x <- as.character(x)
  out <- substr(x, 1, n)
  out[is.na(x) | x == "" | nchar(x) < n] <- NA_character_
  out
}

per_query_metrics <- function(dist_mat, inchikeys, inchikey_chars = 14) {
  ids <- rownames(dist_mat)
  names(inchikeys) <- names(inchikeys) %||% ids
  inchikeys <- inchikeys[ids]

  full_counts <- table(inchikeys[!is.na(inchikeys) & inchikeys != ""])
  multi_full <- names(full_counts[full_counts > 1])
  q_idx <- which(inchikeys %in% multi_full)

  ik_prefix <- prefix14(inchikeys, inchikey_chars)
  prefix_counts <- table(ik_prefix[!is.na(ik_prefix)])
  multi_prefix <- names(prefix_counts[prefix_counts > 1])
  q_prefix_idx <- which(ik_prefix %in% multi_prefix)

  rep_rows <- vector("list", length(q_idx))
  for (ii in seq_along(q_idx)) {
    qi <- q_idx[ii]
    same_full <- which(inchikeys == inchikeys[qi])
    other_full <- setdiff(same_full, qi)
    dists <- dist_mat[qi, ]
    dists[qi] <- Inf
    ranks <- rank(dists, ties.method = "min", na.last = "keep")
    best_rank <- suppressWarnings(min(ranks[other_full], na.rm = TRUE))
    if (!is.finite(best_rank)) best_rank <- NA_real_
    rep_rows[[ii]] <- data.frame(
      query_id = ids[qi],
      inchikey = inchikeys[qi],
      top1 = as.integer(is.finite(best_rank) && best_rank == 1),
      top5 = as.integer(is.finite(best_rank) && best_rank <= 5),
      rr = if (is.finite(best_rank)) 1 / best_rank else NA_real_,
      n_full_relevant = length(other_full),
      stringsAsFactors = FALSE
    )
  }
  rep_df <- if (length(rep_rows)) do.call(rbind, rep_rows) else data.frame()

  pk_rows <- vector("list", length(q_prefix_idx))
  for (ii in seq_along(q_prefix_idx)) {
    qi <- q_prefix_idx[ii]
    dists <- dist_mat[qi, ]
    dists[qi] <- Inf
    ord <- order(dists, na.last = TRUE)
    ord <- ord[ord != qi]
    rel <- ik_prefix[ord] == ik_prefix[qi]
    rel[is.na(rel)] <- FALSE
    n_rel <- sum(!is.na(ik_prefix) & ik_prefix == ik_prefix[qi]) - 1
    p1 <- if (length(rel)) as.numeric(rel[[1]]) else NA_real_
    ap <- 0
    if (n_rel > 0 && any(rel)) {
      cum_rel <- cumsum(rel)
      precision_at_rank <- cum_rel / seq_along(cum_rel)
      ap <- sum(precision_at_rank[rel]) / n_rel
    }
    pk_rows[[ii]] <- data.frame(
      query_id = ids[qi],
      inchikey_prefix = ik_prefix[qi],
      p_at_1 = p1,
      ap = ap,
      n_prefix_relevant = n_rel,
      stringsAsFactors = FALSE
    )
  }
  pk_df <- if (length(pk_rows)) do.call(rbind, pk_rows) else data.frame()

  merge(rep_df, pk_df, by = "query_id", all = TRUE)
}

bootstrap_ci <- function(x, R = 1000, seed = 1L) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n = 0))
  }
  set.seed(seed)
  vals <- replicate(R, mean(sample(x, length(x), replace = TRUE), na.rm = TRUE))
  c(
    estimate = mean(x, na.rm = TRUE),
    ci_low = unname(stats::quantile(vals, 0.025, na.rm = TRUE)),
    ci_high = unname(stats::quantile(vals, 0.975, na.rm = TRUE)),
    n = length(x)
  )
}

summarize_per_query <- function(per_query, dataset, method, R = 1000, seed = 1L) {
  metrics <- c("top1", "top5", "rr", "p_at_1", "ap")
  rows <- lapply(seq_along(metrics), function(i) {
    metric <- metrics[[i]]
    ci <- bootstrap_ci(per_query[[metric]], R = R, seed = seed + i)
    data.frame(
      dataset = dataset,
      method = method,
      metric = metric,
      estimate = ci[["estimate"]],
      ci_low = ci[["ci_low"]],
      ci_high = ci[["ci_high"]],
      n_queries = as.integer(ci[["n"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

mcnemar_dataset <- function(dataset, per_method, reference = "ppm_wasserstein") {
  if (!reference %in% names(per_method)) {
    return(data.frame(
      dataset = character(),
      reference = character(),
      comparator = character(),
      n = integer(),
      ref_only = integer(),
      comparator_only = integer(),
      both_correct = integer(),
      both_wrong = integer(),
      mcnemar_p = numeric(),
      exact_discordant_p = numeric(),
      analysis_level = character(),
      status = character(),
      inference_method = character(),
      stringsAsFactors = FALSE
    ))
  }
  ref <- per_method[[reference]][, c("query_id", "top1")]
  names(ref)[2] <- "top1_ref"
  rows <- list()
  for (method in setdiff(names(per_method), reference)) {
    cur <- per_method[[method]][, c("query_id", "top1")]
    names(cur)[2] <- "top1_method"
    x <- merge(ref, cur, by = "query_id")
    x <- x[is.finite(x$top1_ref) & is.finite(x$top1_method), , drop = FALSE]
    tab <- table(
      reference = factor(x$top1_ref, levels = c(0, 1)),
      comparator = factor(x$top1_method, levels = c(0, 1))
    )
    p <- tryCatch(stats::mcnemar.test(tab, correct = TRUE)$p.value,
                  error = function(e) NA_real_)
    exact_p <- tryCatch({
      discordant <- as.integer(tab["1", "0"] + tab["0", "1"])
      if (!discordant) {
        NA_real_
      } else {
        stats::binom.test(as.integer(tab["1", "0"]), discordant, p = 0.5)$p.value
      }
    }, error = function(e) NA_real_)
    rows[[method]] <- data.frame(
      dataset = dataset,
      reference = reference,
      comparator = method,
      n = nrow(x),
      ref_only = as.integer(tab["1", "0"]),
      comparator_only = as.integer(tab["0", "1"]),
      both_correct = as.integer(tab["1", "1"]),
      both_wrong = as.integer(tab["0", "0"]),
      mcnemar_p = p,
      exact_discordant_p = exact_p,
      analysis_level = "query",
      status = "exploratory",
      inference_method = "mcnemar_optimistic_top1",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

run_dataset <- function(dataset) {
  params <- base_params()
  dataset_dir <- file.path(output_dir, gsub("[^A-Za-z0-9]+", "_", dataset))
  dir.create(dataset_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n=== Dataset: ", dataset, " ===")
  spectra <- if (nzchar(spectra_rds)) {
    message("Loading prepared spectra: ", spectra_rds)
    readRDS(normalizePath(spectra_rds, mustWork = TRUE))
  } else {
    load_dataset(dataset, params)
  }
  saveRDS(spectra, file.path(dataset_dir, tagged_name("spectra_filtered", ".rds")))

  parser_diag <- spectra$msp_parser_diagnostics
  if (!is.null(parser_diag) && !is.null(parser_diag$embedded_record_repairs)) {
    utils::write.csv(
      parser_diag$embedded_record_repairs,
      file.path(dataset_dir, tagged_name("msp_embedded_record_repairs", ".csv")),
      row.names = FALSE
    )
  }
  if (!is.null(parser_diag) && !is.null(parser_diag$record_issues)) {
    utils::write.csv(
      parser_diag$record_issues,
      file.path(dataset_dir, tagged_name("msp_record_issues", ".csv")),
      row.names = FALSE
    )
  }

  inchikeys <- spectra$df_spec$inchikey
  names(inchikeys) <- spectra$df_spec$id
  qc <- data.frame(
    dataset = dataset,
    metric = c("raw_msp_records", "filtered_spectra", "ri_available", "inchikey_available",
               "duplicate_full_inchikey_spectra", "duplicate_prefix_spectra",
               "hard_match_tolerance_ppm", "ppmWass_base_cost_ppm",
               "transition_multiplier", "saturation_width_ppm",
               "sinkhorn_epsilon", "sinkhorn_iterations", "n_cores"),
    value = c(
      if (!is.null(spectra$raw_msp_metadata)) nrow(spectra$raw_msp_metadata) else nrow(spectra$df_spec),
      nrow(spectra$df_spec),
      sum(!is.na(spectra$df_spec$RI)),
      sum(!is.na(spectra$df_spec$inchikey) & nchar(spectra$df_spec$inchikey) >= 14),
      sum(table(inchikeys)[table(inchikeys) >= 2]),
      sum(table(prefix14(inchikeys))[table(prefix14(inchikeys)) >= 2]),
      tol_ppm,
      tol_ppm,
      params$wasserstein_transition_mult,
      tol_ppm * params$wasserstein_transition_mult,
      if (identical(params$ot_method, "exact")) NA_real_ else params$sinkhorn_epsilon,
      if (identical(params$ot_method, "exact")) NA_real_ else params$sinkhorn_niter,
      params$n_cores
    )
  )
  utils::write.csv(
    qc, file.path(dataset_dir, tagged_name("qc_summary", ".csv")),
    row.names = FALSE
  )
  print(qc)

  per_method <- list()
  summary_rows <- list()
  query_bootstrap_rows <- list()
  tie_sensitivity_rows <- list()
  rounded_tie_sensitivity_rows <- list()
  tie_diagnostic_rows <- list()
  tie_diagnostic_summary_rows <- list()
  sanitize_rows <- list()
  pairwise_auc_rows <- list()
  fdr_rows <- list()
  distance_provenance_rows <- list()

  for (mi in seq_along(methods_all)) {
    method <- methods_all[[mi]]
    message("[", dataset, "] ", method)
    p <- params
    p$distance_method <- method
    p <- validate_params(p)

    destination_rds <- file.path(
      dataset_dir, tagged_name(paste0("dist_", method), ".rds")
    )
    destination_csv <- file.path(
      dataset_dir, tagged_name(paste0("dist_", method), ".csv")
    )
    source_rds <- NA_character_
    source_csv <- NA_character_
    if (nzchar(reuse_distance_dir)) {
      source_dataset_dir <- file.path(
        reuse_distance_dir, gsub("[^A-Za-z0-9]+", "_", dataset)
      )
      source_rds <- file.path(
        source_dataset_dir, tagged_name(paste0("dist_", method), ".rds")
      )
      source_csv <- file.path(
        source_dataset_dir, tagged_name(paste0("dist_", method), ".csv")
      )
      if (!file.exists(source_rds) || !file.exists(source_csv)) {
        stop("Missing source distance pair for ", dataset, "/", method)
      }
      dm <- readRDS(source_rds)
      expected_ids <- as.character(spectra$df_spec$id)
      valid_reuse <- is.matrix(dm) && is.numeric(dm) &&
        identical(dim(dm), c(length(expected_ids), length(expected_ids))) &&
        identical(rownames(dm), expected_ids) &&
        identical(colnames(dm), expected_ids) && all(is.finite(dm)) &&
        all(diag(dm) == 0)
      if (!valid_reuse) {
        stop("Reused distance matrix failed structure/ID/finiteness validation: ",
             source_rds)
      }
      copied <- c(
        file.copy(source_rds, destination_rds, overwrite = FALSE,
                  copy.mode = TRUE, copy.date = FALSE),
        file.copy(source_csv, destination_csv, overwrite = FALSE,
                  copy.mode = TRUE, copy.date = FALSE)
      )
      if (!all(copied)) stop("Failed to copy verified reused distance files.")
      source_md5 <- unname(tools::md5sum(c(source_rds, source_csv)))
      destination_md5 <- unname(tools::md5sum(c(destination_rds, destination_csv)))
      if (!identical(source_md5, destination_md5)) {
        stop("Reused distance copy is not bit-identical for ", dataset, "/", method)
      }
    } else {
      dm <- compute_distance_matrix(
        spectra$frag_list,
        spectra$loss_list,
        p,
        progress = FALSE,
        loss_typ_list = spectra$loss_typ_list,
        mref_conf = spectra$mref_confidence,
        loss_anchor_list = spectra$loss_anchor_list,
        loss_pair_list = spectra$loss_pair_list,
        loss_anchor_typ_list = spectra$loss_anchor_typ_list,
        loss_pair_typ_list = spectra$loss_pair_typ_list
      )
      source_md5 <- c(NA_character_, NA_character_)
    }
    sane <- sanitize_dist_matrix_local(dm, dataset, method)
    dm <- sane$mat
    if (!nzchar(reuse_distance_dir)) {
      saveRDS(dm, destination_rds)
      write_dist_csv(dm, destination_csv)
    }
    destination_md5 <- unname(tools::md5sum(c(destination_rds, destination_csv)))
    distance_provenance_rows[[method]] <- data.frame(
      dataset = dataset, method = method,
      ot_method = if (identical(method, "ppm_wasserstein")) p$ot_method else NA_character_,
      ot_estimand = if (identical(method, "ppm_wasserstein")) {
        if (identical(p$ot_method, "exact"))
          "unregularized_exact_transport_cost" else
          "finite_iteration_regularized_plan_transport_cost"
      } else NA_character_,
      matrix_mode = if (nzchar(reuse_distance_dir))
        "verified_reuse" else "computed_current_run",
      distance_generation_commit = distance_generation_commit,
      statistics_recompute_commit = git_commit,
      source_rds = provenance_path(source_rds),
      source_csv = provenance_path(source_csv),
      destination_rds = provenance_path(destination_rds),
      destination_csv = provenance_path(destination_csv),
      source_rds_md5 = source_md5[[1L]], source_csv_md5 = source_md5[[2L]],
      destination_rds_md5 = destination_md5[[1L]],
      destination_csv_md5 = destination_md5[[2L]],
      source_destination_bit_identical = if (nzchar(reuse_distance_dir))
        identical(source_md5, destination_md5) else NA,
      n_rows = nrow(dm), n_columns = ncol(dm), ids_exact = TRUE,
      all_finite = all(is.finite(dm)), diagonal_zero = all(diag(dm) == 0),
      stringsAsFactors = FALSE
    )

    pq <- per_query_metrics_tie_aware(
      dm,
      inchikeys,
      tie_tolerance = 0,
      random_seed = seed + mi * 1000L
    )
    # Binary optimistic Top-1 is retained only for the exploratory legacy
    # McNemar comparison. Fractional expected outcomes are primary elsewhere.
    pq$top1 <- pq$top1_optimistic
    per_method[[method]] <- pq
    utils::write.csv(
      pq, file.path(dataset_dir, tagged_name(paste0("per_query_", method), ".csv")),
      row.names = FALSE
    )
    summary_rows[[method]] <- summarize_tie_aware_metrics(
      pq, dataset, method,
      query_boot_R = n_boot,
      cluster_boot_R = cluster_boot,
      seed = seed + mi * 1000L
    )
    query_bootstrap_rows[[method]] <- query_bootstrap_ci_table(
      summary_rows[[method]]
    )
    tie_sensitivity_rows[[method]] <- summarize_tie_sensitivity(
      pq, dataset, method
    )
    tie_diagnostic_rows[[method]] <- tie_diagnostics(pq, dataset, method)
    tie_diagnostic_summary_rows[[method]] <- summarize_tie_diagnostics(
      pq, dataset, method
    )

    pq_round12 <- per_query_metrics_tie_aware(
      round(dm, 12),
      inchikeys,
      tie_tolerance = 0,
      random_seed = seed + 500000L + mi * 1000L
    )
    rounded <- summarize_tie_sensitivity(pq_round12, dataset, method)
    rounded$tie_definition <- "distance_rounded_12_digits"
    rounded_tie_sensitivity_rows[[method]] <- rounded
    pair_metrics <- ordered_pair_discrimination_metrics(dm, inchikeys)
    pairwise_auc_rows[[method]] <- cbind(
      dataset = dataset, method = method, pair_metrics$auc
    )
    fdr_rows[[method]] <- cbind(
      dataset = dataset, method = method,
      similarity_transform = "ordered_pair_minmax_distance", pair_metrics$fdr
    )
    sanitize_rows[[method]] <- data.frame(
      dataset = dataset,
      method = method,
      n_total = length(dm),
      n_nonfinite = sane$n_nonfinite,
      replacement = sane$replacement,
      stringsAsFactors = FALSE
    )
  }

  per_query_all <- do.call(rbind, lapply(names(per_method), function(m) {
    cbind(dataset = dataset, method = m, per_method[[m]])
  }))
  summary <- do.call(rbind, summary_rows)
  query_bootstrap <- do.call(rbind, query_bootstrap_rows)
  tie_sensitivity <- do.call(rbind, tie_sensitivity_rows)
  tie_sensitivity$tie_definition <- "full_precision_exact_equality"
  rounded_tie_sensitivity <- do.call(rbind, rounded_tie_sensitivity_rows)
  tie_diagnostic <- do.call(rbind, tie_diagnostic_rows)
  tie_diagnostic_summary <- do.call(rbind, tie_diagnostic_summary_rows)
  sanitize_log <- do.call(rbind, sanitize_rows)
  pairwise_auc <- do.call(rbind, pairwise_auc_rows)
  fdr <- do.call(rbind, fdr_rows)
  mcnemar <- mcnemar_dataset(dataset, per_method)

  paired_rows <- list()
  if ("ppm_wasserstein" %in% names(per_method)) {
    specs <- data.frame(
      metric = c("top1", "mrr", "p_at_1"),
      value_col = c("top1_fractional", "rr_fractional", "p_at_1_fractional"),
      cluster_col = c("inchikey", "inchikey", "inchikey_prefix"),
      cluster_unit = c("full_inchikey", "full_inchikey", "inchikey_prefix_14"),
      stringsAsFactors = FALSE
    )
    row_idx <- 0L
    for (comparator in setdiff(names(per_method), "ppm_wasserstein")) {
      for (si in seq_len(nrow(specs))) {
        row_idx <- row_idx + 1L
        diff <- paired_cluster_bootstrap_difference(
          per_method[["ppm_wasserstein"]],
          per_method[[comparator]],
          value_col = specs$value_col[[si]],
          cluster_col = specs$cluster_col[[si]],
          R = cluster_boot,
          seed = seed + 700000L + row_idx
        )
        paired_rows[[row_idx]] <- data.frame(
          dataset = dataset,
          reference = "ppm_wasserstein",
          comparator = comparator,
          metric = specs$metric[[si]],
          tie_policy = "fractional_expected",
          difference = diff[["estimate"]],
          ci_low = diff[["ci_low"]],
          ci_high = diff[["ci_high"]],
          n_queries = as.integer(diff[["n_queries"]]),
          n_clusters = as.integer(diff[["n_clusters"]]),
          analysis_level = "cluster",
          status = "primary",
          inference_method = "paired_cluster_bootstrap",
          cluster_unit = specs$cluster_unit[[si]],
          estimand = paste0(
            "query_weighted_mean_paired_difference_reference_minus_comparator_",
            "conditional_on_fixed_candidate_library; ",
            "clusters_resampled_with_replacement_by_",
            specs$cluster_unit[[si]]
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  paired_cluster <- if (length(paired_rows)) do.call(rbind, paired_rows) else data.frame()
  distance_provenance <- do.call(rbind, distance_provenance_rows)
  rownames(distance_provenance) <- NULL
  utils::write.csv(
    distance_provenance,
    file.path(dataset_dir, tagged_name("distance_matrix_provenance", ".csv")),
    row.names = FALSE, na = ""
  )

  dataset_outputs <- list(
    retrieval_per_query_metrics = per_query_all,
    retrieval_cluster_bootstrap_ci = summary,
    retrieval_query_bootstrap_ci_exploratory = query_bootstrap,
    retrieval_tie_sensitivity_exact = tie_sensitivity,
    retrieval_tie_sensitivity_round12 = rounded_tie_sensitivity,
    retrieval_tie_diagnostics = tie_diagnostic,
    retrieval_tie_affected_summary = tie_diagnostic_summary,
    retrieval_paired_cluster_bootstrap = paired_cluster,
    retrieval_mcnemar = mcnemar,
    distance_sanitize_log = sanitize_log,
    retrieval_ordered_pair_auc = pairwise_auc,
    retrieval_ordered_pair_fdr = fdr
  )
  for (output_stem in names(dataset_outputs)) {
    utils::write.csv(
      dataset_outputs[[output_stem]],
      file.path(dataset_dir, tagged_name(output_stem, ".csv")), row.names = FALSE
    )
  }

  list(
    dataset = dataset,
    spectra = spectra,
    qc = qc,
    per_query = per_query_all,
    summary = summary,
    query_bootstrap = query_bootstrap,
    tie_sensitivity = tie_sensitivity,
    rounded_tie_sensitivity = rounded_tie_sensitivity,
    tie_diagnostic = tie_diagnostic,
    tie_diagnostic_summary = tie_diagnostic_summary,
    paired_cluster = paired_cluster,
    mcnemar = mcnemar,
    sanitize = sanitize_log,
    pairwise_auc = pairwise_auc,
    fdr = fdr,
    distance_provenance = distance_provenance
  )
}

results <- lapply(datasets, run_dataset)
names(results) <- datasets
attr(results, "file_tag") <- file_tag
attr(results, "tol_ppm") <- tol_ppm
attr(results, "package_commit") <- git_commit
attr(results, "package_tree_dirty") <- package_tree_dirty
attr(results, "distance_generation_commit") <- distance_generation_commit
combined_distance_provenance <- do.call(
  rbind, lapply(results, function(value) value$distance_provenance)
)
rownames(combined_distance_provenance) <- NULL
utils::write.csv(
  combined_distance_provenance,
  file.path(output_dir, tagged_name("combined_distance_matrix_provenance", ".csv")),
  row.names = FALSE, na = ""
)
results_path <- file.path(output_dir, tagged_name("main_retrieval_results", ".rds"))
saveRDS(results, results_path)

# Finalization is a separate process so that an interrupted figure device or
# logging connection cannot invalidate already completed distance/statistics
# work. It is also independently rerunnable from the saved result object.
finalizer <- file.path(repo_dir, "inst", "scripts",
                       "finalize_main_retrieval_publication.R")
status <- system2(
  command = file.path(R.home("bin"), "Rscript"),
  args = c(
    finalizer, paste0("--output-dir=", normalizePath(output_dir)),
    paste0("--results-rds=", normalizePath(results_path)),
    paste0("--file-tag=", file_tag)
  ),
  stdout = "",
  stderr = ""
)
if (length(status) != 1L || is.na(status) || status != 0L) {
  stop("Main result finalization failed with status ", status,
       ". Distance and dataset-level outputs remain saved in ", output_dir)
}
message("All done.")
