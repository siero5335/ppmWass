#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", full_args, value = TRUE)
runtime_script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
} else {
  NA_character_
}

parse_csv <- function(x, type = c("character", "integer")) {
  type <- match.arg(type)
  vals <- strsplit(x, ",", fixed = TRUE)[[1]]
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  if (identical(type, "integer")) {
    vals <- as.integer(vals)
    if (any(!is.finite(vals))) {
      stop("Expected comma-separated integers, got: ", x)
    }
  }
  vals
}

parse_args <- function(args) {
  opts <- list(
    repo = getwd(),
    input = "RECETOX_merged.msp",
    output_dir = "benchmark_output_runtime_micro",
    methods = c("cosine", "weighted_cosine", "composite",
                "entropy_weighted", "entropy_unweighted",
                "hellinger", "ppm_wasserstein"),
    sizes = c(20L, 40L, 80L),
    n_replicates = 3L,
    parallel_mode = "both",
    n_cores = NULL,
    seed = 20260429L,
    min_mz = 35,
    max_mz = 650,
    tol_ppm = 15,
    ot_method = "exact",
    sinkhorn_niter = 100L,
    sinkhorn_epsilon = 0.05,
    transition_mult = 3
  )

  for (arg in args) {
    if (startsWith(arg, "--repo=")) {
      opts$repo <- sub("^--repo=", "", arg)
    } else if (startsWith(arg, "--input=")) {
      opts$input <- sub("^--input=", "", arg)
    } else if (startsWith(arg, "--output-dir=")) {
      opts$output_dir <- sub("^--output-dir=", "", arg)
    } else if (startsWith(arg, "--methods=")) {
      opts$methods <- parse_csv(sub("^--methods=", "", arg), "character")
    } else if (startsWith(arg, "--sizes=")) {
      opts$sizes <- parse_csv(sub("^--sizes=", "", arg), "integer")
    } else if (startsWith(arg, "--replicates=")) {
      opts$n_replicates <- as.integer(sub("^--replicates=", "", arg))
    } else if (startsWith(arg, "--parallel=")) {
      opts$parallel_mode <- sub("^--parallel=", "", arg)
    } else if (startsWith(arg, "--n-cores=")) {
      opts$n_cores <- as.integer(sub("^--n-cores=", "", arg))
    } else if (startsWith(arg, "--seed=")) {
      opts$seed <- as.integer(sub("^--seed=", "", arg))
    } else if (startsWith(arg, "--tol-ppm=")) {
      opts$tol_ppm <- as.numeric(sub("^--tol-ppm=", "", arg))
    } else if (startsWith(arg, "--ot-method=")) {
      opts$ot_method <- tolower(sub("^--ot-method=", "", arg))
    } else if (startsWith(arg, "--sinkhorn-niter=")) {
      opts$sinkhorn_niter <- as.integer(sub("^--sinkhorn-niter=", "", arg))
    } else if (startsWith(arg, "--sinkhorn-epsilon=")) {
      opts$sinkhorn_epsilon <- as.numeric(sub("^--sinkhorn-epsilon=", "", arg))
    } else if (startsWith(arg, "--transition-mult=")) {
      opts$transition_mult <- as.numeric(sub("^--transition-mult=", "", arg))
    } else if (identical(arg, "--help")) {
      cat(
        "Usage:\n",
        "  Rscript inst/scripts/run_runtime_microbenchmark.R [options]\n\n",
        "Options:\n",
        "  --repo=DIR              Development-tree package root.\n",
        "  --input=FILE             MSP file. Default: RECETOX_merged.msp\n",
        "  --output-dir=DIR         Output directory. Default: benchmark_output_runtime_micro\n",
        "  --methods=a,b,c          Methods to benchmark.\n",
        "  --sizes=20,40,80         Subset sizes. Keep small for quick runtime checks.\n",
        "  --replicates=3           Random subsets per size.\n",
        "  --parallel=both          serial, parallel, or both.\n",
        "  --n-cores=4              Parallel cores. Default: detectCores()-1, fallback 4.\n",
        "  --tol-ppm=15             Hard-match/base-cost ppm value.\n",
        "  --ot-method=exact        exact, sinkhorn, or greenkhorn.\n",
        "  --sinkhorn-niter=100     Sinkhorn iterations for ppm_wasserstein.\n",
        "  --sinkhorn-epsilon=0.05  Sinkhorn regularization.\n",
        "  --transition-mult=3      ppmWass saturation multiplier.\n",
        "  --seed=20260429          Random seed.\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  if (!opts$parallel_mode %in% c("serial", "parallel", "both")) {
    stop("--parallel must be serial, parallel, or both.")
  }
  if (!is.finite(opts$n_replicates) || opts$n_replicates < 1) {
    stop("--replicates must be >= 1.")
  }
  if (any(!is.finite(opts$sizes)) || any(opts$sizes < 2)) {
    stop("--sizes must contain integers >= 2.")
  }
  opts$sizes <- sort(unique(as.integer(opts$sizes)))
  opts$n_replicates <- as.integer(opts$n_replicates)
  if (!is.finite(opts$sinkhorn_epsilon) || opts$sinkhorn_epsilon <= 0 ||
      !is.finite(opts$transition_mult) || opts$transition_mult <= 0) {
    stop("--sinkhorn-epsilon and --transition-mult must be positive.")
  }
  if (!opts$ot_method %in% c("exact", "sinkhorn", "greenkhorn")) {
    stop("--ot-method must be exact, sinkhorn, or greenkhorn.")
  }
  opts
}

resolve_runtime_cores <- function(default = 4L, reserve = 1L) {
  n <- parallel::detectCores()
  if (!is.finite(n) || is.na(n)) {
    return(as.integer(default))
  }
  as.integer(max(1L, n - reserve))
}

load_ppmwass <- function(repo) {
  repo <- normalizePath(repo, mustWork = TRUE)
  if (!file.exists(file.path(repo, "DESCRIPTION"))) {
    stop("--repo is not an R package root: ", repo)
  }
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to load the audited development tree.")
  }
  pkgload::load_all(repo, quiet = TRUE)
  invisible(repo)
}

subset_by_position <- function(spectra, pos) {
  out <- spectra
  list_names <- c("frag_list", "loss_list", "loss_typ_list", "mref_confidence",
                  "loss_anchor_list", "loss_pair_list",
                  "loss_anchor_typ_list", "loss_pair_typ_list")
  for (nm in list_names) {
    if (!is.null(out[[nm]])) {
      out[[nm]] <- out[[nm]][pos]
    }
  }
  if (!is.null(out$df_spec)) {
    out$df_spec <- out$df_spec[pos, , drop = FALSE]
  }
  out
}

filter_runtime_spectra <- function(spectra) {
  has_peaks <- function(x) {
    if (is.null(x)) {
      return(FALSE)
    }
    if (is.data.frame(x)) {
      return(nrow(x) > 0)
    }
    if (is.matrix(x)) {
      return(nrow(x) > 0)
    }
    length(x) > 0
  }

  keep <- vapply(spectra$frag_list, has_peaks, logical(1))
  if (!is.null(spectra$loss_list)) {
    keep <- keep & vapply(spectra$loss_list, has_peaks, logical(1))
  }

  if (!any(keep)) {
    stop("No spectra with valid fragment/loss peaks were found.")
  }
  subset_by_position(spectra, which(keep))
}

benchmark_one <- function(spectra, method, size, rep_idx, mode, params, seed) {
  # Keep the sampled spectra identical across methods and execution modes so
  # runtime comparisons refer to exactly the same query/library workload.
  subset_seed <- as.integer(seed + size * 1000L + rep_idx * 17L)
  set.seed(subset_seed)
  n_total <- length(spectra$frag_list)
  pos <- if (size == n_total) seq_len(n_total) else sample(seq_len(n_total), size)
  spectra_sub <- subset_by_position(spectra, pos)

  p <- params
  p$methods_key <- NULL
  p$distance_method <- method
  p$use_parallel <- identical(mode, "parallel")

  gc(verbose = FALSE)
  elapsed <- system.time({
    D <- compute_distance_matrix_search(
      query_frag_list = spectra_sub$frag_list,
      query_loss_list = spectra_sub$loss_list,
      lib_frag_list = spectra_sub$frag_list,
      lib_loss_list = spectra_sub$loss_list,
      params = p,
      progress = FALSE
    )
  })[["elapsed"]]

  data.frame(
    method = method,
    size = size,
    n_query = size,
    n_library = size,
    n_cells = size * size,
    replicate = rep_idx,
    subset_seed = subset_seed,
    mode = mode,
    n_cores = if (identical(mode, "parallel")) p$n_cores else 1L,
    elapsed_sec = unname(elapsed),
    sec_per_1k_cells = unname(elapsed) / (size * size) * 1000,
    all_finite = all(is.finite(D)),
    stringsAsFactors = FALSE
  )
}

summarize_runtime <- function(results) {
  split_key <- interaction(results$method, results$size, results$mode, drop = TRUE)
  rows <- lapply(split(results, split_key), function(df) {
    data.frame(
      method = df$method[1],
      size = df$size[1],
      n_cells = df$n_cells[1],
      mode = df$mode[1],
      n_cores = df$n_cores[1],
      mean_sec = mean(df$elapsed_sec),
      median_sec = stats::median(df$elapsed_sec),
      sd_sec = stats::sd(df$elapsed_sec),
      min_sec = min(df$elapsed_sec),
      max_sec = max(df$elapsed_sec),
      mean_sec_per_1k_cells = mean(df$sec_per_1k_cells),
      all_finite = all(df$all_finite),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  fastest <- stats::aggregate(mean_sec ~ size + mode, data = out, min)
  names(fastest)[names(fastest) == "mean_sec"] <- "fastest_mean_sec"
  out <- merge(out, fastest, by = c("size", "mode"), all.x = TRUE)
  out$relative_to_fastest <- out$mean_sec / out$fastest_mean_sec
  out[order(out$size, out$mode, out$mean_sec), ]
}

summarize_parallel_speedup <- function(summary_df) {
  serial <- summary_df[summary_df$mode == "serial",
                       c("method", "size", "mean_sec")]
  parallel <- summary_df[summary_df$mode == "parallel",
                         c("method", "size", "mean_sec", "n_cores")]
  if (nrow(serial) == 0 || nrow(parallel) == 0) {
    return(data.frame())
  }
  names(serial)[3] <- "serial_mean_sec"
  names(parallel)[3] <- "parallel_mean_sec"
  out <- merge(serial, parallel, by = c("method", "size"))
  out$speedup <- out$serial_mean_sec / out$parallel_mean_sec
  out[order(out$size, -out$speedup), ]
}

write_markdown_report <- function(path, opts, spectra, summary_df, speedup_df) {
  fmt <- function(x, digits = 3) {
    ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))
  }

  lines <- c(
    "# ppmWass runtime microbenchmark",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Settings",
    "",
    paste0("- input: `", normalizePath(opts$input, winslash = "/", mustWork = FALSE), "`"),
    paste0("- audited repo: `", normalizePath(opts$repo, winslash = "/", mustWork = TRUE), "`"),
    paste0("- spectra loaded: ", length(spectra$frag_list)),
    paste0("- methods: ", paste(opts$methods, collapse = ", ")),
    paste0("- sizes: ", paste(opts$sizes, collapse = ", ")),
    paste0("- replicates: ", opts$n_replicates),
    paste0("- parallel mode: ", opts$parallel_mode),
    paste0("- n_cores: ", opts$n_cores),
    paste0("- tol_ppm: ", opts$tol_ppm),
    paste0("- sinkhorn_niter: ", opts$sinkhorn_niter),
    paste0("- sinkhorn_epsilon: ", opts$sinkhorn_epsilon),
    paste0("- transition_multiplier: ", opts$transition_mult),
    paste0("- saturation_width_ppm: ", opts$tol_ppm * opts$transition_mult),
    "",
    "## Mean elapsed time",
    "",
    "| size | mode | method | mean sec | sec / 1k cells | relative to fastest | finite |",
    "|---:|---|---|---:|---:|---:|---|"
  )

  for (i in seq_len(nrow(summary_df))) {
    r <- summary_df[i, ]
    lines <- c(lines, paste0(
      "| ", r$size,
      " | ", r$mode,
      " | ", r$method,
      " | ", fmt(r$mean_sec),
      " | ", fmt(r$mean_sec_per_1k_cells, 4),
      " | ", fmt(r$relative_to_fastest, 2), "x",
      " | ", if (isTRUE(r$all_finite)) "yes" else "no",
      " |"
    ))
  }

  if (nrow(speedup_df) > 0) {
    lines <- c(
      lines,
      "",
      "## Parallel speedup",
      "",
      "| size | method | serial sec | parallel sec | n cores | speedup |",
      "|---:|---|---:|---:|---:|---:|"
    )
    for (i in seq_len(nrow(speedup_df))) {
      r <- speedup_df[i, ]
      lines <- c(lines, paste0(
        "| ", r$size,
        " | ", r$method,
        " | ", fmt(r$serial_mean_sec),
        " | ", fmt(r$parallel_mean_sec),
        " | ", r$n_cores,
        " | ", fmt(r$speedup, 2), "x",
        " |"
      ))
    }
  }

  lines <- c(
    lines,
    "",
    "## Interpretation note",
    "",
    "This is a small runtime-only benchmark. It is intended to document the computational cost of each distance method under the same query-vs-library workload, not to re-estimate retrieval accuracy.",
    "",
    "Because `ppm_wasserstein` solves an optimal-transport problem for each pair, it is expected to be substantially slower than vector-style distances such as cosine, weighted cosine, Hellinger, or entropy-based methods. The relevant manuscript framing is therefore high-resolution reranking / accuracy-oriented retrieval rather than fastest all-vs-all screening."
  )

  writeLines(lines, con = path)
}

opts <- parse_args(args)
opts$repo <- load_ppmwass(opts$repo)
if (identical(opts$ot_method, "exact") &&
    !requireNamespace("transport", quietly = TRUE)) {
  stop("--ot-method=exact requires transport; publication runs fail closed.")
}

if (!file.exists(opts$input)) {
  stop("Input MSP not found: ", opts$input)
}

if (is.null(opts$n_cores) || !is.finite(opts$n_cores) || opts$n_cores < 1) {
  opts$n_cores <- resolve_runtime_cores(default = 4L)
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

params <- eihrms_default_params()
params$min_mz <- opts$min_mz
params$max_mz <- opts$max_mz
params$tol_ppm <- opts$tol_ppm
params$ot_method <- opts$ot_method
params$sinkhorn_niter <- opts$sinkhorn_niter
params$sinkhorn_epsilon <- opts$sinkhorn_epsilon
params$wasserstein_transition_mult <- opts$transition_mult
params$n_cores <- opts$n_cores
params$methods_key <- opts$methods

message("Loading spectra from: ", opts$input)
spectra <- build_spectra_from_msp(opts$input, params, require_ri = FALSE, progress = TRUE)
spectra <- filter_runtime_spectra(spectra)

n_total <- length(spectra$frag_list)
opts$sizes <- opts$sizes[opts$sizes <= n_total]
if (length(opts$sizes) == 0) {
  stop("No requested sizes are <= number of spectra loaded: ", n_total)
}

modes <- switch(
  opts$parallel_mode,
  serial = "serial",
  parallel = "parallel",
  both = c("serial", "parallel")
)

message("Spectra loaded: ", n_total)
message("Methods: ", paste(opts$methods, collapse = ", "))
message("Sizes: ", paste(opts$sizes, collapse = ", "))
message("Modes: ", paste(modes, collapse = ", "))
message("Replicates: ", opts$n_replicates)

rows <- list()
k <- 0L
for (size in opts$sizes) {
  for (mode in modes) {
    for (method in opts$methods) {
      for (rep_idx in seq_len(opts$n_replicates)) {
        message(sprintf(
          "[runtime] size=%d mode=%s method=%s rep=%d/%d",
          size, mode, method, rep_idx, opts$n_replicates
        ))
        k <- k + 1L
        rows[[k]] <- benchmark_one(
          spectra = spectra,
          method = method,
          size = size,
          rep_idx = rep_idx,
          mode = mode,
          params = params,
          seed = opts$seed
        )
      }
    }
  }
}

results <- do.call(rbind, rows)
summary_df <- summarize_runtime(results)
speedup_df <- summarize_parallel_speedup(summary_df)

write.csv(results, file.path(opts$output_dir, "runtime_microbenchmark_results.csv"),
          row.names = FALSE)
write.csv(summary_df, file.path(opts$output_dir, "runtime_microbenchmark_summary.csv"),
          row.names = FALSE)
write.csv(speedup_df, file.path(opts$output_dir, "runtime_microbenchmark_parallel_speedup.csv"),
          row.names = FALSE)

write_markdown_report(
  path = file.path(opts$output_dir, "runtime_microbenchmark_summary.md"),
  opts = opts,
  spectra = spectra,
  summary_df = summary_df,
  speedup_df = speedup_df
)

writeLines(capture.output(sessionInfo()),
           con = file.path(opts$output_dir, "sessionInfo.txt"))

git_commit <- tryCatch(
  suppressWarnings(trimws(system2(
    "git", c("-C", opts$repo, "rev-parse", "HEAD"),
    stdout = TRUE, stderr = FALSE
  ))),
  error = function(e) NA_character_
)
if (length(git_commit) != 1L || is.na(git_commit) || !nzchar(git_commit)) {
  git_commit <- NA_character_
}
git_status <- tryCatch(
  system2(
    "git", c("-C", opts$repo, "status", "--porcelain"),
    stdout = TRUE, stderr = FALSE
  ),
  error = function(e) character()
)
package_tree_dirty <- length(git_status) > 0L
exact_command <- paste(
  c(file.path(R.home("bin"), "Rscript"), commandArgs()), collapse = " "
)
run_parameters <- data.frame(
  parameter = c(
    "timestamp_utc", "package_commit", "package_tree_dirty", "command",
    "repo", "input", "methods",
    "sizes", "replicates", "parallel_mode", "n_cores", "seed",
    "hard_match_tolerance_ppm", "ppmWass_base_cost_ppm",
    "transition_multiplier", "saturation_width_ppm",
    "ot_method", "ot_estimand", "sinkhorn_epsilon", "sinkhorn_iterations",
    "ot_marginal_residual_tolerance", "solver_backend"
  ),
  value = c(
    format(Sys.time(), tz = "UTC"), git_commit, package_tree_dirty,
    exact_command, opts$repo,
    normalizePath(opts$input), paste(opts$methods, collapse = ","),
    paste(opts$sizes, collapse = ","), opts$n_replicates,
    opts$parallel_mode, opts$n_cores, opts$seed, opts$tol_ppm,
    opts$tol_ppm, opts$transition_mult, opts$tol_ppm * opts$transition_mult,
    opts$ot_method,
    if (identical(opts$ot_method, "exact"))
      "unregularized_exact_transport_cost" else
      "finite_iteration_regularized_plan_transport_cost",
    if (identical(opts$ot_method, "exact")) NA else opts$sinkhorn_epsilon,
    if (identical(opts$ot_method, "exact")) NA else opts$sinkhorn_niter,
    1e-8,
    if (identical(opts$ot_method, "exact"))
      "transport exact unregularized OT" else
      paste0("approxOT finite-iteration ", opts$ot_method,
             " with exact fallback on invalid plans")
  ),
  stringsAsFactors = FALSE
)
write.csv(
  run_parameters, file.path(opts$output_dir, "runtime_run_parameters.csv"),
  row.names = FALSE
)
flatten_parameter <- function(value) {
  if (is.null(value)) return("<NULL>")
  if (is.function(value)) return("<function>")
  if (is.atomic(value)) return(paste(as.character(value), collapse = ";"))
  paste(capture.output(str(value, give.attr = FALSE)), collapse = " ")
}
effective_params <- params
effective_params$methods_key <- NULL
write.csv(
  data.frame(
    parameter = names(effective_params),
    value = vapply(effective_params, flatten_parameter, character(1)),
    stringsAsFactors = FALSE
  ),
  file.path(opts$output_dir, "runtime_effective_parameters.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    file = c(normalizePath(opts$input), runtime_script_path),
    md5 = unname(tools::md5sum(c(normalizePath(opts$input), runtime_script_path))),
    stringsAsFactors = FALSE
  ),
  file.path(opts$output_dir, "runtime_input_checksums.csv"), row.names = FALSE
)

message("Wrote runtime benchmark artifacts to: ",
        normalizePath(opts$output_dir, winslash = "/"))
