#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(
    input = NULL,
    output_dir = "benchmark_output",
    methods = c("cosine", "entropy", "wasserstein", "weighted_cosine", "composite"),
    seed = 1,
    demo = TRUE
  )

  for (arg in args) {
    if (identical(arg, "--no-demo")) {
      opts$demo <- FALSE
    } else if (startsWith(arg, "--input=")) {
      opts$input <- sub("^--input=", "", arg)
      opts$demo <- FALSE
    } else if (startsWith(arg, "--output-dir=")) {
      opts$output_dir <- sub("^--output-dir=", "", arg)
    } else if (startsWith(arg, "--methods=")) {
      opts$methods <- strsplit(sub("^--methods=", "", arg), ",", fixed = TRUE)[[1]]
    } else if (startsWith(arg, "--seed=")) {
      opts$seed <- as.integer(sub("^--seed=", "", arg))
    } else if (identical(arg, "--help")) {
      cat(
        "Usage:\n",
        "  Rscript inst/scripts/run_benchmark_reproducible.R [--input=FILE] [--output-dir=DIR] [--methods=cosine,entropy,...] [--seed=1]\n",
        "  If --input is omitted, a synthetic demo dataset is used.\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  opts
}

opts <- parse_args(args)

if (!requireNamespace("ppmWass", quietly = TRUE)) {
  stop("Package 'ppmWass' is not installed. Install it first with remotes::install_local('.').")
}

suppressPackageStartupMessages(library(ppmWass))

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

params <- eihrms_default_params()

build_input_spectra <- function() {
  if (isTRUE(opts$demo) || is.null(opts$input)) {
    demo_df <- generate_demo_msdial(
      n_known = 12,
      n_unknown = 8,
      n_samples = 8,
      seed = opts$seed
    )
    quant <- prepare_quant_data(demo_df, params, write_cleaned_csv = FALSE)
    return(build_spectra(demo_df, quant$final, params, progress = FALSE))
  }

  df <- read_ms_dial(opts$input)
  quant <- prepare_quant_data(df, params, write_cleaned_csv = FALSE)
  build_spectra(df, quant$final, params, progress = FALSE)
}

make_sweep_grid <- function() {
  weight_profiles <- data.frame(
    w_frag = c(0.5, 0.7, 0.9),
    w_loss = c(0.5, 0.3, 0.1),
    stringsAsFactors = FALSE
  )
  tol_beta <- expand.grid(
    tol_ppm = c(10, 20, 30),
    loss_raw_beta = c(0.5, 0.8),
    stringsAsFactors = FALSE
  )

  do.call(rbind, lapply(seq_len(nrow(weight_profiles)), function(i) {
    cbind(weight_profiles[rep(i, nrow(tol_beta)), , drop = FALSE], tol_beta)
  }))
}

choose_runtime_sizes <- function(n_total) {
  sizes <- unique(pmin(c(10, 25, 50, 100, 250, 500, 1000), n_total))
  sizes <- sizes[sizes >= 2]
  sort(unique(as.integer(sizes)))
}

spectra <- build_input_spectra()

benchmark_res <- run_benchmark(
  spectra_result = spectra,
  methods = opts$methods,
  params = params
)
readr::write_csv(benchmark_res$summary, file.path(opts$output_dir, "benchmark_summary.csv"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_cmp <- plot_benchmark_comparison(benchmark_res)
  ggplot2::ggsave(file.path(opts$output_dir, "benchmark_comparison.png"), p_cmp, width = 10, height = 8)

  p_fdr <- plot_fdr_curves(benchmark_res)
  if (!is.null(p_fdr)) {
    ggplot2::ggsave(file.path(opts$output_dir, "benchmark_fdr.png"), p_fdr, width = 8, height = 5)
  }
}

sweep_res <- sweep_distance_params(
  spectra_result = spectra,
  grid = make_sweep_grid(),
  methods = opts$methods,
  params = params,
  metrics = c("replicate_top1", "replicate_mrr", "auc_roc", "map", "fdr_at_07", "silhouette", "homolog_score"),
  progress = TRUE
)
readr::write_csv(sweep_res$summary, file.path(opts$output_dir, "distance_param_sweep.csv"))

runtime_sizes <- choose_runtime_sizes(length(spectra$frag_list))
runtime_res <- benchmark_distance_runtime(
  spectra_result = spectra,
  methods = opts$methods,
  sizes = runtime_sizes,
  n_replicates = 3,
  params = params,
  seed = opts$seed,
  progress = TRUE
)
readr::write_csv(runtime_res$results, file.path(opts$output_dir, "distance_runtime_results.csv"))
readr::write_csv(runtime_res$summary, file.path(opts$output_dir, "distance_runtime_summary.csv"))

writeLines(capture.output(sessionInfo()), con = file.path(opts$output_dir, "sessionInfo.txt"))

message("Wrote benchmark artifacts to: ", normalizePath(opts$output_dir, winslash = "/"))
