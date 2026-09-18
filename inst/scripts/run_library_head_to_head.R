#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(
    query = NULL,
    library = NULL,
    output_dir = "library_benchmark_output",
    methods = c("cosine", "entropy", "wasserstein", "weighted_cosine", "composite", "ppm_wasserstein"),
    ri_tolerance = NULL
  )

  for (arg in args) {
    if (startsWith(arg, "--query=")) {
      opts$query <- sub("^--query=", "", arg)
    } else if (startsWith(arg, "--library=")) {
      opts$library <- sub("^--library=", "", arg)
    } else if (startsWith(arg, "--output-dir=")) {
      opts$output_dir <- sub("^--output-dir=", "", arg)
    } else if (startsWith(arg, "--methods=")) {
      opts$methods <- strsplit(sub("^--methods=", "", arg), ",", fixed = TRUE)[[1]]
    } else if (startsWith(arg, "--ri-tolerance=")) {
      opts$ri_tolerance <- as.numeric(sub("^--ri-tolerance=", "", arg))
    } else if (identical(arg, "--help")) {
      cat(
        "Usage:\n",
        "  Rscript inst/scripts/run_library_head_to_head.R --query=query.msp --library=library.msp [--output-dir=DIR] [--methods=cosine,entropy,...] [--ri-tolerance=50]\n",
        sep = ""
      )
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  if (is.null(opts$query) || is.null(opts$library)) {
    stop("--query and --library are required.")
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

query_spectra <- build_spectra_from_msp(opts$query, params = params, require_ri = !is.null(opts$ri_tolerance), progress = TRUE)
library_spectra <- build_spectra_from_msp(opts$library, params = params, require_ri = !is.null(opts$ri_tolerance), progress = TRUE)

bench <- benchmark_library_search(
  query_spectra = query_spectra,
  library_spectra = library_spectra,
  methods = opts$methods,
  params = params,
  ri_tolerance = opts$ri_tolerance,
  keep_search_results = TRUE,
  search_top_n = 20,
  progress = TRUE
)

readr::write_csv(bench$summary, file.path(opts$output_dir, "library_benchmark_summary.csv"))

for (method in opts$methods) {
  if (!is.null(bench[[method]]$hits)) {
    readr::write_csv(
      bench[[method]]$hits,
      file.path(opts$output_dir, paste0("library_hits_", method, ".csv"))
    )
  }
}

writeLines(capture.output(sessionInfo()), con = file.path(opts$output_dir, "sessionInfo.txt"))

message("Wrote library benchmark artifacts to: ", normalizePath(opts$output_dir, winslash = "/"))
