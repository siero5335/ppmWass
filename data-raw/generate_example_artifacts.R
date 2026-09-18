load_package_for_examples <- function() {
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }

  if (requireNamespace("ppmWass", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ppmWass))
    return(invisible(TRUE))
  }

  stop("Install ppmWass or pkgload before running this script.")
}

load_package_for_examples()

dir.create("inst/extdata/benchmark_demo", recursive = TRUE, showWarnings = FALSE)
dir.create("inst/extdata/library_demo", recursive = TRUE, showWarnings = FALSE)

params <- eihrms_default_params()

demo_df <- generate_demo_msdial(
  n_known = 10,
  n_unknown = 4,
  n_samples = 6,
  seed = 1
)
quant <- prepare_quant_data(demo_df, params, write_cleaned_csv = FALSE)
spectra <- build_spectra(demo_df, quant$final, params, progress = FALSE)

n_demo <- nrow(spectra$df_spec)
replicate_group <- rep(seq_len(ceiling(n_demo / 2)), each = 2)[seq_len(n_demo)]
spectra$df_spec$inchikey <- sprintf("SIMKEY%08d-DEMO", replicate_group)

benchmark_methods <- c("cosine", "hellinger", "weighted_cosine", "composite")
benchmark_res <- run_benchmark(
  spectra_result = spectra,
  methods = benchmark_methods,
  params = params
)

readr::write_csv(
  benchmark_res$summary,
  "inst/extdata/benchmark_demo/benchmark_summary.csv"
)

sweep_res <- sweep_distance_params(
  spectra_result = spectra,
  grid = data.frame(
    tol_ppm = c(10, 10, 20, 20),
    w_frag = c(0.7, 0.5, 0.9, 0.6),
    w_loss = c(0.3, 0.5, 0.1, 0.4),
    stringsAsFactors = FALSE
  ),
  methods = c("cosine", "weighted_cosine"),
  params = params,
  n_replicates = 3,
  sample_frac = 0.8,
  progress = FALSE
)

readr::write_csv(
  sweep_res$summary,
  "inst/extdata/benchmark_demo/distance_param_sweep.csv"
)

runtime_sizes <- sort(unique(pmin(c(5, 10, length(spectra$frag_list)), length(spectra$frag_list))))
runtime_res <- benchmark_distance_runtime(
  spectra_result = spectra,
  methods = c("cosine", "weighted_cosine"),
  sizes = runtime_sizes,
  n_replicates = 2,
  params = params,
  seed = 1,
  progress = FALSE
)

readr::write_csv(
  runtime_res$summary,
  "inst/extdata/benchmark_demo/distance_runtime_summary.csv"
)

writeLines(
  c(
    "Synthetic benchmark artifacts generated from data-raw/generate_example_artifacts.R",
    "Seed: 1",
    paste("Methods:", paste(benchmark_methods, collapse = ", ")),
    "InChIKeys were assigned synthetically in replicate pairs for tabular examples."
  ),
  "inst/extdata/benchmark_demo/provenance.txt"
)

make_spec <- function(mz, intensity) {
  cbind(mz = mz, intensity = intensity)
}

library_df <- data.frame(
  id = c("LIB_A", "LIB_B", "LIB_C", "LIB_D"),
  RI = c(1000, 1120, 1260, 1390),
  inchikey = c(
    "AAAAAAAAAAAAAA-TESTAA",
    "BBBBBBBBBBBBBB-TESTBB",
    "CCCCCCCCCCCCCC-TESTCC",
    "DDDDDDDDDDDDDD-TESTDD"
  ),
  formula = c("C6H6", "C7H8", "C8H10O", "C5H10O2"),
  compound_class = c("aromatic", "aromatic", "oxygenated", "ester"),
  stringsAsFactors = FALSE
)

library_frag <- list(
  LIB_A = make_spec(c(50, 77, 91, 105), c(0.30, 1.00, 0.80, 0.35)),
  LIB_B = make_spec(c(55, 79, 93, 107), c(0.28, 1.00, 0.76, 0.32)),
  LIB_C = make_spec(c(60, 88, 103, 121), c(0.25, 1.00, 0.74, 0.40)),
  LIB_D = make_spec(c(43, 61, 74, 87), c(0.42, 1.00, 0.66, 0.36))
)

query_df <- data.frame(
  id = c("Q_A", "Q_B", "Q_C"),
  RI = c(1005, 1115, 1268),
  inchikey = c(
    "AAAAAAAAAAAAAA-QUERY1",
    "BBBBBBBBBBBBBB-QUERY2",
    "CCCCCCCCCCCCCC-QUERY3"
  ),
  formula = c("C6H6", "C7H8", "C8H10O"),
  compound_class = c("aromatic", "aromatic", "oxygenated"),
  stringsAsFactors = FALSE
)

query_frag <- list(
  Q_A = make_spec(c(50, 77, 91, 105), c(0.33, 1.00, 0.77, 0.31)),
  Q_B = make_spec(c(55, 79, 93, 107), c(0.30, 1.00, 0.73, 0.29)),
  Q_C = make_spec(c(60, 88, 103, 121), c(0.27, 1.00, 0.70, 0.38))
)

write_msp(
  df_spec = query_df,
  frag_list = query_frag,
  file = "inst/extdata/library_demo/query_demo.msp"
)
write_msp(
  df_spec = library_df,
  frag_list = library_frag,
  file = "inst/extdata/library_demo/library_demo.msp"
)

query_spectra <- build_spectra_from_msp(
  "inst/extdata/library_demo/query_demo.msp",
  params = params,
  require_ri = FALSE,
  progress = FALSE
)
library_spectra <- build_spectra_from_msp(
  "inst/extdata/library_demo/library_demo.msp",
  params = params,
  require_ri = FALSE,
  progress = FALSE
)

library_res <- benchmark_library_search(
  query_spectra = query_spectra,
  library_spectra = library_spectra,
  methods = c("cosine", "weighted_cosine", "composite"),
  params = params,
  keep_search_results = TRUE,
  search_top_n = 3,
  progress = FALSE
)

readr::write_csv(
  library_res$summary,
  "inst/extdata/library_demo/library_benchmark_summary.csv"
)

for (method in c("cosine", "weighted_cosine", "composite")) {
  readr::write_csv(
    library_res[[method]]$hits,
    file.path("inst/extdata/library_demo", paste0("library_hits_", method, ".csv"))
  )
}

writeLines(
  c(
    "Synthetic library benchmark artifacts generated from data-raw/generate_example_artifacts.R",
    "Query MSP: query_demo.msp",
    "Library MSP: library_demo.msp"
  ),
  "inst/extdata/library_demo/provenance.txt"
)
