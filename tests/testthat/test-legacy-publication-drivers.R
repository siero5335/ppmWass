legacy_publication_script_path <- function(file) {
  path <- testthat::test_path("..", "..", "inst", "scripts", file)
  if (!file.exists(path)) {
    path <- system.file("scripts", file, package = "ppmWass")
  }
  path
}

test_that("superseded publication drivers fail before reading inputs", {
  legacy <- c(
    "run_main_retrieval_tol15.R",
    "run_hrei_msdb_benchmark.R",
    "run_additional_experiments.R",
    "run_tier26_auroc_7methods.R",
    "run_tol15_posthoc_additional.R",
    "run_retrieval_bootstrap_stats.R",
    "run_controlled_bootstrap_stats.R"
  )

  for (file in legacy) {
    path <- legacy_publication_script_path(file)
    expect_true(file.exists(path), info = file)
    expect_error(
      sys.source(path, envir = new.env(parent = baseenv())),
      "PPMWASS_LEGACY_DRIVER_DISABLED",
      info = file
    )
  }
})

test_that("superseded full-run wrapper is a fail-closed compatibility stub", {
  path <- legacy_publication_script_path("run_additional_experiments_full.sh")
  expect_true(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "PPMWASS_LEGACY_DRIVER_DISABLED", fixed = TRUE)
  expect_match(text, "exit 2", fixed = TRUE)
  expect_false(grepl("Rscript", text, fixed = TRUE))
})

test_that("distributed benchmark scripts contain no maximum-finite imputation", {
  scripts_dir <- dirname(legacy_publication_script_path(
    "run_main_retrieval_publication.R"
  ))
  scripts <- list.files(
    scripts_dir,
    pattern = "[.](R|r)$",
    full.names = TRUE,
    recursive = TRUE
  )
  text <- paste(
    vapply(scripts, function(path) {
      paste(readLines(path, warn = FALSE), collapse = "\n")
    }, character(1L)),
    collapse = "\n"
  )

  prohibited <- c(
    "max[[:space:]]*\\([[:space:]]*finite_vals",
    "dm[[:space:]]*\\[[[:space:]]*bad[[:space:]]*\\][[:space:]]*<-[[:space:]]*(replacement|replace_val)",
    "assignInNamespace[[:space:]]*\\([[:space:]]*[\"']compute_distance_matrix"
  )
  for (pattern in prohibited) {
    expect_false(grepl(pattern, text, perl = TRUE), info = pattern)
  }
})
