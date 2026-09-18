toy_lc_spectra <- function() {
  frag <- list(
    q1 = matrix(c(100, 2, 150, 1), ncol = 2, byrow = TRUE),
    q2 = matrix(c(120, 2, 170, 1), ncol = 2, byrow = TRUE),
    q3 = matrix(c(140, 2, 190, 1), ncol = 2, byrow = TRUE)
  )
  metadata <- data.frame(
    id = names(frag),
    precursor_mz = c(300, 300.002, 450),
    rt = c(5.0, 5.1, 9.0),
    ion_mode = c("POSITIVE", "pos", "negative"),
    precursor_type = c("[M+H]+", "[M+H]+", "[M-H]-"),
    collision_energy = c("20 eV", "22", "20"),
    stringsAsFactors = FALSE
  )
  as_lc_spectra(frag, metadata)
}

test_that("LC defaults are fragment-only and preserve core parameters", {
  params <- lc_default_params()
  expect_equal(params$w_frag, 1)
  expect_equal(params$w_loss, 0)
  expect_equal(params$max_mz, 2000)
  expect_false(params$use_mref_confidence)
})

test_that("LC spectra canonicalize metadata and normalize TIC", {
  x <- toy_lc_spectra()
  expect_s3_class(x, "LCSpectra")
  expect_equal(x$metadata$ion_mode, c("positive", "positive", "negative"))
  expect_equal(x$metadata$adduct, c("[M+H]+", "[M+H]+", "[M-H]-"))
  expect_equal(x$metadata$collision_energy, c(20, 22, 20))
  expect_equal(unname(vapply(x$frag_list, function(s) sum(s[, 2]), numeric(1))), rep(1, 3))
})

test_that("LC candidate gates retain only compatible entries", {
  x <- toy_lc_spectra()
  query <- as_lc_spectra(x$frag_list[1], x$metadata[1, , drop = FALSE])
  plan <- lc_search_plan(
    precursor_ppm = 10,
    rt_window = 0.25,
    require_ion_mode = TRUE,
    require_adduct = TRUE,
    collision_energy_tolerance = 3,
    first_pass_method = "hellinger",
    rerank_method = "hellinger"
  )
  candidates <- lc_candidate_sets(query, x, plan)
  expect_equal(candidates$indices$q1, c(1L, 2L))
  expect_equal(candidates$audit$initial, 3)
  expect_equal(candidates$audit$final, 2)
})

test_that("missing metadata policy is explicit", {
  x <- toy_lc_spectra()
  query_meta <- x$metadata[1, , drop = FALSE]
  query_meta$precursor_mz <- NA_real_
  query <- as_lc_spectra(x$frag_list[1], query_meta)

  allow <- lc_search_plan(
    precursor_ppm = 10,
    require_ion_mode = FALSE,
    missing_metadata = "allow",
    first_pass_method = "hellinger",
    rerank_method = "hellinger"
  )
  exclude <- lc_search_plan(
    precursor_ppm = 10,
    require_ion_mode = FALSE,
    missing_metadata = "exclude",
    first_pass_method = "hellinger",
    rerank_method = "hellinger"
  )
  error_plan <- lc_search_plan(
    precursor_ppm = 10,
    require_ion_mode = FALSE,
    missing_metadata = "error",
    first_pass_method = "hellinger",
    rerank_method = "hellinger"
  )

  expect_equal(length(lc_candidate_sets(query, x, allow)$indices[[1]]), 3L)
  expect_equal(length(lc_candidate_sets(query, x, exclude)$indices[[1]]), 0L)
  expect_error(lc_candidate_sets(query, x, error_plan), "precursor_mz")
})

test_that("LC search applies gates before spectral ranking", {
  library <- toy_lc_spectra()
  query <- as_lc_spectra(library$frag_list[1], library$metadata[1, , drop = FALSE])
  plan <- lc_search_plan(
    precursor_ppm = 10,
    rt_window = 0.25,
    require_ion_mode = TRUE,
    require_adduct = TRUE,
    first_pass_method = "hellinger",
    rerank_method = "hellinger",
    first_pass_top_k = 2,
    final_top_k = 1
  )
  result <- search_lc_library(query, library, plan, progress = FALSE)
  expect_equal(result$results$lib_id, "q1")
  expect_equal(result$candidate_summary$final, 2)
  expect_equal(result$candidate_summary$n_returned, 1)
})

test_that("LC MSP import preserves precursor and acquisition metadata", {
  path <- tempfile(fileext = ".msp")
  writeLines(c(
    "Name: compound",
    "PrecursorMZ: 321.1234",
    "PrecursorType: [M+H]+",
    "IonMode: Positive",
    "RetentionTime: 4.2",
    "CollisionEnergy: 25 eV",
    "Num Peaks: 2",
    "100.1 50",
    "150.2 100",
    ""
  ), path)
  x <- suppressMessages(read_lc_msp(path, progress = 0))
  expect_s3_class(x, "LCSpectra")
  expect_equal(x$metadata$precursor_mz, 321.1234)
  expect_equal(x$metadata$adduct, "[M+H]+")
  expect_equal(x$metadata$ion_mode, "positive")
  expect_equal(x$metadata$rt, 4.2)
  expect_equal(x$metadata$collision_energy, 25)
})
