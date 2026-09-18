test_that("MSP search attaches the retained record's metadata in ranked order", {
  path <- tempfile(fileext = ".msp")
  on.exit(unlink(path), add = TRUE)
  writeLines(c(
    "NAME: A", "CAS: excluded", "Num Peaks: 1", "100 100", "",
    "NAME: B", "RI: 1000", "CAS: cas-b",
    "INCHIKEY: BBBBBBBBBBBBBB-BBBBBBBBBB-B", "Num Peaks: 1", "120 100", "",
    "NAME: A", "RI: 1100", "CAS: cas-a",
    "INCHIKEY: AAAAAAAAAAAAAA-AAAAAAAAAA-A", "Num Peaks: 1", "100 100", "",
    "NAME: A", "RI: 1200", "CAS: duplicate", "Num Peaks: 1", "150 100", "",
    "NAME: empty", "RI: 1300", "Num Peaks: 0", ""
  ), path)
  params <- eihrms_default_params()
  params$use_parallel <- FALSE
  params$distance_method <- "cosine"
  params$w_frag <- 1
  params$w_loss <- 0
  spectra <- build_spectra_from_msp(path, params, require_ri = TRUE, progress = FALSE)
  expect_equal(spectra$df_spec$id, c("B", "A"))
  expect_equal(nrow(spectra$msp_metadata), 5L)

  hits <- search_library(spectra, spectra, params, top_n = 2)
  expect_equal(nrow(hits), 4L)
  expect_equal(hits$cas, c("cas-b", "cas-a", "cas-a", "cas-b"))
  expect_equal(hits$ri, c(1000, 1100, 1100, 1000))
  expect_equal(hits$name, hits$library_id)
  expect_equal(hits$inchikey, ifelse(hits$library_id == "A",
    "AAAAAAAAAAAAAA-AAAAAAAAAA-A", "BBBBBBBBBBBBBB-BBBBBBBBBB-B"))
  expect_equal(hits$rank, rep(1:2, 2))

  # Older serialized objects do not contain df_spec$name.
  spectra$df_spec$name <- NULL
  listed <- search_library(spectra, spectra, params, top_n = 2, as_list = TRUE)
  expect_equal(listed$A$library_id, c("A", "B"))
  expect_equal(listed$A$name, c("A", "B"))
  expect_equal(listed$A$cas, c("cas-a", "cas-b"))
})

test_that("MSP conversion preserves the original name with a custom ID column", {
  raw <- tibble::tibble(name = "Original name", record_id = "ID1", spectrum = "100:1")
  converted <- convert_msp_to_internal(raw, id_column = "record_id")
  expect_equal(converted$id, "ID1")
  expect_equal(converted$name, "Original name")
})

test_that("one-entry libraries retain IDs and RI values for every query", {
  sp <- matrix(c(100, 1), ncol = 2)
  empty <- matrix(numeric(0), ncol = 2)
  query <- list(frag_list = list(q1 = sp, q2 = sp),
                loss_list = list(q1 = empty, q2 = empty), ri = c(q1 = 1000, q2 = 1020))
  library <- list(frag_list = list(lib = sp), loss_list = list(lib = empty), ri = c(lib = 1010))
  params <- eihrms_default_params()
  params$use_parallel <- FALSE
  params$distance_method <- "cosine"
  params$w_frag <- 1
  params$w_loss <- 0
  hits <- search_library(query, library, params)
  expect_equal(hits$query_id, c("q1", "q2"))
  expect_equal(hits$library_id, c("lib", "lib"))
  expect_equal(hits$rank, c(1L, 1L))
  expect_equal(hits$similarity, c(1, 1))
  expect_equal(hits$delta_ri, c(10, 10))
})

test_that("library search does not turn infinite distances into finite hits", {
  local_mocked_bindings(compute_distance_matrix_search = function(...) {
    matrix(c(Inf, 0.2, NA_real_, -Inf), nrow = 1,
           dimnames = list("q", c("infinite", "valid", "missing", "negative_inf")))
  })
  sp <- matrix(c(100, 1), ncol = 2)
  query <- list(frag_list = list(q = sp), loss_list = list(q = sp))
  library <- list(frag_list = list(infinite = sp, valid = sp, missing = sp, negative_inf = sp))
  hits <- search_library(query, library, eihrms_default_params())
  expect_equal(hits$library_id, "valid")
  expect_equal(hits$similarity, 0.8)
})
