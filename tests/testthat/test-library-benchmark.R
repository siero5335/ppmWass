make_library_benchmark_inputs <- function() {
  empty_loss <- matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity")))

  query <- list(
    df_spec = data.frame(
      id = c("q1", "q2"),
      inchikey = c("AAAAAAAAAAAAAA-111", "BBBBBBBBBBBBBB-222"),
      stringsAsFactors = FALSE
    ),
    frag_list = list(
      q1 = matrix(c(50, 0.6, 100, 0.4), ncol = 2, byrow = TRUE),
      q2 = matrix(c(70, 0.7, 150, 0.3), ncol = 2, byrow = TRUE)
    ),
    loss_list = list(q1 = empty_loss, q2 = empty_loss),
    loss_typ_list = NULL,
    loss_anchor_list = NULL,
    loss_pair_list = NULL,
    loss_anchor_typ_list = NULL,
    loss_pair_typ_list = NULL,
    mref_confidence = c(q1 = 1, q2 = 1),
    ri = c(q1 = 100, q2 = 220)
  )

  library <- list(
    df_spec = data.frame(
      id = c("l1", "l2", "l3"),
      inchikey = c("AAAAAAAAAAAAAA-XYZ", "BBBBBBBBBBBBBB-XYZ", "CCCCCCCCCCCCCC-XYZ"),
      stringsAsFactors = FALSE
    ),
    frag_list = list(
      l1 = matrix(c(50, 0.62, 100, 0.38), ncol = 2, byrow = TRUE),
      l2 = matrix(c(70, 0.68, 150, 0.32), ncol = 2, byrow = TRUE),
      l3 = matrix(c(200, 1.0), ncol = 2, byrow = TRUE)
    ),
    loss_list = list(l1 = empty_loss, l2 = empty_loss, l3 = empty_loss),
    loss_typ_list = NULL,
    loss_anchor_list = NULL,
    loss_pair_list = NULL,
    loss_anchor_typ_list = NULL,
    loss_pair_typ_list = NULL,
    mref_confidence = c(l1 = 1, l2 = 1, l3 = 1),
    ri = c(l1 = 110, l2 = 215, l3 = 600)
  )

  list(query = query, library = library)
}

test_that("benchmark_library_search summarizes query-library retrieval", {
  inputs <- make_library_benchmark_inputs()
  params <- ppmWass::eihrms_default_params()

  res <- ppmWass::benchmark_library_search(
    query_spectra = inputs$query,
    library_spectra = inputs$library,
    methods = "cosine",
    params = params,
    top_k = c(1, 2),
    keep_search_results = TRUE,
    search_top_n = 2,
    progress = FALSE
  )

  expect_true(is.list(res))
  expect_equal(nrow(res$summary), 1)
  expect_equal(res$summary$method, "cosine")
  expect_equal(res$summary$n_queries, 2)
  expect_true(res$summary$p_at_1 >= 0.5)
  expect_true("hits" %in% names(res$cosine))
})

test_that("queries with every candidate RI-excluded are scored as misses", {
  inputs <- make_library_benchmark_inputs()
  inputs$library$ri[] <- 2000
  params <- eihrms_default_params()
  params$use_parallel <- FALSE
  result <- benchmark_library_search(inputs$query, inputs$library,
    methods = "cosine", params = params, ri_tolerance = 50,
    top_k = c(1, 5), progress = FALSE)
  expect_equal(result$summary$n_queries, 2L)
  expect_equal(result$summary$map, 0)
  expect_equal(result$summary$mrr, 0)
  expect_equal(result$summary$p_at_1, 0)
  expect_equal(result$summary$p_at_5, 0)
  expect_equal(nrow(result$cosine$hits), 0L)
  expect_true(all(is.infinite(result$cosine$dist_matrix)))
})

test_that("missing RI candidates are excluded consistently from metrics and hits", {
  inputs <- make_library_benchmark_inputs()
  inputs$library$ri[] <- NA_real_
  params <- eihrms_default_params()
  params$use_parallel <- FALSE
  result <- benchmark_library_search(inputs$query, inputs$library,
    methods = "cosine", params = params, ri_tolerance = 50, progress = FALSE)
  expect_equal(result$summary$map, 0)
  expect_equal(result$summary$mrr, 0)
  expect_equal(nrow(result$cosine$hits), 0L)
})

test_that("nonfinite relevant candidates remain misses in average precision", {
  distance <- matrix(c(0.1, Inf, NA_real_, -Inf, 0.2), nrow = 1)
  query_keys <- "AAAAAAAAAAAAAA-111"
  library_keys <- c(rep(query_keys, 4), "BBBBBBBBBBBBBB-222")
  result <- ppmWass:::evaluate_library_search_retrieval(
    distance, query_keys, library_keys, top_k = c(1, 5))
  expect_equal(result$n_queries, 1L)
  expect_equal(result$mrr, 1)
  expect_equal(result$map, 1 / 4)
  expect_equal(result[["P@1"]], 1)
  expect_equal(result[["P@5"]], 1 / 2)

  distance[] <- Inf
  result <- ppmWass:::evaluate_library_search_retrieval(
    distance, query_keys, library_keys)
  expect_equal(result$map, 0)
  expect_equal(result$mrr, 0)
  expect_equal(result[["P@1"]], 0)
})
