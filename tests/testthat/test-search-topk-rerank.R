toy_topk_frag <- function() {
  list(
    q1 = matrix(c(100, 1, 101, 0.5), ncol = 2, byrow = TRUE),
    q2 = matrix(c(130, 1, 131, 0.5), ncol = 2, byrow = TRUE),
    q3 = matrix(c(170, 1, 171, 0.5), ncol = 2, byrow = TRUE)
  )
}

test_that("first_pass_top_k = n_library agrees with full rerank search ranking", {
  query <- toy_topk_frag()[1:2]
  library <- toy_topk_frag()
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0
  params$use_parallel <- FALSE

  res <- search_topk_rerank(
    query_frag_list = query,
    lib_frag_list = library,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "weighted_cosine",
    first_pass_top_k = length(library),
    final_top_k = 2,
    progress = FALSE
  )

  params_full <- params
  params_full$distance_method <- "weighted_cosine"
  empty_q <- rep(list(matrix(numeric(0), ncol = 2)), length(query))
  empty_l <- rep(list(matrix(numeric(0), ncol = 2)), length(library))
  D_full <- compute_distance_matrix_search(
    query, empty_q, library, empty_l, params_full, progress = FALSE
  )

  for (qid in names(query)) {
    expected <- names(library)[ppmWass:::top_k_order(D_full[qid, ], 2)]
    observed <- res$results$lib_id[res$results$query_id == qid]
    expect_equal(observed, expected)
  }
})

test_that("matching first-pass and rerank methods short-circuit distances", {
  frag <- toy_topk_frag()
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  res <- search_topk_rerank(
    query_frag_list = frag,
    lib_frag_list = frag,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "hellinger",
    first_pass_top_k = 3,
    final_top_k = 2,
    exclude_self = FALSE,
    progress = FALSE
  )

  expect_equal(res$results$rerank_distance, res$results$first_pass_distance)
})

test_that("exclude_self removes identical query-library IDs", {
  frag <- toy_topk_frag()
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  res <- search_topk_rerank(
    query_frag_list = frag,
    lib_frag_list = frag,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "hellinger",
    first_pass_top_k = 3,
    final_top_k = 1,
    exclude_self = TRUE,
    progress = FALSE
  )

  expect_false(any(res$results$query_id == res$results$lib_id))
})

test_that("identical query and library IDs auto-enable library by library Top-K", {
  frag <- toy_topk_frag()
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  res <- search_topk_rerank(
    query_frag_list = frag,
    lib_frag_list = frag,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "hellinger",
    first_pass_top_k = 2,
    final_top_k = 1,
    progress = FALSE
  )

  expect_equal(sort(unique(res$results$query_id)), sort(names(frag)))
  expect_false(any(res$results$query_id == res$results$lib_id))
})

test_that("fragment-only search works with NULL loss lists", {
  frag <- toy_topk_frag()
  params <- eihrms_default_params()

  res <- search_topk_rerank(
    query_frag_list = frag[1],
    query_loss_list = NULL,
    lib_frag_list = frag,
    lib_loss_list = NULL,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "weighted_cosine",
    first_pass_top_k = 3,
    final_top_k = 2,
    progress = FALSE
  )

  expect_equal(nrow(res$results), 2L)
  expect_equal(res$params_first_pass$w_loss, 0)
  expect_equal(res$params_rerank$w_loss, 0)
})

test_that("one-sided NULL loss input is rejected", {
  frag <- toy_topk_frag()
  loss <- lapply(frag, function(x) matrix(c(10, 1), ncol = 2, byrow = TRUE))
  params <- eihrms_default_params()

  expect_error(
    search_topk_rerank(
      query_frag_list = frag,
      query_loss_list = loss,
      lib_frag_list = frag,
      lib_loss_list = NULL,
      params = params,
      first_pass_method = "hellinger",
      rerank_method = "hellinger",
      first_pass_top_k = 2,
      final_top_k = 1,
      progress = FALSE
    ),
    "both be NULL or both be non-NULL"
  )
})

test_that("split-loss reranking is explicit unsupported in v1", {
  frag <- toy_topk_frag()
  loss <- lapply(frag, function(x) matrix(c(10, 1), ncol = 2, byrow = TRUE))
  params <- eihrms_default_params()
  params$use_split_loss <- TRUE

  expect_error(
    search_topk_rerank(
      query_frag_list = frag,
      query_loss_list = loss,
      lib_frag_list = frag,
      lib_loss_list = loss,
      params = params,
      first_pass_method = "hellinger",
      rerank_method = "hellinger",
      first_pass_top_k = 2,
      final_top_k = 1,
      progress = FALSE
    ),
    "does not support split-loss reranking"
  )
})

test_that("renamed typical-loss arguments are accepted", {
  frag <- toy_topk_frag()
  loss <- lapply(frag, function(x) matrix(c(10, 1), ncol = 2, byrow = TRUE))
  typ <- lapply(frag, function(x) matrix(c(12, 1), ncol = 2, byrow = TRUE))
  params <- eihrms_default_params()
  params$distance_method <- "hellinger"
  params$use_typical_loss <- TRUE

  res <- search_topk_rerank(
    query_frag_list = frag[1],
    query_loss_list = loss[1],
    lib_frag_list = frag,
    lib_loss_list = loss,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "hellinger",
    first_pass_top_k = 2,
    final_top_k = 1,
    query_loss_typ_list = typ[1],
    lib_loss_typ_list = typ,
    progress = FALSE
  )

  expect_equal(nrow(res$results), 1L)
  expect_true(isTRUE(res$params_rerank$use_typical_loss))
})

test_that("mref confidence arguments match lower-level naming", {
  frag <- toy_topk_frag()
  loss <- lapply(frag, function(x) matrix(c(10, 1), ncol = 2, byrow = TRUE))
  params <- eihrms_default_params()
  params$use_mref_confidence <- TRUE

  res <- search_topk_rerank(
    query_frag_list = frag[1],
    query_loss_list = loss[1],
    lib_frag_list = frag,
    lib_loss_list = loss,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "hellinger",
    first_pass_top_k = 2,
    final_top_k = 1,
    query_mref_conf = 1,
    lib_mref_conf = c(1, 0.5, 0.25),
    progress = FALSE
  )

  expect_equal(nrow(res$results), 1L)
  expect_true(isTRUE(res$params_rerank$use_mref_confidence))
})

test_that("parallel and serial reranking agree", {
  frag <- toy_topk_frag()
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  serial <- search_topk_rerank(
    query_frag_list = frag,
    lib_frag_list = frag,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "weighted_cosine",
    first_pass_top_k = 3,
    final_top_k = 2,
    exclude_self = FALSE,
    use_parallel = FALSE,
    progress = FALSE
  )

  parallel <- search_topk_rerank(
    query_frag_list = frag,
    lib_frag_list = frag,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "weighted_cosine",
    first_pass_top_k = 3,
    final_top_k = 2,
    exclude_self = FALSE,
    use_parallel = TRUE,
    n_cores = 2,
    progress = FALSE
  )

  expect_equal(parallel$results, serial$results)
  expect_equal(
    parallel$first_pass_distance_matrix,
    serial$first_pass_distance_matrix
  )
})

test_that("finite candidate shortage returns short results instead of erroring", {
  frag <- list(a = matrix(c(100, 1), ncol = 2, byrow = TRUE))
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  expect_warning(
    res <- search_topk_rerank(
      query_frag_list = frag,
      lib_frag_list = frag,
      params = params,
      first_pass_method = "hellinger",
      rerank_method = "hellinger",
      first_pass_top_k = 1,
      final_top_k = 1,
      exclude_self = TRUE,
      progress = FALSE
    ),
    "Fewer than final_top_k"
  )

  expect_equal(nrow(res$results), 0L)
  expect_equal(res$candidate_summary$n_returned, 0L)
})

test_that("ties are determined by library order", {
  query <- list(q = matrix(c(100, 1), ncol = 2, byrow = TRUE))
  library <- list(
    b = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    a = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    c = matrix(c(150, 1), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  res <- search_topk_rerank(
    query_frag_list = query,
    lib_frag_list = library,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "hellinger",
    first_pass_top_k = 3,
    final_top_k = 2,
    progress = FALSE
  )

  expect_equal(res$results$lib_id, c("b", "a"))
})

test_that("return_first_pass = FALSE drops first-pass columns", {
  frag <- toy_topk_frag()
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  res <- search_topk_rerank(
    query_frag_list = frag[1],
    lib_frag_list = frag,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "hellinger",
    first_pass_top_k = 3,
    final_top_k = 1,
    return_first_pass = FALSE,
    progress = FALSE
  )

  expect_false("first_pass_rank" %in% names(res$results))
  expect_false("first_pass_distance" %in% names(res$results))
  expect_null(res$first_pass_distance_matrix)
})

test_that("return_first_pass = TRUE returns a fully aligned distance matrix", {
  frag <- toy_topk_frag()
  query <- frag[c("q2", "q1")]
  library <- frag[c("q3", "q1", "q2")]
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  res <- search_topk_rerank(
    query_frag_list = query,
    lib_frag_list = library,
    params = params,
    first_pass_method = "hellinger",
    rerank_method = "weighted_cosine",
    first_pass_top_k = 3,
    final_top_k = 1,
    exclude_self = FALSE,
    return_first_pass = TRUE,
    progress = FALSE
  )

  params_first <- params
  params_first$distance_method <- "hellinger"
  empty_query <- rep(list(matrix(numeric(0), ncol = 2)), length(query))
  empty_library <- rep(list(matrix(numeric(0), ncol = 2)), length(library))
  expected <- compute_distance_matrix_search(
    query, empty_query, library, empty_library, params_first, progress = FALSE
  )

  expect_true(is.matrix(res$first_pass_distance_matrix))
  expect_identical(rownames(res$first_pass_distance_matrix), names(query))
  expect_identical(colnames(res$first_pass_distance_matrix), names(library))
  expect_equal(res$first_pass_distance_matrix, expected)
})

test_that("first_pass_top_k larger than library size is clamped", {
  frag <- toy_topk_frag()[1:2]
  params <- eihrms_default_params()
  params$w_frag <- 1
  params$w_loss <- 0

  expect_message(
    res <- search_topk_rerank(
      query_frag_list = frag[1],
      lib_frag_list = frag,
      params = params,
      first_pass_method = "hellinger",
      rerank_method = "hellinger",
      first_pass_top_k = 10,
      final_top_k = 2,
      exclude_self = FALSE,
      progress = FALSE
    ),
    "clamping"
  )

  expect_equal(res$candidate_summary$n_candidates, 2L)
})

test_that("final_top_k larger than first_pass_top_k is an error", {
  frag <- toy_topk_frag()
  params <- eihrms_default_params()

  expect_error(
    search_topk_rerank(
      query_frag_list = frag,
      lib_frag_list = frag,
      params = params,
      first_pass_method = "hellinger",
      rerank_method = "hellinger",
      first_pass_top_k = 1,
      final_top_k = 2,
      progress = FALSE
    ),
    "final_top_k"
  )
})
