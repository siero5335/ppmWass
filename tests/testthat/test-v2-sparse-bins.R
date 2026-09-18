test_that("v2 Hellinger keeps matrix shape", {
  frag <- list(
    a = matrix(c(100, 1, 101, 2), ncol = 2, byrow = TRUE),
    b = matrix(c(100.001, 1, 120, 1), ncol = 2, byrow = TRUE)
  )
  Q <- as_spectra_set(frag, ids = names(frag), bin_ppm = 2, smear_ppm = 20)
  D <- compute_distance_matrix_search_v2(
    Q, Q,
    params = list(distance_method = "hellinger", w_frag = 1, w_loss = 0),
    progress = FALSE
  )
  expect_true(is.matrix(D))
  expect_equal(dim(D), c(2L, 2L))
  expect_equal(rownames(D), names(frag))
  expect_equal(colnames(D), names(frag))
})

test_that("default search backend remains exact pair-loop", {
  frag <- list(
    q1 = matrix(c(100, 1, 101, 2), ncol = 2, byrow = TRUE),
    q2 = matrix(c(120, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    q1 = matrix(c(10, 1), ncol = 2, byrow = TRUE),
    q2 = matrix(c(20, 1), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$distance_method <- "cosine"

  D_default <- compute_distance_matrix_search(frag, loss, frag, loss, params, progress = FALSE)
  D_pair <- ppmWass:::compute_distance_matrix_search_pairloop(
    frag, loss, frag, loss, params, progress = FALSE
  )

  expect_equal(D_default, D_pair)
})

test_that("pair-loop parallel search matches serial search", {
  frag <- list(
    a = matrix(c(100, 1, 101, 2), ncol = 2, byrow = TRUE),
    b = matrix(c(100.001, 1, 120, 1), ncol = 2, byrow = TRUE),
    c = matrix(c(300, 5), ncol = 2, byrow = TRUE),
    d = matrix(c(150, 1, 151, 3), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    a = matrix(c(10, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(20, 1), ncol = 2, byrow = TRUE),
    c = matrix(c(30, 1), ncol = 2, byrow = TRUE),
    d = matrix(c(40, 1), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$distance_method <- "cosine"
  params$use_parallel <- FALSE
  serial <- compute_distance_matrix_search(frag, loss, frag, loss, params, progress = FALSE)

  params$use_parallel <- TRUE
  params$n_cores <- 2
  parallel <- compute_distance_matrix_search(frag, loss, frag, loss, params, progress = FALSE)

  expect_equal(parallel, serial)
})

test_that("prepared pair-loop search matches per-pair combined_distance", {
  frag <- list(
    a = matrix(c(101, 2, 100, 1, 140, 3), ncol = 2, byrow = TRUE),
    b = matrix(c(100.001, 1, 120, 2), ncol = 2, byrow = TRUE),
    c = matrix(c(300, 5, 150, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    a = matrix(c(10, 1, 11, 3), ncol = 2, byrow = TRUE),
    b = matrix(c(10.0001, 2, 25, 1), ncol = 2, byrow = TRUE),
    c = matrix(c(30, 1, 31, 2), ncol = 2, byrow = TRUE)
  )

  for (method in c("hellinger", "cosine", "weighted_cosine", "composite")) {
    params <- eihrms_default_params()
    params$distance_method <- method
    params$use_parallel <- FALSE
    D <- compute_distance_matrix_search(frag, loss, frag, loss, params, progress = FALSE)

    D_ref <- matrix(0, length(frag), length(frag), dimnames = list(names(frag), names(frag)))
    for (i in seq_along(frag)) {
      for (j in seq_along(frag)) {
        D_ref[i, j] <- ppmWass:::combined_distance(
          frag[[i]], frag[[j]], loss[[i]], loss[[j]], params
        )
      }
    }

    expect_equal(D, D_ref, tolerance = 1e-14, info = method)
  }
})

test_that("prepared symmetric pair-loop matches per-pair combined_distance", {
  frag <- list(
    a = matrix(c(101, 2, 100, 1, 140, 3), ncol = 2, byrow = TRUE),
    b = matrix(c(100.001, 1, 120, 2), ncol = 2, byrow = TRUE),
    c = matrix(c(300, 5, 150, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    a = matrix(c(10, 1, 11, 3), ncol = 2, byrow = TRUE),
    b = matrix(c(10.0001, 2, 25, 1), ncol = 2, byrow = TRUE),
    c = matrix(c(30, 1, 31, 2), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$distance_method <- "weighted_cosine"
  params$use_parallel <- FALSE

  D <- ppmWass:::compute_distance_matrix(frag, loss, params, progress = FALSE)
  D_ref <- matrix(0, length(frag), length(frag), dimnames = list(names(frag), names(frag)))
  for (i in seq_along(frag)) {
    for (j in seq_along(frag)) {
      if (i != j) {
        D_ref[i, j] <- ppmWass:::combined_distance(
          frag[[i]], frag[[j]], loss[[i]], loss[[j]], params
        )
      }
    }
  }

  expect_equal(D, D_ref, tolerance = 1e-14)
})

test_that("parallel core resolver handles requested core counts", {
  expect_gte(ppmWass:::resolve_parallel_cores(1), 1)
  expect_gte(ppmWass:::resolve_parallel_cores(2), 1)
})

test_that("sparse-bin shim opt-in matches direct v2 dense API", {
  frag <- list(
    q1 = matrix(c(100, 1, 101, 2), ncol = 2, byrow = TRUE),
    q2 = matrix(c(100.001, 1, 120, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    q1 = matrix(c(10, 1), ncol = 2, byrow = TRUE),
    q2 = matrix(c(20, 1), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$backend <- "sparse_bins"
  params$distance_method <- "cosine"
  params$min_mz <- 90
  params$max_mz <- 130
  params$bin_ppm <- 2
  params$smear_ppm <- 20

  D_shim <- compute_distance_matrix_search(frag, loss, frag, loss, params, progress = FALSE)
  Q <- as_spectra_set(
    frag, loss, ids = names(frag),
    mz_min = params$min_mz, mz_max = params$max_mz,
    bin_ppm = params$bin_ppm, smear_ppm = params$smear_ppm
  )
  D_v2 <- compute_distance_matrix_search_v2(Q, Q, params, progress = FALSE)

  expect_equal(D_shim, D_v2)
})

test_that("sparse-bin shim rejects common backend typo", {
  frag <- list(
    a = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(101, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    a = matrix(c(10, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(11, 1), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$backedn <- "sparse_bins"

  expect_error(
    compute_distance_matrix_search(frag, loss, frag, loss, params, progress = FALSE),
    "did you mean 'backend'"
  )
})

test_that("direct v2 APIs reject common backend typo", {
  frag <- list(
    a = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(101, 1), ncol = 2, byrow = TRUE)
  )
  Q <- as_spectra_set(frag, ids = names(frag))
  params <- list(distance_method = "cosine", backedn = "sparse_bins")

  expect_error(
    compute_distance_matrix_search_v2(Q, Q, params, progress = FALSE),
    "did you mean 'backend'"
  )
  expect_error(
    compute_distance_topk_v2(Q, Q, params, progress = FALSE),
    "did you mean 'backend'"
  )
})

test_that("v2 direct API stops for unsupported legacy-only options", {
  frag <- list(
    a = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(101, 1), ncol = 2, byrow = TRUE)
  )
  Q <- as_spectra_set(frag, ids = names(frag))
  params <- eihrms_default_params()
  params$distance_method <- "cosine"
  params$use_typical_loss <- TRUE

  expect_error(
    compute_distance_matrix_search_v2(Q, Q, params, progress = FALSE),
    "cannot safely reproduce"
  )
})

test_that("sparse-bin shim falls back instead of silently ignoring unsupported options", {
  frag <- list(
    a = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(101, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    a = matrix(c(10, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(11, 1), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$backend <- "sparse_bins"
  params$distance_method <- "cosine"
  params$use_typical_loss <- TRUE

  expect_message(
    D <- compute_distance_matrix_search(frag, loss, frag, loss, params, progress = TRUE),
    "falling back to pair_loop"
  )
  expect_true(is.matrix(D))
  expect_equal(dim(D), c(2L, 2L))
})

test_that("sparse-bin shim falls back for non-batch distance methods", {
  frag <- list(
    a = matrix(c(100, 1, 101, 2), ncol = 2, byrow = TRUE),
    b = matrix(c(100.001, 1, 120, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    a = matrix(c(10, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(20, 1), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$backend <- "sparse_bins"
  params$distance_method <- "composite"

  expect_message(
    D <- compute_distance_matrix_search(frag, loss, frag, loss, params, progress = TRUE),
    "distance_method 'composite'"
  )
  expect_true(is.matrix(D))
  expect_equal(dim(D), c(2L, 2L))
})

test_that("sparse-bin shim forwards dense-size guard", {
  frag <- list(
    a = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(101, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    a = matrix(c(10, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(11, 1), ncol = 2, byrow = TRUE)
  )
  params <- eihrms_default_params()
  params$backend <- "sparse_bins"
  params$distance_method <- "cosine"
  params$max_dense_cells <- 3
  params$min_mz <- 90
  params$max_mz <- 110

  expect_error(
    compute_distance_matrix_search(frag, loss, frag, loss, params, progress = FALSE),
    "Refusing to materialise"
  )
})

test_that("top-k refinement keeps prefilter distances aligned to returned candidates", {
  frag <- list(
    a = matrix(c(100, 1, 101, 2), ncol = 2, byrow = TRUE),
    b = matrix(c(100.001, 1, 120, 1), ncol = 2, byrow = TRUE),
    c = matrix(c(300, 5), ncol = 2, byrow = TRUE)
  )
  Q <- as_spectra_set(frag, ids = names(frag), bin_ppm = 2, smear_ppm = 20)
  res <- compute_distance_topk_v2(
    Q, Q,
    params = list(distance_method = "cosine", w_frag = 1, w_loss = 0),
    top_k = 3,
    expensive_method = "hellinger",
    expensive_top_k = 2,
    return_prefilter = TRUE,
    progress = FALSE
  )

  expect_equal(dim(res$topk_idx), c(3L, 2L))
  expect_equal(dim(res$topk_dist), c(3L, 2L))
  expect_equal(dim(res$topk_dist_pre), c(3L, 2L))
  expect_equal(rownames(res$topk_dist_pre), names(frag))
})

test_that("top-k refinement parallel path matches serial path", {
  frag <- list(
    a = matrix(c(100, 1, 101, 2), ncol = 2, byrow = TRUE),
    b = matrix(c(100.001, 1, 120, 1), ncol = 2, byrow = TRUE),
    c = matrix(c(300, 5), ncol = 2, byrow = TRUE),
    d = matrix(c(150, 1, 151, 3), ncol = 2, byrow = TRUE)
  )
  Q <- as_spectra_set(frag, ids = names(frag), bin_ppm = 2, smear_ppm = 20)
  params <- list(distance_method = "cosine", w_frag = 1, w_loss = 0)

  serial <- compute_distance_topk_v2(
    Q, Q, params,
    top_k = 3,
    expensive_method = "hellinger",
    expensive_top_k = 2,
    progress = FALSE
  )

  params$use_parallel <- TRUE
  params$n_cores <- 2
  parallel <- compute_distance_topk_v2(
    Q, Q, params,
    top_k = 3,
    expensive_method = "hellinger",
    expensive_top_k = 2,
    progress = FALSE
  )

  expect_equal(parallel$topk_idx, serial$topk_idx)
  expect_equal(parallel$topk_dist, serial$topk_dist)
})

test_that("top-k ordering pads unavailable finite candidates instead of stopping", {
  ord <- ppmWass:::top_k_order(c(0.2, Inf, NA, 0.1), 3)
  expect_equal(ord, c(4L, 1L, NA_integer_))

  ord_all_bad <- ppmWass:::top_k_order(c(Inf, NA), 2)
  expect_equal(ord_all_bad, c(NA_integer_, NA_integer_))
})

test_that("top-k ordering is stable for ties by library index", {
  ord <- ppmWass:::top_k_order(c(0.1, 0.2, 0.1, 0.1), 3)
  expect_equal(ord, c(1L, 3L, 4L))
})

test_that("top-k handles rectangular query-library search and clamps k", {
  query <- list(
    q1 = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    q2 = matrix(c(130, 1), ncol = 2, byrow = TRUE)
  )
  library <- list(
    l1 = matrix(c(100.001, 1), ncol = 2, byrow = TRUE),
    l2 = matrix(c(131, 1), ncol = 2, byrow = TRUE),
    l3 = matrix(c(200, 1), ncol = 2, byrow = TRUE)
  )
  Q <- as_spectra_set(query, ids = names(query), mz_min = 90, mz_max = 210)
  L <- as_spectra_set(library, ids = names(library), mz_min = 90, mz_max = 210)

  res <- compute_distance_topk_v2(
    Q, L,
    params = list(distance_method = "cosine", w_frag = 1, w_loss = 0),
    top_k = 10,
    block_size = 1,
    progress = FALSE
  )

  expect_equal(dim(res$topk_idx), c(2L, 3L))
  expect_equal(res$query_ids, names(query))
  expect_equal(res$library_ids, names(library))
  expect_true(all(res$topk_idx >= 1L & res$topk_idx <= 3L))
})

test_that("top-k requires a batchable stage-1 method", {
  frag <- list(
    a = matrix(c(100, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(101, 1), ncol = 2, byrow = TRUE)
  )
  Q <- as_spectra_set(frag, ids = names(frag))

  expect_error(
    compute_distance_topk_v2(
      Q, Q,
      params = list(distance_method = "ppm_wasserstein"),
      progress = FALSE
    ),
    "requires a batchable"
  )
})

test_that("topk_accuracy can exclude trivial self hits", {
  result <- list(
    query_ids = c("a", "b"),
    library_ids = c("a", "x", "b"),
    topk_idx = matrix(c(1L, 2L, 3L, 2L), nrow = 2, byrow = TRUE),
    topk_dist = matrix(c(0, 0.2, 0, 0.3), nrow = 2, byrow = TRUE)
  )
  class(result) <- c("TopKResult", "list")

  expect_equal(topk_accuracy(result, top_k = 1), c(top_1 = 1))
  expect_equal(
    topk_accuracy(result, top_k = 1, ground_truth = list(a = "x", b = "x"), exclude_self = TRUE),
    c(top_1 = 1)
  )
})
