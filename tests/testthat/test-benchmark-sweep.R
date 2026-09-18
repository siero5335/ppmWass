make_benchmark_spectra_result <- function() {
  empty_loss <- matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity")))
  frag_list <- list(
    A = matrix(c(50, 0.5, 100, 0.5), ncol = 2, byrow = TRUE),
    B = matrix(c(50, 0.52, 100, 0.48), ncol = 2, byrow = TRUE),
    C = matrix(c(70, 0.5, 150, 0.5), ncol = 2, byrow = TRUE)
  )
  loss_list <- list(A = empty_loss, B = empty_loss, C = empty_loss)
  df_spec <- data.frame(
    id = c("A", "B", "C"),
    inchikey = c("AAAAAAAAAAAAAA-BBB", "AAAAAAAAAAAAAA-CCC", "ZZZZZZZZZZZZZZ-DDD"),
    compound_class = c("alkane", "alkane", "aromatic"),
    compound_class_source = c("external", "external", "external"),
    stringsAsFactors = FALSE
  )
  ri <- c(A = 100, B = 180, C = 420)

  list(
    df_spec = df_spec,
    frag_list = frag_list,
    loss_list = loss_list,
    loss_typ_list = NULL,
    loss_anchor_list = NULL,
    loss_pair_list = NULL,
    loss_anchor_typ_list = NULL,
    loss_pair_typ_list = NULL,
    mref_confidence = c(A = 1, B = 1, C = 1),
    ri = ri
  )
}

test_that("sweep_distance_params expands grids and returns summary rows", {
  spectra <- make_benchmark_spectra_result()
  params <- ppmWass::eihrms_default_params()
  grid <- data.frame(
    tol_ppm = c(10, 20),
    w_frag = c(0.7, 0.5),
    w_loss = c(0.3, 0.5),
    stringsAsFactors = FALSE
  )

  res <- ppmWass::sweep_distance_params(
    spectra_result = spectra,
    grid = grid,
    methods = c("cosine", "hellinger"),
    params = params,
    metrics = c("auc_roc", "map"),
    progress = FALSE
  )

  expect_true(is.list(res))
  expect_equal(nrow(res$summary), 4)
  expect_true(all(c("sweep_id", "replicate", "n_spectra", "tol_ppm", "w_frag", "w_loss", "method", "auc_roc", "map") %in% names(res$summary)))
})

test_that("sweep_distance_params rejects unknown parameter names", {
  spectra <- make_benchmark_spectra_result()

  expect_error(
    ppmWass::sweep_distance_params(
      spectra_result = spectra,
      grid = data.frame(not_a_param = 1),
      methods = "cosine",
      progress = FALSE
    ),
    "Unknown parameter name"
  )
})

test_that("sweep_distance_params supports repeated subsampling", {
  spectra <- make_benchmark_spectra_result()

  res <- ppmWass::sweep_distance_params(
    spectra_result = spectra,
    grid = data.frame(tol_ppm = 15),
    methods = "cosine",
    n_replicates = 3,
    sample_frac = 2 / 3,
    progress = FALSE
  )

  expect_equal(nrow(res$summary), 3)
  expect_equal(sort(unique(res$summary$replicate)), 1:3)
  expect_true(all(res$summary$n_spectra == 2))
})

test_that("benchmark_distance_runtime returns timing summaries", {
  spectra <- make_benchmark_spectra_result()
  params <- ppmWass::eihrms_default_params()

  res <- ppmWass::benchmark_distance_runtime(
    spectra_result = spectra,
    methods = "cosine",
    sizes = c(2, 3),
    n_replicates = 2,
    params = params,
    progress = FALSE
  )

  expect_true(is.list(res))
  expect_equal(nrow(res$results), 4)
  expect_equal(nrow(res$summary), 2)
  expect_true(all(res$results$elapsed_sec >= 0))
  expect_true(all(c("method", "size", "n_pairs", "mean", "sd", "median", "min", "max") %in% names(res$summary)))
})

test_that("package rank AUC remains finite for large pair counts", {
  rank_auc <- getFromNamespace("rank_binary_auc", "ppmWass")
  n_positive <- 50000L
  n_negative <- 50000L
  score <- c(rep(1, n_positive), rep(0, n_negative))
  is_positive <- c(rep(TRUE, n_positive), rep(FALSE, n_negative))

  expect_gt(as.double(n_positive) * as.double(n_negative),
            .Machine$integer.max)
  expect_equal(rank_auc(score, is_positive), 1, tolerance = 1e-12)
})
