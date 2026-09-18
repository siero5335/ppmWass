test_that("default parameters validate and allow ppm_wasserstein", {
  params <- ppmWass::eihrms_default_params()
  params$distance_method <- "ppm_wasserstein"

  expect_equal(ppmWass:::validate_params(params), params)
})

test_that("validate_params rejects invalid numeric ranges", {
  params <- ppmWass::eihrms_default_params()
  params$tol_ppm <- 0
  expect_error(ppmWass:::validate_params(params), "tol_ppm")

  params <- ppmWass::eihrms_default_params()
  params$loss_min <- 50
  params$loss_max <- 20
  expect_error(ppmWass:::validate_params(params), "loss_max must be greater than loss_min")

  params <- ppmWass::eihrms_default_params()
  params$analog_sim_threshold <- 1.5
  expect_error(ppmWass:::validate_params(params), "analog_sim_threshold")

  params <- ppmWass::eihrms_default_params()
  params$n_cores <- 0
  expect_error(ppmWass:::validate_params(params), "n_cores")

  params <- ppmWass::eihrms_default_params()
  params$w_frag <- 0
  params$w_loss <- 0
  expect_error(ppmWass:::validate_params(params), "w_frag \\+ w_loss")
})

test_that("validate_params checks sparse-bin backend parameters", {
  params <- ppmWass::eihrms_default_params()
  params$backend <- "sparse_bins"
  params$bin_ppm <- 2
  params$smear_ppm <- 15
  params$smear_kind <- "tri"
  params$max_dense_cells <- 1000
  params$prefilter_method <- "cosine"
  params$prefilter_top_k <- 50
  expect_equal(ppmWass:::validate_params(params), params)

  params_bad <- params
  params_bad$backend <- "fast"
  expect_error(ppmWass:::validate_params(params_bad), "backend")

  params_bad <- params
  params_bad$bin_ppm <- 0
  expect_error(ppmWass:::validate_params(params_bad), "bin_ppm")

  params_bad <- params
  params_bad$smear_kind <- "gaussian"
  expect_error(ppmWass:::validate_params(params_bad), "smear_kind")

  params_bad <- params
  params_bad$max_dense_cells <- 1.5
  expect_error(ppmWass:::validate_params(params_bad), "max_dense_cells")

  params_bad <- params
  params_bad$prefilter_top_k <- 0
  expect_error(ppmWass:::validate_params(params_bad), "prefilter_top_k")

  params_bad <- params
  params_bad$backedn <- "sparse_bins"
  expect_error(ppmWass:::validate_params(params_bad), "did you mean 'backend'")
})

test_that("validate_params warns for non-normalized weights and conflicting heuristics", {
  params <- ppmWass::eihrms_default_params()
  params$w_frag <- 0.8
  params$w_loss <- 0.5
  expect_warning(ppmWass:::validate_params(params), "w_frag \\+ w_loss is not 1.0")

  params <- ppmWass::eihrms_default_params()
  params$mref_conf_w_nloss <- 0.5
  params$mref_conf_w_intensity <- 0.5
  params$mref_conf_w_hi <- 0.5
  expect_warning(ppmWass:::validate_params(params), "Mref confidence component weights do not sum to 1")

  params <- ppmWass::eihrms_default_params()
  params$loss_top_peaks <- 80
  params$loss_max_peaks <- 40
  expect_warning(ppmWass:::validate_params(params), "loss_max_peaks is smaller than loss_top_peaks")

  params <- ppmWass::eihrms_default_params()
  params$derivatization_auto_min_int <- 0.2
  params$derivatization_auto_strong_int <- 0.1
  expect_warning(ppmWass:::validate_params(params), "derivatization_auto_strong_int is smaller")
})

test_that("validate_params checks normalize_cols type", {
  params <- ppmWass::eihrms_default_params()
  params$normalize_cols <- c("Sample_1", "")

  expect_error(ppmWass:::validate_params(params), "normalize_cols")
})
