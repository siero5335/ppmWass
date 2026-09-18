test_that("signal-noise mixture preserves the requested TIC weighting", {
  signal <- matrix(c(100, 3, 110, 1), ncol = 2, byrow = TRUE)
  noise <- matrix(c(200, 1, 210, 1), ncol = 2, byrow = TRUE)
  mixed <- mix_signal_noise_spectrum(signal, noise, signal_fraction = 0.2, merge_ppm = 0)

  expect_equal(sum(mixed[, 2]), 1)
  expect_equal(sum(mixed[mixed[, 1] < 150, 2]), 0.2)
  expect_equal(sum(mixed[mixed[, 1] > 150, 2]), 0.8)
})

test_that("scaling occurs before final normalization", {
  signal <- matrix(c(100, 100), ncol = 2)
  noise <- matrix(c(200, 1), ncol = 2)
  low <- mix_signal_noise_spectrum(signal, noise, signal_fraction = 0.05, merge_ppm = 0)
  expect_equal(unname(low[low[, 1] == 100, 2]), 0.05)
  expect_equal(unname(low[low[, 1] == 200, 2]), 0.95)
})

test_that("list mixtures avoid same-ID interferents and report diagnostics", {
  spectra <- list(
    a = matrix(c(100, 1), ncol = 2),
    b = matrix(c(200, 1), ncol = 2),
    c = matrix(c(300, 1), ncol = 2)
  )
  out <- mix_signal_noise_spectra(
    spectra,
    signal_fraction = 0.1,
    noise_components = 2,
    merge_ppm = 0,
    seed = 42
  )
  expect_equal(names(out$spectra), names(spectra))
  expect_equal(nrow(out$diagnostics), 3L)
  expect_true(all(out$diagnostics$realized_signal_fraction == 0.1))
  expect_false(any(vapply(seq_len(nrow(out$diagnostics)), function(i) {
    out$diagnostics$signal_id[[i]] %in%
      strsplit(out$diagnostics$noise_ids[[i]], ";", fixed = TRUE)[[1]]
  }, logical(1))))
})

test_that("mixture generation is reproducible and restores RNG state", {
  spectra <- list(
    a = matrix(c(100, 1, 101, 1), ncol = 2, byrow = TRUE),
    b = matrix(c(200, 1, 201, 1), ncol = 2, byrow = TRUE)
  )
  set.seed(99)
  before <- .Random.seed
  one <- mix_signal_noise_spectra(spectra, signal_fraction = 0.2,
                                  dropout_prob = 0.5, seed = 7)
  expect_equal(.Random.seed, before)
  two <- mix_signal_noise_spectra(spectra, signal_fraction = 0.2,
                                  dropout_prob = 0.5, seed = 7)
  expect_equal(one, two)
})

test_that("noise is required below full signal", {
  signal <- matrix(c(100, 1), ncol = 2)
  expect_error(
    mix_signal_noise_spectrum(signal, matrix(numeric(0), ncol = 2), signal_fraction = 0.5),
    "noise must contain"
  )
})
