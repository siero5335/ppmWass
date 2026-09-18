test_that("perturb_spectrum preserves structure", {
  spec <- data.frame(mz = c(73, 147, 207), intensity = c(100, 50, 30))

  # No noise -> identical
  p0 <- perturb_spectrum(spec, seed = 1)
  expect_equal(p0, spec)

  # m/z noise -> m/z changes but intensity unchanged
  p1 <- perturb_spectrum(spec, mz_noise_ppm = 10, seed = 1)
  expect_equal(nrow(p1), 3)
  expect_false(all(p1[, 1] == spec[, 1]))
  expect_equal(p1[, 2], spec[, 2])

  # Intensity noise -> m/z unchanged, intensity changes
  p2 <- perturb_spectrum(spec, intensity_noise_sd = 0.3, seed = 1)
  expect_equal(p2[, 1], spec[, 1])
  expect_false(all(p2[, 2] == spec[, 2]))

  # Dropout -> fewer or equal peaks
  p3 <- perturb_spectrum(spec, dropout_prob = 0.5, seed = 1)
  expect_true(nrow(p3) <= 3)

  # Reproducibility
  p4a <- perturb_spectrum(spec, mz_noise_ppm = 10, seed = 42)
  p4b <- perturb_spectrum(spec, mz_noise_ppm = 10, seed = 42)
  expect_equal(p4a, p4b)
})

test_that("perturb_spectrum handles matrix input", {
  spec <- matrix(c(73, 147, 207, 100, 50, 30), ncol = 2,
                 dimnames = list(NULL, c("mz", "intensity")))
  p <- perturb_spectrum(spec, mz_noise_ppm = 10, seed = 1)
  expect_true(is.matrix(p))
  expect_equal(ncol(p), 2)
  expect_equal(nrow(p), 3)
})

test_that("perturb_spectrum handles empty spectrum", {
  spec <- data.frame(mz = numeric(0), intensity = numeric(0))
  p <- perturb_spectrum(spec, mz_noise_ppm = 10, seed = 1)
  expect_equal(nrow(p), 0)
})

test_that("perturb_spectrum handles full dropout", {
  spec <- data.frame(mz = c(73), intensity = c(100))
  # With prob = 0.99 and seed that causes dropout
  # Run multiple seeds to find one that drops the single peak
  any_empty <- FALSE
  for (s in 1:50) {
    p <- perturb_spectrum(spec, dropout_prob = 0.99, seed = s)
    if (nrow(p) == 0) { any_empty <- TRUE; break }
  }
  expect_true(any_empty)
})

test_that("perturb_spectra_list works on list", {
  specs <- list(
    a = data.frame(mz = c(73, 147), intensity = c(100, 50)),
    b = data.frame(mz = c(207, 281), intensity = c(80, 40))
  )
  result <- perturb_spectra_list(specs, mz_noise_ppm = 10, seed = 1)
  expect_equal(names(result), c("a", "b"))
  expect_equal(nrow(result$a), 2)
  expect_equal(nrow(result$b), 2)

  # Reproducibility
  r2 <- perturb_spectra_list(specs, mz_noise_ppm = 10, seed = 1)
  expect_equal(result, r2)
})
