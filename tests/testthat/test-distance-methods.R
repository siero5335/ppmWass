mk_spec_distance <- function(mz, intensity) {
  cbind(mz = as.numeric(mz), intensity = as.numeric(intensity))
}

normalize_distance_it <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x) | x < 0] <- 0
  s <- sum(x)
  if (s > 0) x / s else x
}

test_that("Wasserstein distances are symmetric and zero for identical spectra", {
  skip_if_not_installed("transport")

  A <- mk_spec_distance(c(50, 100, 150), normalize_distance_it(c(1, 2, 3)))
  B <- mk_spec_distance(c(52, 101, 148), normalize_distance_it(c(3, 2, 1)))

  dAA <- ppmWass:::wasserstein_distance(A, A, ppm = 20, align = FALSE)
  dAB <- ppmWass:::wasserstein_distance(A, B, ppm = 20, align = FALSE)
  dBA <- ppmWass:::wasserstein_distance(B, A, ppm = 20, align = FALSE)
  dAB_align <- ppmWass:::wasserstein_distance(A, B, ppm = 20, align = TRUE)

  expect_equal(dAA, 0, tolerance = 1e-12)
  expect_equal(dAB, dBA, tolerance = 1e-12)
  expect_true(dAB >= 0 && dAB <= 1)
  expect_true(dAB_align >= 0 && dAB_align <= 1)
})

test_that("ppm-aware Wasserstein distance stays within bounds", {
  skip_if_not_installed("transport")

  A <- mk_spec_distance(c(73.0474, 147.0654, 221.1172), normalize_distance_it(c(0.5, 0.3, 0.2)))
  B <- mk_spec_distance(c(73.0478, 147.0660, 235.1328), normalize_distance_it(c(0.45, 0.35, 0.2)))

  d <- ppmWass:::compute_distance(
    A, B,
    method = "ppm_wasserstein",
    ppm = 20,
    wasserstein_transition_mult = 3
  )

  expect_true(is.finite(d))
  expect_true(d >= 0 && d <= 1)
})

test_that("Composite distance is symmetric and bounded", {
  A <- mk_spec_distance(c(45, 73.0474, 147.0654), normalize_distance_it(c(0.2, 0.5, 0.3)))
  B <- mk_spec_distance(c(45, 73.0474, 147.0654), normalize_distance_it(c(0.3, 0.4, 0.3)))

  dAA <- ppmWass:::composite_distance(A, A, ppm = 20)
  dAB <- ppmWass:::composite_distance(A, B, ppm = 20)
  dBA <- ppmWass:::composite_distance(B, A, ppm = 20)

  expect_equal(dAA, 0, tolerance = 1e-12)
  expect_equal(dAB, dBA, tolerance = 1e-12)
  expect_true(dAB >= 0 && dAB <= 1)
})

test_that("Entropy distance wrappers agree with msentropy", {
  skip_if_not_installed("msentropy")

  A <- mk_spec_distance(c(43.0184, 58.0419, 91.0548), normalize_distance_it(c(0.4, 0.35, 0.25)))
  B <- mk_spec_distance(c(43.0185, 58.0417, 105.0704), normalize_distance_it(c(0.45, 0.30, 0.25)))

  expected_weighted <- 1 - msentropy::msentropy_similarity(
    A, B,
    ms2_tolerance_in_da = -1,
    ms2_tolerance_in_ppm = 20,
    clean_spectra = FALSE,
    weighted = TRUE
  )
  expected_unweighted <- 1 - msentropy::msentropy_similarity(
    A, B,
    ms2_tolerance_in_da = -1,
    ms2_tolerance_in_ppm = 20,
    clean_spectra = FALSE,
    weighted = FALSE
  )

  expect_equal(ppmWass:::entropy_distance(A, B, ppm = 20), expected_weighted, tolerance = 1e-12)
  expect_equal(ppmWass:::entropy_unweighted_distance(A, B, ppm = 20), expected_unweighted, tolerance = 1e-12)
  expect_equal(ppmWass:::entropy_distance(A, A, ppm = 20), 0, tolerance = 1e-12)
})
