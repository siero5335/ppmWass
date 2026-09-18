mk_spec <- function(mz, it) {
  m <- cbind(mz = as.numeric(mz), intensity = as.numeric(it))
  m
}

empty_spec <- function() {
  matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity")))
}

normalize_it <- function(x) {
  x <- as.numeric(x)
  x[x < 0 | !is.finite(x)] <- 0
  s <- sum(x)
  if (s > 0) x / s else x
}

test_that("Hellinger distance is symmetric and handles basic edge cases", {
  A <- mk_spec(c(50, 100, 150), normalize_it(c(1, 2, 3)))
  B <- mk_spec(c(50, 100, 150), normalize_it(c(3, 2, 1)))

  dAB <- ppmWass:::hellinger_distance(A, B, ppm = 20)
  dBA <- ppmWass:::hellinger_distance(B, A, ppm = 20)
  expect_equal(dAB, dBA, tolerance = 1e-12)

  dAA <- ppmWass:::hellinger_distance(A, A, ppm = 20)
  expect_true(abs(dAA) < 1e-12)

  d_empty <- ppmWass:::hellinger_distance(empty_spec(), A, ppm = 20)
  expect_equal(d_empty, 1)

  # Completely disjoint support should give Hellinger distance ~ 1
  C <- mk_spec(c(10, 20), c(0.6, 0.4))
  D <- mk_spec(c(200, 220), c(0.5, 0.5))
  dCD <- ppmWass:::hellinger_distance(C, D, ppm = 5)
  expect_true(abs(dCD - 1) < 1e-12)
})

test_that("parse_ei treats NULL and non-scalar inputs as empty spectra", {
  expect_equal(nrow(ppmWass:::parse_ei(NULL)), 0)
  expect_equal(nrow(ppmWass:::parse_ei(NA_character_)), 0)
  expect_equal(nrow(ppmWass:::parse_ei(c("50:100", "60:200"))), 0)
})

test_that("Hellinger distance obeys triangle inequality when aligned on a common m/z grid", {
  # IMPORTANT: Triangle inequality is guaranteed for Hellinger on a fixed probability space.
  # This test uses spectra sharing the *same* m/z values so pairwise alignment produces the same grid.
  mz <- c(45, 73.0474, 91.0548, 147.0654)
  A <- mk_spec(mz, normalize_it(c(1, 2, 3, 4)))
  B <- mk_spec(mz, normalize_it(c(4, 3, 2, 1)))
  C <- mk_spec(mz, normalize_it(c(2, 2, 2, 2)))

  dAC <- ppmWass:::hellinger_distance(A, C, ppm = 20)
  dAB <- ppmWass:::hellinger_distance(A, B, ppm = 20)
  dBC <- ppmWass:::hellinger_distance(B, C, ppm = 20)

  expect_true(dAC <= dAB + dBC + 1e-12)
})

test_that("Cosine distances are symmetric and handle identity/empty cases", {
  A <- mk_spec(c(50, 100, 150), normalize_it(c(1, 2, 3)))
  B <- mk_spec(c(50, 100, 150), normalize_it(c(3, 2, 1)))

  dAB <- ppmWass:::cosine_distance(A, B, ppm = 20)
  dBA <- ppmWass:::cosine_distance(B, A, ppm = 20)
  expect_equal(dAB, dBA, tolerance = 1e-12)

  dAA <- ppmWass:::cosine_distance(A, A, ppm = 20)
  expect_true(abs(dAA) < 1e-12)

  d_empty <- ppmWass:::cosine_distance(empty_spec(), A, ppm = 20)
  expect_equal(d_empty, 1)

  dwAB <- ppmWass:::weighted_cosine_distance(A, B, ppm = 20)
  dwBA <- ppmWass:::weighted_cosine_distance(B, A, ppm = 20)
  expect_equal(dwAB, dwBA, tolerance = 1e-12)

  dwAA <- ppmWass:::weighted_cosine_distance(A, A, ppm = 20)
  expect_true(abs(dwAA) < 1e-12)
})
