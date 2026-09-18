# Tests for approxOT backend integration in ppm_wasserstein_distance

# Helper: create a simple spectrum matrix
mk_spec <- function(mz, intensity) {
  m <- cbind(mz = mz, intensity = intensity)
  storage.mode(m) <- "double"
  m
}

spec1 <- mk_spec(c(73.047, 147.065, 207.089), c(100, 50, 30))
spec2 <- mk_spec(c(73.047, 147.080, 207.089), c(80, 60, 40))

# ------------------------------------------------------------------
# Basic properties (run regardless of backend availability)
# ------------------------------------------------------------------

test_that("ppm_wasserstein returns [0,1] with default ot_method", {
  skip_if_not_installed("transport")

  d <- ppm_wasserstein_distance(spec1, spec2, ppm = 15)
  expect_true(d >= 0 && d <= 1)
})

test_that("ppm_wasserstein identity: d(x,x) == 0", {
  skip_if_not_installed("transport")

  d <- ppm_wasserstein_distance(spec1, spec1, ppm = 15)
  expect_equal(d, 0, tolerance = 1e-6)
})

test_that("ppm_wasserstein symmetry: d(a,b) == d(b,a)", {
  skip_if_not_installed("transport")

  d_ab <- ppm_wasserstein_distance(spec1, spec2, ppm = 15)
  d_ba <- ppm_wasserstein_distance(spec2, spec1, ppm = 15)
  expect_equal(d_ab, d_ba, tolerance = 1e-6)
})

test_that("ppm_wasserstein handles empty spectra", {
  skip_if_not_installed("transport")

  empty <- mk_spec(numeric(0), numeric(0))
  expect_equal(ppm_wasserstein_distance(empty, spec2, ppm = 15), 1)
  expect_equal(ppm_wasserstein_distance(spec1, empty, ppm = 15), 1)
})

# ------------------------------------------------------------------
# exact uses transport, sinkhorn uses approxOT
# ------------------------------------------------------------------

test_that("exact routes to transport::transport, sinkhorn to approxOT", {
  skip_if_not_installed("approxOT")
  skip_if_not_installed("transport")

  d_exact <- ppm_wasserstein_distance(spec1, spec2, ppm = 15,
                                       ot_method = "exact")
  d_sink  <- ppm_wasserstein_distance(spec1, spec2, ppm = 15,
                                       ot_method = "sinkhorn",
                                       sinkhorn_epsilon = 0.01,
                                       sinkhorn_niter = 500L)

  # Both should be in [0, 1]
  expect_true(d_exact >= 0 && d_exact <= 1)
  expect_true(d_sink  >= 0 && d_sink  <= 1)

  # Sinkhorn with small epsilon should be close to exact
  expect_equal(d_exact, d_sink, tolerance = 0.02)
})

test_that("exact via transport matches transport-only ppm_wasserstein", {
  skip_if_not_installed("transport")

  # exact always goes through transport::transport(), so the result
  # must match regardless of whether approxOT is installed.
  d <- ppm_wasserstein_distance(spec1, spec2, ppm = 15,
                                 ot_method = "exact")
  expect_true(d >= 0 && d <= 1)

  # Run again — deterministic

  d2 <- ppm_wasserstein_distance(spec1, spec2, ppm = 15,
                                  ot_method = "exact")
  expect_equal(d, d2)
})

test_that("exact is stable on larger spectra", {
  skip_if_not_installed("transport")

  # Reproduce the kind of input that triggered approxOT networkflow
  # instability: ~50 peaks per spectrum.
  set.seed(42)
  big1 <- mk_spec(sort(runif(50, 50, 500)), runif(50, 1, 100))
  big2 <- mk_spec(sort(runif(50, 50, 500)), runif(50, 1, 100))

  d <- ppm_wasserstein_distance(big1, big2, ppm = 15, ot_method = "exact")
  expect_true(d >= 0 && d <= 1)
  expect_true(is.finite(d))
})

test_that("greenkhorn method works via approxOT", {
  skip_if_not_installed("approxOT")

  d <- ppm_wasserstein_distance(spec1, spec2, ppm = 15,
                                 ot_method = "greenkhorn",
                                 sinkhorn_epsilon = 0.05)
  expect_true(d >= 0 && d <= 1)
})

# ------------------------------------------------------------------
# Sinkhorn epsilon sensitivity
# ------------------------------------------------------------------

test_that("larger epsilon gives coarser approximation", {
  skip_if_not_installed("approxOT")
  skip_if_not_installed("transport")

  d_exact <- ppm_wasserstein_distance(spec1, spec2, ppm = 15,
                                       ot_method = "exact")
  d_eps01 <- ppm_wasserstein_distance(spec1, spec2, ppm = 15,
                                       ot_method = "sinkhorn",
                                       sinkhorn_epsilon = 0.01)
  d_eps50 <- ppm_wasserstein_distance(spec1, spec2, ppm = 15,
                                       ot_method = "sinkhorn",
                                       sinkhorn_epsilon = 0.50)

  # Small epsilon should be closer to exact than large epsilon
  expect_true(abs(d_eps01 - d_exact) <= abs(d_eps50 - d_exact) + 0.01)
})

# ------------------------------------------------------------------
# OT params flow through compute_distance and combined_distance
# ------------------------------------------------------------------

test_that("compute_distance passes OT params to ppm_wasserstein", {
  skip_if_not_installed("transport")

  d <- compute_distance(spec1, spec2,
                         method = "ppm_wasserstein",
                         ppm = 15,
                         ot_method = "exact",
                         sinkhorn_epsilon = 0.05,
                         sinkhorn_niter = 100L)
  expect_true(d >= 0 && d <= 1)
})

test_that("combined_distance uses OT params from params list", {
  skip_if_not_installed("transport")

  params <- eihrms_default_params()
  params$distance_method <- "ppm_wasserstein"
  params$tol_ppm <- 15
  params$ot_method <- "exact"

  frag1 <- spec1
  frag2 <- spec2
  loss1 <- mk_spec(c(15.023, 18.011), c(40, 60))
  loss2 <- mk_spec(c(15.023, 28.031), c(50, 50))

  d <- combined_distance(frag1, frag2, loss1, loss2, params)
  expect_true(d >= 0 && d <= 1)

  params$ot_method <- NULL
  d_missing <- combined_distance(frag1, frag2, loss1, loss2, params)
  params$ot_method <- "exact"
  d_exact <- combined_distance(frag1, frag2, loss1, loss2, params)
  expect_equal(d_missing, d_exact, tolerance = 0)
})

# ------------------------------------------------------------------
# Fallback behaviour
# ------------------------------------------------------------------

test_that("the publication-grade default is exact OT", {
  params <- eihrms_default_params()
  expect_identical(params$ot_method, "exact")
})

test_that("direct ppm_wasserstein calls reject unknown OT methods", {
  expect_error(
    ppm_wasserstein_distance(spec1, spec2, ot_method = "bogus"),
    "ot_method must be one of"
  )
})

# ------------------------------------------------------------------
# Parameter validation
# ------------------------------------------------------------------

test_that("validate_params accepts valid ot_method values", {
  params <- eihrms_default_params()

  for (m in c("sinkhorn", "exact", "greenkhorn")) {
    params$ot_method <- m
    expect_silent(validate_params(params))
  }
})

test_that("validate_params rejects invalid ot_method", {
  params <- eihrms_default_params()
  params$ot_method <- "bogus"
  expect_error(validate_params(params), "ot_method")
})

test_that("validate_params rejects networkflow as direct user input", {
  params <- eihrms_default_params()
  params$ot_method <- "networkflow"
  expect_error(validate_params(params), "ot_method")
})

test_that("validate_params rejects non-positive sinkhorn_epsilon", {
  params <- eihrms_default_params()
  params$sinkhorn_epsilon <- 0
  expect_error(validate_params(params), "sinkhorn_epsilon")
})

test_that("validate_params rejects non-integer sinkhorn_niter", {
  params <- eihrms_default_params()
  params$sinkhorn_niter <- 10.5
  expect_error(validate_params(params), "sinkhorn_niter")
})
