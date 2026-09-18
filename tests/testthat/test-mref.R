test_that("estimate_Mref_HRMS recovers Mref when loss-supported", {
  # Synthetic spectrum: Mref at 200 with a supporting H2O loss peak
  mz <- c(181.9894, 200.0000, 91.0548)
  it <- c(0.06, 0.03, 0.20)
  it <- it / sum(it)

  res <- ppmWass:::estimate_Mref_HRMS(
    mz = mz,
    it = it,
    ppm = 20,
    min_rel_int = 0.005,
    typical_losses = ppmWass:::TYPICAL_LOSSES_CORE
  )

  expect_true(is.finite(res$Mref))
  expect_true(abs(res$Mref - 200.0) < 0.01)
})
