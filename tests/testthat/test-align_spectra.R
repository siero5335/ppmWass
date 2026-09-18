test_that("align_spectra optimized matches reference implementation", {
  set.seed(1)

  for (rep in 1:20) {
    nA <- sample(c(1, 5, 20, 80), 1)
    nB <- sample(c(1, 5, 20, 80), 1)

    mzA <- sort(stats::runif(nA, 35, 300))
    mzB <- sort(stats::runif(nB, 35, 300))
    intA <- stats::runif(nA)
    intB <- stats::runif(nB)

    specA <- cbind(mzA, intA)
    specB <- cbind(mzB, intB)

    ref <- ppmWass:::align_spectra_reference(specA, specB, ppm = 20)
    opt <- ppmWass:::align_spectra(specA, specB, ppm = 20)

    expect_equal(opt$mz, ref$mz, tolerance = 1e-12)
    expect_equal(opt$p, ref$p, tolerance = 1e-12)
    expect_equal(opt$q, ref$q, tolerance = 1e-12)
    expect_true(abs(sum(opt$p) - 1) < 1e-12)
    expect_true(abs(sum(opt$q) - 1) < 1e-12)
  }
})
