test_that("Derivatization raw-loss downweighting affects truncation (loss_max_peaks)", {
  params <- ppmWass::eihrms_default_params()
  params$tol_ppm <- 20
  params$loss_min <- 5
  params$loss_max <- 200
  params$loss_top_peaks <- 10
  params$loss_max_peaks <- 1

  # Construct peaks such that Mref=max(mz)=300, and one anchored loss matches TMS nominal 90 (exact 90.0501).
  # Without downweighting, the 90.0501 loss dominates and will be the only kept loss peak.
  # With downweighting, another loss (60) should survive truncation.
  peaks <- rbind(
    c(300.0000, 0.2),
    c(209.9499, 1.0),  # 300 - 209.9499 = 90.0501
    c(240.0000, 0.9),  # loss 60
    c(250.0000, 0.8)   # loss 50
  )
  colnames(peaks) <- c("mz", "intensity")

  params$derivatization_raw_loss_weight <- 1
  out0 <- ppmWass:::build_loss_peaks(peaks, params, return_info = TRUE, deriv_type = "TMS")
  mz0 <- out0$lossA_peaks[1, 1]

  params$derivatization_raw_loss_weight <- 0.01
  out1 <- ppmWass:::build_loss_peaks(peaks, params, return_info = TRUE, deriv_type = "TMS")
  mz1 <- out1$lossA_peaks[1, 1]

  expect_true(abs(mz0 - 90.05) < 0.01)
  expect_true(abs(mz1 - 60.0) < 1e-6)
})
