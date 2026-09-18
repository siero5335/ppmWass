test_that("detect_derivatization_from_frag flags TMS when diagnostic ions are present", {
  params <- ppmWass::eihrms_default_params()
  params$derivatization_type <- "BOTH"
  params$derivatization_auto_ppm <- 20
  params$derivatization_auto_min_int <- 0.02
  params$derivatization_auto_strong_int <- 0.08

  # Construct a simple spectrum with a strong 73.0474 ion
  peaks_raw <- rbind(
    c(73.0474, 100),
    c(147.0654, 30),
    c(91.0548, 500),
    c(43.0548, 370)
  )
  colnames(peaks_raw) <- c("mz", "intensity")

  params_detect <- params
  params_detect$use_sqrt <- FALSE
  frag <- ppmWass:::prep_peaks(peaks_raw, params_detect)

  info <- ppmWass:::detect_derivatization_from_frag(frag, params)
  expect_true(info$type %in% c("TMS", "TBDMS", "none"))
  expect_equal(info$type, "TMS")
  expect_true(is.finite(info$i73) && info$i73 > 0)
})

test_that("confusion_summary metrics do not inherit recall names", {
  cs <- ppmWass::confusion_summary(
    truth = c("none", "TMS", "TBDMS"),
    pred = c("none", "none", "TBDMS")
  )

  expect_false("none" %in% rownames(cs$metrics))
  expect_equal(names(cs$metrics), c("accuracy", "recall_none", "recall_TMS", "recall_TBDMS"))
})

test_that("roc_curve_binary custom thresholds include ROC endpoints", {
  roc <- ppmWass::roc_curve_binary(
    scores = c(0.1, 0.4, 0.8, 0.9),
    truth_positive = c(FALSE, TRUE, TRUE, FALSE),
    thresholds = c(0.5)
  )

  expect_equal(roc$threshold[1], Inf)
  expect_equal(roc$threshold[nrow(roc)], -Inf)
  expect_equal(roc$tpr[1], 0)
  expect_equal(roc$fpr[1], 0)
  expect_equal(roc$tpr[nrow(roc)], 1)
  expect_equal(roc$fpr[nrow(roc)], 1)
})
