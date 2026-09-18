test_that("evaluate_homolog_detection reports inferred label usage", {
  dist_mat <- matrix(
    c(
      0.0, 0.1, 0.8,
      0.1, 0.0, 0.7,
      0.8, 0.7, 0.0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  df_spec <- data.frame(
    id = c("a", "b", "c"),
    compound_class = c("alkane", "alkane", "aromatic"),
    compound_class_source = c("inferred", "inferred", "inferred"),
    stringsAsFactors = FALSE
  )
  ri <- c(a = 100, b = 220, c = 450)

  res <- ppmWass::evaluate_homolog_detection(dist_mat, df_spec, ri, ri_diff_threshold = 50, sim_threshold = 0.5)

  expect_true(isTRUE(res$uses_inferred_classes))
  expect_match(res$ground_truth_note, "internal consistency")
  expect_equal(res$n_homolog_pairs, 1)
  expect_equal(res$n_true_homologs, 1)
})

test_that("benchmark similarity helpers separate ranking from calibrated thresholding", {
  dist_mat <- matrix(
    c(
      0.0, 0.4, 1.7,
      0.4, 0.0, -0.2,
      1.7, -0.2, 0.0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  ik <- c("AAAAAAAAAAAAAA", "AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB")

  raw_pairs <- ppmWass:::extract_similarity_pairs(dist_mat, ik)
  expect_equal(sort(raw_pairs$similarity), sort(c(0.6, -0.7, 1.2)))

  calibrated_pairs <- ppmWass:::extract_similarity_pairs(
    dist_mat,
    ik,
    similarity_transform = "minmax"
  )
  expect_true(all(calibrated_pairs$similarity >= 0 & calibrated_pairs$similarity <= 1))
  expect_equal(sort(calibrated_pairs$similarity), sort(c(1, 0, 1 - (0.4 + 0.2) / 1.9)))

  fdr <- ppmWass::evaluate_fdr(dist_mat, ik, thresholds = 0.5)
  expect_equal(fdr$n_matches, 2)
  expect_equal(fdr$n_true_pos, 1)
  expect_equal(fdr$n_false_pos, 1)

  df_spec <- data.frame(
    id = c("a", "b", "c"),
    compound_class = c("alkane", "alkane", "alkane"),
    stringsAsFactors = FALSE
  )
  ri <- c(a = 100, b = 220, c = 450)
  res <- ppmWass::evaluate_homolog_detection(dist_mat, df_spec, ri, sim_threshold = 0.5)
  expect_equal(res$n_similar_pairs, 2)
})

test_that("AUC uses raw monotone distance ranks by default", {
  dist_mat <- matrix(
    c(
      0.0, 1.5, 2.0,
      1.5, 0.0, 3.0,
      2.0, 3.0, 0.0
    ),
    nrow = 3,
    byrow = TRUE
  )
  ik <- c("AAAAAAAAAAAAAA", "AAAAAAAAAAAAAA", "BBBBBBBBBBBBBB")

  auc <- ppmWass::evaluate_auc_roc(dist_mat, ik)
  expect_equal(auc$auc, 1)
})

test_that("create_homolog_ground_truth annotates connectivity limitation", {
  df_spec <- data.frame(
    id = c("x1", "x2", "x3", "y1"),
    inchikey = c("ABCDEFGHIJKLMN-AAAAAA", "ABCDEFGHIJKLMN-BBBBBB", "ABCDEFGHIJKLMN-CCCCCC", "ZZZZZZZZZZZZZZ-DDDDDD"),
    stringsAsFactors = FALSE
  )

  gt <- ppmWass::create_homolog_ground_truth(df_spec, min_series_length = 3)

  expect_equal(nrow(gt), 3)
  expect_match(attr(gt, "ground_truth_note"), "not curated homolog-series labels")
})
