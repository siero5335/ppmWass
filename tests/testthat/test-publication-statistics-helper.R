stats_helper_path <- system.file(
  "scripts", "lib", "publication_retrieval_statistics.R",
  package = "ppmWass"
)
if (!nzchar(stats_helper_path) || !file.exists(stats_helper_path)) {
  stats_helper_path <- testthat::test_path(
    "..", "..", "inst", "scripts", "lib",
    "publication_retrieval_statistics.R"
  )
}
source(stats_helper_path, local = TRUE)

test_that("fractional MAP is the exact expectation within tie blocks", {
  result <- prefix_ranking_metrics_fractional(
    dists = c(0, 0, 0),
    relevant = c(TRUE, TRUE, FALSE),
    ks = c(1L, 2L, 3L)
  )

  # Average over the three equally likely placements of the irrelevant item.
  expected_ap <- mean(c(1, (1 + 2 / 3) / 2, (1 / 2 + 2 / 3) / 2))
  expect_equal(result$average_precision_fractional, expected_ap, tolerance = 1e-12)
  expect_equal(result$precision_at_1_fractional, 2 / 3, tolerance = 1e-12)
  expect_equal(result$precision_at_2_fractional, 2 / 3, tolerance = 1e-12)
  expect_equal(result$precision_at_3_fractional, 2 / 3, tolerance = 1e-12)
})

test_that("ordered-pair discrimination retains both directional orientations", {
  ids <- c("a", "b", "c")
  distance <- matrix(
    c(
      0, 0.1, 0.9,
      0.8, 0, 0.2,
      0.3, 0.7, 0
    ),
    nrow = 3L, byrow = TRUE,
    dimnames = list(ids, ids)
  )
  inchikey <- c(
    a = "AAAAAAAAAAAAAA-ONE",
    b = "AAAAAAAAAAAAAA-TWO",
    c = "CCCCCCCCCCCCCC-ONE"
  )

  result <- ordered_pair_discrimination_metrics(
    distance, inchikey, thresholds = c(0.5)
  )
  expect_equal(result$auc$n_ordered_pairs, 6L)
  expect_equal(result$auc$n_positive, 2L)
  expect_equal(result$auc$n_negative, 4L)
  expect_identical(
    result$auc$pair_orientation,
    "ordered_query_to_library"
  )
  expect_equal(nrow(result$fdr), 1L)
})

test_that("rank AUC remains finite above the integer multiplication limit", {
  n_positive <- 50000L
  n_negative <- 50000L
  distance <- c(rep(0, n_positive), rep(1, n_negative))
  is_positive <- c(rep(TRUE, n_positive), rep(FALSE, n_negative))

  expect_gt(as.double(n_positive) * as.double(n_negative),
            .Machine$integer.max)
  expect_equal(binary_rank_auc(distance, is_positive), 1, tolerance = 1e-12)
})

test_that("ordered-pair metrics reject matrices with ambiguous ID orientation", {
  distance <- diag(2)
  rownames(distance) <- c("q1", "q2")
  colnames(distance) <- c("q2", "q1")
  expect_error(
    ordered_pair_discrimination_metrics(
      distance,
      c(q1 = "AAAAAAAAAAAAAA-X", q2 = "AAAAAAAAAAAAAA-Y")
    ),
    "identical row/column IDs"
  )
})

test_that("ordered-pair FDR distinguishes no valid pairs from no calls", {
  ids <- c("a", "b")
  distance <- matrix(c(0, 1, 1, 0), nrow = 2L, byrow = TRUE,
                     dimnames = list(ids, ids))
  result <- ordered_pair_discrimination_metrics(
    distance,
    c(a = "AAAAAAAAAAAAAA-X", b = "BBBBBBBBBBBBBB-X"),
    thresholds = 1.1
  )
  expect_equal(result$fdr$n_matches, 0L)
  expect_equal(result$fdr$fdr, 0)

  singleton <- matrix(0, nrow = 1L, dimnames = list("a", "a"))
  empty_result <- ordered_pair_discrimination_metrics(
    singleton, c(a = "AAAAAAAAAAAAAA-X"), thresholds = 0.5
  )
  expect_equal(empty_result$auc$n_ordered_pairs, 0L)
  expect_true(is.na(empty_result$auc$auc))
  expect_equal(empty_result$fdr$n_matches, 0L)
  expect_true(is.na(empty_result$fdr$fdr))
})

test_that("prefix-only query sets retain an explicit empty full-key schema", {
  ids <- c("a", "b", "c")
  distance <- matrix(
    c(0, 0.1, 0.8, 0.2, 0, 0.7, 0.9, 0.6, 0),
    nrow = 3L, byrow = TRUE, dimnames = list(ids, ids)
  )
  inchikey <- c(
    a = "AAAAAAAAAAAAAA-ONE",
    b = "AAAAAAAAAAAAAA-TWO",
    c = "CCCCCCCCCCCCCC-ONE"
  )

  result <- per_query_metrics_tie_aware(distance, inchikey)
  expect_setequal(result$query_id, c("a", "b"))
  expect_true(all(is.na(result$top1_fractional)))
  expect_true(all(is.finite(result$p_at_1_fractional)))
  expect_type(result$full_has_relevant_tie, "logical")
  expect_type(result$full_n_better, "integer")
  expect_silent(
    summarize_tie_aware_metrics(
      result, "prefix-only", "cosine", query_boot_R = 5L,
      cluster_boot_R = 5L
    )
  )
})

test_that("full-key-only query sets retain an explicit empty prefix schema", {
  ids <- c("a", "b", "c")
  distance <- matrix(
    c(0, 0.1, 0.8, 0.2, 0, 0.7, 0.9, 0.6, 0),
    nrow = 3L, byrow = TRUE, dimnames = list(ids, ids)
  )
  # Short identifiers are valid for the full-key grouping used by this helper,
  # but intentionally do not yield a 14-character connectivity prefix.
  inchikey <- c(a = "KEY-A", b = "KEY-A", c = "KEY-C")

  result <- per_query_metrics_tie_aware(distance, inchikey)
  expect_setequal(result$query_id, c("a", "b"))
  expect_true(all(is.finite(result$top1_fractional)))
  expect_true(all(is.na(result$p_at_1_fractional)))
  expect_type(result$prefix_has_relevant_tie, "logical")
  expect_type(result$prefix_n_better, "integer")
})

test_that("no duplicate query set returns a stable zero-row schema", {
  ids <- c("a", "b")
  distance <- matrix(c(0, 1, 1, 0), nrow = 2L, byrow = TRUE,
                     dimnames = list(ids, ids))
  result <- per_query_metrics_tie_aware(
    distance,
    c(a = "AAAAAAAAAAAAAA-X", b = "BBBBBBBBBBBBBB-X")
  )
  expect_equal(nrow(result), 0L)
  expect_true(all(c(
    "query_id", "top1_fractional", "p_at_1_fractional",
    "average_precision_fractional"
  ) %in% names(result)))
})
