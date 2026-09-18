publication_stats_helper_test_path <- function() {
  path <- system.file(
    "scripts", "lib", "publication_retrieval_statistics.R",
    package = "ppmWass"
  )
  if (!nzchar(path) || !file.exists(path)) {
    path <- testthat::test_path(
      "..", "..", "inst", "scripts", "lib",
      "publication_retrieval_statistics.R"
    )
  }
  path
}

test_that("MSP peak parsing preserves input precision and decimal intensity", {
  parse <- getFromNamespace("parse_msp_peaks", "ppmWass")
  out <- parse(c(
    "100.123456789 0.125",
    "2e2:3.75e-1",
    "300 4; 400 5",
    "145.10437 2983677 \"Theoretical m/z 145.104324, Mass diff 0\""
  ))

  expect_identical(
    out,
    paste(
      "100.123456789:0.125",
      "2e2:3.75e-1",
      "300:4",
      "400:5",
      "145.10437:2983677"
    )
  )
  expect_match(out, "0.125", fixed = TRUE)
})

test_that("MSP peak parsing reports rather than silently truncates odd tokens", {
  parse <- getFromNamespace("parse_msp_peaks", "ppmWass")
  expect_warning(parse("100 1 200"), "Unpaired numeric token")
  expect_warning(parse("200"), "Unpaired numeric token")
  expect_warning(parse("100 NA"), "missing value")
  expect_warning(parse("NA 1"), "missing value")
  expect_warning(parse("100 0"), "Invalid non-positive")
})

test_that("MSP parser supports mixed encodings and a final record without blank line", {
  path <- tempfile(fileext = ".msp")
  writeChar(paste(c(
    "Name: mixed",
    "InChIKey: AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    "Num Peaks: 5",
    "100.1:1.25 110.2:2.5",
    "120.3 3.75; 130.4 4.5",
    "1.405e2 5e-1"
  ), collapse = "\n"), path, eos = NULL)

  x <- read_msp(path, progress = 0)
  expect_equal(nrow(x), 1)
  expect_identical(
    x$spectrum,
    "100.1:1.25 110.2:2.5 120.3:3.75 130.4:4.5 1.405e2:5e-1"
  )
})

test_that("blank MSP records are ignored without creating phantom spectra", {
  path <- tempfile(fileext = ".msp")
  writeLines(c(
    "", "",
    "Name: first", "Num Peaks: 1", "100 1",
    "", "", "",
    "Name: second", "Num Peaks: 1", "200 2",
    "", ""
  ), path)

  x <- read_msp(path, progress = 0)
  expect_identical(x$name, c("first", "second"))
  expect_identical(x$spectrum, c("100:1", "200:2"))
  expect_equal(attr(x, "msp_parser_diagnostics")$n_records, 2L)
  expect_identical(
    getFromNamespace("parse_msp_peaks", "ppmWass")(character()),
    ""
  )
})

test_that("duplicate metadata and declared peak mismatches are diagnosed", {
  path <- tempfile(fileext = ".msp")
  writeLines(c(
    "Name: duplicate metadata",
    "RI: 1000",
    "RI: 1001",
    "Num Peaks: 3",
    "100 1",
    "Num Peaks: 99",
    "200 2"
  ), path)
  x <- read_msp(path, progress = 0)
  issues <- attr(x, "msp_parser_diagnostics")$record_issues
  expect_equal(x$ri, 1001)
  expect_equal(x$num_peaks, 3)
  expect_equal(x$spectrum, "100:1 200:2")
  expect_true(any(
    issues$issue_type == "duplicate_metadata_field" & issues$field == "ri"
  ))
  expect_true(any(
    issues$issue_type == "duplicate_metadata_field" &
      issues$field == "num_peaks"
  ))
  expect_true(any(
    issues$issue_type == "declared_peak_count_mismatch" &
      issues$declared_peak_count == 3L & issues$parsed_peak_count == 2L
  ))
})

test_that("embedded MSP record starts are repaired without losing the peak", {
  path <- tempfile(fileext = ".msp")
  writeLines(c(
    "Name: first",
    "InChIKey: AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    "Num Peaks: 1",
    "429.08865 1533679NAME: second",
    "InChIKey: CCCCCCCCCCCCCC-DDDDDDDDDD-E",
    "Num Peaks: 1",
    "100.123456789 0.125"
  ), path)

  x <- read_msp(path, progress = 0)
  diag <- attr(x, "msp_parser_diagnostics")
  expect_equal(nrow(x), 2)
  expect_identical(x$name, c("first", "second"))
  expect_identical(x$spectrum, c("429.08865:1533679", "100.123456789:0.125"))
  expect_equal(diag$n_embedded_record_repairs, 1)
})

test_that("directional composite matrices are not mirrored", {
  A <- cbind(mz = c(50, 75, 100), intensity = c(0.6, 0.3, 0.1))
  B <- cbind(mz = c(50.0001, 100.0001), intensity = c(0.7, 0.3))
  frag <- list(A = A, B = B)
  derived <- list(A = A, B = B)

  p <- eihrms_default_params()
  p$distance_method <- "composite"
  p$w_frag <- 1
  p$w_loss <- 0
  p$use_parallel <- FALSE
  p <- validate_params(p)

  direct_ab <- getFromNamespace("composite_distance", "ppmWass")(
    A, B, ppm = p$tol_ppm,
    mass_power = p$mass_power,
    intensity_power = p$intensity_power
  )
  direct_ba <- getFromNamespace("composite_distance", "ppmWass")(
    B, A, ppm = p$tol_ppm,
    mass_power = p$mass_power,
    intensity_power = p$intensity_power
  )
  dm <- compute_distance_matrix(frag, derived, p, progress = FALSE)

  expect_false(isTRUE(all.equal(direct_ab, direct_ba, tolerance = 1e-14)))
  expect_equal(dm["A", "B"], direct_ab)
  expect_equal(dm["B", "A"], direct_ba)
})

test_that("named Mref confidence is aligned by spectrum ID for directional matrices", {
  one_peak <- function(mz) {
    matrix(c(mz, 1), ncol = 2,
           dimnames = list(NULL, c("mz", "intensity")))
  }
  ids <- c("A", "B", "C")
  frag <- setNames(rep(list(one_peak(100)), 3), ids)
  anchor <- setNames(lapply(c(50, 60, 70), one_peak), ids)
  pair <- setNames(rep(list(one_peak(20)), 3), ids)

  p <- eihrms_default_params()
  p$distance_method <- "composite"
  p$use_split_loss <- TRUE
  p$use_mref_confidence <- TRUE
  p$w_frag <- 0
  p$w_loss <- 1
  p$use_parallel <- FALSE
  p <- validate_params(p)

  confidence <- c(A = 0.04, B = 0.36, C = 1)
  ordered <- compute_distance_matrix(
    frag, anchor, p, progress = FALSE,
    mref_conf = confidence,
    loss_anchor_list = anchor,
    loss_pair_list = pair
  )
  permuted <- compute_distance_matrix(
    frag, anchor, p, progress = FALSE,
    mref_conf = confidence[c("C", "A", "B")],
    loss_anchor_list = anchor,
    loss_pair_list = pair
  )
  positional_control <- compute_distance_matrix(
    frag, anchor, p, progress = FALSE,
    mref_conf = unname(confidence[c("C", "A", "B")]),
    loss_anchor_list = anchor,
    loss_pair_list = pair
  )

  expect_equal(permuted, ordered, tolerance = 1e-14)
  expect_gt(max(abs(positional_control - ordered)), 0.1)
})

test_that("directional component matrices preserve the exported component API", {
  A <- cbind(mz = c(50, 75, 100), intensity = c(0.6, 0.3, 0.1))
  B <- cbind(mz = c(50.0001, 100.0001), intensity = c(0.7, 0.3))
  frag <- list(A = A, B = B)
  derived <- list(A = A, B = B)

  p <- eihrms_default_params()
  p$distance_method <- "composite"
  p$w_frag <- 1
  p$w_loss <- 0
  p$use_parallel <- FALSE
  p <- validate_params(p)

  base <- compute_distance_matrix(frag, derived, p, progress = FALSE)
  with_components <- compute_distance_matrices(
    frag, derived, p, progress = FALSE
  )
  expect_equal(with_components$dist, base, tolerance = 1e-14)
  expect_equal(with_components$components$d_frag["A", "B"], base["A", "B"])
  expect_equal(with_components$components$d_frag["B", "A"], base["B", "A"])
  expect_false(isTRUE(all.equal(
    with_components$components$d_frag["A", "B"],
    with_components$components$d_frag["B", "A"],
    tolerance = 1e-14
  )))

  p$return_distance_components <- TRUE
  similarity <- compute_similarity_matrices(
    frag, derived, ri = c(A = NA_real_, B = NA_real_), params = p
  )
  expect_equal(similarity$dist_raw, base, tolerance = 1e-14)
  expect_true(is.list(similarity$dist_components))
})

test_that("directional square matrices retain a zero diagonal with empty channels", {
  one_peak <- function(mz) matrix(c(mz, 1), ncol = 2)
  empty <- matrix(numeric(0), nrow = 0, ncol = 2)
  frag <- list(A = one_peak(100), B = one_peak(110))
  derived <- list(A = empty, B = empty)

  p <- eihrms_default_params()
  p$distance_method <- "composite"
  p$use_parallel <- FALSE
  p <- validate_params(p)

  dm <- compute_distance_matrix(frag, derived, p, progress = FALSE)
  expect_equal(unname(diag(dm)), c(0, 0), tolerance = 0)
})

test_that("publication method symmetry classification matches actual distances", {
  classify <- getFromNamespace("distance_method_is_symmetric", "ppmWass")
  distance <- getFromNamespace("compute_distance", "ppmWass")
  symmetric_methods <- c(
    "ppm_wasserstein", "wasserstein", "cosine", "weighted_cosine",
    "hellinger", "entropy", "entropy_weighted", "entropy_unweighted"
  )
  expect_false(classify("composite"))
  expect_true(all(vapply(symmetric_methods, classify, logical(1))))
  expect_true(classify("ppm_wasserstein", "exact"))
  expect_false(classify("ppm_wasserstein", "sinkhorn"))
  expect_false(classify("ppm_wasserstein", "greenkhorn"))

  spectra <- list(
    A = cbind(mz = c(50, 75, 100), intensity = c(0.6, 0.3, 0.1)),
    B = cbind(mz = c(50.0001, 100.0001), intensity = c(0.7, 0.3)),
    C = cbind(mz = c(50.0002, 75.0002, 125), intensity = c(0.4, 0.4, 0.2))
  )
  pairs <- list(c("A", "B"), c("A", "C"), c("B", "C"))
  evaluate <- function(method, pair) {
    suppressWarnings(c(
      ab = distance(
        spectra[[pair[[1]]]], spectra[[pair[[2]]]], method = method,
        ppm = 15, ot_method = "exact", sinkhorn_niter = 100L
      ),
      ba = distance(
        spectra[[pair[[2]]]], spectra[[pair[[1]]]], method = method,
        ppm = 15, ot_method = "exact", sinkhorn_niter = 100L
      )
    ))
  }

  for (method in symmetric_methods) {
    for (pair in pairs) {
      values <- evaluate(method, pair)
      expect_true(all(is.finite(values)), info = paste(method, pair, collapse = ":"))
      expect_equal(
        unname(values[["ab"]]), unname(values[["ba"]]),
        tolerance = 1e-6,
        info = paste(method, paste(pair, collapse = "->"))
      )
    }
  }
  composite_differences <- vapply(pairs, function(pair) {
    values <- evaluate("composite", pair)
    abs(values[["ab"]] - values[["ba"]])
  }, numeric(1))
  expect_true(any(composite_differences > 1e-12))
})

test_that("exact ppmW is order invariant and square/rectangular paths agree", {
  skip_if_not_installed("transport")
  ppmw <- getFromNamespace("ppm_wasserstein_distance", "ppmWass")
  search_matrix <- getFromNamespace("compute_distance_matrix_search", "ppmWass")
  A <- cbind(
    mz = c(91.05425, 105.069885, 119.08556, 161.13245),
    intensity = c(0.0640527962151771, 0.180230435344229,
                  0.12268667663643, 0.0738341491041677)
  )
  B <- cbind(
    mz = c(81.069855, 91.05417, 105.06977, 119.08542, 161.13237),
    intensity = c(0.0405897497350223, 0.0732282724366475,
                  0.164381906032416, 0.176824758042847,
                  0.0993508930752914)
  )
  d_ab <- ppmw(A, B, ppm = 15, ot_method = "exact")
  d_ba <- ppmw(B, A, ppm = 15, ot_method = "exact")
  expect_equal(d_ab, d_ba, tolerance = 1e-14)

  empty <- matrix(numeric(), nrow = 0L, ncol = 2L)
  frag <- list(A = A, B = B)
  loss <- list(A = empty, B = empty)
  p <- eihrms_default_params()
  p$distance_method <- "ppm_wasserstein"
  p$tol_ppm <- 15
  p$w_frag <- 1
  p$w_loss <- 0
  p$ot_method <- "exact"
  p$use_parallel <- FALSE
  square <- compute_distance_matrix(frag, loss, p, progress = FALSE)
  rectangular <- search_matrix(
    query_frag_list = frag, query_loss_list = loss,
    lib_frag_list = frag, lib_loss_list = loss,
    params = p, progress = FALSE
  )
  expect_equal(square["A", "B"], d_ab, tolerance = 1e-14)
  expect_equal(square["B", "A"], rectangular["B", "A"], tolerance = 1e-14)
  expect_equal(rectangular["A", "B"], rectangular["B", "A"], tolerance = 1e-14)
})

test_that("finite-iteration Sinkhorn is never mirrored as a symmetric matrix", {
  skip_if_not_installed("approxOT")
  ppmw <- getFromNamespace("ppm_wasserstein_distance", "ppmWass")
  A <- cbind(
    mz = c(91.05425, 105.069885, 119.08556, 161.13245),
    intensity = c(0.0640527962151771, 0.180230435344229,
                  0.12268667663643, 0.0738341491041677)
  )
  B <- cbind(
    mz = c(81.069855, 91.05417, 105.06977, 119.08542, 161.13237),
    intensity = c(0.0405897497350223, 0.0732282724366475,
                  0.164381906032416, 0.176824758042847,
                  0.0993508930752914)
  )
  d_ab <- ppmw(
    A, B, ppm = 15, ot_method = "sinkhorn",
    sinkhorn_epsilon = 0.05, sinkhorn_niter = 100L
  )
  d_ba <- ppmw(
    B, A, ppm = 15, ot_method = "sinkhorn",
    sinkhorn_epsilon = 0.05, sinkhorn_niter = 100L
  )
  expect_gt(abs(d_ab - d_ba), 1e-6)

  empty <- matrix(numeric(), nrow = 0L, ncol = 2L)
  frag <- list(A = A, B = B)
  loss <- list(A = empty, B = empty)
  p <- eihrms_default_params()
  p$distance_method <- "ppm_wasserstein"
  p$tol_ppm <- 15
  p$w_frag <- 1
  p$w_loss <- 0
  p$ot_method <- "sinkhorn"
  p$sinkhorn_epsilon <- 0.05
  p$sinkhorn_niter <- 100L
  p$use_parallel <- FALSE
  square <- compute_distance_matrix(frag, loss, p, progress = FALSE)
  expect_equal(square["A", "B"], d_ab, tolerance = 1e-14)
  expect_equal(square["B", "A"], d_ba, tolerance = 1e-14)
  expect_gt(abs(square["A", "B"] - square["B", "A"]), 1e-6)
})

test_that("ppmWass cost scale saturates and constant-cost fallback is finite", {
  ppmw <- getFromNamespace("ppm_wasserstein_distance", "ppmWass")
  A <- cbind(mz = c(50, 60), intensity = c(0.5, 0.5))
  B <- cbind(mz = c(100, 110), intensity = c(0.5, 0.5))
  expect_identical(
    ppmw(A, B, ppm = 15, transition_mult = 3,
         ot_method = "sinkhorn", sinkhorn_epsilon = 0.05,
         sinkhorn_niter = 100L),
    1
  )

  delta <- 15e-6
  center <- 100
  a <- center * (1 - delta / 2)
  b <- center * (1 + delta / 2)
  one_a <- cbind(mz = a, intensity = 1)
  one_b <- cbind(mz = b, intensity = 1)
  expect_equal(
    ppmw(one_a, one_b, ppm = 15, transition_mult = 3,
         ot_method = "sinkhorn", sinkhorn_epsilon = 0.05,
         sinkhorn_niter = 100L),
    1 / 3,
    tolerance = 1e-10
  )
})

test_that("OT marginal checks fix from/to orientation for both backends", {
  ppmw <- getFromNamespace("ppm_wasserstein_distance", "ppmWass")
  tolerance <- getFromNamespace("PPMWASS_OT_MARGINAL_TOLERANCE", "ppmWass")
  A <- cbind(mz = c(100, 150), intensity = c(0.6, 0.4))
  B <- cbind(mz = c(100.001, 140, 200), intensity = c(0.2, 0.3, 0.5))
  valid_plan <- list(
    from = c(1L, 1L, 1L, 2L, 2L),
    to = c(1L, 2L, 3L, 2L, 3L),
    mass = c(0.2, 0.2, 0.2, 0.1, 0.3)
  )
  approximate <- function(...) valid_plan
  exact <- function(...) valid_plan

  approx_log <- new.env(parent = emptyenv())
  exact_log <- new.env(parent = emptyenv())
  d_approx <- ppmw(
    A, B, ot_method = "sinkhorn", .approx_solver = approximate,
    .exact_solver = exact, .diagnostics = approx_log
  )
  d_exact <- ppmw(
    A, B, ot_method = "exact", .approx_solver = approximate,
    .exact_solver = exact, .diagnostics = exact_log
  )
  expect_true(is.finite(d_approx))
  expect_true(is.finite(d_exact))
  expect_identical(approx_log$records[[1]]$selected_path, "approx_primary")
  expect_identical(exact_log$records[[1]]$selected_path, "exact_primary")
  expect_lte(
    approx_log$records[[1]]$attempts$row_residual_linf[[1]], tolerance
  )
  expect_lte(
    exact_log$records[[1]]$attempts$col_residual_linf[[1]], tolerance
  )

  swapped <- valid_plan
  swapped$from <- valid_plan$to
  swapped$to <- valid_plan$from
  assessed <- getFromNamespace("assess_ot_plan", "ppmWass")(
    swapped, matrix(0, 2, 3), c(0.6, 0.4), c(0.2, 0.3, 0.5),
    tolerance = tolerance
  )
  expect_false(assessed$accepted)
  expect_identical(assessed$validation_error, "invalid_plan_index")
})

test_that("installed approxOT and transport plans use the assessed index orientation", {
  skip_if_not_installed("approxOT")
  skip_if_not_installed("transport")
  ppmw <- getFromNamespace("ppm_wasserstein_distance", "ppmWass")
  A <- cbind(mz = c(100, 150), intensity = c(0.6, 0.4))
  B <- cbind(mz = c(100.001, 140, 200), intensity = c(0.2, 0.3, 0.5))
  approx_log <- new.env(parent = emptyenv())
  exact_log <- new.env(parent = emptyenv())

  ppmw(A, B, ot_method = "sinkhorn", .diagnostics = approx_log)
  ppmw(A, B, ot_method = "exact", .diagnostics = exact_log)
  approx_attempt <- approx_log$records[[1]]$attempts[1, ]
  exact_attempt <- exact_log$records[[1]]$attempts[1, ]
  expect_true(approx_attempt$accepted)
  expect_true(exact_attempt$accepted)
  expect_equal(approx_attempt$invalid_index, 0L)
  expect_equal(exact_attempt$invalid_index, 0L)
  expect_lte(approx_attempt$row_residual_linf, 1e-8)
  expect_lte(exact_attempt$col_residual_linf, 1e-8)
})

test_that("finite OT plans with bad marginals use logged exact fallback", {
  ppmw <- getFromNamespace("ppm_wasserstein_distance", "ppmWass")
  A <- cbind(mz = c(100, 150), intensity = c(0.6, 0.4))
  B <- cbind(mz = c(100.001, 140, 200), intensity = c(0.2, 0.3, 0.5))
  valid_plan <- list(
    from = c(1L, 1L, 1L, 2L, 2L),
    to = c(1L, 2L, 3L, 2L, 3L),
    mass = c(0.2, 0.2, 0.2, 0.1, 0.3)
  )
  bad_plan <- valid_plan
  bad_plan$mass <- c(0.2, 0.2, 0.2, 0.1, 0.1)
  collector <- new.env(parent = emptyenv())

  value <- ppmw(
    A, B, ot_method = "sinkhorn",
    .approx_solver = function(...) bad_plan,
    .exact_solver = function(...) valid_plan,
    .diagnostics = collector,
    .context = list(query_id = "query-A", library_id = "library-B")
  )
  record <- collector$records[[1]]
  expect_true(is.finite(value))
  expect_identical(record$selected_path, "exact_fallback")
  expect_true(record$fallback_used)
  expect_identical(record$context$query_id, "query-A")
  expect_identical(record$context$library_id, "library-B")
  expect_equal(nrow(record$attempts), 4L)
  expect_true(is.finite(record$attempts$total_cost[[1]]))
  expect_false(record$attempts$accepted[[1]])
  expect_identical(
    record$attempts$validation_error[[1]],
    "marginal_residual_exceeds_tolerance"
  )
  expect_true(record$attempts$accepted[[4]])
})

test_that("OT residual threshold and original retry error are auditable", {
  ppmw <- getFromNamespace("ppm_wasserstein_distance", "ppmWass")
  A <- cbind(mz = c(100, 150), intensity = c(0.6, 0.4))
  B <- cbind(mz = c(100.001, 140, 200), intensity = c(0.2, 0.3, 0.5))
  valid_plan <- list(
    from = c(1L, 1L, 1L, 2L, 2L),
    to = c(1L, 2L, 3L, 2L, 3L),
    mass = c(0.2, 0.2, 0.2, 0.1, 0.3)
  )
  near_plan <- valid_plan
  near_plan$mass[1:2] <- near_plan$mass[1:2] + c(5e-8, -5e-8)

  loose <- new.env(parent = emptyenv())
  strict <- new.env(parent = emptyenv())
  ppmw(
    A, B, ot_method = "sinkhorn", .approx_solver = function(...) near_plan,
    .exact_solver = function(...) valid_plan,
    .marginal_tolerance = 1e-6, .diagnostics = loose
  )
  ppmw(
    A, B, ot_method = "sinkhorn", .approx_solver = function(...) near_plan,
    .exact_solver = function(...) valid_plan,
    .marginal_tolerance = 1e-8, .diagnostics = strict
  )
  expect_identical(loose$records[[1]]$selected_path, "approx_primary")
  expect_identical(strict$records[[1]]$selected_path, "exact_fallback")
  expect_identical(strict$records[[1]]$marginal_tolerance, 1e-8)

  calls <- 0L
  retry_log <- new.env(parent = emptyenv())
  retry_solver <- function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop("forced primary failure")
    valid_plan
  }
  ppmw(
    A, B, ot_method = "sinkhorn", .approx_solver = retry_solver,
    .exact_solver = function(...) valid_plan,
    .diagnostics = retry_log
  )
  retry_record <- retry_log$records[[1]]
  expect_identical(retry_record$selected_path, "approx_retry_500")
  expect_match(retry_record$original_solver_error, "forced primary failure")
  expect_identical(retry_record$attempts$solver_error[[1]],
                   "forced primary failure")
})

test_that("requested-only approximate OT executes exactly the requested call", {
  requested <- getFromNamespace(
    "ppm_wasserstein_requested_approx", "ppmWass"
  )
  A <- cbind(mz = c(100, 150), intensity = c(0.6, 0.4))
  B <- cbind(mz = c(100.001, 140, 200), intensity = c(0.2, 0.3, 0.5))
  near_plan <- list(
    from = c(1L, 1L, 1L, 2L, 2L),
    to = c(1L, 2L, 3L, 2L, 3L),
    mass = c(0.2 + 5e-8, 0.2 - 5e-8, 0.2, 0.1, 0.3)
  )
  calls <- list()
  solver <- function(...) {
    calls[[length(calls) + 1L]] <<- list(...)
    near_plan
  }
  collector <- new.env(parent = emptyenv())

  value <- requested(
    A, B, ot_method = "sinkhorn", sinkhorn_epsilon = 0.037,
    sinkhorn_niter = 37L, .marginal_tolerance = 1e-8,
    .approx_solver = solver, .diagnostics = collector
  )
  record <- collector$records[[1L]]

  expect_true(is.finite(value))
  expect_length(calls, 1L)
  expect_identical(calls[[1L]]$method, "sinkhorn")
  expect_identical(calls[[1L]]$niter, 37L)
  expect_equal(calls[[1L]]$epsilon, 0.037, tolerance = 0)
  expect_identical(record$selected_path, "approx_requested_only")
  expect_false(record$fallback_used)
  expect_true(record$requested_only)
  expect_identical(record$requested_niter, 37L)
  expect_equal(record$requested_epsilon, 0.037, tolerance = 0)
  expect_equal(nrow(record$attempts), 1L)
  expect_identical(record$attempts$step, "approx_requested_only")
  expect_false(record$attempts$accepted)
  expect_identical(
    record$attempts$validation_error,
    "marginal_residual_exceeds_tolerance"
  )
})

test_that("requested-only approximate OT never retries or falls back on error", {
  requested <- getFromNamespace(
    "ppm_wasserstein_requested_approx", "ppmWass"
  )
  A <- cbind(mz = c(100, 150), intensity = c(0.6, 0.4))
  B <- cbind(mz = c(100.001, 140, 200), intensity = c(0.2, 0.3, 0.5))
  calls <- 0L
  collector <- new.env(parent = emptyenv())
  value <- requested(
    A, B, ot_method = "sinkhorn", sinkhorn_epsilon = 0.05,
    sinkhorn_niter = 29L,
    .approx_solver = function(...) {
      calls <<- calls + 1L
      stop("forced requested-only failure")
    },
    .diagnostics = collector
  )
  record <- collector$records[[1L]]

  expect_true(is.na(value))
  expect_identical(calls, 1L)
  expect_identical(record$selected_path, "approx_requested_only_invalid")
  expect_false(record$fallback_used)
  expect_identical(record$status, "raw_invalid")
  expect_equal(nrow(record$attempts), 1L)
  expect_identical(record$attempts$niter, 29L)
  expect_match(record$attempts$solver_error, "forced requested-only failure")
  expect_false(any(grepl("retry|exact", record$attempts$step)))
})

test_that("constant ground cost keeps raw failure audit and analytic distance", {
  requested <- getFromNamespace(
    "ppm_wasserstein_requested_approx", "ppmWass"
  )
  A <- cbind(mz = c(50, 60), intensity = c(0.5, 0.5))
  B <- cbind(mz = c(500, 600), intensity = c(0.5, 0.5))
  calls <- 0L
  nonfinite_plan <- list(
    from = c(1L, 2L, 1L, 2L), to = c(1L, 1L, 2L, 2L),
    mass = rep(NaN, 4L)
  )
  collector <- new.env(parent = emptyenv())
  value <- requested(
    A, B, ppm = 15, transition_mult = 3, ot_method = "sinkhorn",
    sinkhorn_epsilon = 0.05, sinkhorn_niter = 30L,
    .approx_solver = function(...) {
      calls <<- calls + 1L
      nonfinite_plan
    },
    .diagnostics = collector
  )
  record <- collector$records[[1L]]

  expect_equal(value, 1, tolerance = 0)
  expect_identical(calls, 1L)
  expect_identical(
    record$selected_path,
    "analytic_constant_cost_with_raw_requested_diagnostic"
  )
  expect_identical(
    record$status, "analytic_constant_cost_raw_diagnostic_invalid"
  )
  expect_false(record$fallback_used)
  expect_equal(nrow(record$attempts), 1L)
  expect_identical(record$attempts$step, "approx_requested_only")
  expect_identical(record$attempts$niter, 30L)
  expect_identical(record$attempts$validation_error, "nonfinite_plan_mass")
  expect_identical(record$attempts$nonfinite_mass, 4L)
})

test_that("NA-safe publication grouping retains exact-OT inactive settings", {
  group_key <- getFromNamespace("na_safe_group_key", "ppmWass")
  data <- data.frame(
    solver_backend = c("exact", "exact", "sinkhorn"),
    sinkhorn_epsilon_effective = c(NA_real_, NA_real_, 0.05),
    sinkhorn_iterations_effective = c(NA_integer_, NA_integer_, 100L),
    transition_multiplier = c(3, 3, 3),
    stringsAsFactors = FALSE
  )
  key <- group_key(data, names(data))
  groups <- split(seq_len(nrow(data)), key)

  expect_length(groups, 2L)
  expect_true(any(vapply(groups, identical, logical(1), c(1L, 2L))))
  expect_identical(sort(unlist(groups, use.names = FALSE)), 1:3)
})

test_that("Sinkhorn publication summaries require every replicate to be finite", {
  summarize_metric <- getFromNamespace(
    "summarize_publication_metric", "ppmWass"
  )
  incomplete <- summarize_metric(c(1, NA_real_, 0), require_complete = TRUE)
  descriptive <- summarize_metric(c(1, NA_real_, 0), require_complete = FALSE)

  expect_identical(incomplete$n_finite, 2L)
  expect_false(incomplete$complete)
  expect_true(is.na(incomplete$mean))
  expect_true(is.na(incomplete$sd))
  expect_equal(descriptive$mean, 0.5)
  expect_equal(descriptive$sd, stats::sd(c(1, 0)))
})

test_that("generic distance calls forward OT diagnostics and pair context", {
  skip_if_not_installed("transport")
  distance <- getFromNamespace("compute_distance", "ppmWass")
  collector <- new.env(parent = emptyenv())
  A <- cbind(mz = c(100, 150), intensity = c(0.6, 0.4))
  B <- cbind(mz = c(100.001, 140, 200), intensity = c(0.2, 0.3, 0.5))

  value <- distance(
    A, B, method = "ppm_wasserstein", ppm = 15, ot_method = "exact",
    ot_diagnostics = collector,
    ot_context = list(query_id = "query-A", library_id = "library-B"),
    ot_marginal_tolerance = 1e-8
  )

  expect_true(is.finite(value))
  expect_length(collector$records, 1L)
  expect_identical(collector$records[[1]]$selected_path, "exact_primary")
  expect_identical(collector$records[[1]]$context$query_id, "query-A")
  expect_identical(collector$records[[1]]$context$library_id, "library-B")
  expect_identical(collector$records[[1]]$marginal_tolerance, 1e-8)
})

test_that("derived representation exposes accurate compatibility aliases", {
  path <- tempfile(fileext = ".msp")
  writeLines(c(
    "Name: example",
    "InChIKey: AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    "RI: 1000",
    "Num Peaks: 4",
    "50.1 10",
    "75.2 20",
    "100.3 30",
    "125.4 40"
  ), path)
  p <- eihrms_default_params()
  p$use_split_loss <- TRUE
  x <- build_spectra_from_msp(path, p, progress = FALSE)
  expect_identical(x$derived_list, x$loss_list)
  expect_identical(x$derived_anchor_list, x$loss_anchor_list)
  expect_identical(x$derived_pair_list, x$loss_pair_list)
})

test_that("anchored and pairwise derived channels obey their definitions", {
  build <- getFromNamespace("build_loss_peaks", "ppmWass")
  peaks <- cbind(
    mz = c(50, 75, 100, 140),
    intensity = c(0.40, 0.30, 0.20, 0.10)
  )
  p <- eihrms_default_params()
  p$loss_min <- 1
  p$loss_max <- 500
  p$loss_top_peaks <- 4L
  p$loss_max_peaks <- 100L
  p$derivatization_raw_loss_weight <- 1
  info <- build(peaks, p, return_info = TRUE)

  expect_true(nrow(info$lossA_peaks) > 0)
  expect_true(nrow(info$lossB_peaks) > 0)
  expect_equal(sum(info$lossA_peaks[, 2]), 1, tolerance = 1e-14)
  expect_equal(sum(info$lossB_peaks[, 2]), 1, tolerance = 1e-14)
  expect_lte(nrow(info$lossA_peaks), p$loss_max_peaks)
  expect_lte(nrow(info$lossB_peaks), p$loss_max_peaks)

  # Exercise the package implementation with a controlled reference override:
  # anchored values change, while pairwise construction never sees Mref.
  anchor_error_da <- 0.25
  shifted <- build(
    peaks, p, return_info = TRUE,
    reference_mz = info$Mref + anchor_error_da
  )
  expect_identical(shifted$lossB_peaks, info$lossB_peaks)
  expect_false(identical(shifted$lossA_peaks, info$lossA_peaks))
  expect_equal(
    shifted$lossA_peaks[, 1] - info$lossA_peaks[, 1],
    rep(anchor_error_da, nrow(info$lossA_peaks)),
    tolerance = 1e-14
  )
  expect_false(identical(shifted$loss_peaks, info$loss_peaks))
})

test_that("derived channels enforce range, cap, and one-stage normalization", {
  build <- getFromNamespace("build_loss_peaks", "ppmWass")
  peaks <- cbind(
    mz = c(50, 75, 100, 140),
    intensity = c(0.40, 0.30, 0.20, 0.10)
  )
  p <- eihrms_default_params()
  p$loss_min <- 30
  p$loss_max <- 90
  p$loss_top_peaks <- 4L
  p$loss_max_peaks <- 2L
  p$derivatization_raw_loss_weight <- 1

  out <- build(peaks, p, return_info = TRUE, reference_mz = 140)
  scaled_peaks <- peaks
  scaled_peaks[, 2] <- scaled_peaks[, 2] * 10
  scaled <- build(scaled_peaks, p, return_info = TRUE, reference_mz = 140)
  for (channel in c("loss_peaks", "lossA_peaks", "lossB_peaks")) {
    x <- out[[channel]]
    expect_lte(nrow(x), p$loss_max_peaks)
    expect_true(all(x[, 1] >= p$loss_min & x[, 1] <= p$loss_max))
    expect_equal(sum(x[, 2]), 1, tolerance = 1e-14)
    expect_equal(scaled[[channel]], x, tolerance = 1e-14)
  }

  # After range filtering the anchored raw weights are 0.4, 0.3, 0.2;
  # the cap retains 0.4 and 0.3, then normalizes once at summarization.
  expect_equal(out$lossA_peaks[, 1], c(65, 90), tolerance = 0)
  expect_equal(out$lossA_peaks[, 2], c(3 / 7, 4 / 7), tolerance = 1e-14)
})

test_that("derived construction rejects NA and malformed inputs explicitly", {
  build <- getFromNamespace("build_loss_peaks", "ppmWass")
  p <- eihrms_default_params()
  peaks <- cbind(mz = c(50, 100), intensity = c(0.6, 0.4))

  bad_mz <- peaks
  bad_mz[1, 1] <- NA_real_
  expect_error(build(bad_mz, p), "NA removal is not implicit")
  bad_intensity <- peaks
  bad_intensity[2, 2] <- Inf
  expect_error(build(bad_intensity, p), "finite numeric")
  expect_error(build(c(50, 1, 100, 1), p), "matrix-like")
  expect_error(build(matrix(c(50, 100), ncol = 1), p), "matrix-like")
  expect_error(build(peaks, p, reference_mz = NA_real_), "reference_mz")
})

test_that("empty derived branches are safe", {
  build <- getFromNamespace("build_loss_peaks", "ppmWass")
  empty <- matrix(numeric(0), nrow = 0, ncol = 2,
                  dimnames = list(NULL, c("mz", "intensity")))
  result <- build(empty, eihrms_default_params(), return_info = TRUE)
  expect_equal(nrow(result$loss_peaks), 0)
  expect_equal(nrow(result$lossA_peaks), 0)
  expect_equal(nrow(result$lossB_peaks), 0)
  expect_true(is.na(result$Mref))
})

test_that("explicit-TIC diagnostics use accurate terminology", {
  signal <- cbind(mz = c(100, 150), intensity = c(1, 1))
  noise <- cbind(mz = c(100, 200), intensity = c(1, 3))
  result <- getFromNamespace("mix_signal_noise_impl", "ppmWass")(
    signal, noise, signal_fraction = 0.2, dropout_prob = 0,
    merge_ppm = 10, seed = 11
  )
  d <- result$diagnostics
  expect_equal(unname(d[["requested_premerge_signal_tic_weight"]]), 0.2)
  expect_equal(unname(d[["assigned_premerge_signal_tic_weight"]]), 0.2)
  expect_equal(
    unname(d[["source_labelled_postmerge_signal_tic_fraction"]]), 0.2
  )
  expect_equal(sum(result$spectrum[, 2]), 1, tolerance = 1e-14)
})

test_that("tie-aware statistics use fractional expected outcomes", {
  helper <- publication_stats_helper_test_path()
  env <- new.env(parent = baseenv())
  sys.source(helper, envir = env)
  ids <- c("q", "relevant", "irrelevant")
  keys <- c(q = "AAAAAAAAAAAAAA-ONE", relevant = "AAAAAAAAAAAAAA-ONE",
            irrelevant = "BBBBBBBBBBBBBB-TWO")
  dm <- matrix(c(
    0, 0.1, 0.1,
    0.1, 0, 0.5,
    0.1, 0.5, 0
  ), nrow = 3, byrow = TRUE, dimnames = list(ids, ids))
  pq <- env$per_query_metrics_tie_aware(dm, keys, tie_tolerance = 0,
                                        random_seed = 1)
  q <- pq[pq$query_id == "q", ]
  expect_equal(q$top1_optimistic, 1)
  expect_equal(q$top1_fractional, 0.5)
  expect_equal(q$top1_pessimistic, 0)
  expect_equal(q$rr_fractional, 0.75)
  expect_true(q$full_top1_tie_block_crossing)
  expect_true(q$full_top1_optimistic_pessimistic_crossing)
  expect_true(q$full_top1_affected)
  expect_false(q$full_top5_affected)
})

test_that("Top-5 and Top-10 tie sensitivity and crossing flags are reported", {
  helper <- publication_stats_helper_test_path()
  env <- new.env(parent = baseenv())
  sys.source(helper, envir = env)

  at_five <- env$first_relevant_tie_metrics(
    dists = c(1:4, 5, 5, 6),
    relevant = c(rep(FALSE, 4), TRUE, FALSE, FALSE)
  )
  at_ten <- env$first_relevant_tie_metrics(
    dists = c(1:9, 10, 10, 11),
    relevant = c(rep(FALSE, 9), TRUE, FALSE, FALSE)
  )
  expect_equal(at_five$top5_optimistic, 1)
  expect_equal(at_five$top5_pessimistic, 0)
  expect_equal(at_five$top10_optimistic, at_five$top10_pessimistic)
  expect_equal(at_ten$top10_optimistic, 1)
  expect_equal(at_ten$top10_pessimistic, 0)

  ids <- c("a", "b", "c")
  keys <- stats::setNames(rep("AAAAAAAAAAAAAA-ONE", 3), ids)
  dm <- matrix(
    c(0, 0.1, 0.1, 0.1, 0, 0.1, 0.1, 0.1, 0),
    nrow = 3, byrow = TRUE, dimnames = list(ids, ids)
  )
  pq <- env$per_query_metrics_tie_aware(
    dm, keys, tie_tolerance = 0, random_seed = 1
  )
  # Each first-relevant block has size two, but both tied candidates are
  # relevant. Thus tie_size > 1 is descriptive and does not imply that the
  # optimistic and pessimistic Top-K outcomes differ.
  expect_true(all(pq$full_has_relevant_tie))
  expect_true(all(pq$full_top1_tie_block_crossing))
  expect_false(any(pq$full_top1_affected))
  expect_false(any(pq$full_top5_affected))
  expect_false(any(pq$full_top10_affected))

  sensitivity <- env$summarize_tie_sensitivity(pq, "TEST", "method")
  expect_true(all(c("top5", "top10") %in% sensitivity$metric))

  diagnostic <- env$summarize_tie_diagnostics(pq, "TEST", "method")
  top1 <- diagnostic[diagnostic$metric == "top1", ]
  expect_equal(top1$n_with_relevant_tie, 3)
  expect_equal(top1$n_tie_blocks_crossing_cutoff, 3)
  expect_equal(top1$n_affected, 0)
  expect_equal(top1$proportion_affected, 0)
  expect_equal(top1$n_tied_not_affected, 3)
})

test_that("cluster and exploratory query bootstrap tables state their estimands", {
  helper <- publication_stats_helper_test_path()
  env <- new.env(parent = baseenv())
  sys.source(helper, envir = env)

  ids <- c("a", "b", "c")
  keys <- stats::setNames(rep("AAAAAAAAAAAAAA-ONE", 3), ids)
  dm <- matrix(
    c(0, 0.1, 0.1, 0.1, 0, 0.1, 0.1, 0.1, 0),
    nrow = 3, byrow = TRUE, dimnames = list(ids, ids)
  )
  pq <- env$per_query_metrics_tie_aware(dm, keys, random_seed = 1)
  cluster <- env$summarize_tie_aware_metrics(
    pq, "TEST", "method", query_boot_R = 20, cluster_boot_R = 20, seed = 1
  )
  expect_true(all(cluster$analysis_level == "cluster"))
  expect_true(all(cluster$status == "primary"))
  expect_true(all(nzchar(cluster$cluster_unit)))
  expect_true(all(grepl("query_weighted_mean", cluster$estimand)))
  expect_true(all(cluster$query_ci_analysis_level == "query"))
  expect_true(all(cluster$query_ci_status == "exploratory"))

  query <- env$query_bootstrap_ci_table(cluster)
  expect_true(all(query$analysis_level == "query"))
  expect_true(all(query$status == "exploratory"))
  expect_true(all(query$inference_method == "iid_query_bootstrap"))
})
