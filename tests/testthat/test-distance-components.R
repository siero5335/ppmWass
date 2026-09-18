expected_pairwise_frag_dist <- function(frag, params) {
  ids <- names(frag)
  n <- length(frag)
  out <- matrix(0, n, n, dimnames = list(ids, ids))
  mz_range <- if (!is.null(params$min_mz) && !is.null(params$max_mz)) {
    params$max_mz - params$min_mz
  } else {
    NULL
  }

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      d <- switch(
        params$distance_method,
        hellinger = ppmWass:::hellinger_distance(frag[[i]], frag[[j]], ppm = params$tol_ppm),
        wasserstein = ppmWass:::wasserstein_distance(
          frag[[i]], frag[[j]],
          ppm = params$tol_ppm,
          align = params$wasserstein_align
        ),
        stop("Unsupported test method: ", params$distance_method)
      )
      out[i, j] <- d
      out[j, i] <- d
    }
  }

  out
}

test_that("component distance path matches base matrix for hellinger", {
  params <- ppmWass::eihrms_default_params()
  params$distance_method <- "hellinger"
  params$use_parallel <- FALSE

  frag <- list(
    A = matrix(c(50, 1, 150, 1), ncol = 2, byrow = TRUE),
    B = matrix(c(55, 1, 250, 1), ncol = 2, byrow = TRUE),
    C = matrix(c(60, 1, 350, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    A = matrix(numeric(0), ncol = 2),
    B = matrix(numeric(0), ncol = 2),
    C = matrix(numeric(0), ncol = 2)
  )

  base <- ppmWass:::compute_distance_matrix(frag, loss, params, progress = FALSE)
  with_components <- ppmWass:::compute_distance_matrices(frag, loss, params, progress = FALSE)
  expected_frag <- expected_pairwise_frag_dist(frag, params)

  expect_equal(with_components$dist, base, tolerance = 1e-12)
  expect_equal(with_components$components$d_frag, expected_frag, tolerance = 1e-12)
})

test_that("component distance path matches base matrix for wasserstein", {
  skip_if_not_installed("transport")

  params <- ppmWass::eihrms_default_params()
  params$distance_method <- "wasserstein"
  params$use_parallel <- FALSE
  params$min_mz <- 35
  params$max_mz <- 650

  frag <- list(
    A = matrix(c(50, 1, 150, 1), ncol = 2, byrow = TRUE),
    B = matrix(c(55, 1, 250, 1), ncol = 2, byrow = TRUE),
    C = matrix(c(60, 1, 350, 1), ncol = 2, byrow = TRUE)
  )
  loss <- list(
    A = matrix(numeric(0), ncol = 2),
    B = matrix(numeric(0), ncol = 2),
    C = matrix(numeric(0), ncol = 2)
  )

  base <- ppmWass:::compute_distance_matrix(frag, loss, params, progress = FALSE)
  with_components <- ppmWass:::compute_distance_matrices(frag, loss, params, progress = FALSE)
  expected_frag <- expected_pairwise_frag_dist(frag, params)

  expect_equal(with_components$dist, base, tolerance = 1e-12)
  expect_equal(with_components$components$d_frag, expected_frag, tolerance = 1e-12)
})
