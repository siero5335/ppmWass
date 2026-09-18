test_that("run_umap validates sample size and neighbor range", {
  skip_if_not_installed("uwot")

  sim <- diag(2)
  expect_error(
    ppmWass::run_umap(sim),
    "at least 3 spectra"
  )

  sim3 <- diag(3)
  expect_error(
    ppmWass::run_umap(sim3, n_neighbors = 1),
    "between 2 and nrow"
  )
})

test_that("cluster_optics returns clusters named by similarity matrix rownames", {
  skip_if_not_installed("dbscan")

  sim <- matrix(
    c(
      1.0, 0.9, 0.2,
      0.9, 1.0, 0.2,
      0.2, 0.2, 1.0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )

  res <- ppmWass::cluster_optics(sim, min_pts = 2)
  expect_equal(names(res$clusters), rownames(sim))
})
