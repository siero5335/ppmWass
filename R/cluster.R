#' Cluster with OPTICS on Distance Matrix
#'
#' @param sim_cluster Similarity matrix for clustering.
#' @param min_pts Minimum points for OPTICS (defaults to max(3, floor(n/20))).
#' @param eps_cl Epsilon for extractDBSCAN.
#' @return A list with optics result and cluster labels.
#' @export
cluster_optics <- function(sim_cluster, min_pts = NULL, eps_cl = 0.73) {
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required for OPTICS clustering.")
  }

  n <- nrow(sim_cluster)
  if (n < 2) {
    stop("cluster_optics requires at least 2 spectra.")
  }
  if (is.null(min_pts)) {
    min_pts <- min(n, max(2, floor(n / 20)))
  } else {
    min_pts <- as.integer(min_pts)
    if (!is.finite(min_pts) || min_pts < 2 || min_pts > n) {
      stop("min_pts must be an integer between 2 and nrow(sim_cluster).")
    }
  }

  dist_for_cluster <- stats::as.dist(1 - sim_cluster)

  res_h2 <- dbscan::optics(dist_for_cluster, minPts = min_pts)
  res_d2 <- dbscan::extractDBSCAN(res_h2, eps_cl = eps_cl)
  clusters <- as.factor(res_d2$cluster)
  names(clusters) <- rownames(sim_cluster)

  list(
    optics = res_h2,
    clusters = clusters,
    dist = dist_for_cluster
  )
}

#' Run UMAP from Similarity Matrix
#'
#' @param sim_cluster Similarity matrix for clustering.
#' @param n_neighbors Number of neighbors (default min(10, floor(n/3))).
#' @param min_dist UMAP min_dist.
#' @param n_components Output dimensions.
#' @param seed Random seed.
#' @return A data frame with UMAP coordinates.
#' @export
run_umap <- function(sim_cluster, n_neighbors = NULL, min_dist = 0.05, n_components = 2, seed = 71) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for UMAP.")
  }

  n <- nrow(sim_cluster)
  if (n < 3) {
    stop("run_umap requires at least 3 spectra.")
  }
  if (is.null(n_neighbors)) {
    n_neighbors <- min(10, max(2, floor(n / 3)))
  } else {
    n_neighbors <- as.integer(n_neighbors)
    if (!is.finite(n_neighbors) || n_neighbors < 2 || n_neighbors >= n) {
      stop("n_neighbors must be an integer between 2 and nrow(sim_cluster) - 1.")
    }
  }

  dist_for_cluster <- stats::as.dist(1 - sim_cluster)
  um <- with_seed(seed, {
    uwot::umap(
      dist_for_cluster,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      n_components = n_components
    )
  })

  res_umap <- as.data.frame(um)
  if (n_components == 2) {
    colnames(res_umap) <- c("V1", "V2")
  }

  res_umap
}
