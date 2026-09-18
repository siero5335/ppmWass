#' Benchmark Spectral Similarity Methods
#'
#' Functions for evaluating and comparing spectral similarity metrics
#' using reference libraries (NIST, MoNA, MassBank).
#'
#' @name benchmark
NULL

#' Run Benchmark Suite
#'
#' Comprehensive benchmark comparing multiple distance methods.
#'
#' @param spectra_result Output from build_spectra() or build_spectra_from_msp().
#' @param methods Character vector of distance methods to compare.
#' @param ground_truth Optional data frame with true groupings (e.g., homolog series).
#' @param params Parameter list.
#' @param cluster_method Clustering method ("kmeans" or "dbscan") for evaluation.
#' @param kmeans_nstart Number of random starts for k-means.
#' @param dbscan_eps Epsilon for DBSCAN (only used if cluster_method = "dbscan").
#' @param dbscan_minPts Minimum points for DBSCAN (only used if cluster_method = "dbscan").
#' @param umap_neighbors Number of UMAP neighbors (NULL uses min(15, n - 1)).
#' @param apply_weights If TRUE, apply spectral weighting before distance calculation.
#' @param weight_intensity_power Intensity power for preprocessing weights.
#' @param weight_freq_power m/z frequency power for preprocessing weights.
#' @param weight_bin_width Bin width for m/z frequency calculation.
#' @param apply_loss_weights If TRUE, apply intensity weighting to loss spectra.
#' @param threshold_similarity_transform Distance-to-similarity transform used
#'   for threshold-based metrics (FDR and homolog detection).
#' @param threshold_similarity_scale Optional scale for the exponential transform.
#' @return A list containing benchmark results for each method.
#' @export
run_benchmark <- function(spectra_result, 
                          methods = c("entropy", "hellinger", "cosine"),
                          ground_truth = NULL,
                          params = eihrms_default_params(),
                          cluster_method = c("kmeans", "dbscan"),
                          kmeans_nstart = 10,
                          dbscan_eps = 0.5,
                          dbscan_minPts = 3,
                          umap_neighbors = NULL,
                          apply_weights = FALSE,
                          weight_intensity_power = 0.5,
                          weight_freq_power = 0.5,
                          weight_bin_width = 1.0,
                          apply_loss_weights = TRUE,
                          threshold_similarity_transform = c("minmax", "linear_clamp", "reciprocal", "exponential", "linear"),
                          threshold_similarity_scale = NULL) {
  cluster_method <- match.arg(cluster_method)
  threshold_similarity_transform <- match.arg(threshold_similarity_transform)
  
  message("=== Running Benchmark Suite ===")
  message("Methods: ", paste(methods, collapse = ", "))
  message("Spectra: ", length(spectra_result$frag_list))
  
  frag_list <- spectra_result$frag_list
  loss_list <- spectra_result$loss_list
  if (isTRUE(apply_weights)) {
    message("Applying spectral weights before benchmarking...")
    mz_freq <- compute_mz_frequency(frag_list, bin_width = weight_bin_width)
    weighted <- weight_spectra(
      frag_list = frag_list,
      loss_list = loss_list,
      mz_freq = mz_freq,
      intensity_power = weight_intensity_power,
      freq_power = weight_freq_power,
      bin_width = weight_bin_width,
      apply_loss_weights = apply_loss_weights
    )
    frag_list <- weighted$frag_list
    loss_list <- weighted$loss_list
  }

  results <- list()
  homolog_note_printed <- FALSE
  
  for (method in methods) {
    message("\n--- Evaluating: ", method, " ---")
    
    # Update params with current method
    params_method <- params
    params_method$distance_method <- method
    
    # Compute distance matrix
    dist_mat <- compute_distance_matrix(
      frag_list,
      loss_list,
      params_method,
      progress = FALSE,
      loss_typ_list = spectra_result$loss_typ_list,
      mref_conf = spectra_result$mref_confidence,
      loss_anchor_list = spectra_result$loss_anchor_list,
      loss_pair_list = spectra_result$loss_pair_list,
      loss_anchor_typ_list = spectra_result$loss_anchor_typ_list,
      loss_pair_typ_list = spectra_result$loss_pair_typ_list
    )
    
    # Evaluate
    results[[method]] <- list(
      method = method,
      dist_matrix = dist_mat,
      metrics = list()
    )
    
    # 1. Replicate consistency (if InChIKey available)
    if ("inchikey" %in% names(spectra_result$df_spec)) {
      rep_metrics <- evaluate_replicate_consistency(
        dist_mat, 
        spectra_result$df_spec$inchikey,
        spectra_result$df_spec$id
      )
      results[[method]]$metrics$replicate <- rep_metrics
      message("  Replicate Top-1 Accuracy: ", round(rep_metrics$top1_accuracy, 3))
    }
    
    # 2. Clustering quality
    if (!is.null(ground_truth)) {
      cluster_metrics <- evaluate_clustering_quality(
        dist_mat,
        ground_truth,
        spectra_result$df_spec$id,
        cluster_method = cluster_method,
        kmeans_nstart = kmeans_nstart,
        dbscan_eps = dbscan_eps,
        dbscan_minPts = dbscan_minPts,
        umap_neighbors = umap_neighbors
      )
      results[[method]]$metrics$clustering <- cluster_metrics
      message("  NMI: ", round(cluster_metrics$nmi, 3))
      message("  ARI: ", round(cluster_metrics$ari, 3))
    }
    
    # 3. Silhouette analysis with compound class
    if ("compound_class" %in% names(spectra_result$df_spec)) {
      sil_metrics <- evaluate_silhouette(
        dist_mat,
        spectra_result$df_spec$compound_class
      )
      results[[method]]$metrics$silhouette <- sil_metrics
      message("  Mean Silhouette: ", round(sil_metrics$mean_silhouette, 3))
    }
    
    # 4. Homolog detection (if RI available)
    if (!all(is.na(spectra_result$ri))) {
      homolog_metrics <- evaluate_homolog_detection(
        dist_mat,
        spectra_result$df_spec,
        spectra_result$ri,
        similarity_transform = threshold_similarity_transform,
        similarity_scale = threshold_similarity_scale
      )
      results[[method]]$metrics$homolog <- homolog_metrics
      message("  Homolog Detection Score: ", round(homolog_metrics$score, 3))
      if (!homolog_note_printed && isTRUE(homolog_metrics$uses_inferred_classes)) {
        message("  Note: homolog detection used inferred compound_class labels and should be treated as an internal consistency check.")
        homolog_note_printed <- TRUE
      }
    }

    # 5. AUC-ROC
    if ("inchikey" %in% names(spectra_result$df_spec)) {
      auc_metrics <- evaluate_auc_roc(dist_mat, spectra_result$df_spec$inchikey)
      results[[method]]$metrics$auc_roc <- auc_metrics
      if (is.finite(auc_metrics$auc)) {
        message("  AUC-ROC: ", round(auc_metrics$auc, 3))
      }
    }

    # 6. Precision@K and MAP
    if ("inchikey" %in% names(spectra_result$df_spec)) {
      pk_metrics <- evaluate_precision_at_k(dist_mat, spectra_result$df_spec$inchikey)
      results[[method]]$metrics$precision_at_k <- pk_metrics
      if (is.finite(pk_metrics$mean_average_precision)) {
        message("  MAP: ", round(pk_metrics$mean_average_precision, 3))
      }
    }

    # 7. FDR by threshold
    if ("inchikey" %in% names(spectra_result$df_spec)) {
      fdr_metrics <- evaluate_fdr(
        dist_mat,
        spectra_result$df_spec$inchikey,
        similarity_transform = threshold_similarity_transform,
        similarity_scale = threshold_similarity_scale
      )
      results[[method]]$metrics$fdr <- fdr_metrics
      fdr_idx <- which.min(abs(fdr_metrics$threshold - 0.7))
      if (length(fdr_idx) == 1 && is.finite(fdr_metrics$fdr[fdr_idx])) {
        message("  FDR@0.7: ", round(fdr_metrics$fdr[fdr_idx], 3))
      }
    }
  }
  
  # Summary comparison
  results$summary <- summarize_benchmark(results, methods)
  
  results
}

#' Evaluate Replicate Consistency
#'
#' Tests if spectra from the same compound (same InChIKey) are most similar.
#'
#' @param dist_mat Distance matrix.
#' @param inchikeys Vector of InChIKeys for each spectrum.
#' @param ids Vector of spectrum IDs.
#' @return List with accuracy metrics.
#' @export
evaluate_replicate_consistency <- function(dist_mat, inchikeys, ids) {
  
  # Find compounds with multiple spectra
  inchikey_counts <- table(inchikeys[!is.na(inchikeys) & inchikeys != ""])
  multi_inchikeys <- names(inchikey_counts[inchikey_counts > 1])
  
  if (length(multi_inchikeys) == 0) {
    message("  No replicate spectra found (same InChIKey)")
    return(list(
      top1_accuracy = NA,
      top5_accuracy = NA,
      mrr = NA,
      n_queries = 0
    ))
  }
  
  # For each spectrum of multi-replicate compounds, check if nearest neighbor
  # is another spectrum of the same compound
  top1_correct <- 0
  top5_correct <- 0
  query_idx <- which(inchikeys %in% multi_inchikeys)
  n_queries <- length(query_idx)
  reciprocal_ranks <- numeric(n_queries)
  rr_idx <- 0
  
  for (ik in multi_inchikeys) {
    idx_same <- which(inchikeys == ik)
    
    for (i in idx_same) {
      # Find distances to all other spectra of same compound
      other_same <- setdiff(idx_same, i)
      
      # Get distances (excluding self)
      dists <- dist_mat[i, ]
      dists[i] <- Inf  # Exclude self
      
      # Rank all spectra by distance
      ranks <- rank(dists, ties.method = "min")
      
      # Check if any same-compound spectrum is in top-k
      ranks_same <- ranks[other_same]
      best_rank <- min(ranks_same)
      
      rr_idx <- rr_idx + 1

      if (best_rank == 1) top1_correct <- top1_correct + 1
      if (best_rank <= 5) top5_correct <- top5_correct + 1
      
      reciprocal_ranks[rr_idx] <- 1 / best_rank
    }
  }

  list(
    top1_accuracy = top1_correct / n_queries,
    top5_accuracy = top5_correct / n_queries,
    mrr = mean(reciprocal_ranks),  # Mean Reciprocal Rank
    n_queries = n_queries,
    n_compounds = length(multi_inchikeys)
  )
}

#' @keywords internal
#' @noRd
normalize_inchikey_prefix <- function(inchikeys, inchikey_chars = 14) {
  prefix <- substr(as.character(inchikeys), 1, inchikey_chars)
  invalid <- is.na(inchikeys) | inchikeys == ""
  prefix[invalid] <- NA_character_
  prefix
}

#' @keywords internal
#' @noRd
distance_to_similarity <- function(dist_mat,
                                   transform = c("linear", "linear_clamp", "minmax", "reciprocal", "exponential"),
                                   scale = NULL) {
  transform <- match.arg(transform)
  d <- as.matrix(dist_mat)

  if (transform == "linear") {
    return(1 - d)
  }
  if (transform == "linear_clamp") {
    return(1 - pmin(pmax(d, 0), 1))
  }
  if (transform == "reciprocal") {
    return(1 / (1 + pmax(d, 0)))
  }
  if (transform == "exponential") {
    if (is.null(scale)) {
      vals <- d[upper.tri(d)]
      vals <- vals[is.finite(vals) & vals > 0]
      scale <- if (length(vals) > 0) stats::median(vals) else 1
    }
    if (length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("scale must be a positive finite number for exponential similarity.")
    }
    return(exp(-pmax(d, 0) / scale))
  }

  vals <- d[upper.tri(d)]
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return(matrix(NA_real_, nrow = nrow(d), ncol = ncol(d), dimnames = dimnames(d)))
  }
  lo <- min(vals)
  hi <- max(vals)
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    sim <- matrix(0.5, nrow = nrow(d), ncol = ncol(d), dimnames = dimnames(d))
  } else {
    sim <- 1 - (d - lo) / (hi - lo)
    sim <- pmin(pmax(sim, 0), 1)
  }
  diag(sim) <- 1
  sim
}

#' @keywords internal
#' @noRd
extract_similarity_pairs <- function(dist_mat,
                                     ik_prefix,
                                     similarity_transform = "linear",
                                     similarity_scale = NULL) {
  idx <- which(upper.tri(dist_mat), arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(list(similarity = numeric(0), same_compound = logical(0)))
  }

  valid <- !is.na(ik_prefix[idx[, 1]]) & !is.na(ik_prefix[idx[, 2]])
  if (!any(valid)) {
    return(list(similarity = numeric(0), same_compound = logical(0)))
  }

  sim_mat <- distance_to_similarity(
    dist_mat,
    transform = similarity_transform,
    scale = similarity_scale
  )
  similarity <- sim_mat[idx][valid]
  same_compound <- ik_prefix[idx[, 1]][valid] == ik_prefix[idx[, 2]][valid]

  list(similarity = similarity, same_compound = same_compound)
}

rank_binary_auc <- function(score, is_positive) {
  n_pos <- sum(is_positive)
  n_neg <- sum(!is_positive)
  if (!n_pos || !n_neg) return(NA_real_)

  rank_score <- rank(score, ties.method = "average")
  n_pos_double <- as.double(n_pos)
  n_neg_double <- as.double(n_neg)
  (sum(rank_score[is_positive]) -
     n_pos_double * (n_pos_double + 1) / 2) /
    (n_pos_double * n_neg_double)
}

#' Evaluate AUC-ROC for Spectral Similarity
#'
#' Calculates AUC for distinguishing same-compound vs different-compound pairs.
#'
#' @param dist_mat Distance matrix.
#' @param inchikeys Vector of InChIKeys for each spectrum.
#' @param inchikey_chars Number of InChIKey characters for grouping.
#' @param similarity_transform Monotone transform used to score pairs.
#'   The default keeps raw distance ranks for AUC.
#' @param similarity_scale Optional scale for the exponential transform.
#' @return List with AUC, n_positive, n_negative.
#' @export
evaluate_auc_roc <- function(dist_mat,
                             inchikeys,
                             inchikey_chars = 14,
                             similarity_transform = "linear",
                             similarity_scale = NULL) {
  if (nrow(dist_mat) < 3) {
    return(list(auc = NA_real_, n_positive = 0L, n_negative = 0L))
  }

  ik_prefix <- normalize_inchikey_prefix(inchikeys, inchikey_chars)
  pairs <- extract_similarity_pairs(
    dist_mat,
    ik_prefix,
    similarity_transform = similarity_transform,
    similarity_scale = similarity_scale
  )
  similarity <- pairs$similarity
  same_compound <- pairs$same_compound

  n_pos <- sum(same_compound)
  n_neg <- sum(!same_compound)
  if (n_pos == 0 || n_neg == 0) {
    return(list(auc = NA_real_, n_positive = n_pos, n_negative = n_neg))
  }

  auc <- rank_binary_auc(similarity, same_compound)

  list(auc = auc, n_positive = n_pos, n_negative = n_neg)
}

#' Evaluate Precision at K
#'
#' Calculates mean precision at K and mean average precision for replicate queries.
#'
#' @param dist_mat Distance matrix.
#' @param inchikeys Vector of InChIKeys for each spectrum.
#' @param k_values K values to evaluate.
#' @param inchikey_chars Number of InChIKey characters for grouping.
#' @return List with precision_at_k, mean_average_precision, and n_queries.
#' @export
evaluate_precision_at_k <- function(dist_mat, inchikeys, k_values = c(1, 5, 10, 20), inchikey_chars = 14) {
  ik_prefix <- normalize_inchikey_prefix(inchikeys, inchikey_chars)

  ik_counts <- table(ik_prefix[!is.na(ik_prefix)])
  multi_ik <- names(ik_counts[ik_counts > 1])
  query_idx <- which(ik_prefix %in% multi_ik)
  n_queries <- length(query_idx)

  if (n_queries == 0) {
    return(list(
      precision_at_k = stats::setNames(rep(NA_real_, length(k_values)), paste0("P@", k_values)),
      mean_average_precision = NA_real_,
      n_queries = 0L
    ))
  }

  precision_sums <- numeric(length(k_values))
  ap_sum <- 0

  for (qi in query_idx) {
    dists <- dist_mat[qi, ]
    sorted_idx <- setdiff(order(dists, na.last = TRUE), qi)
    if (length(sorted_idx) == 0) next

    is_relevant <- ik_prefix[sorted_idx] == ik_prefix[qi]
    is_relevant[is.na(is_relevant)] <- FALSE

    for (ki in seq_along(k_values)) {
      k <- min(k_values[ki], length(sorted_idx))
      precision_sums[ki] <- precision_sums[ki] + sum(is_relevant[seq_len(k)]) / k
    }

    n_relevant_total <- sum(ik_prefix == ik_prefix[qi], na.rm = TRUE) - 1
    if (n_relevant_total > 0) {
      cum_rel <- cumsum(is_relevant)
      precision_at_rank <- cum_rel / seq_along(cum_rel)
      ap_sum <- ap_sum + sum(precision_at_rank[is_relevant]) / n_relevant_total
    }
  }

  precision_at_k <- precision_sums / n_queries
  names(precision_at_k) <- paste0("P@", k_values)

  list(
    precision_at_k = precision_at_k,
    mean_average_precision = ap_sum / n_queries,
    n_queries = n_queries
  )
}

#' Evaluate False Discovery Rate at Similarity Thresholds
#'
#' @param dist_mat Distance matrix.
#' @param inchikeys Vector of InChIKeys for each spectrum.
#' @param thresholds Similarity thresholds to evaluate.
#' @param inchikey_chars Number of InChIKey characters for grouping.
#' @param similarity_transform Distance-to-similarity transform. The default
#'   min-max calibration makes threshold grids usable for unbounded distances.
#' @param similarity_scale Optional scale for the exponential transform.
#' @return Data frame with threshold, match counts, and FDR.
#' @export
evaluate_fdr <- function(dist_mat,
                         inchikeys,
                         thresholds = seq(0.5, 0.95, by = 0.05),
                         inchikey_chars = 14,
                         similarity_transform = c("minmax", "linear_clamp", "reciprocal", "exponential", "linear"),
                         similarity_scale = NULL) {
  similarity_transform <- match.arg(similarity_transform)
  ik_prefix <- normalize_inchikey_prefix(inchikeys, inchikey_chars)
  pairs <- extract_similarity_pairs(
    dist_mat,
    ik_prefix,
    similarity_transform = similarity_transform,
    similarity_scale = similarity_scale
  )
  similarity <- pairs$similarity
  same_compound <- pairs$same_compound

  out <- data.frame(
    threshold = thresholds,
    n_matches = integer(length(thresholds)),
    n_true_pos = integer(length(thresholds)),
    n_false_pos = integer(length(thresholds)),
    fdr = numeric(length(thresholds)),
    stringsAsFactors = FALSE
  )

  if (length(similarity) == 0) {
    out$fdr[] <- NA_real_
    return(out)
  }

  for (i in seq_along(thresholds)) {
    above <- similarity >= thresholds[i]
    n_above <- sum(above)
    n_tp <- sum(above & same_compound)
    n_fp <- sum(above & !same_compound)

    out$n_matches[i] <- n_above
    out$n_true_pos[i] <- n_tp
    out$n_false_pos[i] <- n_fp
    out$fdr[i] <- if (n_above > 0) n_fp / n_above else 0
  }

  out
}

#' Evaluate Clustering Quality
#'
#' Compares clustering results with ground truth labels.
#'
#' @param dist_mat Distance matrix.
#' @param ground_truth Data frame with 'id' and 'group' columns.
#' @param ids Vector of spectrum IDs.
#' @param cluster_method Clustering method ("kmeans" or "dbscan").
#' @param kmeans_nstart Number of random starts for k-means.
#' @param dbscan_eps Epsilon for DBSCAN (only used if cluster_method = "dbscan").
#' @param dbscan_minPts Minimum points for DBSCAN (only used if cluster_method = "dbscan").
#' @param umap_neighbors Number of UMAP neighbors (NULL uses min(15, n - 1)).
#' @return List with NMI, ARI, and other metrics.
#' @export
evaluate_clustering_quality <- function(dist_mat,
                                        ground_truth,
                                        ids,
                                        cluster_method = c("kmeans", "dbscan"),
                                        kmeans_nstart = 10,
                                        dbscan_eps = 0.5,
                                        dbscan_minPts = 3,
                                        umap_neighbors = NULL) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for clustering evaluation.")
  }
  cluster_method <- match.arg(cluster_method)
  
  # Match IDs
  common_ids <- intersect(ids, ground_truth$id)
  if (length(common_ids) < 10) {
    warning("Too few common IDs for clustering evaluation")
    return(list(nmi = NA, ari = NA))
  }
  
  idx <- match(common_ids, ids)
  true_labels <- ground_truth$group[match(common_ids, ground_truth$id)]
  
  # Subset distance matrix
  dist_sub <- dist_mat[idx, idx]
  
  # Perform clustering (DBSCAN)
  # First do UMAP to get coordinates
  if (is.null(umap_neighbors)) {
    umap_neighbors <- min(15, nrow(dist_sub) - 1)
  }

  umap_res <- with_seed(42, {
    uwot::umap(stats::as.dist(dist_sub), n_neighbors = umap_neighbors)
  })

  if (cluster_method == "kmeans") {
    k <- length(unique(true_labels))
    if (k < 2) {
      return(list(nmi = NA, ari = NA))
    }
    km <- stats::kmeans(umap_res, centers = k, nstart = kmeans_nstart)
    pred_labels <- km$cluster
  } else {
    if (!requireNamespace("dbscan", quietly = TRUE)) {
      stop("Package 'dbscan' is required for clustering evaluation.")
    }
    db_res <- dbscan::dbscan(umap_res, eps = dbscan_eps, minPts = dbscan_minPts)
    pred_labels <- db_res$cluster
  }
  
  # Calculate NMI and ARI
  nmi <- calculate_nmi(true_labels, pred_labels)
  ari <- calculate_ari(true_labels, pred_labels)
  
  list(
    nmi = nmi,
    ari = ari,
    n_true_clusters = length(unique(true_labels)),
    n_pred_clusters = length(unique(pred_labels[pred_labels > 0])),
    n_samples = length(common_ids)
  )
}

#' Evaluate Silhouette Score
#'
#' Calculates silhouette score using compound class as labels.
#'
#' @param dist_mat Distance matrix.
#' @param labels Vector of class labels.
#' @return List with silhouette metrics.
#' @export
evaluate_silhouette <- function(dist_mat, labels) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required for silhouette evaluation.")
  }
  
  # Remove unclassified and NA
  valid_idx <- !is.na(labels) & labels != "" & labels != "unclassified"
  
  if (sum(valid_idx) < 10) {
    return(list(mean_silhouette = NA, per_class = list()))
  }
  
  dist_sub <- dist_mat[valid_idx, valid_idx]
  labels_sub <- labels[valid_idx]
  
  # Need at least 2 classes
  if (length(unique(labels_sub)) < 2) {
    return(list(mean_silhouette = NA, per_class = list()))
  }
  
  # Calculate silhouette
  sil <- cluster::silhouette(as.integer(factor(labels_sub)), stats::as.dist(dist_sub))
  
  # Per-class average
  per_class <- tapply(sil[, 3], labels_sub, mean)
  
  list(
    mean_silhouette = mean(sil[, 3]),
    per_class = as.list(per_class),
    n_samples = sum(valid_idx),
    n_classes = length(unique(labels_sub))
  )
}

#' Evaluate Homolog Detection
#'
#' Tests if homologous series (same class, different RI) are detected.
#'
#' Note: When compound_class is inferred by this package, the score reflects
#' internal consistency between the similarity metric and inferred classes,
#' not an external ground truth taxonomy.
#'
#' @param dist_mat Distance matrix.
#' @param df_spec Data frame with spectral info.
#' @param ri Named vector of retention indices.
#' @param ri_diff_threshold Minimum RI difference for homolog pairs.
#' @param sim_threshold Similarity threshold for considering pairs.
#' @param similarity_transform Distance-to-similarity transform. The default
#'   min-max calibration makes the threshold meaningful for unbounded distances.
#' @param similarity_scale Optional scale for the exponential transform.
#' @return List with homolog detection metrics.
#' @export
evaluate_homolog_detection <- function(dist_mat, 
                                       df_spec,
                                       ri,
                                       ri_diff_threshold = 50,
                                       sim_threshold = 0.5,
                                       similarity_transform = c("minmax", "linear_clamp", "reciprocal", "exponential", "linear"),
                                       similarity_scale = NULL) {
  similarity_transform <- match.arg(similarity_transform)
  
  ids <- df_spec$id
  n <- length(ids)
  
  # Convert distance to similarity
  sim_mat <- distance_to_similarity(
    dist_mat,
    transform = similarity_transform,
    scale = similarity_scale
  )
  
  # Find pairs that are similar but have different RI
  similar_pairs <- which(sim_mat > sim_threshold & upper.tri(sim_mat), arr.ind = TRUE)
  uses_inferred_classes <- FALSE
  ground_truth_note <- "No compound_class labels were available; the score only reflects the fraction of similar pairs with RI separation above the threshold."
  
  if (nrow(similar_pairs) == 0) {
    return(list(
      score = 0,
      n_homolog_pairs = 0,
      n_true_homologs = NA,
      n_similar_pairs = 0,
      uses_inferred_classes = uses_inferred_classes,
      ground_truth_note = ground_truth_note
    ))
  }
  
  # Check RI difference
  ri_diffs <- abs(ri[ids[similar_pairs[, 1]]] - ri[ids[similar_pairs[, 2]]])
  
  # Homolog candidates: similar spectra with large RI difference
  is_homolog_candidate <- ri_diffs > ri_diff_threshold
  n_homolog_pairs <- sum(is_homolog_candidate, na.rm = TRUE)
  
  # If compound class is available, check if same class
  if ("compound_class" %in% names(df_spec)) {
    classes <- df_spec$compound_class
    class_sources <- if ("compound_class_source" %in% names(df_spec)) unique(stats::na.omit(df_spec$compound_class_source)) else character(0)
    uses_inferred_classes <- length(class_sources) > 0 && all(class_sources == "inferred")
    ground_truth_note <- if (uses_inferred_classes) {
      paste(
        "compound_class was inferred internally by this package.",
        "Interpret homolog detection as an internal consistency score, not external ground truth."
      )
    } else {
      "compound_class labels from df_spec were used as homolog-class labels."
    }
    same_class <- classes[similar_pairs[, 1]] == classes[similar_pairs[, 2]]
    
    # True homologs: similar, large RI diff, same class
    true_homologs <- is_homolog_candidate & same_class
    n_true_homologs <- sum(true_homologs, na.rm = TRUE)
    
    # Score: ratio of true homologs among candidates
    score <- if (n_homolog_pairs > 0) n_true_homologs / n_homolog_pairs else 0
  } else {
    score <- n_homolog_pairs / nrow(similar_pairs)
    n_true_homologs <- NA
  }
  
  list(
    score = score,
    n_homolog_pairs = n_homolog_pairs,
    n_true_homologs = if (exists("n_true_homologs")) n_true_homologs else NA,
    n_similar_pairs = nrow(similar_pairs),
    uses_inferred_classes = uses_inferred_classes,
    ground_truth_note = ground_truth_note
  )
}

#' Calculate Normalized Mutual Information
#' @keywords internal
#' @noRd
calculate_nmi <- function(true_labels, pred_labels) {
  # Remove noise cluster (0) from prediction
  valid <- pred_labels != 0
  if (sum(valid) < 2) return(NA)
  
  true_sub <- true_labels[valid]
  pred_sub <- pred_labels[valid]
  
  # Contingency table
  cont <- table(true_sub, pred_sub)
  
  # Marginal probabilities
  n <- sum(cont)
  p_true <- rowSums(cont) / n
  p_pred <- colSums(cont) / n
  
  # Joint probability
  p_joint <- cont / n
  
  # Mutual information
  mi <- 0
  for (i in seq_len(nrow(cont))) {
    for (j in seq_len(ncol(cont))) {
      if (p_joint[i, j] > 0) {
        mi <- mi + p_joint[i, j] * log(p_joint[i, j] / (p_true[i] * p_pred[j]))
      }
    }
  }
  
  # Entropies
  h_true <- -sum(p_true[p_true > 0] * log(p_true[p_true > 0]))
  h_pred <- -sum(p_pred[p_pred > 0] * log(p_pred[p_pred > 0]))
  
  # NMI
  if (h_true + h_pred == 0) return(0)
  2 * mi / (h_true + h_pred)
}

#' Calculate Adjusted Rand Index
#' @keywords internal
#' @noRd
calculate_ari <- function(true_labels, pred_labels) {
  valid <- pred_labels != 0
  if (sum(valid) < 2) return(NA)
  
  true_sub <- true_labels[valid]
  pred_sub <- pred_labels[valid]
  
  # Contingency table
  cont <- table(true_sub, pred_sub)
  
  # Sum of combinations
  sum_comb_rows <- sum(choose(rowSums(cont), 2))
  sum_comb_cols <- sum(choose(colSums(cont), 2))
  sum_comb_cells <- sum(choose(cont, 2))
  n <- sum(cont)
  sum_comb_n <- choose(n, 2)
  
  # Expected index
  expected <- sum_comb_rows * sum_comb_cols / sum_comb_n
  
  # Max index
  max_index <- (sum_comb_rows + sum_comb_cols) / 2
  
  # ARI
  if (max_index - expected == 0) return(0)
  (sum_comb_cells - expected) / (max_index - expected)
}

#' Summarize Benchmark Results
#' @keywords internal
#' @noRd
summarize_benchmark <- function(results, methods) {
  
  summary_df <- data.frame(method = methods, stringsAsFactors = FALSE)
  
  for (m in methods) {
    res <- results[[m]]$metrics
    
    if (!is.null(res$replicate)) {
      summary_df$replicate_top1[summary_df$method == m] <- res$replicate$top1_accuracy
      summary_df$replicate_mrr[summary_df$method == m] <- res$replicate$mrr
    }
    
    if (!is.null(res$clustering)) {
      summary_df$nmi[summary_df$method == m] <- res$clustering$nmi
      summary_df$ari[summary_df$method == m] <- res$clustering$ari
    }
    
    if (!is.null(res$silhouette)) {
      summary_df$silhouette[summary_df$method == m] <- res$silhouette$mean_silhouette
    }
    
    if (!is.null(res$homolog)) {
      summary_df$homolog_score[summary_df$method == m] <- res$homolog$score
    }

    if (!is.null(res$auc_roc)) {
      summary_df$auc_roc[summary_df$method == m] <- res$auc_roc$auc
    }

    if (!is.null(res$precision_at_k)) {
      summary_df$map[summary_df$method == m] <- res$precision_at_k$mean_average_precision
      summary_df$p_at_1[summary_df$method == m] <- res$precision_at_k$precision_at_k[["P@1"]]
    }

    if (!is.null(res$fdr)) {
      fdr_idx <- which.min(abs(res$fdr$threshold - 0.7))
      summary_df$fdr_at_07[summary_df$method == m] <- res$fdr$fdr[fdr_idx]
    }
  }
  
  summary_df
}

#' Create Ground Truth by InChIKey Connectivity Layer
#'
#' Groups spectra by the first 14 characters of InChIKey (connectivity layer).
#' This effectively groups identical structures (or closely related stereoisomers),
#' and does NOT represent homologous series (e.g., CnH2n+2 progression).
#'
#' @param df_spec Data frame with spectral info including InChIKey.
#' @param min_series_length Minimum compounds in a series.
#' @return Data frame with 'id' and 'group' columns.
#' @export
create_homolog_ground_truth <- function(df_spec, min_series_length = 3) {
  
  if (!"inchikey" %in% names(df_spec)) {
    stop("InChIKey column required for ground truth generation")
  }
  
  # Extract first 14 characters of InChIKey (connectivity layer)
  df_spec$inchikey_prefix <- substr(df_spec$inchikey, 1, 14)
  
  # Group by InChIKey prefix
  groups <- split(df_spec$id, df_spec$inchikey_prefix)
  

  # Filter to groups with multiple members
  valid_groups <- groups[sapply(groups, length) >= min_series_length]
  
  if (length(valid_groups) == 0) {
    message("No InChIKey connectivity groups found with minimum length ", min_series_length)
    return(NULL)
  }
  
  # Create ground truth data frame
  result <- do.call(rbind, lapply(names(valid_groups), function(g) {
    data.frame(
      id = valid_groups[[g]],
      group = g,
      stringsAsFactors = FALSE
    )
  }))
  
  message("Found ", length(valid_groups), " InChIKey groups with ",
          nrow(result), " total compounds (connectivity groups, not curated homolog series)")
  attr(result, "ground_truth_note") <- paste(
    "Groups are defined by the first 14 characters of InChIKey.",
    "They approximate structural identity/connectivity and are not curated homolog-series labels."
  )
  
  result
}

#' Plot Benchmark Comparison
#'
#' Creates a visual comparison of benchmark results.
#'
#' @param benchmark_results Output from run_benchmark().
#' @param output_file Optional path to save the plot.
#' @return A ggplot object.
#' @export
plot_benchmark_comparison <- function(benchmark_results, output_file = NULL) {
  
  summary_df <- benchmark_results$summary
  
  # Reshape for plotting
  summary_long <- tidyr::pivot_longer(
    summary_df,
    cols = -method,
    names_to = "metric",
    values_to = "value"
  )
  
  summary_long <- summary_long[!is.na(summary_long$value), ]
  
  p <- ggplot2::ggplot(summary_long, ggplot2::aes(x = method, y = value, fill = method)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Benchmark Comparison of Distance Methods",
      x = "Method",
      y = "Score"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = 10, height = 8)
    message("Benchmark plot saved to: ", output_file)
  }
  
  p
}

#' Plot FDR Curves for Multiple Methods
#'
#' @param benchmark_results Output from run_benchmark().
#' @param methods Methods to include (default: all).
#' @return A ggplot object, or NULL if no FDR metrics are available.
#' @export
plot_fdr_curves <- function(benchmark_results, methods = NULL) {
  if (is.null(methods)) {
    methods <- setdiff(names(benchmark_results), "summary")
  }

  fdr_list <- list()
  for (m in methods) {
    fdr_df <- benchmark_results[[m]]$metrics$fdr
    if (!is.null(fdr_df)) {
      fdr_df$method <- m
      fdr_list[[m]] <- fdr_df
    }
  }

  if (length(fdr_list) == 0) {
    message("No FDR data available.")
    return(NULL)
  }

  fdr_all <- do.call(rbind, fdr_list)
  rownames(fdr_all) <- NULL

  ggplot2::ggplot(fdr_all, ggplot2::aes(x = threshold, y = fdr, color = method)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      title = "False Discovery Rate by Similarity Threshold",
      x = "Similarity Threshold",
      y = "FDR",
      color = "Method"
    ) +
    ggplot2::theme_minimal()
}


#' Evaluate Noise Robustness of Distance Methods
#'
#' For each noise level, generates perturbed query spectra from the library,
#' computes distances between perturbed queries and original library spectra,
#' and evaluates matching accuracy.
#'
#' @param spectra_result Output of \code{build_spectra_from_msp()} or
#'   equivalent list containing frag_list, loss_list, df_spec.
#' @param methods Character vector of distance methods to evaluate.
#' @param mz_noise_levels Numeric vector of ppm noise levels to test.
#' @param intensity_noise_sd Numeric. Fixed intensity noise for all levels.
#' @param dropout_prob Numeric. Fixed dropout probability for all levels.
#' @param n_replicates Integer. Number of independent noise replicates.
#' @param top_k Integer vector. Evaluate Top-K accuracy for each K.
#' @param params Parameter list from \code{eihrms_default_params()}.
#' @param seed Base random seed for reproducibility.
#' @return A list with \code{results} (per-replicate) and \code{summary}
#'   (averaged over replicates) data frames.
#' @export
evaluate_noise_robustness <- function(spectra_result,
                                       methods = c("cosine", "entropy_weighted",
                                                    "ppm_wasserstein"),
                                       mz_noise_levels = c(0, 5, 10, 15, 20, 30, 50),
                                       intensity_noise_sd = 0.1,
                                       dropout_prob = 0,
                                       n_replicates = 5,
                                       top_k = c(1, 3, 5, 10),
                                       params = eihrms_default_params(),
                                       seed = 42) {
  frag_list <- spectra_result$frag_list
  loss_list <- spectra_result$loss_list
  n <- length(frag_list)
  ids <- names(frag_list)

  message("=== Noise Robustness Evaluation ===")
  message("Spectra: ", n, " | Methods: ", paste(methods, collapse = ", "))
  message("Noise levels (ppm): ", paste(mz_noise_levels, collapse = ", "))
  message("Replicates: ", n_replicates)

  # Collect all results
  all_rows <- list()
  row_idx <- 0

  for (method in methods) {
    params_m <- params
    params_m$distance_method <- method

    for (noise_ppm in mz_noise_levels) {
      for (rep in seq_len(n_replicates)) {
        rep_seed <- seed + (match(method, methods) - 1) * 10000 +
          match(noise_ppm, mz_noise_levels) * 100 + rep

        message(sprintf("  %s | noise=%d ppm | rep %d/%d",
                        method, noise_ppm, rep, n_replicates))

        # Perturb spectra
        perturbed_frag <- perturb_spectra_list(
          frag_list,
          mz_noise_ppm = noise_ppm,
          intensity_noise_sd = intensity_noise_sd,
          dropout_prob = dropout_prob,
          seed = rep_seed
        )
        perturbed_loss <- perturb_spectra_list(
          loss_list,
          mz_noise_ppm = noise_ppm,
          intensity_noise_sd = intensity_noise_sd,
          dropout_prob = dropout_prob,
          seed = rep_seed + 50000
        )

        # Compute query (perturbed) vs library (original) distance matrix
        dist_mat <- compute_distance_matrix_search(
          query_frag_list = perturbed_frag,
          query_loss_list = perturbed_loss,
          lib_frag_list = frag_list,
          lib_loss_list = loss_list,
          params = params_m,
          progress = FALSE
        )

        # Evaluate Top-K retrieval accuracy
        # For each query i, check if original i is in top-K nearest
        topk_acc <- numeric(length(top_k))
        names(topk_acc) <- paste0("top_", top_k)

        for (ki in seq_along(top_k)) {
          k <- top_k[ki]
          hits <- 0
          for (i in seq_len(n)) {
            dists_i <- dist_mat[i, ]
            # Rank library entries by distance (ascending)
            ranked <- order(dists_i)
            top_ids <- ids[ranked[seq_len(min(k, length(ranked)))]]
            if (ids[i] %in% top_ids) hits <- hits + 1
          }
          topk_acc[ki] <- hits / n
        }

        row_idx <- row_idx + 1
        row <- data.frame(
          method = method,
          mz_noise_ppm = noise_ppm,
          replicate = rep,
          stringsAsFactors = FALSE
        )
        for (nm in names(topk_acc)) row[[nm]] <- topk_acc[nm]
        all_rows[[row_idx]] <- row
      }
    }
  }

  results <- do.call(rbind, all_rows)
  rownames(results) <- NULL

  # Compute summary (mean and sd over replicates)
  topk_cols <- paste0("top_", top_k)
  summary_rows <- list()
  si <- 0

  for (method in methods) {
    for (noise_ppm in mz_noise_levels) {
      subset <- results[results$method == method &
                          results$mz_noise_ppm == noise_ppm, ]
      si <- si + 1
      srow <- data.frame(
        method = method,
        mz_noise_ppm = noise_ppm,
        stringsAsFactors = FALSE
      )
      for (tc in topk_cols) {
        srow[[paste0("mean_", tc)]] <- mean(subset[[tc]], na.rm = TRUE)
        srow[[paste0("sd_", tc)]] <- stats::sd(subset[[tc]], na.rm = TRUE)
      }
      summary_rows[[si]] <- srow
    }
  }

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  list(results = results, summary = summary_df)
}
