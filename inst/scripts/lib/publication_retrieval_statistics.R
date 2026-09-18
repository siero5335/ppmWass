# Tie-aware and cluster-aware retrieval statistics used for the publication
# rerun. These helpers intentionally keep optimistic and pessimistic sensitivity
# columns alongside the fractional (random-order expectation) primary values.

prefix_inchikey <- function(x, n = 14L) {
  x <- as.character(x)
  out <- substr(x, 1L, n)
  out[is.na(x) | x == "" | nchar(x) < n] <- NA_character_
  out
}

tie_equal <- function(x, y, tolerance = 0) {
  if (isTRUE(tolerance <= 0)) {
    return(is.finite(x) & is.finite(y) & x == y)
  }
  is.finite(x) & is.finite(y) & abs(x - y) <= tolerance
}

first_relevant_tie_metrics <- function(dists, relevant, ks = c(1L, 5L, 10L),
                                       tolerance = 0) {
  keep <- is.finite(dists)
  dists <- dists[keep]
  relevant <- as.logical(relevant[keep])
  relevant[is.na(relevant)] <- FALSE

  if (!length(dists) || !any(relevant)) {
    out <- list(
      best_distance = NA_real_, n_better = NA_integer_, tie_size = NA_integer_,
      relevant_in_tie = NA_integer_, optimistic_rank = NA_real_,
      pessimistic_rank = NA_real_, expected_rr = NA_real_
    )
    for (k in ks) {
      out[[paste0("top", k, "_optimistic")]] <- NA_real_
      out[[paste0("top", k, "_pessimistic")]] <- NA_real_
      out[[paste0("top", k, "_fractional")]] <- NA_real_
    }
    return(out)
  }

  best_distance <- min(dists[relevant])
  better <- if (tolerance <= 0) {
    dists < best_distance
  } else {
    dists < (best_distance - tolerance)
  }
  tied <- tie_equal(dists, best_distance, tolerance)
  n_better <- sum(better)
  tie_size <- sum(tied)
  relevant_in_tie <- sum(tied & relevant)

  optimistic_rank <- n_better + 1L
  pessimistic_rank <- n_better + (tie_size - relevant_in_tie) + 1L

  # If r relevant items occupy a uniformly random subset of m tied positions,
  # P(first relevant position = j) = choose(m-j, r-1) / choose(m, r).
  max_j <- tie_size - relevant_in_tie + 1L
  j <- seq_len(max_j)
  log_prob <- lchoose(tie_size - j, relevant_in_tie - 1L) -
    lchoose(tie_size, relevant_in_tie)
  prob <- exp(log_prob)
  prob <- prob / sum(prob)
  absolute_rank <- n_better + j
  expected_rr <- sum(prob / absolute_rank)

  out <- list(
    best_distance = best_distance,
    n_better = as.integer(n_better),
    tie_size = as.integer(tie_size),
    relevant_in_tie = as.integer(relevant_in_tie),
    optimistic_rank = as.numeric(optimistic_rank),
    pessimistic_rank = as.numeric(pessimistic_rank),
    expected_rr = expected_rr
  )
  for (k in ks) {
    out[[paste0("top", k, "_optimistic")]] <- as.numeric(optimistic_rank <= k)
    out[[paste0("top", k, "_pessimistic")]] <- as.numeric(pessimistic_rank <= k)
    out[[paste0("top", k, "_fractional")]] <- sum(prob[absolute_rank <= k])
  }
  out
}

random_first_relevant_rank <- function(n_better, tie_size, relevant_in_tie) {
  if (!all(is.finite(c(n_better, tie_size, relevant_in_tie)))) return(NA_real_)
  labels <- c(rep(TRUE, relevant_in_tie), rep(FALSE, tie_size - relevant_in_tie))
  order <- sample(labels, length(labels), replace = FALSE)
  n_better + which(order)[[1]]
}

# Exact expectation of precision and average precision under a uniformly
# random ordering within every exact distance-tie block. This avoids using the
# input column order as an undocumented tie breaker for MAP/P@K.
prefix_ranking_metrics_fractional <- function(dists, relevant,
                                               ks = c(1L, 5L, 10L)) {
  keep <- is.finite(dists)
  dists <- as.numeric(dists[keep])
  relevant <- as.logical(relevant[keep])
  relevant[is.na(relevant)] <- FALSE
  n_relevant <- sum(relevant)
  out <- list(average_precision_fractional = NA_real_)
  for (k in ks) out[[paste0("precision_at_", k, "_fractional")]] <- NA_real_
  if (!length(dists) || !n_relevant) return(out)

  ordered_distances <- sort(unique(dists))
  n_before <- 0L
  relevant_before <- 0L
  expected_ap_numerator <- 0
  expected_relevant_at_k <- stats::setNames(numeric(length(ks)), as.character(ks))

  for (distance_value in ordered_distances) {
    in_block <- dists == distance_value
    block_size <- sum(in_block)
    relevant_in_block <- sum(relevant[in_block])
    if (relevant_in_block > 0L) {
      positions <- seq_len(block_size)
      expected_prior_relevant <- if (block_size == 1L) {
        0
      } else {
        (positions - 1) * (relevant_in_block - 1) / (block_size - 1)
      }
      probability_relevant <- relevant_in_block / block_size
      expected_ap_numerator <- expected_ap_numerator + sum(
        probability_relevant *
          (relevant_before + 1 + expected_prior_relevant) /
          (n_before + positions)
      )
    }
    for (k in ks) {
      slots <- max(0L, min(block_size, as.integer(k) - n_before))
      if (slots > 0L) {
        expected_relevant_at_k[[as.character(k)]] <-
          expected_relevant_at_k[[as.character(k)]] +
          slots * relevant_in_block / block_size
      }
    }
    n_before <- n_before + block_size
    relevant_before <- relevant_before + relevant_in_block
  }

  out$average_precision_fractional <- expected_ap_numerator / n_relevant
  for (k in ks) {
    denominator <- min(as.integer(k), length(dists))
    out[[paste0("precision_at_", k, "_fractional")]] <-
      expected_relevant_at_k[[as.character(k)]] / denominator
  }
  out
}

# Descriptive pair-level metrics use every ordered query-to-library pair. For
# symmetric methods this merely duplicates each unordered pair and leaves AUC
# and FDR proportions unchanged; for directional Composite it preserves the
# implemented query orientation instead of silently selecting one triangle.
binary_rank_auc <- function(distance, is_positive) {
  n_positive <- sum(is_positive)
  n_negative <- sum(!is_positive)
  if (!n_positive || !n_negative) return(NA_real_)

  score_rank <- rank(-distance, ties.method = "average")
  # Coerce counts before multiplication. Large nominal libraries can have
  # n_positive * n_negative > .Machine$integer.max even though the rank-based
  # AUC itself is well within double precision.
  n_positive_double <- as.double(n_positive)
  n_negative_double <- as.double(n_negative)
  (sum(score_rank[is_positive]) -
     n_positive_double * (n_positive_double + 1) / 2) /
    (n_positive_double * n_negative_double)
}

ordered_pair_discrimination_metrics <- function(
    dist_mat, inchikeys, thresholds = seq(0.5, 0.95, by = 0.05)) {
  query_ids <- rownames(dist_mat)
  library_ids <- colnames(dist_mat)
  if (is.null(query_ids) || is.null(library_ids) ||
      !identical(query_ids, library_ids)) {
    stop("ordered pair metrics require a square matrix with identical row/column IDs.")
  }
  if (is.null(names(inchikeys))) names(inchikeys) <- query_ids
  prefixes <- prefix_inchikey(inchikeys[query_ids])
  pair_index <- which(row(dist_mat) != col(dist_mat), arr.ind = TRUE)
  distance <- dist_mat[pair_index]
  query_prefix <- prefixes[pair_index[, 1L]]
  library_prefix <- prefixes[pair_index[, 2L]]
  valid <- is.finite(distance) & !is.na(query_prefix) & !is.na(library_prefix)
  distance <- distance[valid]
  same_prefix <- query_prefix[valid] == library_prefix[valid]
  n_positive <- sum(same_prefix)
  n_negative <- sum(!same_prefix)
  auc <- binary_rank_auc(distance, same_prefix)

  if (!length(distance)) {
    similarity <- numeric()
  } else {
    distance_min <- min(distance)
    distance_max <- max(distance)
    similarity <- if (distance_max > distance_min) {
      1 - (distance - distance_min) / (distance_max - distance_min)
    } else {
      rep(0.5, length(distance))
    }
  }
  fdr <- do.call(rbind, lapply(thresholds, function(threshold) {
    selected <- similarity >= threshold
    n_matches <- sum(selected)
    n_true <- sum(selected & same_prefix)
    n_false <- sum(selected & !same_prefix)
    data.frame(
      threshold = threshold, n_matches = n_matches,
      n_true_positive = n_true, n_false_positive = n_false,
      # With no valid ordered pairs the estimand is undefined. If valid pairs
      # exist but this threshold makes no calls, retain the legacy convention
      # FDR = 0 for an empty selected set.
      fdr = if (!length(distance)) {
        NA_real_
      } else if (n_matches) {
        n_false / n_matches
      } else {
        0
      },
      stringsAsFactors = FALSE
    )
  }))
  list(
    auc = data.frame(
      auc = auc, n_ordered_pairs = length(distance),
      n_positive = n_positive, n_negative = n_negative,
      pair_orientation = "ordered_query_to_library",
      similarity_for_auc = "negative_raw_distance_rank",
      stringsAsFactors = FALSE
    ),
    fdr = fdr
  )
}

per_query_metrics_tie_aware <- function(dist_mat, inchikeys,
                                        tie_tolerance = 0,
                                        random_seed = 20260718L) {
  ids <- rownames(dist_mat)
  if (is.null(names(inchikeys))) names(inchikeys) <- ids
  inchikeys <- as.character(inchikeys[ids])
  prefixes <- prefix_inchikey(inchikeys)

  full_counts <- table(inchikeys[!is.na(inchikeys) & inchikeys != ""])
  prefix_counts <- table(prefixes[!is.na(prefixes) & prefixes != ""])
  full_queries <- which(inchikeys %in% names(full_counts[full_counts > 1L]))
  prefix_queries <- which(prefixes %in% names(prefix_counts[prefix_counts > 1L]))

  set.seed(random_seed)
  full_rows <- lapply(full_queries, function(qi) {
    d <- dist_mat[qi, ]
    d[qi] <- Inf
    relevant <- inchikeys == inchikeys[[qi]]
    relevant[qi] <- FALSE
    metrics <- first_relevant_tie_metrics(
      d, relevant, tolerance = tie_tolerance
    )
    full_top1_tie_block_crossing <-
      metrics$n_better < 1L & metrics$n_better + metrics$tie_size > 1L
    full_top5_tie_block_crossing <-
      metrics$n_better < 5L & metrics$n_better + metrics$tie_size > 5L
    full_top10_tie_block_crossing <-
      metrics$n_better < 10L & metrics$n_better + metrics$tie_size > 10L
    full_top1_affected <-
      metrics$top1_optimistic != metrics$top1_pessimistic
    full_top5_affected <-
      metrics$top5_optimistic != metrics$top5_pessimistic
    full_top10_affected <-
      metrics$top10_optimistic != metrics$top10_pessimistic
    random_rank <- random_first_relevant_rank(
      metrics$n_better, metrics$tie_size, metrics$relevant_in_tie
    )
    data.frame(
      query_id = ids[[qi]],
      top1_optimistic = metrics$top1_optimistic,
      top1_fractional = metrics$top1_fractional,
      top1_pessimistic = metrics$top1_pessimistic,
      top1_random = as.numeric(random_rank == 1),
      top5_optimistic = metrics$top5_optimistic,
      top5_fractional = metrics$top5_fractional,
      top5_pessimistic = metrics$top5_pessimistic,
      top5_random = as.numeric(random_rank <= 5),
      top10_optimistic = metrics$top10_optimistic,
      top10_fractional = metrics$top10_fractional,
      top10_pessimistic = metrics$top10_pessimistic,
      top10_random = as.numeric(random_rank <= 10),
      rr_optimistic = 1 / metrics$optimistic_rank,
      rr_fractional = metrics$expected_rr,
      rr_pessimistic = 1 / metrics$pessimistic_rank,
      rr_random = 1 / random_rank,
      full_best_distance = metrics$best_distance,
      full_n_better = metrics$n_better,
      full_tie_size = metrics$tie_size,
      full_relevant_in_tie = metrics$relevant_in_tie,
      full_has_relevant_tie = metrics$tie_size > 1L,
      full_top1_tie_block_crossing = full_top1_tie_block_crossing,
      full_top1_optimistic_pessimistic_crossing = full_top1_affected,
      full_top1_affected = full_top1_affected,
      full_top5_tie_block_crossing = full_top5_tie_block_crossing,
      full_top5_optimistic_pessimistic_crossing = full_top5_affected,
      full_top5_affected = full_top5_affected,
      full_top10_tie_block_crossing = full_top10_tie_block_crossing,
      full_top10_optimistic_pessimistic_crossing = full_top10_affected,
      full_top10_affected = full_top10_affected,
      n_full_relevant = sum(relevant),
      stringsAsFactors = FALSE
    )
  })
  if (length(full_rows)) {
    full_df <- do.call(rbind, full_rows)
  } else {
    full_columns <- c(
      "top1_optimistic", "top1_fractional", "top1_pessimistic", "top1_random",
      "top5_optimistic", "top5_fractional", "top5_pessimistic", "top5_random",
      "top10_optimistic", "top10_fractional", "top10_pessimistic", "top10_random",
      "rr_optimistic", "rr_fractional", "rr_pessimistic", "rr_random",
      "full_best_distance", "full_n_better", "full_tie_size",
      "full_relevant_in_tie", "full_has_relevant_tie",
      "full_top1_tie_block_crossing",
      "full_top1_optimistic_pessimistic_crossing", "full_top1_affected",
      "full_top5_tie_block_crossing",
      "full_top5_optimistic_pessimistic_crossing", "full_top5_affected",
      "full_top10_tie_block_crossing",
      "full_top10_optimistic_pessimistic_crossing", "full_top10_affected",
      "n_full_relevant"
    )
    full_df <- data.frame(query_id = character(), stringsAsFactors = FALSE)
    for (column in full_columns) full_df[[column]] <- numeric()
    for (column in c(
      "full_has_relevant_tie", "full_top1_tie_block_crossing",
      "full_top1_optimistic_pessimistic_crossing", "full_top1_affected",
      "full_top5_tie_block_crossing",
      "full_top5_optimistic_pessimistic_crossing", "full_top5_affected",
      "full_top10_tie_block_crossing",
      "full_top10_optimistic_pessimistic_crossing", "full_top10_affected"
    )) full_df[[column]] <- logical()
    for (column in c(
      "full_n_better", "full_tie_size", "full_relevant_in_tie",
      "n_full_relevant"
    )) full_df[[column]] <- integer()
  }

  prefix_rows <- lapply(prefix_queries, function(qi) {
    d <- dist_mat[qi, ]
    d[qi] <- Inf
    relevant <- prefixes == prefixes[[qi]]
    relevant[qi] <- FALSE
    metrics <- first_relevant_tie_metrics(
      d, relevant, tolerance = tie_tolerance
    )
    ranking_metrics <- prefix_ranking_metrics_fractional(
      d, relevant, ks = c(1L, 5L, 10L)
    )
    prefix_p_at_1_tie_block_crossing <-
      metrics$n_better < 1L & metrics$n_better + metrics$tie_size > 1L
    prefix_p_at_1_affected <-
      metrics$top1_optimistic != metrics$top1_pessimistic
    random_rank <- random_first_relevant_rank(
      metrics$n_better, metrics$tie_size, metrics$relevant_in_tie
    )
    data.frame(
      query_id = ids[[qi]],
      p_at_1_optimistic = metrics$top1_optimistic,
      p_at_1_fractional = metrics$top1_fractional,
      p_at_1_pessimistic = metrics$top1_pessimistic,
      p_at_1_random = as.numeric(random_rank == 1),
      precision_at_5_fractional =
        ranking_metrics$precision_at_5_fractional,
      precision_at_10_fractional =
        ranking_metrics$precision_at_10_fractional,
      average_precision_fractional =
        ranking_metrics$average_precision_fractional,
      prefix_best_distance = metrics$best_distance,
      prefix_n_better = metrics$n_better,
      prefix_tie_size = metrics$tie_size,
      prefix_relevant_in_tie = metrics$relevant_in_tie,
      prefix_has_relevant_tie = metrics$tie_size > 1L,
      prefix_p_at_1_tie_block_crossing = prefix_p_at_1_tie_block_crossing,
      prefix_p_at_1_optimistic_pessimistic_crossing =
        prefix_p_at_1_affected,
      prefix_p_at_1_affected = prefix_p_at_1_affected,
      n_prefix_relevant = sum(relevant),
      stringsAsFactors = FALSE
    )
  })
  if (length(prefix_rows)) {
    prefix_df <- do.call(rbind, prefix_rows)
  } else {
    prefix_columns <- c(
      "p_at_1_optimistic", "p_at_1_fractional", "p_at_1_pessimistic",
      "p_at_1_random", "precision_at_5_fractional",
      "precision_at_10_fractional", "average_precision_fractional",
      "prefix_best_distance", "prefix_n_better", "prefix_tie_size",
      "prefix_relevant_in_tie", "prefix_has_relevant_tie",
      "prefix_p_at_1_tie_block_crossing",
      "prefix_p_at_1_optimistic_pessimistic_crossing",
      "prefix_p_at_1_affected", "n_prefix_relevant"
    )
    prefix_df <- data.frame(query_id = character(), stringsAsFactors = FALSE)
    for (column in prefix_columns) prefix_df[[column]] <- numeric()
    for (column in c(
      "prefix_has_relevant_tie", "prefix_p_at_1_tie_block_crossing",
      "prefix_p_at_1_optimistic_pessimistic_crossing",
      "prefix_p_at_1_affected"
    )) prefix_df[[column]] <- logical()
    for (column in c(
      "prefix_n_better", "prefix_tie_size", "prefix_relevant_in_tie",
      "n_prefix_relevant"
    )) prefix_df[[column]] <- integer()
  }

  metric_df <- merge(full_df, prefix_df, by = "query_id", all = TRUE)
  metadata <- data.frame(
    query_id = ids,
    inchikey = inchikeys,
    inchikey_prefix = prefixes,
    stringsAsFactors = FALSE
  )
  merge(metadata, metric_df, by = "query_id", all.y = TRUE)
}

bootstrap_mean <- function(values, R = 1000L, seed = 1L) {
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(c(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, n = 0))
  }
  set.seed(seed)
  draws <- replicate(R, mean(sample(values, length(values), replace = TRUE)))
  c(
    estimate = mean(values),
    ci_low = unname(stats::quantile(draws, 0.025, type = 8)),
    ci_high = unname(stats::quantile(draws, 0.975, type = 8)),
    n = length(values)
  )
}

cluster_bootstrap_mean <- function(data, value_col, cluster_col,
                                   R = 5000L, seed = 1L) {
  keep <- is.finite(data[[value_col]]) & !is.na(data[[cluster_col]]) &
    data[[cluster_col]] != ""
  x <- data[keep, c(value_col, cluster_col), drop = FALSE]
  names(x) <- c("value", "cluster")
  clusters <- unique(x$cluster)
  if (!nrow(x) || !length(clusters)) {
    return(c(
      estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      n_queries = 0, n_clusters = 0
    ))
  }

  by_cluster <- split(x$value, x$cluster)
  set.seed(seed)
  draws <- replicate(R, {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    mean(unlist(by_cluster[sampled], use.names = FALSE))
  })
  c(
    estimate = mean(x$value),
    ci_low = unname(stats::quantile(draws, 0.025, type = 8)),
    ci_high = unname(stats::quantile(draws, 0.975, type = 8)),
    n_queries = nrow(x),
    n_clusters = length(clusters)
  )
}

summarize_tie_aware_metrics <- function(per_query, dataset, method,
                                        query_boot_R = 1000L,
                                        cluster_boot_R = 20000L,
                                        seed = 1L) {
  specifications <- data.frame(
    metric = c("top1", "mrr", "p_at_1", "map", "p_at_5", "p_at_10"),
    value_col = c(
      "top1_fractional", "rr_fractional", "p_at_1_fractional",
      "average_precision_fractional", "precision_at_5_fractional",
      "precision_at_10_fractional"
    ),
    cluster_col = c(
      "inchikey", "inchikey", "inchikey_prefix", "inchikey_prefix",
      "inchikey_prefix", "inchikey_prefix"
    ),
    cluster_unit = c(
      "full_inchikey", "full_inchikey", "inchikey_prefix_14",
      "inchikey_prefix_14", "inchikey_prefix_14", "inchikey_prefix_14"
    ),
    stringsAsFactors = FALSE
  )

  rows <- lapply(seq_len(nrow(specifications)), function(i) {
    spec <- specifications[i, ]
    q <- bootstrap_mean(
      per_query[[spec$value_col]], R = query_boot_R, seed = seed + i
    )
    cl <- cluster_bootstrap_mean(
      per_query, spec$value_col, spec$cluster_col,
      R = cluster_boot_R, seed = seed + 100L + i
    )
    data.frame(
      dataset = dataset,
      method = method,
      metric = spec$metric,
      tie_policy = "fractional_expected",
      estimate = cl[["estimate"]],
      cluster_ci_low = cl[["ci_low"]],
      cluster_ci_high = cl[["ci_high"]],
      query_ci_low = q[["ci_low"]],
      query_ci_high = q[["ci_high"]],
      n_queries = as.integer(cl[["n_queries"]]),
      n_clusters = as.integer(cl[["n_clusters"]]),
      analysis_level = "cluster",
      status = "primary",
      inference_method = "cluster_bootstrap",
      cluster_unit = spec$cluster_unit,
      estimand = paste0(
        "query_weighted_mean_conditional_on_fixed_candidate_library; ",
        "clusters_resampled_with_replacement_by_", spec$cluster_unit
      ),
      query_ci_analysis_level = "query",
      query_ci_status = "exploratory",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_tie_sensitivity <- function(per_query, dataset, method) {
  specs <- list(
    top1 = c(optimistic = "top1_optimistic", fractional = "top1_fractional",
             pessimistic = "top1_pessimistic", random = "top1_random"),
    top5 = c(optimistic = "top5_optimistic", fractional = "top5_fractional",
             pessimistic = "top5_pessimistic", random = "top5_random"),
    top10 = c(optimistic = "top10_optimistic", fractional = "top10_fractional",
              pessimistic = "top10_pessimistic", random = "top10_random"),
    mrr = c(optimistic = "rr_optimistic", fractional = "rr_fractional",
            pessimistic = "rr_pessimistic", random = "rr_random"),
    p_at_1 = c(optimistic = "p_at_1_optimistic", fractional = "p_at_1_fractional",
               pessimistic = "p_at_1_pessimistic", random = "p_at_1_random")
  )
  rows <- list()
  k <- 0L
  for (metric in names(specs)) {
    for (policy in names(specs[[metric]])) {
      k <- k + 1L
      value <- per_query[[specs[[metric]][[policy]]]]
      rows[[k]] <- data.frame(
        dataset = dataset,
        method = method,
        metric = metric,
        tie_policy = policy,
        estimate = mean(value, na.rm = TRUE),
        n_queries = sum(is.finite(value)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

tie_diagnostics <- function(per_query, dataset, method) {
  data.frame(
    dataset = dataset,
    method = method,
    query_id = per_query$query_id,
    inchikey = per_query$inchikey,
    inchikey_prefix = per_query$inchikey_prefix,
    full_n_better = per_query$full_n_better,
    full_has_relevant_tie = per_query$full_has_relevant_tie,
    full_rank1_tie = per_query$full_n_better == 0 & per_query$full_tie_size > 1,
    full_top1_tie_block_crossing = per_query$full_top1_tie_block_crossing,
    full_top1_optimistic_pessimistic_crossing =
      per_query$full_top1_optimistic_pessimistic_crossing,
    full_top1_affected = per_query$full_top1_affected,
    full_top5_tie_block_crossing = per_query$full_top5_tie_block_crossing,
    full_top5_optimistic_pessimistic_crossing =
      per_query$full_top5_optimistic_pessimistic_crossing,
    full_top5_affected = per_query$full_top5_affected,
    full_top10_tie_block_crossing = per_query$full_top10_tie_block_crossing,
    full_top10_optimistic_pessimistic_crossing =
      per_query$full_top10_optimistic_pessimistic_crossing,
    full_top10_affected = per_query$full_top10_affected,
    full_tie_size = per_query$full_tie_size,
    full_relevant_in_tie = per_query$full_relevant_in_tie,
    prefix_n_better = per_query$prefix_n_better,
    prefix_has_relevant_tie = per_query$prefix_has_relevant_tie,
    prefix_rank1_tie = per_query$prefix_n_better == 0 & per_query$prefix_tie_size > 1,
    prefix_p_at_1_tie_block_crossing =
      per_query$prefix_p_at_1_tie_block_crossing,
    prefix_p_at_1_optimistic_pessimistic_crossing =
      per_query$prefix_p_at_1_optimistic_pessimistic_crossing,
    prefix_p_at_1_affected = per_query$prefix_p_at_1_affected,
    prefix_tie_size = per_query$prefix_tie_size,
    prefix_relevant_in_tie = per_query$prefix_relevant_in_tie,
    stringsAsFactors = FALSE
  )
}

summarize_tie_diagnostics <- function(per_query, dataset, method) {
  specifications <- data.frame(
    identity_level = c(
      "full_inchikey", "full_inchikey", "full_inchikey",
      "inchikey_prefix_14"
    ),
    metric = c("top1", "top5", "top10", "p_at_1"),
    cutoff = c(1L, 5L, 10L, 1L),
    tie_size_col = c(
      "full_tie_size", "full_tie_size", "full_tie_size",
      "prefix_tie_size"
    ),
    tie_block_crossing_col = c(
      "full_top1_tie_block_crossing",
      "full_top5_tie_block_crossing",
      "full_top10_tie_block_crossing",
      "prefix_p_at_1_tie_block_crossing"
    ),
    affected_col = c(
      "full_top1_affected", "full_top5_affected", "full_top10_affected",
      "prefix_p_at_1_affected"
    ),
    stringsAsFactors = FALSE
  )

  rows <- lapply(seq_len(nrow(specifications)), function(i) {
    spec <- specifications[i, ]
    affected <- as.logical(per_query[[spec$affected_col]])
    tie_block_crossing <-
      as.logical(per_query[[spec$tie_block_crossing_col]])
    tie_size <- per_query[[spec$tie_size_col]]
    evaluated <- !is.na(affected)
    has_relevant_tie <- evaluated & is.finite(tie_size) & tie_size > 1L
    affected_evaluated <- evaluated & affected
    tie_crossing_evaluated <- evaluated & !is.na(tie_block_crossing) &
      tie_block_crossing
    n_queries <- sum(evaluated)
    n_affected <- sum(affected_evaluated)
    n_with_relevant_tie <- sum(has_relevant_tie)
    data.frame(
      dataset = dataset,
      method = method,
      identity_level = spec$identity_level,
      metric = spec$metric,
      cutoff = spec$cutoff,
      n_queries_evaluated = n_queries,
      n_with_relevant_tie = n_with_relevant_tie,
      proportion_with_relevant_tie = if (n_queries) {
        n_with_relevant_tie / n_queries
      } else {
        NA_real_
      },
      n_tie_blocks_crossing_cutoff = sum(tie_crossing_evaluated),
      n_optimistic_pessimistic_crossings = n_affected,
      n_affected = n_affected,
      proportion_affected = if (n_queries) n_affected / n_queries else NA_real_,
      n_tied_not_affected = sum(has_relevant_tie & !affected_evaluated),
      tie_block_crossing_definition =
        "n_better < cutoff and n_better + tie_size > cutoff",
      affected_definition = paste0(
        "optimistic_hit_differs_from_pessimistic_hit_at_cutoff; ",
        "tie_size_greater_than_1_alone_is_not_affected"
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

query_bootstrap_ci_table <- function(cluster_summary) {
  data.frame(
    dataset = cluster_summary$dataset,
    method = cluster_summary$method,
    metric = cluster_summary$metric,
    tie_policy = cluster_summary$tie_policy,
    estimate = cluster_summary$estimate,
    ci_low = cluster_summary$query_ci_low,
    ci_high = cluster_summary$query_ci_high,
    n_queries = cluster_summary$n_queries,
    analysis_level = "query",
    status = "exploratory",
    inference_method = "iid_query_bootstrap",
    resampling_unit = "query",
    estimand = paste0(
      "query_weighted_mean_conditional_on_fixed_candidate_library; ",
      "iid_queries_resampled_with_replacement"
    ),
    stringsAsFactors = FALSE
  )
}

paired_cluster_bootstrap_difference <- function(reference, comparator,
                                                value_col, cluster_col,
                                                R = 20000L, seed = 1L) {
  a <- reference[, c("query_id", value_col, cluster_col), drop = FALSE]
  b <- comparator[, c("query_id", value_col), drop = FALSE]
  names(a) <- c("query_id", "reference", "cluster")
  names(b) <- c("query_id", "comparator")
  x <- merge(a, b, by = "query_id")
  x <- x[is.finite(x$reference) & is.finite(x$comparator) &
           !is.na(x$cluster) & x$cluster != "", , drop = FALSE]
  x$difference <- x$reference - x$comparator
  clusters <- unique(x$cluster)
  if (!nrow(x) || !length(clusters)) {
    return(c(
      estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      n_queries = 0, n_clusters = 0
    ))
  }
  by_cluster <- split(x$difference, x$cluster)
  set.seed(seed)
  draws <- replicate(R, {
    sampled <- sample(clusters, length(clusters), replace = TRUE)
    mean(unlist(by_cluster[sampled], use.names = FALSE))
  })
  c(
    estimate = mean(x$difference),
    ci_low = unname(stats::quantile(draws, 0.025, type = 8)),
    ci_high = unname(stats::quantile(draws, 0.975, type = 8)),
    n_queries = nrow(x),
    n_clusters = length(clusters)
  )
}
