# Cluster-level interpretation and substructure-oriented summaries.
# These helpers turn per-pair distance components into cluster-level labels
# and aggregate typical neutral-loss patterns within each cluster.

#' Extract upper-triangle (i<j) values from a square matrix for a given index set
#'
#' @keywords internal
#' @noRd
extract_pair_values <- function(mat, idx) {
  if (is.null(mat) || length(idx) < 2) return(numeric(0))
  sub <- mat[idx, idx, drop = FALSE]
  sub[upper.tri(sub, diag = FALSE)]
}

#' Compute mean of pairwise minima without enumerating all pairs
#'
#' For values x_1..x_n, average_{i<j} min(x_i, x_j).
#'
#' @keywords internal
#' @noRd
pairwise_min_mean <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x) | x < 0] <- 0
  n <- length(x)
  if (n < 2) return(NA_real_)
  xs <- sort(x, decreasing = FALSE)
  # sum_{i=1}^{n-1} x_i * (n - i)
  w <- rev(seq_len(n) - 1)  # (n-1, n-2, ..., 0)
  sum_min <- sum(xs * w)
  sum_min / (n * (n - 1) / 2)
}

#' Convert a typical-loss spectrum matrix to a named numeric vector over all typical losses
#'
#' @keywords internal
#' @noRd
typical_spec_to_named_vector <- function(typ_spec, formulas = NULL) {
  # Create a full vector over the specified typical-loss set.
  if (is.null(formulas)) {
    formulas <- vapply(TYPICAL_LOSSES, function(z) z$formula, character(1))
  }
  v <- stats::setNames(numeric(length(formulas)), formulas)
  if (is.null(typ_spec) || nrow(typ_spec) == 0) return(v)

  rn <- rownames(typ_spec)
  if (is.null(rn) || length(rn) != nrow(typ_spec)) return(v)
  ints <- as.numeric(typ_spec[, 2])
  ints[!is.finite(ints) | ints < 0] <- 0
  m <- tapply(ints, rn, sum)
  nm <- intersect(names(m), names(v))
  v[nm] <- as.numeric(m[nm])
  v
}

#' Summarize cluster drivers from distance-component matrices
#'
#' Requires `dist_components` as returned by `compute_similarity_matrices()` when
#' `params$return_distance_components = TRUE`.
#'
#' @param clusters Cluster labels (factor/character/numeric). Names should be compound IDs.
#' @param df_spec Data frame with at least columns `id`, `RI`, `known`, `compound_class`.
#' @param dist_components List of component matrices (e.g., c_frag_total, c_loss_total, ...).
#' @param dist_raw Optional full distance matrix (for medoid/cohesion).
#' @param sim_cluster Optional similarity matrix used for clustering (for cohesion summaries).
#' @param exclude_noise Exclude cluster label "0".
#' @param driver_threshold Threshold for declaring a dominant driver.
#' @param subdriver_threshold Threshold for declaring anchored vs pairwise dominance (within loss).
#' @param typical_threshold Threshold for declaring typical-loss dominance (within loss).
#' @return A tibble with one row per cluster.
#' @export
summarize_cluster_drivers <- function(clusters,
                                     df_spec,
                                     dist_components,
                                     dist_raw = NULL,
                                     sim_cluster = NULL,
                                     exclude_noise = TRUE,
                                     driver_threshold = 0.60,
                                     subdriver_threshold = 0.60,
                                     typical_threshold = 0.30) {
  if (is.null(dist_components) || !is.list(dist_components)) {
    stop("dist_components must be a list of matrices (set params$return_distance_components=TRUE).")
  }
  if (is.null(dist_components$c_frag_total) || is.null(dist_components$c_loss_total)) {
    stop("dist_components must include c_frag_total and c_loss_total.")
  }

  # Normalize clusters to a named character vector keyed by compound id.
  if (is.null(names(clusters))) {
    if (!is.null(df_spec$id) && length(clusters) == nrow(df_spec)) {
      names(clusters) <- df_spec$id
    }
  }
  if (is.null(names(clusters))) stop("clusters must be named by compound id (or match df_spec order).")

  cl <- as.character(clusters)
  ids <- names(cl)

  # Exclude noise cluster label "0" if requested.
  keep_ids <- rep(TRUE, length(ids))
  if (isTRUE(exclude_noise)) {
    keep_ids <- keep_ids & !is.na(cl) & cl != "0"
  }
  ids <- ids[keep_ids]
  cl <- cl[keep_ids]

  # Align df_spec
  df_spec <- df_spec %>% dplyr::filter(.data$id %in% ids)
  df_spec <- df_spec %>% dplyr::mutate(.cluster = cl[.data$id])

  # Convenience: component matrices (must be square, with dimnames)
  comp <- dist_components
  mat_ids <- rownames(comp$c_frag_total)
  if (is.null(mat_ids)) stop("dist_components matrices must have rownames (compound ids).")

  # Helper to get indices for ids in the component matrices
  id_to_idx <- stats::setNames(seq_along(mat_ids), mat_ids)
  get_idx <- function(members) unname(id_to_idx[members])

  split_mode <- !is.null(comp$c_anchor_in_loss) && !is.null(comp$c_pair_in_loss)
  typical_mode_combined <- !is.null(comp$c_loss_typ_in_loss)
  typical_mode_split <- !is.null(comp$c_anchor_typ_in_anchor) || !is.null(comp$c_pair_typ_in_pair)

  # Cluster-level summaries
  out <- df_spec %>%
    dplyr::group_by(.data$.cluster) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_known = sum(.data$known == "T", na.rm = TRUE),
      frac_known = ifelse(n > 0, n_known / n, NA_real_),
      RI_min = {
        x <- .data$RI
        x <- x[is.finite(x)]
        if (length(x) == 0) NA_real_ else min(x)
      },
      RI_median = {
        x <- .data$RI
        x <- x[is.finite(x)]
        if (length(x) == 0) NA_real_ else stats::median(x)
      },
      RI_max = {
        x <- .data$RI
        x <- x[is.finite(x)]
        if (length(x) == 0) NA_real_ else max(x)
      },
      top_class = {
        cc <- .data$compound_class
        cc <- cc[!is.na(cc) & nzchar(cc)]
        if (length(cc) == 0) "" else names(sort(table(cc), decreasing = TRUE))[1]
      },
      .groups = "drop"
    )

  # Add component-derived summaries cluster-by-cluster
  out <- out %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pairs = ifelse(.data$n >= 2, .data$n * (.data$n - 1) / 2, 0),
      mean_sim_cluster = {
        if (is.null(sim_cluster) || .data$n < 2) {
          NA_real_
        } else {
          members <- df_spec$id[df_spec$.cluster == .data$.cluster]
          idx <- get_idx(members)
          vals <- extract_pair_values(sim_cluster, idx)
          if (length(vals) == 0) NA_real_ else mean(vals, na.rm = TRUE)
        }
      },
      medoid_id = {
        if (is.null(dist_raw) || .data$n < 2) {
          NA_character_
        } else {
          members <- df_spec$id[df_spec$.cluster == .data$.cluster]
          if (!is.null(rownames(dist_raw))) {
            members <- intersect(members, rownames(dist_raw))
            if (length(members) < 2) {
              NA_character_
            } else {
              dsub <- dist_raw[members, members, drop = FALSE]
              diag(dsub) <- NA_real_
              md <- suppressWarnings(rowMeans(dsub, na.rm = TRUE))
              if (!any(is.finite(md))) NA_character_ else members[which.min(md)]
            }
          } else {
            idx <- get_idx(members)
            ok <- !is.na(idx)
            idx <- idx[ok]
            members <- members[ok]
            if (length(idx) < 2) {
              NA_character_
            } else {
              dsub <- dist_raw[idx, idx, drop = FALSE]
              diag(dsub) <- NA_real_
              md <- suppressWarnings(rowMeans(dsub, na.rm = TRUE))
              if (!any(is.finite(md))) NA_character_ else members[which.min(md)]
            }
          }
        }
      },
      c_frag_mean = {
        members <- df_spec$id[df_spec$.cluster == .data$.cluster]
        idx <- get_idx(members)
        v <- extract_pair_values(comp$c_frag_total, idx)
        if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
      },
      c_loss_mean = {
        members <- df_spec$id[df_spec$.cluster == .data$.cluster]
        idx <- get_idx(members)
        v <- extract_pair_values(comp$c_loss_total, idx)
        if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
      },
      conf_pair_mean = {
        if (is.null(comp$conf_pair)) {
          NA_real_
        } else {
          members <- df_spec$id[df_spec$.cluster == .data$.cluster]
          idx <- get_idx(members)
          v <- extract_pair_values(comp$conf_pair, idx)
          if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
        }
      },
      c_anchor_in_loss_mean = {
        if (!split_mode) {
          NA_real_
        } else {
          members <- df_spec$id[df_spec$.cluster == .data$.cluster]
          idx <- get_idx(members)
          v <- extract_pair_values(comp$c_anchor_in_loss, idx)
          if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
        }
      },
      c_pair_in_loss_mean = {
        if (!split_mode) {
          NA_real_
        } else {
          members <- df_spec$id[df_spec$.cluster == .data$.cluster]
          idx <- get_idx(members)
          v <- extract_pair_values(comp$c_pair_in_loss, idx)
          if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
        }
      },
      c_loss_typ_in_loss_mean = {
        if (!typical_mode_combined) {
          NA_real_
        } else {
          members <- df_spec$id[df_spec$.cluster == .data$.cluster]
          idx <- get_idx(members)
          v <- extract_pair_values(comp$c_loss_typ_in_loss, idx)
          if (length(v) == 0) NA_real_ else mean(v, na.rm = TRUE)
        }
      },
      driver_main = {
        if (!is.finite(.data$c_frag_mean) || !is.finite(.data$c_loss_mean)) {
          NA_character_
        } else if (.data$c_frag_mean >= driver_threshold) {
          "frag"
        } else if (.data$c_loss_mean >= driver_threshold) {
          "loss"
        } else {
          "mixed"
        }
      },
      driver_loss_sub = {
        if (is.na(.data$driver_main) || .data$driver_main != "loss" || !split_mode) {
          NA_character_
        } else {
          a <- .data$c_anchor_in_loss_mean
          b <- .data$c_pair_in_loss_mean
          if (!is.finite(a) || !is.finite(b)) {
            NA_character_
          } else if (a >= subdriver_threshold) {
            "anchored"
          } else if (b >= subdriver_threshold) {
            "pairwise"
          } else {
            "mixed"
          }
        }
      },
      driver_typical_flag = {
        if (!is.finite(.data$c_loss_typ_in_loss_mean)) {
          NA_character_
        } else if (.data$c_loss_typ_in_loss_mean >= typical_threshold) {
          "typical-rich"
        } else {
          "typical-minor"
        }
      },
      label = {
        if (is.na(.data$driver_main)) {
          NA_character_
        } else if (.data$driver_main == "frag") {
          "frag-driven"
        } else if (.data$driver_main == "loss") {
          if (split_mode && is.finite(.data$c_anchor_in_loss_mean) && is.finite(.data$c_pair_in_loss_mean)) {
            base <- paste0("loss-", .data$driver_loss_sub)
          } else {
            base <- "loss"
          }
          if (!is.na(.data$driver_typical_flag) && .data$driver_typical_flag == "typical-rich") {
            paste0(base, "+typical")
          } else {
            base
          }
        } else {
          "mixed"
        }
      }
    ) %>%
    dplyr::ungroup()

  out %>% dplyr::rename(cluster = .data$.cluster)
}

#' Summarize typical neutral losses within clusters
#'
#' Computes presence, mean intensity, and an efficient estimate of average pairwise overlap
#' (mean of min(int_i, int_j) over i<j) for each typical loss within each cluster.
#'
#' @param typ_list Named list of typical-loss spectra (as from build_spectra()$loss_typ_list, etc.).
#' @param clusters Named vector of cluster labels keyed by compound id.
#' @param exclude_noise Exclude cluster label "0".
#' @param use_idf If TRUE, compute a simple IDF-style weight using global presence.
#' @param eps Minimum intensity treated as present.
#' @param typical_losses Typical-loss definition list used to map formulas, classes, and exact masses.
#' @return A tibble with cluster × formula rows.
#' @export
summarize_cluster_typical_losses <- function(typ_list,
                                            clusters,
                                            exclude_noise = TRUE,
                                            use_idf = TRUE,
                                            eps = 0,
                                            typical_losses = TYPICAL_LOSSES) {
  if (is.null(typ_list) || length(typ_list) == 0) {
    return(tibble::tibble())
  }
  if (is.null(names(typ_list))) stop("typ_list must be a named list keyed by compound id.")
  if (is.null(names(clusters))) stop("clusters must be a named vector keyed by compound id.")

  ids <- intersect(names(typ_list), names(clusters))
  if (length(ids) == 0) return(tibble::tibble())
  cl <- as.character(clusters[ids])
  if (isTRUE(exclude_noise)) {
    keep <- !is.na(cl) & cl != "0"
    ids <- ids[keep]
    cl <- cl[keep]
  }
  if (length(ids) == 0) return(tibble::tibble())

  # Precompute mappings
  formulas <- vapply(typical_losses, function(z) z$formula, character(1))
  class_map <- stats::setNames(vapply(typical_losses, function(z) z$class, character(1)), formulas)
  exact_map <- stats::setNames(vapply(typical_losses, function(z) z$exact, numeric(1)), formulas)

  # Build matrix (n_compounds × n_formulas) of typical-loss intensities
  V <- matrix(0, nrow = length(ids), ncol = length(formulas),
              dimnames = list(ids, formulas))
  for (i in seq_along(ids)) {
    v <- typical_spec_to_named_vector(typ_list[[ids[i]]], formulas = formulas)
    V[i, ] <- v[formulas]
  }

  present <- V > eps
  global_presence <- colMeans(present, na.rm = TRUE)
  idf <- rep(1, length(formulas))
  names(idf) <- formulas
  if (isTRUE(use_idf)) {
    # Simple IDF-like weight: higher when global presence is low.
    idf <- -log(pmax(global_presence, 1e-6))
    names(idf) <- formulas
  }

  # Cluster-wise summaries
  df_ids <- tibble::tibble(id = ids, cluster = cl)
  clusters_unique <- sort(unique(cl))

  out_list <- vector("list", length(clusters_unique))
  for (k in seq_along(clusters_unique)) {
    cc <- clusters_unique[k]
    members <- df_ids$id[df_ids$cluster == cc]
    if (length(members) == 0) next
    X <- V[members, , drop = FALSE]
    n <- nrow(X)
presence <- colMeans(X > eps, na.rm = TRUE)
mean_int <- colMeans(X, na.rm = TRUE)

# Dispersion / consistency metrics (include zeros for absent losses)
sd_int <- apply(X, 2, stats::sd)
sd_int[!is.finite(sd_int)] <- NA_real_
cv_int <- sd_int / (mean_int + 1e-12)
cv_int[!is.finite(cv_int)] <- NA_real_
stability <- 1 / (1 + pmax(cv_int, 0))
stability[!is.finite(stability)] <- NA_real_

# "Uniformity" (participation ratio): 1 when evenly distributed across members, small if concentrated in few.
s1 <- colSums(X, na.rm = TRUE)
s2 <- colSums(X^2, na.rm = TRUE)
uniformity <- ifelse(s2 > 0, (s1^2) / (n * s2), NA_real_)
uniformity[!is.finite(uniformity)] <- NA_real_

# Efficient overlap: mean_{i<j} min(x_i, x_j)
overlap <- apply(X, 2, pairwise_min_mean)
overlap[!is.finite(overlap)] <- NA_real_

score <- overlap
score_idf <- overlap * idf

# Consistency-weighted scores (separate axis from IDF)
score_consistency <- overlap * stability
score_consistency_idf <- score_consistency * idf

out_list[[k]] <- tibble::tibble(
  cluster = cc,
  n_members = n,
  formula = formulas,
  class = unname(class_map[formulas]),
  exact = unname(exact_map[formulas]),
  presence = as.numeric(presence),
  mean_intensity = as.numeric(mean_int),
  sd_intensity = as.numeric(sd_int),
  cv_intensity = as.numeric(cv_int),
  stability = as.numeric(stability),
  uniformity = as.numeric(uniformity),
  pairwise_overlap = as.numeric(overlap),
  idf = as.numeric(idf[formulas]),
  score = as.numeric(score),
  score_idf = as.numeric(score_idf),
  score_consistency = as.numeric(score_consistency),
  score_consistency_idf = as.numeric(score_consistency_idf)
)
  }

  dplyr::bind_rows(out_list) %>%
    dplyr::arrange(.data$cluster, dplyr::desc(.data$score_consistency_idf), dplyr::desc(.data$score_idf), dplyr::desc(.data$score))
}

#' Annotate clusters with driver labels and typical-loss motifs
#'
#' This is a convenience wrapper that produces (i) membership table,
#' (ii) driver summary per cluster, and (iii) typical-loss motif tables.
#'
#' @param spectra Result of build_spectra().
#' @param clusters Cluster vector (e.g., result of cluster_optics()$clusters). Names should be IDs.
#' @param similarity Result of compute_similarity_matrices().
#' @param params Parameter list.
#' @param top_n_losses Number of top typical losses to include as a semicolon-separated label.
#' @param exclude_noise Exclude cluster label "0".
#' @param use_idf Use IDF weighting for typical-loss motif ranking.
#' @return A list with `membership`, `cluster_summary`, `typical_loss_long`.
#' @export
annotate_clusters <- function(spectra,
                             clusters,
                             similarity,
                             params = eihrms_default_params(),
                             top_n_losses = 5,
                             exclude_noise = TRUE,
                             use_idf = TRUE) {
  params <- validate_params(params)
  typical_losses_universe <- get_typical_losses_universe(params)
  if (is.null(spectra$df_spec) || is.null(spectra$df_spec$id)) {
    stop("spectra must be the result of build_spectra() (needs df_spec with id).")
  }
  if (is.null(names(clusters))) {
    # Try aligning by df_spec order
    if (length(clusters) == nrow(spectra$df_spec)) {
      names(clusters) <- spectra$df_spec$id
    }
  }
  if (is.null(names(clusters))) stop("clusters must be a named vector keyed by compound id.")

  membership <- spectra$df_spec %>%
    dplyr::mutate(cluster = as.character(clusters[.data$id]))

  # Driver summary requires dist_components
  cluster_summary <- NULL
  if (!is.null(similarity$dist_components)) {
    cluster_summary <- summarize_cluster_drivers(
      clusters = clusters,
      df_spec = spectra$df_spec,
      dist_components = similarity$dist_components,
      dist_raw = similarity$dist_raw,
      sim_cluster = similarity$sim_cluster,
      exclude_noise = exclude_noise
    )
  }

  # Typical-loss motifs (combined / anchor / pair)
  typ_long <- tibble::tibble()
  add_channel <- function(df, channel) {
    if (nrow(df) == 0) return(df)
    df$channel <- channel
    df
  }

  if (!is.null(spectra$loss_typ_list)) {
    df <- summarize_cluster_typical_losses(
      spectra$loss_typ_list,
      clusters = clusters,
      exclude_noise = exclude_noise,
      use_idf = use_idf,
      typical_losses = typical_losses_universe
    )
    typ_long <- dplyr::bind_rows(typ_long, add_channel(df, "combined"))
  }
  if (!is.null(spectra$loss_anchor_typ_list)) {
    df <- summarize_cluster_typical_losses(
      spectra$loss_anchor_typ_list,
      clusters = clusters,
      exclude_noise = exclude_noise,
      use_idf = use_idf,
      typical_losses = typical_losses_universe
    )
    typ_long <- dplyr::bind_rows(typ_long, add_channel(df, "anchor"))
  }
  if (!is.null(spectra$loss_pair_typ_list)) {
    df <- summarize_cluster_typical_losses(
      spectra$loss_pair_typ_list,
      clusters = clusters,
      exclude_noise = exclude_noise,
      use_idf = use_idf,
      typical_losses = typical_losses_universe
    )
    typ_long <- dplyr::bind_rows(typ_long, add_channel(df, "pair"))
  }

  # Add compact top-loss labels into cluster_summary (if available)
  if (!is.null(cluster_summary) && nrow(typ_long) > 0 && is.finite(top_n_losses) && top_n_losses > 0) {
      top_labels_idf <- typ_long %>%
    dplyr::group_by(.data$channel, .data$cluster) %>%
    dplyr::slice_max(order_by = .data$score_idf, n = top_n_losses, with_ties = FALSE) %>%
    dplyr::summarise(top_typical_losses = paste(.data$formula, collapse = ";"), .groups = "drop")

  top_labels_cons <- typ_long %>%
    dplyr::group_by(.data$channel, .data$cluster) %>%
    dplyr::slice_max(order_by = .data$score_consistency_idf, n = top_n_losses, with_ties = FALSE) %>%
    dplyr::summarise(top_typical_losses_consistent = paste(.data$formula, collapse = ";"), .groups = "drop")

  # Prefer anchor labels if split-loss is enabled; otherwise combined.
  pref_channel <- if (isTRUE(params$use_split_loss)) "anchor" else "combined"
  lab_pref <- top_labels_idf %>% dplyr::filter(.data$channel == pref_channel) %>%
    dplyr::select(.data$cluster, top_typical_losses)

  lab_pref_cons <- top_labels_cons %>% dplyr::filter(.data$channel == pref_channel) %>%
    dplyr::select(.data$cluster, top_typical_losses_consistent)

  cluster_summary <- cluster_summary %>%
    dplyr::left_join(lab_pref, by = "cluster") %>%
    dplyr::left_join(lab_pref_cons, by = "cluster")
}

  list(
    membership = membership,
    cluster_summary = cluster_summary,
    typical_loss_long = typ_long
  )
}
