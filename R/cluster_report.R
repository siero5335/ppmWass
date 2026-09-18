#' Cluster explanation report generation
#'
#' Automatically extracts "medoid vs top-similar within cluster" pairs and ties
#' together `distance_breakdown_pair()` and `typical_loss_match_table()` to generate
#' human-readable cluster explanation reports (Markdown) plus machine-readable CSV.
#'
#' @keywords internal
#' @noRd

#' @keywords internal
#' @noRd
sanitize_filename <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

#' @keywords internal
#' @noRd
compute_medoid_id <- function(dist_raw, members) {
  if (is.null(dist_raw) || length(members) < 2) return(NA_character_)
  members <- intersect(members, rownames(dist_raw))
  if (length(members) < 2) return(NA_character_)
  dsub <- dist_raw[members, members, drop = FALSE]
  diag(dsub) <- NA_real_
  md <- suppressWarnings(rowMeans(dsub, na.rm = TRUE))
  members[which.min(md)]
}

#' @keywords internal
#' @noRd
markdown_table <- function(df, digits = 3, max_rows = 10) {
  if (is.null(df) || nrow(df) == 0) return("_none_")
  if (nrow(df) > max_rows) df <- df[seq_len(max_rows), , drop = FALSE]

  fmt <- function(x) {
    if (is.numeric(x)) return(format(round(x, digits), trim = TRUE))
    as.character(x)
  }

  cols <- colnames(df)
  header <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", length(cols)), collapse = " | "), " |")
  rows <- apply(df, 1, function(r) paste0("| ", paste(mapply(fmt, r), collapse = " | "), " |"))

  paste(c(header, sep, rows), collapse = "\n")
}

#' @keywords internal
#' @noRd
select_medoid_partners <- function(sim_mat, members, medoid_id, top_k = 5, min_similarity = 0) {
  if (is.null(sim_mat) || is.na(medoid_id) || !medoid_id %in% rownames(sim_mat)) return(character(0))
  members <- intersect(members, colnames(sim_mat))
  members <- setdiff(members, medoid_id)
  if (length(members) == 0) return(character(0))

  s <- sim_mat[medoid_id, members]
  s[!is.finite(s)] <- -Inf
  ord <- order(s, decreasing = TRUE)
  members <- members[ord]
  s <- s[ord]

  keep <- which(s >= min_similarity)
  if (length(keep) == 0) return(character(0))

  members[keep][seq_len(min(top_k, length(keep)))]
}

#' @keywords internal
#' @noRd
select_medoid_worst_within <- function(sim_mat, members, medoid_id) {
  if (is.null(sim_mat) || is.na(medoid_id) || !medoid_id %in% rownames(sim_mat)) return(NA_character_)
  members <- intersect(members, colnames(sim_mat))
  members <- setdiff(members, medoid_id)
  if (length(members) == 0) return(NA_character_)
  s <- sim_mat[medoid_id, members]
  s[!is.finite(s)] <- Inf
  members[which.min(s)]
}

#' @keywords internal
#' @noRd
select_medoid_top_outside <- function(sim_mat, members, medoid_id, top_k = 3, min_similarity = 0) {
  if (is.null(sim_mat) || is.na(medoid_id) || !medoid_id %in% rownames(sim_mat)) return(character(0))
  all_ids <- colnames(sim_mat)
  if (is.null(all_ids)) return(character(0))
  outside <- setdiff(all_ids, members)
  outside <- setdiff(outside, medoid_id)
  if (length(outside) == 0) return(character(0))
  s <- sim_mat[medoid_id, outside]
  s[!is.finite(s)] <- -Inf
  ord <- order(s, decreasing = TRUE)
  outside <- outside[ord]
  s <- s[ord]
  keep <- which(s >= min_similarity)
  if (length(keep) == 0) return(character(0))
  outside[keep][seq_len(min(top_k, length(keep)))]
}

#' @keywords internal
#' @noRd
build_medoid_partner_table <- function(sim_mat,
                                       members,
                                       medoid_id,
                                       clusters = NULL,
                                       top_pairs = 5,
                                       min_similarity = 0,
                                       include_within_worst = TRUE,
                                       outside_top_k = 3,
                                       outside_min_similarity = NULL,
                                       include_outside = TRUE) {
  if (is.null(outside_min_similarity)) outside_min_similarity <- min_similarity

  within_top <- select_medoid_partners(sim_mat, members, medoid_id, top_k = top_pairs, min_similarity = min_similarity)
  out <- data.frame(
    partner = within_top,
    pair_type = rep("within_top", length(within_top)),
    similarity = if (length(within_top) > 0) as.numeric(sim_mat[medoid_id, within_top]) else numeric(0),
    stringsAsFactors = FALSE
  )

  if (isTRUE(include_within_worst)) {
    worst <- select_medoid_worst_within(sim_mat, members, medoid_id)
    if (!is.na(worst) && nzchar(worst) && !(worst %in% out$partner)) {
      out <- rbind(out, data.frame(
        partner = worst,
        pair_type = "within_worst",
        similarity = as.numeric(sim_mat[medoid_id, worst]),
        stringsAsFactors = FALSE
      ))
    }
  }

  if (isTRUE(include_outside) && is.finite(outside_top_k) && outside_top_k > 0) {
    outs <- select_medoid_top_outside(sim_mat, members, medoid_id, top_k = outside_top_k, min_similarity = outside_min_similarity)
    if (length(outs) > 0) {
      add <- data.frame(
        partner = outs,
        pair_type = rep("outside_top", length(outs)),
        similarity = as.numeric(sim_mat[medoid_id, outs]),
        stringsAsFactors = FALSE
      )
      # do not duplicate partners (shouldn't happen, but safe)
      add <- add[!(add$partner %in% out$partner), , drop = FALSE]
      out <- rbind(out, add)
    }
  }

  # Add partner cluster labels if available
  if (!is.null(clusters)) {
    if (!is.null(names(clusters))) {
      out$partner_cluster <- as.character(clusters[out$partner])
    } else {
      out$partner_cluster <- NA_character_
    }
  } else {
    out$partner_cluster <- NA_character_
  }

  out
}


#' Generate cluster explanation reports (Markdown + CSV)
#'
#' @param spectra Output of build_spectra().
#' @param similarity Output of compute_similarity_matrices().
#' @param clusters Cluster labels (named by compound id).
#' @param params Parameter list.
#' @param cluster_annotation Optional output from annotate_clusters(). If provided, medoid and labels are reused.
#' @param output_dir Base output directory.
#' @param report_subdir Subdirectory for markdown reports.
#' @param top_pairs Top medoid-partner pairs per cluster.
#' @param top_losses Top typical-loss matches per pair.
#' @param min_similarity Minimum similarity for selecting partners.
#' @param exclude_noise Exclude cluster label "0".
#' @param overwrite Overwrite existing markdown files.
#' @param include_boundary_pairs Include high-similarity neighbors outside the cluster in the CSV summary.
#' @param include_within_worst Include low-similarity within-cluster pairs in the CSV summary.
#' @param outside_top_k Number of outside-cluster neighbors to keep per medoid.
#' @param outside_min_similarity Minimum similarity for outside-cluster neighbors.
#' @param include_cluster_motifs If TRUE, include cluster-level typical-loss motifs in reports.
#' @param motif_top_n Number of cluster motifs to show per report.
#' @param motif_rank_metric Ranking metric used for motif selection.
#' @param motif_channel Optional typical-loss channel to use for motifs.
#' @param include_pair_plots If TRUE, render per-pair spectrum plots.
#' @param plots_subdir Subdirectory for generated pair plots.
#' @param plot_channels Channels to include in pair plots.
#' @param plot_as_panel If TRUE, combine selected channels into a panel plot.
#' @param plot_include_contribution If TRUE, include contribution bars in pair plots.
#' @param plot_panel_height Optional panel height override.
#' @param plot_use_aligned_bins If TRUE, plot aligned spectra bins.
#' @param plot_normalize Intensity normalization used for plots.
#' @param plot_scale Intensity scaling factor for plots.
#' @param plot_label_top_n Number of top peaks to label in plots.
#' @param plot_label_digits Number of digits for peak labels.
#' @param plot_width Plot width in inches.
#' @param plot_height Plot height in inches.
#' @param plot_dpi Plot DPI.
#' @return A list with `pairs`, `matches`, and `index_file`.
#' @export
generate_cluster_explanation_reports <- function(spectra,
                                                similarity,
                                                clusters,
                                                params = eihrms_default_params(),
                                                cluster_annotation = NULL,
                                                output_dir = ".",
                                                report_subdir = "cluster_reports",
                                                top_pairs = 5,
                                                top_losses = 10,
                                                min_similarity = 0,
                                                exclude_noise = TRUE,
                                                overwrite = TRUE,
                                                include_boundary_pairs = TRUE,
                                                include_within_worst = TRUE,
                                                outside_top_k = 3,
                                                outside_min_similarity = NULL,
                                                include_cluster_motifs = TRUE,
                                                motif_top_n = 10,
                                                motif_rank_metric = c("score_consistency_idf", "score_idf"),
                                                motif_channel = NULL,
                                                include_pair_plots = TRUE,
                                                plots_subdir = "plots",
                                                plot_channels = c("frag", "loss"),
                                                plot_as_panel = TRUE,
                                                plot_include_contribution = TRUE,
                                                plot_panel_height = NULL,
                                                plot_use_aligned_bins = TRUE,
                                                plot_normalize = c("max", "sum", "none"),
                                                plot_scale = 100,
                                                plot_label_top_n = 8,
                                                plot_label_digits = 4,
                                                plot_width = 10,
                                                plot_height = 4,
                                                plot_dpi = 180) {
  params <- validate_params(params)
  if (is.null(spectra$df_spec)) stop("spectra must be the result of build_spectra().")
  if (is.null(similarity$dist_raw) || is.null(similarity$sim_cluster)) {
    stop("similarity must be the result of compute_similarity_matrices() (needs dist_raw and sim_cluster).")
  }

  # Ensure clusters are named by id.
  if (is.null(names(clusters))) {
    if (!is.null(spectra$df_spec$id) && length(clusters) == nrow(spectra$df_spec)) {
      names(clusters) <- spectra$df_spec$id
    }
  }
  if (is.null(names(clusters))) stop("clusters must be named by compound id (or match df_spec order).")

  cl <- as.character(clusters)
  ids <- names(cl)

  sim_mat <- similarity$sim_cluster
  dist_raw <- similarity$dist_raw

  report_dir <- file.path(output_dir, report_subdir)
  if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  plot_normalize <- match.arg(plot_normalize)
  if (is.null(plot_channels) || length(plot_channels) == 0) plot_channels <- c("frag", "loss")
  plot_channels <- unique(as.character(plot_channels))

  # Allow automatic channel selection based on params
  if (length(plot_channels) == 1 && tolower(plot_channels[1]) == "auto" || any(tolower(plot_channels) == "auto")) {
    plot_channels <- c("frag", if (isTRUE(params$use_split_loss)) c("loss_anchor", "loss_pair") else "loss")
  }
  plot_channels <- unique(as.character(plot_channels))

  plots_dir <- NULL
  if (isTRUE(include_pair_plots)) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plots_dir <- file.path(report_dir, plots_subdir)
      if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

      if (isTRUE(plot_as_panel) && !requireNamespace("gridExtra", quietly = TRUE)) {
        warning("plot_as_panel=TRUE but package 'gridExtra' is not installed; falling back to per-channel plots.")
        plot_as_panel <- FALSE
      }
    } else {
      warning("include_pair_plots=TRUE but package 'ggplot2' is not installed; skipping pair plot images.")
    }
  }


  cluster_summary <- NULL
  if (!is.null(cluster_annotation) && !is.null(cluster_annotation$cluster_summary)) {
    cluster_summary <- cluster_annotation$cluster_summary
  }

  clusters_unique <- sort(unique(cl))
  if (isTRUE(exclude_noise)) clusters_unique <- setdiff(clusters_unique, "0")

  pairs_out <- list()
  matches_out <- list()
  report_files <- character(0)

  df_spec <- spectra$df_spec
  rownames(df_spec) <- df_spec$id

  for (cluster_id in clusters_unique) {
    members <- ids[cl == cluster_id]
    members <- members[members %in% rownames(sim_mat)]
    if (length(members) < 2) next

    # Medoid + label
    medoid_id <- NA_character_
    label <- NA_character_
    top_losses_label <- ""

    if (!is.null(cluster_summary)) {
      row <- cluster_summary[as.character(cluster_summary$cluster) == as.character(cluster_id), , drop = FALSE]
      if (nrow(row) > 0) {
        medoid_id <- row$medoid_id[1]
        label <- row$label[1]
        if (!is.null(row$top_typical_losses_consistent) &&
    !is.na(row$top_typical_losses_consistent[1]) &&
    nzchar(as.character(row$top_typical_losses_consistent[1]))) {
  top_losses_label <- as.character(row$top_typical_losses_consistent[1])
} else if (!is.null(row$top_typical_losses)) {
  top_losses_label <- as.character(row$top_typical_losses[1])
}
      }
    }
    if (is.na(medoid_id) || !nzchar(medoid_id)) {
      medoid_id <- compute_medoid_id(dist_raw, members)
    }
    if (is.na(label) || !nzchar(label)) label <- paste0("cluster ", cluster_id)

partner_tbl <- build_medoid_partner_table(
  sim_mat = sim_mat,
  members = members,
  medoid_id = medoid_id,
  clusters = clusters,
  top_pairs = top_pairs,
  min_similarity = min_similarity,
  include_within_worst = include_within_worst,
  outside_top_k = outside_top_k,
  outside_min_similarity = outside_min_similarity,
  include_outside = include_boundary_pairs
)
if (nrow(partner_tbl) == 0) next

partners <- partner_tbl$partner

det_list <- vector("list", nrow(partner_tbl))
names(det_list) <- partner_tbl$partner

# Pair-level breakdowns
pair_rows <- vector("list", nrow(partner_tbl))

    for (k in seq_len(nrow(partner_tbl))) {
      p <- partner_tbl$partner[k]
      pair_type <- partner_tbl$pair_type[k]
      sim_val <- partner_tbl$similarity[k]
      partner_cluster <- partner_tbl$partner_cluster[k]
      det <- distance_breakdown_pair(medoid_id, p, spectra, params, include_typical = TRUE, top_n = top_losses)
      det_list[[p]] <- det

      row <- data.frame(
        cluster = as.character(cluster_id),
        label = as.character(label),
        medoid = medoid_id,
        partner = p,
        pair_type = as.character(pair_type),
        partner_cluster = as.character(partner_cluster),
        similarity = as.numeric(sim_val),
        d_total = as.numeric(det$d_total),
        d_frag = as.numeric(det$d_frag),
        d_loss = as.numeric(det$d_loss),
        c_frag_total = if (!is.null(det$c_frag_total)) as.numeric(det$c_frag_total) else NA_real_,
        c_loss_total = if (!is.null(det$c_loss_total)) as.numeric(det$c_loss_total) else NA_real_,
        deriv_medoid = if (medoid_id %in% rownames(df_spec)) as.character(df_spec[medoid_id, "derivatization_type"]) else NA_character_,
        deriv_partner = if (p %in% rownames(df_spec)) as.character(df_spec[p, "derivatization_type"]) else NA_character_,
        stringsAsFactors = FALSE
      )

      if (!is.null(det$c_anchor_total)) {
        row$c_anchor_total <- as.numeric(det$c_anchor_total)
        row$c_pair_total <- as.numeric(det$c_pair_total)
        row$w_anchor_final <- if (!is.null(det$w_anchor_final)) as.numeric(det$w_anchor_final) else NA_real_
        row$w_pair_final <- if (!is.null(det$w_pair_final)) as.numeric(det$w_pair_final) else NA_real_
      }
      if (!is.null(det$conf_pair)) {
        row$conf_pair <- as.numeric(det$conf_pair)
      } else if (!is.null(det$confA) && !is.null(det$confB)) {
        a <- ifelse(is.finite(det$confA), det$confA, 0)
        b <- ifelse(is.finite(det$confB), det$confB, 0)
        row$conf_pair <- sqrt(a * b)
      }

      pair_rows[[k]] <- row

      # Typical-loss overlap tables
      if (!is.null(det$typical_loss_matches) && nrow(det$typical_loss_matches) > 0) {
        mt <- det$typical_loss_matches
        mt$cluster <- as.character(cluster_id)
        mt$medoid <- medoid_id
        mt$partner <- p
        mt$pair_type <- pair_type
        mt$partner_cluster <- partner_cluster
        mt$channel <- "combined"
        matches_out[[length(matches_out) + 1]] <- mt
      }
      if (!is.null(det$typical_anchor_matches) && nrow(det$typical_anchor_matches) > 0) {
        mt <- det$typical_anchor_matches
        mt$cluster <- as.character(cluster_id)
        mt$medoid <- medoid_id
        mt$partner <- p
        mt$pair_type <- pair_type
        mt$partner_cluster <- partner_cluster
        mt$channel <- "anchor"
        matches_out[[length(matches_out) + 1]] <- mt
      }
      if (!is.null(det$typical_pair_matches) && nrow(det$typical_pair_matches) > 0) {
        mt <- det$typical_pair_matches
        mt$cluster <- as.character(cluster_id)
        mt$medoid <- medoid_id
        mt$partner <- p
        mt$pair_type <- pair_type
        mt$partner_cluster <- partner_cluster
        mt$channel <- "pair"
        matches_out[[length(matches_out) + 1]] <- mt
      }
    }

    pairs_df <- do.call(rbind, pair_rows)
    pairs_out[[length(pairs_out) + 1]] <- pairs_df

    # Markdown report
    fname <- paste0("cluster_", sanitize_filename(cluster_id), "_", sanitize_filename(label), ".md")
    fpath <- file.path(report_dir, fname)
    if (file.exists(fpath) && !isTRUE(overwrite)) {
      report_files <- c(report_files, fpath)
      next
    }

    members_df <- df_spec[members, , drop = FALSE]
    n <- nrow(members_df)
    frac_known <- ifelse(n > 0, mean(members_df$known == "T", na.rm = TRUE), NA_real_)
    ri_min <- suppressWarnings(min(members_df$RI, na.rm = TRUE))
    ri_max <- suppressWarnings(max(members_df$RI, na.rm = TRUE))

    # Split pair summary tables
    pairs_within <- pairs_df[pairs_df$pair_type == "within_top", , drop = FALSE]
    pairs_boundary <- pairs_df[pairs_df$pair_type != "within_top", , drop = FALSE]

    # Cluster motif summary (typical-loss) if available
    motif_tbl <- NULL
    metric <- match.arg(motif_rank_metric)
    motif_channel_use <- motif_channel
    if (is.null(motif_channel_use) || !nzchar(motif_channel_use)) {
      motif_channel_use <- if (isTRUE(params$use_split_loss)) "anchor" else "combined"
    }
    if (isTRUE(include_cluster_motifs) &&
        !is.null(cluster_annotation) &&
        !is.null(cluster_annotation$typical_loss_long) &&
        is.data.frame(cluster_annotation$typical_loss_long) &&
        nrow(cluster_annotation$typical_loss_long) > 0) {

      tl <- cluster_annotation$typical_loss_long
      sub <- tl[as.character(tl$cluster) == as.character(cluster_id) &
                  as.character(tl$channel) == as.character(motif_channel_use), , drop = FALSE]
      if (nrow(sub) > 0) {
        if (!metric %in% colnames(sub)) metric <- "score_idf"
        # Order by selected metric then fall back to score_idf/score for stability
        if ("score_idf" %in% colnames(sub) && "score" %in% colnames(sub)) {
          sub <- sub[order(-sub[[metric]], -sub$score_idf, -sub$score), , drop = FALSE]
        } else {
          sub <- sub[order(-sub[[metric]]), , drop = FALSE]
        }
        keep_cols <- intersect(
          c("formula", "class", "exact", "presence", "mean_intensity", "pairwise_overlap",
            "score_consistency_idf", "score_idf", "idf"),
          colnames(sub)
        )
        motif_tbl <- sub[seq_len(min(motif_top_n, nrow(sub))), keep_cols, drop = FALSE]
      }
    }

    lines <- c(
      paste0("# Cluster ", cluster_id, " - ", label),
      "",
      paste0("- n: ", n),
      paste0("- known fraction: ", ifelse(is.finite(frac_known), round(frac_known, 3), "NA")),
      paste0("- RI range: ", ifelse(is.finite(ri_min), round(ri_min, 1), "NA"), " - ", ifelse(is.finite(ri_max), round(ri_max, 1), "NA")),
      paste0("- medoid: `", medoid_id, "`"),
      if (nzchar(top_losses_label)) paste0("- top typical losses: ", top_losses_label) else NULL,
      if (!is.null(motif_tbl) && nrow(motif_tbl) > 0) paste0("- motif channel: ", motif_channel_use, " (rank: ", metric, ")") else NULL,
      ""
    )

    if (!is.null(motif_tbl) && nrow(motif_tbl) > 0) {
      lines <- c(lines, "## Cluster typical-loss motifs", "", markdown_table(motif_tbl, digits = 4, max_rows = nrow(motif_tbl)), "")
    }

    lines <- c(lines, "## Medoid vs top-similar pairs", "", markdown_table(pairs_within, digits = 3, max_rows = nrow(pairs_within)), "")

    if ((isTRUE(include_boundary_pairs) || isTRUE(include_within_worst)) && nrow(pairs_boundary) > 0) {
      lines <- c(lines, "## Boundary pairs", "", markdown_table(pairs_boundary, digits = 3, max_rows = nrow(pairs_boundary)), "")
    }

    lines <- c(lines, "## Pair explanations", "")

    type_label <- function(x) {
      if (x == "within_top") return("within (top-similar)")
      if (x == "within_worst") return("within (least similar)")
      if (x == "outside_top") return("outside (nearest neighbor)")
      x
    }

    for (k in seq_len(nrow(partner_tbl))) {
      p <- partner_tbl$partner[k]
      pair_type <- partner_tbl$pair_type[k]
      partner_cluster <- partner_tbl$partner_cluster[k]
      det <- det_list[[p]]
      s <- partner_tbl$similarity[k]

      
      plot_md <- character(0)
      if (!is.null(plots_dir) && length(plot_channels) > 0) {
        plot_rel <- character(0)
        # Prefer a single panel plot (mirror plots + contribution bar)
        if (isTRUE(plot_as_panel) && requireNamespace("gridExtra", quietly = TRUE)) {
          ch_tag <- sanitize_filename(paste(plot_channels, collapse = "_"))
          fn_img <- paste0(
            "cluster_", sanitize_filename(cluster_id), "_",
            sanitize_filename(medoid_id), "_vs_", sanitize_filename(p), "_",
            sanitize_filename(pair_type), "_panel_", ch_tag, ".png"
          )
          out_img <- file.path(plots_dir, fn_img)
          rel_img <- paste0(plots_subdir, "/", fn_img)

          ok <- TRUE
          if (!file.exists(out_img) || isTRUE(overwrite)) {
            ok <- FALSE
            tryCatch({
              eff_rows <- length(plot_channels) + ifelse(isTRUE(plot_include_contribution), 0.9, 0)
              h_panel <- if (!is.null(plot_panel_height) && is.finite(plot_panel_height) && plot_panel_height > 0) {
                plot_panel_height
              } else {
                plot_height * eff_rows
              }

              plt <- plot_pair_summary_panel(
                spectra = spectra,
                idA = medoid_id,
                idB = p,
                det = det,
                channels = plot_channels,
                include_contribution = plot_include_contribution,
                params = params,
                use_aligned_bins = plot_use_aligned_bins,
                normalize = plot_normalize,
                scale = plot_scale,
                label_top_n = plot_label_top_n,
                label_digits = plot_label_digits
              )
              ggplot2::ggsave(
                filename = out_img,
                plot = plt,
                width = plot_width,
                height = h_panel,
                dpi = plot_dpi
              )
              ok <- TRUE
            }, error = function(e) {
              ok <- FALSE
            })
          }

          if (isTRUE(ok) && file.exists(out_img)) plot_rel <- c(plot_rel, rel_img)
        } else {
          # Fallback: per-channel plots + separate contribution bar
          for (ch in plot_channels) {
            fn_img <- paste0(
              "cluster_", sanitize_filename(cluster_id), "_",
              sanitize_filename(medoid_id), "_vs_", sanitize_filename(p), "_",
              sanitize_filename(pair_type), "_", sanitize_filename(ch), ".png"
            )
            out_img <- file.path(plots_dir, fn_img)
            rel_img <- paste0(plots_subdir, "/", fn_img)

            ok <- TRUE
            if (!file.exists(out_img) || isTRUE(overwrite)) {
              ok <- FALSE
              tryCatch({
                plt <- plot_pair_spectrum_mirror(
                  spectra = spectra,
                  idA = medoid_id,
                  idB = p,
                  channel = ch,
                  params = params,
                  use_aligned_bins = plot_use_aligned_bins,
                  normalize = plot_normalize,
                  scale = plot_scale,
                  label_top_n = plot_label_top_n,
                  label_digits = plot_label_digits,
                  title = paste0("Mirror spectrum: ", ch),
                  subtitle = paste0(medoid_id, " vs ", p, " (", type_label(pair_type), ", sim=", round(s, 3), ")")
                )
                ggplot2::ggsave(
                  filename = out_img,
                  plot = plt,
                  width = plot_width,
                  height = plot_height,
                  dpi = plot_dpi
                )
                ok <- TRUE
              }, error = function(e) {
                ok <- FALSE
              })
            }
            if (isTRUE(ok) && file.exists(out_img)) plot_rel <- c(plot_rel, rel_img)
          }

          if (isTRUE(plot_include_contribution)) {
            fn_img <- paste0(
              "cluster_", sanitize_filename(cluster_id), "_",
              sanitize_filename(medoid_id), "_vs_", sanitize_filename(p), "_",
              sanitize_filename(pair_type), "_contribution.png"
            )
            out_img <- file.path(plots_dir, fn_img)
            rel_img <- paste0(plots_subdir, "/", fn_img)
            ok <- TRUE
            if (!file.exists(out_img) || isTRUE(overwrite)) {
              ok <- FALSE
              tryCatch({
                plt <- plot_pair_contribution_bar(
                  det,
                  subtitle = paste0(medoid_id, " vs ", p, " (", type_label(pair_type), ", sim=", round(s, 3), ")")
                )
                ggplot2::ggsave(
                  filename = out_img,
                  plot = plt,
                  width = plot_width,
                  height = max(2.2, plot_height * 0.8),
                  dpi = plot_dpi
                )
                ok <- TRUE
              }, error = function(e) {
                ok <- FALSE
              })
            }
            if (isTRUE(ok) && file.exists(out_img)) plot_rel <- c(plot_rel, rel_img)
          }
        }

        if (length(plot_rel) > 0) {
          plot_md <- c(
            "**Spectra (mirror plots)**",
            "",
            vapply(plot_rel, function(x) paste0("![](", x, ")"), character(1)),
            ""
          )
        }
      }
format_frac <- function(x) {
  if (is.null(x) || !is.finite(x)) return("NA")
  format(round(as.numeric(x), 3), nsmall = 3, trim = TRUE)
}

# Contribution summary line (auto: most detailed available)
contrib_line <- NULL
within_loss_line <- NULL

if (!is.null(det$c_anchor_total) && !is.null(det$c_pair_total) &&
    (!is.null(det$c_anchor_raw_total) || !is.null(det$c_pair_raw_total))) {
  a_raw <- det$c_anchor_raw_total
  a_typ <- det$c_anchor_typ_total
  b_raw <- det$c_pair_raw_total
  b_typ <- det$c_pair_typ_total

  if ((is.finite(a_typ) && a_typ > 0) || (is.finite(b_typ) && b_typ > 0)) {
    contrib_line <- paste0(
      "- contribution: frag=", format_frac(det$c_frag_total),
      ", A_raw=", format_frac(a_raw), ", A_typ=", format_frac(a_typ),
      ", B_raw=", format_frac(b_raw), ", B_typ=", format_frac(b_typ)
    )
  } else {
    contrib_line <- paste0(
      "- contribution: frag=", format_frac(det$c_frag_total),
      ", loss_anchor=", format_frac(det$c_anchor_total),
      ", loss_pair=", format_frac(det$c_pair_total)
    )
  }

} else if (!is.null(det$c_loss_raw_total) && !is.null(det$c_loss_typ_total)) {
  if (is.finite(det$c_loss_typ_total) && det$c_loss_typ_total > 0) {
    contrib_line <- paste0(
      "- contribution: frag=", format_frac(det$c_frag_total),
      ", loss_raw=", format_frac(det$c_loss_raw_total),
      ", loss_typical=", format_frac(det$c_loss_typ_total)
    )
  } else {
    contrib_line <- paste0(
      "- contribution: frag=", format_frac(det$c_frag_total),
      ", loss=", format_frac(det$c_loss_total)
    )
  }

} else {
  contrib_line <- paste0(
    "- contribution: frag=", format_frac(det$c_frag_total),
    ", loss=", format_frac(det$c_loss_total)
  )
}



      lines <- c(lines,
                 paste0("### `", medoid_id, "` vs `", p, "` (sim=", round(s, 3), ", type=", type_label(pair_type), ")"),
                 "",
                 plot_md,
                 paste0("- partner cluster: ", ifelse(!is.na(partner_cluster) && nzchar(partner_cluster), partner_cluster, "NA")),
                 paste0("- d_total: ", round(det$d_total, 4), "; d_frag: ", round(det$d_frag, 4), "; d_loss: ", round(det$d_loss, 4)),
                 contrib_line,

                 if (!is.null(det$conf_pair)) paste0("- pair conf: ", round(det$conf_pair, 3)) else NULL,
                 paste0("- derivatization: medoid=", as.character(df_spec[medoid_id, "derivatization_type"]), ", partner=", as.character(df_spec[p, "derivatization_type"])),
                 ""
      )


      if (!is.null(det$typical_loss_matches)) {
        lines <- c(lines, "**Typical-loss overlap (combined)**", "", markdown_table(det$typical_loss_matches, digits = 4, max_rows = top_losses), "")
      }
      if (!is.null(det$typical_anchor_matches)) {
        lines <- c(lines, "**Typical-loss overlap (anchored A)**", "", markdown_table(det$typical_anchor_matches, digits = 4, max_rows = top_losses), "")
      }
      if (!is.null(det$typical_pair_matches)) {
        lines <- c(lines, "**Typical-loss overlap (pairwise B)**", "", markdown_table(det$typical_pair_matches, digits = 4, max_rows = top_losses), "")
      }
    }

writeLines(lines, fpath)
    report_files <- c(report_files, fpath)
  }

  pairs_all <- if (length(pairs_out) > 0) do.call(rbind, pairs_out) else data.frame()
  matches_all <- if (length(matches_out) > 0) {
    df <- do.call(rbind, matches_out)
    cols_front <- c("cluster", "medoid", "partner", "pair_type", "partner_cluster", "channel")
    cols_front <- cols_front[cols_front %in% colnames(df)]
    df <- df[, c(cols_front, setdiff(colnames(df), cols_front)), drop = FALSE]
    df
  } else data.frame()

  if (nrow(pairs_all) > 0) {
    readr::write_csv(pairs_all, file.path(report_dir, "cluster_report_pairs.csv"))
  }
  if (nrow(matches_all) > 0) {
    readr::write_csv(matches_all, file.path(report_dir, "cluster_report_typical_loss_matches.csv"))
  }

  index_file <- file.path(report_dir, "cluster_report_index.md")
  idx_lines <- c(
    "# Cluster report index",
    "",
    paste0("Generated: ", as.character(Sys.time())),
    "",
    paste0("Reports in: `", report_dir, "`"),
    "",
    "## Report files",
    ""
  )

  if (!is.null(plots_dir)) {
    idx_lines <- append(idx_lines, values = c(paste0("Pair plot images in: `", plots_dir, "`"), ""), after = 5)
  }
  if (length(report_files) == 0) {
    idx_lines <- c(idx_lines, "_No reports generated (no clusters with >=2 members or no partners met the criteria)._", "")
  } else {
    for (fp in report_files) {
      idx_lines <- c(idx_lines, paste0("- ", basename(fp)))
    }
    idx_lines <- c(idx_lines, "")
  }
  writeLines(idx_lines, index_file)

  list(
    pairs = pairs_all,
    matches = matches_all,
    report_files = report_files,
    index_file = index_file,
    report_dir = report_dir
  )
}
