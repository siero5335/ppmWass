#' Plot Similarity Heatmap
#'
#' @param sim_cluster Similarity matrix for clustering.
#' @param df_spec Data frame with compound_class and known.
#' @param threshold Similarity threshold for visualization.
#' @param log_transform If TRUE, plot log10(similarity) for contrast.
#' @param palette Vector of three colors for low/mid/high values. If NULL, use palette_name.
#' @param palette_name Named preset palette. Options: bwr, coolwarm, viridis, magma, cividis, gray.
#' @param layout_preset Heatmap layout preset. Options: default, soft, minimal, high_contrast.
#' @param class_palette Vector of colors for compound classes. If NULL, uses class_palette_name.
#' @param class_palette_name Named preset for class colors. Options: bright, pastel, muted, vivid, gray.
#' @param known_colors Named vector for known/unknown colors. Defaults depend on layout_preset.
#' @param legend_title Override legend title (defaults to log10(sim) or sim).
#' @param show_row_names Show row names (NULL uses preset/default).
#' @param show_column_names Show column names (NULL uses preset/default).
#' @param row_fontsize Row label font size.
#' @param column_fontsize Column label font size.
#' @param file Optional PDF output path.
#' @param width PDF width.
#' @param height PDF height.
#' @return Heatmap object (invisibly).
#' @export
plot_similarity_heatmap <- function(
  sim_cluster,
  df_spec,
  threshold = 0.2,
  log_transform = TRUE,
  palette = c("#2c7bb6", "#f7f7f7", "#d7191c"),
  palette_name = NULL,
  layout_preset = c("default", "soft", "minimal", "high_contrast"),
  class_palette = NULL,
  class_palette_name = NULL,
  known_colors = NULL,
  legend_title = NULL,
  show_row_names = NULL,
  show_column_names = NULL,
  row_fontsize = 6,
  column_fontsize = 6,
  file = NULL,
  width = 12,
  height = 10
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required for heatmap plotting.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required for heatmap plotting.")
  }

  safe_log10 <- function(x) {
    nz <- x[x > 0]
    if (length(nz) == 0) return(log10(x + 1e-12))
    log10(x + min(nz) / 2)
  }

  resolve_palette <- function(palette, palette_name) {
    if (!is.null(palette) && length(palette) == 1 && is.null(palette_name)) {
      palette_name <- palette
      palette <- NULL
    }
    if (!is.null(palette)) {
      if (length(palette) < 3) stop("palette must have at least 3 colors.")
      return(palette[1:3])
    }
    presets <- list(
      bwr = c("#2c7bb6", "#f7f7f7", "#d7191c"),
      coolwarm = c("#3b4cc0", "#f7f7f7", "#b40426"),
      viridis = c("#440154", "#21908C", "#FDE725"),
      magma = c("#000004", "#B63679", "#FCFDBF"),
      cividis = c("#00204c", "#7b9f8a", "#ffd966"),
      gray = c("#000000", "#888888", "#ffffff")
    )
    if (is.null(palette_name)) palette_name <- "bwr"
    if (!palette_name %in% names(presets)) {
      stop("Unknown palette_name. Use one of: ", paste(names(presets), collapse = ", "))
    }
    presets[[palette_name]]
  }

  palette <- resolve_palette(palette, palette_name)

  layout_preset <- match.arg(layout_preset)
  if (is.null(class_palette_name)) {
    class_palette_name <- switch(layout_preset,
      default = "bright",
      soft = "pastel",
      minimal = "muted",
      high_contrast = "vivid"
    )
  }
  if (is.null(known_colors)) {
    known_colors <- switch(layout_preset,
      default = c("T" = "darkgreen", "F" = "gray70"),
      soft = c("T" = "#2b8cbe", "F" = "gray80"),
      minimal = c("T" = "black", "F" = "gray70"),
      high_contrast = c("T" = "black", "F" = "gray60")
    )
  }

  mat_for_heat <- sim_cluster
  mat_for_heat[mat_for_heat <= threshold] <- 0
  mat_plot <- if (isTRUE(log_transform)) safe_log10(mat_for_heat) else mat_for_heat
  dist_for_heatmap <- stats::as.dist(1 - sim_cluster)

  class_levels <- unique(df_spec$compound_class)
  class_colors <- resolve_discrete_palette(length(class_levels), class_palette, class_palette_name)
  names(class_colors) <- class_levels

  row_anno <- ComplexHeatmap::rowAnnotation(
    Class = df_spec$compound_class,
    Known = df_spec$known,
    col = list(
      Class = class_colors,
      Known = known_colors
    )
  )

  rng <- range(mat_plot, finite = TRUE)
  mid <- if (isTRUE(log_transform)) -1 else 0
  mid <- min(max(mid, rng[1]), rng[2])
  col_fun <- circlize::colorRamp2(
    c(rng[1], mid, rng[2]),
    palette
  )

  if (is.null(show_row_names)) show_row_names <- (nrow(sim_cluster) <= 100)
  if (is.null(show_column_names)) show_column_names <- (nrow(sim_cluster) <= 100)

  legend_name <- if (!is.null(legend_title)) legend_title else if (isTRUE(log_transform)) "log10(sim)" else "sim"

  ht <- ComplexHeatmap::Heatmap(
    mat_plot,
    col = col_fun,
    name = legend_name,
    row_names_gp = grid::gpar(fontsize = row_fontsize),
    column_names_gp = grid::gpar(fontsize = column_fontsize),
    clustering_distance_columns = dist_for_heatmap,
    clustering_distance_rows = dist_for_heatmap,
    clustering_method_columns = "ward.D2",
    clustering_method_rows = "ward.D2",
    right_annotation = row_anno,
    show_row_names = show_row_names,
    show_column_names = show_column_names
  )

  if (!is.null(file)) {
    grDevices::pdf(file, width = width, height = height)
    ComplexHeatmap::draw(ht)
    grDevices::dev.off()
  }

  invisible(ht)
}

#' Plot UMAP Results
#'
#' @param res_umap UMAP result data frame with V1/V2, cluster, known, compound_class, RI.
#' @param distance_method Label for plot title.
#' @param point_size Point size for scatter.
#' @param point_alpha Alpha for points.
#' @param cluster_palette Vector of colors for cluster levels (optional).
#' @param cluster_palette_name Named preset for clusters (bright, pastel, muted, vivid, gray).
#' @param class_palette Vector of colors for compound classes (optional).
#' @param class_palette_name Named preset for classes (bright, pastel, muted, vivid, gray).
#' @param ri_palette Vector of colors for RI gradient (optional).
#' @param ri_palette_name Named preset for RI gradient (viridis, magma, cividis, inferno, plasma).
#' @param ri_na_color Color for NA RI values.
#' @param file Optional PDF output path.
#' @param width PDF width.
#' @param height PDF height.
#' @return A list of ggplot objects.
#' @export
plot_umap_results <- function(
  res_umap,
  distance_method = "",
  point_size = 2,
  point_alpha = 0.7,
  cluster_palette = NULL,
  cluster_palette_name = NULL,
  class_palette = NULL,
  class_palette_name = NULL,
  ri_palette = NULL,
  ri_palette_name = "plasma",
  ri_na_color = "gray50",
  file = NULL,
  width = 25,
  height = 5
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting UMAP results.")
  }
  if (!all(c("V1", "V2") %in% names(res_umap))) {
    stop("res_umap must contain columns V1 and V2.")
  }
  if (!"cluster" %in% names(res_umap)) res_umap$cluster <- factor("all")
  if (!"known" %in% names(res_umap)) res_umap$known <- "unknown"
  if (!"compound_class" %in% names(res_umap)) res_umap$compound_class <- "unknown"
  if (!"RI" %in% names(res_umap)) res_umap$RI <- NA_real_

  p1 <- ggplot2::ggplot(res_umap, ggplot2::aes(x = .data$V1, y = .data$V2, color = .data$cluster, shape = .data$known)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("UMAP: DBSCAN Cluster (", distance_method, ")"))

  cluster_levels <- levels(as.factor(res_umap$cluster))
  if (length(cluster_levels) > 0) {
    cluster_colors <- resolve_discrete_palette(length(cluster_levels), cluster_palette, cluster_palette_name)
    names(cluster_colors) <- cluster_levels
    p1 <- p1 + ggplot2::scale_color_manual(values = cluster_colors)
  }

  p2 <- ggplot2::ggplot(res_umap, ggplot2::aes(x = .data$V1, y = .data$V2, color = .data$compound_class)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "UMAP: Compound Class")

  class_levels <- unique(res_umap$compound_class)
  if (length(class_levels) > 0) {
    class_colors <- resolve_discrete_palette(length(class_levels), class_palette, class_palette_name)
    names(class_colors) <- class_levels
    p2 <- p2 + ggplot2::scale_color_manual(values = class_colors)
  }

  ri_scale <- if (is.null(ri_palette)) {
    ggplot2::scale_color_viridis_c(option = ri_palette_name, na.value = ri_na_color)
  } else {
    ggplot2::scale_color_gradientn(colors = ri_palette, na.value = ri_na_color)
  }

  p3 <- ggplot2::ggplot(res_umap, ggplot2::aes(x = .data$V1, y = .data$V2, color = .data$RI)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ri_scale +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "UMAP: Retention Index")

  if (!is.null(file)) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Package 'gridExtra' is required for arranging UMAP plots.")
    }
    grDevices::pdf(file, width = width, height = height)
    gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
    grDevices::dev.off()
  }

  list(cluster = p1, class = p2, ri = p3)
}


#' Plot Noise Robustness Results
#'
#' Creates a line plot showing Top-K accuracy as a function of m/z noise level,
#' with one line per distance method.
#'
#' @param robustness_result Output of \code{evaluate_noise_robustness()}.
#' @param k Which Top-K to plot. Default: 1 (Top-1 accuracy).
#' @param file Optional file path to save the plot (PDF).
#' @param width Plot width in inches. Default: 8.
#' @param height Plot height in inches. Default: 5.
#' @return A ggplot object (invisibly if file is specified).
#' @export
plot_noise_robustness <- function(robustness_result, k = 1, file = NULL,
                                   width = 8, height = 5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_noise_robustness()")
  }

  summary_df <- robustness_result$summary
  mean_col <- paste0("mean_top_", k)
  sd_col <- paste0("sd_top_", k)

  if (!mean_col %in% colnames(summary_df)) {
    stop("Top-", k, " accuracy not found in results. Available: ",
         paste(grep("^mean_top_", colnames(summary_df), value = TRUE),
               collapse = ", "))
  }

  plot_df <- summary_df[, c("method", "mz_noise_ppm", mean_col, sd_col)]
  colnames(plot_df) <- c("method", "mz_noise_ppm", "accuracy", "sd")

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = .data$mz_noise_ppm,
                                    y = .data$accuracy,
                                    color = .data$method,
                                    group = .data$method)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(.data$accuracy - .data$sd, 0),
                   ymax = pmin(.data$accuracy + .data$sd, 1)),
      width = 1, linewidth = 0.4
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(
      title = paste0("Noise Robustness: Top-", k, " Retrieval Accuracy"),
      x = "m/z Noise (ppm SD)",
      y = paste0("Top-", k, " Accuracy"),
      color = "Method"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  if (!is.null(file)) {
    ggplot2::ggsave(file, p, width = width, height = height)
    return(invisible(p))
  }

  p
}
