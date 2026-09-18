#' @keywords internal
#' @noRd
resolve_named_options <- function(defaults, overrides = NULL, arg_name = "options") {
  out <- defaults
  if (is.null(overrides)) {
    return(out)
  }
  if (!is.list(overrides)) {
    stop(arg_name, " must be a list when provided.")
  }
  unknown <- setdiff(names(overrides), names(defaults))
  if (length(unknown) > 0) {
    stop(arg_name, " contains unknown names: ", paste(unknown, collapse = ", "))
  }
  out[names(overrides)] <- overrides
  out
}

#' Run EI-HRMS Similarity Pipeline
#'
#' `heatmap_opts`, `umap_opts`, and `cluster_report_opts` provide grouped
#' configuration without removing the existing scalar arguments. When supplied,
#' list entries override the corresponding individual arguments.
#'
#' Supported `heatmap_opts` names:
#' `threshold`, `log_transform`, `palette`, `palette_name`, `layout_preset`,
#' `class_palette`, `class_palette_name`, `known_colors`, `legend_title`,
#' `show_row_names`, `show_column_names`, `row_fontsize`, `column_fontsize`,
#' `width`, `height`.
#'
#' Supported `umap_opts` names:
#' `point_size`, `point_alpha`, `cluster_palette`, `cluster_palette_name`,
#' `class_palette`, `class_palette_name`, `ri_palette`, `ri_palette_name`,
#' `ri_na_color`, `width`, `height`, `seed`.
#'
#' Supported `cluster_report_opts` names:
#' `subdir`, `top_pairs`, `top_losses`, `min_similarity`, `overwrite`,
#' `include_pair_plots`, `plots_subdir`, `plot_channels`, `plot_as_panel`,
#' `plot_include_contribution`, `plot_panel_height`, `plot_use_aligned_bins`,
#' `plot_normalize`, `plot_scale`, `plot_label_top_n`, `plot_label_digits`,
#' `plot_width`, `plot_height`, `plot_dpi`.
#'
#' @param file Path to MS-DIAL export file.
#' @param df Optional MS-DIAL tibble or `msdial_import` object to use instead of reading from file.
#' @param params Parameter list.
#' @param output_dir Output directory for files.
#' @param write_outputs If TRUE, write CSV/PDF outputs.
#' @param progress If TRUE, print progress.
#' @param export_unknown_as_id If TRUE, write unknowns using compound_id in cleaned CSV.
#' @param heatmap_threshold Similarity threshold for heatmap display.
#' @param heatmap_log_transform If TRUE, plot log10(similarity) heatmap.
#' @param heatmap_palette Vector of three colors for low/mid/high values.
#' @param heatmap_palette_name Named preset palette for heatmap. If set, overrides heatmap_palette.
#' @param heatmap_layout_preset Heatmap layout preset (default, soft, minimal, high_contrast).
#' @param heatmap_class_palette Vector of colors for class annotations.
#' @param heatmap_class_palette_name Named preset for class colors (bright, pastel, muted, vivid, gray).
#' @param heatmap_known_colors Named vector for known/unknown colors.
#' @param heatmap_legend_title Override legend title for heatmap.
#' @param heatmap_show_row_names Show row names (NULL uses preset/default).
#' @param heatmap_show_column_names Show column names (NULL uses preset/default).
#' @param heatmap_row_fontsize Row label font size.
#' @param heatmap_column_fontsize Column label font size.
#' @param heatmap_width PDF width for heatmap.
#' @param heatmap_height PDF height for heatmap.
#' @param umap_point_size Point size for UMAP scatter plots.
#' @param umap_point_alpha Alpha for UMAP points.
#' @param umap_cluster_palette Vector of colors for cluster levels.
#' @param umap_cluster_palette_name Named preset for cluster colors (bright, pastel, muted, vivid, gray).
#' @param umap_class_palette Vector of colors for compound classes.
#' @param umap_class_palette_name Named preset for class colors (bright, pastel, muted, vivid, gray).
#' @param umap_ri_palette Vector of colors for RI gradient.
#' @param umap_ri_palette_name Named preset for RI gradient (viridis, magma, cividis, inferno, plasma).
#' @param umap_ri_na_color Color for NA RI values.
#' @param umap_width PDF width for UMAP plot.
#' @param umap_height PDF height for UMAP plot.
#' @param umap_seed Random seed for UMAP.
#' @param optics_eps eps_cl for OPTICS clustering.
#' @param write_cluster_annotations If TRUE, compute cluster-level driver labels and typical-loss motifs
#'   (requires params$return_distance_components = TRUE to label drivers).
#' @param cluster_top_n_losses Number of top typical losses to include as a compact per-cluster label.
#' @param cluster_use_idf If TRUE, apply a simple IDF weighting when ranking typical losses in clusters.
#' @param write_derivatization_evaluation If TRUE, write derivatization auto-detection evaluation outputs (CSV).
#' @param derivatization_truth_col Optional column name in the input `df` containing true derivatization labels
#'   (e.g., "none"/"TMS"/"TBDMS"). If provided, confusion/ROC outputs are generated.
#' @param derivatization_eval_subdir Subdirectory under output_dir for derivatization evaluation outputs.
#' @param write_cluster_reports If TRUE, generate per-cluster explanation reports (Markdown/CSV).
#' @param cluster_report_subdir Subdirectory under output_dir for cluster reports.
#' @param cluster_report_top_pairs Number of medoid-vs-member pairs to include per cluster.
#' @param cluster_report_top_losses Number of typical-loss matches to include per pair.
#' @param cluster_report_min_similarity Minimum similarity for selecting medoid partners.
#' @param cluster_report_overwrite If TRUE, overwrite existing Markdown reports.
#' @param cluster_report_include_pair_plots If TRUE, render per-pair spectrum plots in reports.
#' @param cluster_report_plots_subdir Subdirectory for generated pair plots.
#' @param cluster_report_plot_channels Channels to include in pair plots.
#' @param cluster_report_plot_as_panel If TRUE, combine selected channels into a panel plot.
#' @param cluster_report_plot_include_contribution If TRUE, include contribution bars in pair plots.
#' @param cluster_report_plot_panel_height Optional panel height override.
#' @param cluster_report_plot_use_aligned_bins If TRUE, plot aligned spectra bins.
#' @param cluster_report_plot_normalize Intensity normalization used for plots.
#' @param cluster_report_plot_scale Intensity scaling factor for plots.
#' @param cluster_report_plot_label_top_n Number of top peaks to label in plots.
#' @param cluster_report_plot_label_digits Number of digits for peak labels.
#' @param cluster_report_plot_width Plot width in inches.
#' @param cluster_report_plot_height Plot height in inches.
#' @param cluster_report_plot_dpi Plot DPI.
#' @param heatmap_opts Optional list overriding heatmap-related arguments.
#' @param umap_opts Optional list overriding UMAP-related arguments.
#' @param cluster_report_opts Optional list overriding cluster-report-related arguments.
#' @param deduplicate If TRUE, keep the current prepare-stage deduplication behavior for known compounds.
#' @return A list with all intermediate and final results.
#' @export
run_eihrms_similarity <- function(
  file = NULL,
  df = NULL,
  params = eihrms_default_params(),
  output_dir = ".",
  write_outputs = TRUE,
  progress = TRUE,
  export_unknown_as_id = TRUE,
  heatmap_threshold = 0.2,
  heatmap_log_transform = TRUE,
  heatmap_palette = c("#2c7bb6", "#f7f7f7", "#d7191c"),
  heatmap_palette_name = NULL,
  heatmap_layout_preset = "default",
  heatmap_class_palette = NULL,
  heatmap_class_palette_name = NULL,
  heatmap_known_colors = NULL,
  heatmap_legend_title = NULL,
  heatmap_show_row_names = NULL,
  heatmap_show_column_names = NULL,
  heatmap_row_fontsize = 6,
  heatmap_column_fontsize = 6,
  heatmap_width = 12,
  heatmap_height = 10,
  umap_point_size = 2,
  umap_point_alpha = 0.7,
  umap_cluster_palette = NULL,
  umap_cluster_palette_name = NULL,
  umap_class_palette = NULL,
  umap_class_palette_name = NULL,
  umap_ri_palette = NULL,
  umap_ri_palette_name = "plasma",
  umap_ri_na_color = "gray50",
  umap_width = 25,
  umap_height = 5,
  umap_seed = 71,
  optics_eps = 0.73,
  write_cluster_annotations = TRUE,
  cluster_top_n_losses = 5,
  cluster_use_idf = TRUE,
  write_derivatization_evaluation = FALSE,
  derivatization_truth_col = NULL,
  derivatization_eval_subdir = "derivatization_eval",
  write_cluster_reports = TRUE,
  cluster_report_subdir = "cluster_reports",
  cluster_report_top_pairs = 5,
  cluster_report_top_losses = 10,
  cluster_report_min_similarity = 0,
  cluster_report_overwrite = TRUE,
  cluster_report_include_pair_plots = TRUE,
  cluster_report_plots_subdir = "plots",
  cluster_report_plot_channels = c("frag", "loss"),
  cluster_report_plot_as_panel = TRUE,
  cluster_report_plot_include_contribution = TRUE,
  cluster_report_plot_panel_height = NULL,
  cluster_report_plot_use_aligned_bins = TRUE,
  cluster_report_plot_normalize = "max",
  cluster_report_plot_scale = 100,
  cluster_report_plot_label_top_n = 8,
  cluster_report_plot_label_digits = 4,
  cluster_report_plot_width = 10,
  cluster_report_plot_height = 4,
  cluster_report_plot_dpi = 180,
  heatmap_opts = NULL,
  umap_opts = NULL,
  cluster_report_opts = NULL,
  deduplicate = TRUE
) {

  params <- validate_params(params)
  if (is.null(df)) {
    if (is.null(file)) stop("Either 'file' or 'df' must be provided.")
    if (isTRUE(progress)) message("Reading input: ", file)
    df <- read_ms_dial(file)
  } else {
    if (isTRUE(progress)) message("Using provided data frame input.")
  }
  df_tbl <- as_feature_table(df)

  heatmap_cfg <- resolve_named_options(
    list(
      threshold = heatmap_threshold,
      log_transform = heatmap_log_transform,
      palette = heatmap_palette,
      palette_name = heatmap_palette_name,
      layout_preset = heatmap_layout_preset,
      class_palette = heatmap_class_palette,
      class_palette_name = heatmap_class_palette_name,
      known_colors = heatmap_known_colors,
      legend_title = heatmap_legend_title,
      show_row_names = heatmap_show_row_names,
      show_column_names = heatmap_show_column_names,
      row_fontsize = heatmap_row_fontsize,
      column_fontsize = heatmap_column_fontsize,
      width = heatmap_width,
      height = heatmap_height
    ),
    heatmap_opts,
    arg_name = "heatmap_opts"
  )

  umap_cfg <- resolve_named_options(
    list(
      point_size = umap_point_size,
      point_alpha = umap_point_alpha,
      cluster_palette = umap_cluster_palette,
      cluster_palette_name = umap_cluster_palette_name,
      class_palette = umap_class_palette,
      class_palette_name = umap_class_palette_name,
      ri_palette = umap_ri_palette,
      ri_palette_name = umap_ri_palette_name,
      ri_na_color = umap_ri_na_color,
      width = umap_width,
      height = umap_height,
      seed = umap_seed
    ),
    umap_opts,
    arg_name = "umap_opts"
  )

  cluster_report_cfg <- resolve_named_options(
    list(
      subdir = cluster_report_subdir,
      top_pairs = cluster_report_top_pairs,
      top_losses = cluster_report_top_losses,
      min_similarity = cluster_report_min_similarity,
      overwrite = cluster_report_overwrite,
      include_pair_plots = cluster_report_include_pair_plots,
      plots_subdir = cluster_report_plots_subdir,
      plot_channels = cluster_report_plot_channels,
      plot_as_panel = cluster_report_plot_as_panel,
      plot_include_contribution = cluster_report_plot_include_contribution,
      plot_panel_height = cluster_report_plot_panel_height,
      plot_use_aligned_bins = cluster_report_plot_use_aligned_bins,
      plot_normalize = cluster_report_plot_normalize,
      plot_scale = cluster_report_plot_scale,
      plot_label_top_n = cluster_report_plot_label_top_n,
      plot_label_digits = cluster_report_plot_label_digits,
      plot_width = cluster_report_plot_width,
      plot_height = cluster_report_plot_height,
      plot_dpi = cluster_report_plot_dpi
    ),
    cluster_report_opts,
    arg_name = "cluster_report_opts"
  )

  if (isTRUE(progress)) message("Preparing quant data...")
  quant <- prepare_quant_data(
    df,
    params,
    write_cleaned_csv = isTRUE(write_outputs),
    cleaned_csv_path = file.path(output_dir, "cleaned_voc.csv"),
    export_unknown_as_id = export_unknown_as_id,
    deduplicate = deduplicate
  )

  if (isTRUE(progress)) message("Building spectra...")
  spectra <- build_spectra(df, quant$final, params, progress = progress)

  if (isTRUE(write_outputs)) {
    readr::write_csv(spectra$class_summary, file.path(output_dir, "compound_class_annotation_detailed.csv"))
    readr::write_csv(spectra$class_summary_simple, file.path(output_dir, "compound_class_annotation.csv"))
  }



  # Derivatization evaluation (optional)
  derivatization_evaluation <- NULL
  if (isTRUE(write_derivatization_evaluation)) {
    truth_vec <- NULL
    if (!is.null(derivatization_truth_col) && !is.null(df_tbl) && derivatization_truth_col %in% names(df_tbl)) {
      truth_df <- df_tbl %>%
        dplyr::mutate(compound_id = derive_compound_id(.)) %>%
        dplyr::transmute(id = .data$compound_id, truth = .data[[derivatization_truth_col]]) %>%
        dplyr::filter(.data$id %in% spectra$df_spec$id) %>%
        dplyr::distinct(.data$id, .keep_all = TRUE)
      truth_vec <- stats::setNames(as.character(truth_df$truth), as.character(truth_df$id))
    }

    deriv_eval_dir <- file.path(output_dir, derivatization_eval_subdir)
    derivatization_evaluation <- tryCatch({
      run_derivatization_evaluation(
        spectra = spectra,
        truth = truth_vec,
        output_dir = deriv_eval_dir,
        write_outputs = isTRUE(write_outputs),
        allowed = params$derivatization_type
      )
    }, error = function(e) {
      warning("Derivatization evaluation failed: ", conditionMessage(e))
      NULL
    })
  }

  if (isTRUE(progress)) message("Computing distances and similarities...")
  sim <- compute_similarity_matrices(
    spectra$frag_list,
    spectra$loss_list,
    spectra$ri,
    params,
    loss_typ_list = spectra$loss_typ_list,
    mref_conf = spectra$mref_confidence,
    loss_anchor_list = spectra$loss_anchor_list,
    loss_pair_list = spectra$loss_pair_list,
    loss_anchor_typ_list = spectra$loss_anchor_typ_list,
    loss_pair_typ_list = spectra$loss_pair_typ_list
  )

  if (isTRUE(write_outputs)) {
    readr::write_csv(
      as.data.frame(sim$dist_raw) %>% tibble::rownames_to_column("id"),
      file.path(output_dir, paste0("distance_matrix_", params$distance_method, ".csv"))
    )
    readr::write_csv(
      as.data.frame(sim$sim_raw) %>% tibble::rownames_to_column("id"),
      file.path(output_dir, paste0("similarity_matrix_", params$distance_method, ".csv"))
    )
    readr::write_csv(
      as.data.frame(sim$sim_cluster) %>% tibble::rownames_to_column("id"),
      file.path(output_dir, "similarity_combined_cluster.csv")
    )
    readr::write_csv(
      as.data.frame(sim$sim_analog) %>% tibble::rownames_to_column("id"),
      file.path(output_dir, "similarity_combined_analog_kNN.csv")
    )
  }

  if (isTRUE(progress)) message("Clustering on distance matrix (OPTICS)...")
  cluster <- cluster_optics(sim$sim_cluster, eps_cl = optics_eps)

  if (isTRUE(progress)) message("Running UMAP...")
  res_umap <- run_umap(sim$sim_cluster, seed = umap_cfg$seed)

  ids <- names(spectra$frag_list)
  res_umap$Variable <- ids
  res_umap$known <- spectra$df_spec$known
  res_umap$compound_class <- spectra$df_spec$compound_class
  res_umap$RI <- spectra$df_spec$RI
  res_umap$cluster <- cluster$clusters

  # Cluster-level interpretability (drivers + typical-loss motifs)
  cluster_annotation <- NULL
  if (isTRUE(write_cluster_annotations)) {
    cluster_annotation <- tryCatch({
      annotate_clusters(
        spectra = spectra,
        clusters = cluster$clusters,
        similarity = sim,
        params = params,
        top_n_losses = cluster_top_n_losses,
        exclude_noise = TRUE,
        use_idf = cluster_use_idf
      )
    }, error = function(e) {
      warning("Cluster annotation failed: ", conditionMessage(e))
      NULL
    })

    if (isTRUE(write_outputs) && !is.null(cluster_annotation)) {
      # Always write membership (lightweight)
      readr::write_csv(cluster_annotation$membership, file.path(output_dir, "cluster_membership.csv"))

      # Driver summary requires distance components
      if (!is.null(cluster_annotation$cluster_summary)) {
        readr::write_csv(cluster_annotation$cluster_summary, file.path(output_dir, "cluster_driver_summary.csv"))
      }

      # Typical-loss motifs (long format)
      if (!is.null(cluster_annotation$typical_loss_long) && nrow(cluster_annotation$typical_loss_long) > 0) {
        readr::write_csv(cluster_annotation$typical_loss_long, file.path(output_dir, "cluster_typical_loss_motifs.csv"))
      }
    }
  }

  # Cluster explanation reports (optional)
  cluster_reports <- NULL
  if (isTRUE(write_cluster_reports) && isTRUE(write_outputs)) {
    cluster_reports <- tryCatch({
      generate_cluster_explanation_reports(
        spectra = spectra,
        similarity = sim,
        clusters = cluster$clusters,
        params = params,
        cluster_annotation = cluster_annotation,
        output_dir = output_dir,
        report_subdir = cluster_report_cfg$subdir,
        top_pairs = cluster_report_cfg$top_pairs,
        top_losses = cluster_report_cfg$top_losses,
        min_similarity = cluster_report_cfg$min_similarity,
        exclude_noise = TRUE,
        overwrite = cluster_report_cfg$overwrite,
        include_pair_plots = cluster_report_cfg$include_pair_plots,
        plots_subdir = cluster_report_cfg$plots_subdir,
        plot_channels = cluster_report_cfg$plot_channels,
        plot_as_panel = cluster_report_cfg$plot_as_panel,
        plot_include_contribution = cluster_report_cfg$plot_include_contribution,
        plot_panel_height = cluster_report_cfg$plot_panel_height,
        plot_use_aligned_bins = cluster_report_cfg$plot_use_aligned_bins,
        plot_normalize = cluster_report_cfg$plot_normalize,
        plot_scale = cluster_report_cfg$plot_scale,
        plot_label_top_n = cluster_report_cfg$plot_label_top_n,
        plot_label_digits = cluster_report_cfg$plot_label_digits,
        plot_width = cluster_report_cfg$plot_width,
        plot_height = cluster_report_cfg$plot_height,
        plot_dpi = cluster_report_cfg$plot_dpi
      )
    }, error = function(e) {
      warning("Cluster report generation failed: ", conditionMessage(e))
      NULL
    })
  }

  if (isTRUE(write_outputs)) {
    readr::write_csv(res_umap, file.path(output_dir, "umap_results.csv"))

    plot_similarity_heatmap(
      sim$sim_cluster,
      spectra$df_spec,
      threshold = heatmap_cfg$threshold,
      log_transform = heatmap_cfg$log_transform,
      palette = heatmap_cfg$palette,
      palette_name = heatmap_cfg$palette_name,
      layout_preset = heatmap_cfg$layout_preset,
      class_palette = heatmap_cfg$class_palette,
      class_palette_name = heatmap_cfg$class_palette_name,
      known_colors = heatmap_cfg$known_colors,
      legend_title = heatmap_cfg$legend_title,
      show_row_names = heatmap_cfg$show_row_names,
      show_column_names = heatmap_cfg$show_column_names,
      row_fontsize = heatmap_cfg$row_fontsize,
      column_fontsize = heatmap_cfg$column_fontsize,
      file = file.path(output_dir, paste0("similarity_heatmap_", params$distance_method, ".pdf")),
      width = heatmap_cfg$width,
      height = heatmap_cfg$height
    )

    plot_umap_results(
      res_umap,
      distance_method = params$distance_method,
      point_size = umap_cfg$point_size,
      point_alpha = umap_cfg$point_alpha,
      cluster_palette = umap_cfg$cluster_palette,
      cluster_palette_name = umap_cfg$cluster_palette_name,
      class_palette = umap_cfg$class_palette,
      class_palette_name = umap_cfg$class_palette_name,
      ri_palette = umap_cfg$ri_palette,
      ri_palette_name = umap_cfg$ri_palette_name,
      ri_na_color = umap_cfg$ri_na_color,
      file = file.path(output_dir, paste0("umap_visualization_", params$distance_method, ".pdf")),
      width = umap_cfg$width,
      height = umap_cfg$height
    )
  }

  if (isTRUE(progress)) {
    message("\n========================================")
    message("       Analysis Complete")
    message("========================================")
    message("Total compounds: ", length(ids))
    message("Known: ", sum(spectra$df_spec$known == "T"))
    message("Unknown: ", sum(spectra$df_spec$known == "F"))
    message("\n--- Settings ---")
    message("Distance method: ", params$distance_method)
    message("Wasserstein align: ", params$wasserstein_align)
    message("Fragment weight: ", params$w_frag)
    message("Loss weight: ", params$w_loss)
    message("RI mode (cluster): ", params$ri_mode_cluster)
    message("RI mode (analog): ", params$ri_mode_analog)
    message("\n--- High-Resolution MS Features ---")
    message("Tolerance for class detection: ", params$class_detection_ppm, " ppm")
    if (isTRUE(params$use_mref_confidence) && !all(is.na(spectra$mref_confidence))) {
      mc <- mean(spectra$mref_confidence, na.rm = TRUE)
      message("Mean Mref confidence: ", round(mc, 3))
    }
    message("Compounds with detected Neutral Losses: ",
            sum(spectra$df_spec$detected_neutral_losses != "none" & spectra$df_spec$detected_neutral_losses != ""))
    message("Compounds with inferred functional groups: ",
            sum(spectra$df_spec$inferred_functional_groups != ""))
  }

  list(
    params = params,
    quant = quant,
    spectra = spectra,
    similarity = sim,
    cluster = cluster,
    umap = res_umap,
    cluster_annotation = cluster_annotation,
    derivatization_evaluation = derivatization_evaluation,
    cluster_reports = cluster_reports
  )
}
