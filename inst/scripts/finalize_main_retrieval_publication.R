#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]])
}

output_dir <- normalizePath(get_arg("output-dir"), mustWork = TRUE)
requested_file_tag <- get_arg("file-tag", "")
if (nzchar(requested_file_tag) &&
    !grepl("^[A-Za-z0-9._-]+$", requested_file_tag)) {
  stop("Invalid --file-tag")
}
tagged_name <- function(stem, extension = ".csv") {
  paste0(stem, "_", file_tag, extension)
}
results_argument <- get_arg("results-rds", "")
if (!nzchar(results_argument)) {
  if (!nzchar(requested_file_tag)) {
    stop("Provide --results-rds or --file-tag for standalone finalization.")
  }
  file_tag <- requested_file_tag
  results_argument <- file.path(
    output_dir, tagged_name("main_retrieval_results", ".rds")
  )
}
results_file <- normalizePath(results_argument, mustWork = TRUE)
results <- readRDS(results_file)
stored_file_tag <- attr(results, "file_tag", exact = TRUE)
if (is.null(stored_file_tag) || length(stored_file_tag) != 1L ||
    is.na(stored_file_tag) || !nzchar(stored_file_tag)) {
  stored_file_tag <- ""
}
file_tag <- if (nzchar(requested_file_tag)) requested_file_tag else stored_file_tag
if (!nzchar(file_tag)) {
  stop("Result RDS has no file_tag metadata; provide --file-tag explicitly.")
}
if (!grepl("^[A-Za-z0-9._-]+$", file_tag)) stop("Invalid stored file_tag")
if (nzchar(stored_file_tag) && !identical(file_tag, stored_file_tag)) {
  stop(
    "--file-tag=", file_tag,
    " does not match result RDS file_tag=", stored_file_tag, "."
  )
}

bind_component <- function(name) {
  do.call(rbind, lapply(results, function(x) x[[name]]))
}

combined_summary <- bind_component("summary")
combined_query_bootstrap <- bind_component("query_bootstrap")
combined_per_query <- bind_component("per_query")
combined_tie_sensitivity <- bind_component("tie_sensitivity")
combined_rounded_tie_sensitivity <- bind_component("rounded_tie_sensitivity")
combined_tie_diagnostic <- bind_component("tie_diagnostic")
combined_tie_diagnostic_summary <- bind_component("tie_diagnostic_summary")
combined_paired_cluster <- bind_component("paired_cluster")
combined_mcnemar <- bind_component("mcnemar")
combined_sanitize <- bind_component("sanitize")
combined_qc <- bind_component("qc")
combined_pairwise_auc <- bind_component("pairwise_auc")
combined_fdr <- bind_component("fdr")

combined_outputs <- list(
  combined_retrieval_cluster_bootstrap_ci = combined_summary,
  combined_retrieval_query_bootstrap_ci_exploratory = combined_query_bootstrap,
  combined_retrieval_per_query = combined_per_query,
  combined_tie_sensitivity_exact = combined_tie_sensitivity,
  combined_tie_sensitivity_round12 = combined_rounded_tie_sensitivity,
  combined_tie_diagnostics = combined_tie_diagnostic,
  combined_retrieval_tie_affected_summary = combined_tie_diagnostic_summary,
  combined_paired_cluster_bootstrap = combined_paired_cluster,
  combined_retrieval_mcnemar = combined_mcnemar,
  combined_distance_sanitize_log = combined_sanitize,
  combined_qc_summary = combined_qc,
  combined_retrieval_ordered_pair_auc = combined_pairwise_auc,
  combined_retrieval_ordered_pair_fdr = combined_fdr
)
for (output_stem in names(combined_outputs)) {
  utils::write.csv(
    combined_outputs[[output_stem]],
    file.path(output_dir, tagged_name(output_stem, ".csv")), row.names = FALSE
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
method_labels <- c(
  ppm_wasserstein = "ppm-Wasserstein", composite = "Composite",
  entropy_weighted = "Entropy (weighted)",
  entropy_unweighted = "Entropy (unweighted)", cosine = "Cosine",
  weighted_cosine = "Weighted cosine", hellinger = "Hellinger"
)
method_colors <- c(
  ppm_wasserstein = "#D62728", composite = "#9467BD",
  entropy_weighted = "#1F77B4", entropy_unweighted = "#17BECF",
  cosine = "#7F7F7F", weighted_cosine = "#8C564B",
  hellinger = "#2CA02C"
)

figure_df <- subset(combined_summary, metric %in% c("top1", "mrr", "p_at_1"))
figure_df$metric <- factor(
  figure_df$metric, levels = c("top1", "mrr", "p_at_1"),
  labels = c("Top-1", "MRR", "P@1")
)
figure_df$dataset <- factor(figure_df$dataset, levels = unique(figure_df$dataset))
figure_df$method <- factor(
  figure_df$method, levels = names(method_labels), labels = unname(method_labels)
)
color_map <- stats::setNames(
  method_colors[names(method_labels)], unname(method_labels)
)
utils::write.csv(
  figure_df, file.path(output_dir, tagged_name("figure1_clean_retrieval_data", ".csv")),
  row.names = FALSE
)

p <- ggplot2::ggplot(
  figure_df, ggplot2::aes(x = method, y = estimate, fill = method)
) +
  ggplot2::geom_col(width = 0.72, alpha = 0.92) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = cluster_ci_low, ymax = cluster_ci_high),
    width = 0.18, linewidth = 0.45, color = "grey20"
  ) +
  ggplot2::facet_grid(dataset ~ metric, scales = "fixed") +
  ggplot2::scale_fill_manual(values = color_map, guide = "none") +
  ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
  ggplot2::labs(
    title = "Clean replicate retrieval in two GC-HRMS EI libraries",
    x = NULL, y = "Retrieval metric value"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    strip.text = ggplot2::element_text(face = "bold", size = 11),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 8, r = 12, b = 16, l = 14)
  )

ggplot2::ggsave(
  file.path(output_dir, tagged_name("figure1_clean_retrieval_bootstrap", ".png")),
  p, width = 12, height = 7.2, dpi = 300
)
ggplot2::ggsave(
  file.path(output_dir, tagged_name("figure1_clean_retrieval_bootstrap", ".pdf")),
  p, width = 12, height = 7.2, device = "pdf"
)
message("Finalized main retrieval outputs in: ", output_dir)
