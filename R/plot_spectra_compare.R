#' Spectrum comparison plots (mirror / overlay)
#'
#' ggplot2-based utilities to visually compare two spectra. These functions are designed
#' to work directly with the internal spectrum matrices used by the package
#' (two-column matrices with columns `mz` and `intensity`).
#'
#' @details
#' The mirror plot supports two modes:
#'
#' * `use_aligned_bins = TRUE` (default): uses `align_spectra()` to bin/align peaks
#'   using a ppm-dependent window (consistent with distance calculations).
#' * `use_aligned_bins = FALSE`: plots raw peaks as-is; matched peaks are optionally
#'   identified by ppm windows for alpha shading.
#'
#' @name plot_spectrum_mirror
NULL

#' @keywords internal
#' @noRd
require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it (install.packages('ggplot2')).")
  }
}

#' @keywords internal
#' @noRd
rescale_intensity <- function(int, normalize = c("max", "sum", "none"), scale = 100) {
  normalize <- match.arg(normalize)
  if (length(int) == 0) return(int)
  if (normalize == "none") return(int)
  if (!is.finite(scale) || scale <= 0) stop("scale must be a positive number.")

  if (normalize == "max") {
    m <- max(int)
    if (!is.finite(m) || m <= 0) return(rep(0, length(int)))
    return(int / m * scale)
  }
  # normalize == "sum"
  s <- sum(int)
  if (!is.finite(s) || s <= 0) return(rep(0, length(int)))
  int / s * scale
}

#' @keywords internal
#' @noRd
mark_has_match_ppm <- function(mz_query, mz_ref, ppm = 20) {
  if (length(mz_query) == 0) return(logical(0))
  if (length(mz_ref) == 0) return(rep(FALSE, length(mz_query)))
  mz_ref <- sort(mz_ref)
  tol <- mz_query * ppm * 1e-6
  lower <- mz_query - tol
  upper <- mz_query + tol
  left <- findInterval(lower, mz_ref, left.open = TRUE) + 1
  right <- findInterval(upper, mz_ref)
  left <= right
}

#' Prepare aligned (binned) spectra for plotting
#'
#' @keywords internal
#' @noRd
prep_aligned_spectra_wide <- function(specA, specB, ppm = 20, normalize = "max", scale = 100) {
  aligned <- align_spectra(specA, specB, ppm = ppm)
  mz <- aligned$mz
  intA <- aligned$p
  intB <- aligned$q
  # align_spectra normalizes to sum=1; rescale for display
  intA <- rescale_intensity(intA, normalize = normalize, scale = scale)
  intB <- rescale_intensity(intB, normalize = normalize, scale = scale)
  tibble::tibble(
    mz = mz,
    intA = intA,
    intB = intB,
    shared = (intA > 0) & (intB > 0)
  )
}

#' Prepare raw spectra for plotting (no binning)
#'
#' @keywords internal
#' @noRd
prep_raw_spectra_long <- function(specA, specB, ppm = 20, normalize = "max", scale = 100,
                                  sample_names = c("A", "B"), alpha_shared = 1, alpha_unique = 0.4) {
  mzA <- if (nrow(specA) == 0) numeric(0) else specA[, 1]
  mzB <- if (nrow(specB) == 0) numeric(0) else specB[, 1]
  intA <- if (nrow(specA) == 0) numeric(0) else specA[, 2]
  intB <- if (nrow(specB) == 0) numeric(0) else specB[, 2]
  intA <- rescale_intensity(intA, normalize = normalize, scale = scale)
  intB <- rescale_intensity(intB, normalize = normalize, scale = scale)
  mA <- mark_has_match_ppm(mzA, mzB, ppm = ppm)
  mB <- mark_has_match_ppm(mzB, mzA, ppm = ppm)

  dfA <- tibble::tibble(
    mz = mzA,
    intensity = intA,
    sample = sample_names[1],
    shared = mA,
    alpha = ifelse(mA, alpha_shared, alpha_unique)
  )
  dfB <- tibble::tibble(
    mz = mzB,
    intensity = -intB,
    sample = sample_names[2],
    shared = mB,
    alpha = ifelse(mB, alpha_shared, alpha_unique)
  )
  dplyr::bind_rows(dfA, dfB)
}

#' Plot a mirror spectrum comparison (two spectra)
#'
#' @param specA,specB Two spectra matrices (columns `mz`, `intensity`).
#' @param ppm PPM tolerance used for binning/alignment and match shading.
#' @param use_aligned_bins If TRUE, plot aligned/binned spectra using `align_spectra()`.
#' @param normalize Intensity scaling for display: `max` (default), `sum`, or `none`.
#' @param scale Display scale when `normalize != 'none'` (e.g., 100).
#' @param sample_names Length-2 character vector for legend labels.
#' @param alpha_shared Alpha for peaks/bins present in both spectra (raw mode).
#' @param alpha_unique Alpha for peaks present only in one spectrum (raw mode).
#' @param label_top_n Number of peaks to label per spectrum (0 to disable).
#' @param label_min Minimum absolute displayed intensity required for labeling.
#' @param label_digits Number of digits for m/z labels.
#' @param title,subtitle Plot title/subtitle.
#' @param xlim Optional x-axis limits (numeric length 2).
#' @param file Optional output file path (e.g., .pdf, .png) via `ggplot2::ggsave`.
#' @param width,height Output size for `file`.
#' @return A ggplot object.
#' @export
plot_spectrum_mirror <- function(
  specA,
  specB,
  ppm = 20,
  use_aligned_bins = TRUE,
  normalize = c("max", "sum", "none"),
  scale = 100,
  sample_names = c("A", "B"),
  alpha_shared = 1,
  alpha_unique = 0.4,
  label_top_n = 10,
  label_min = NULL,
  label_digits = 4,
  title = NULL,
  subtitle = NULL,
  xlim = NULL,
  file = NULL,
  width = 9,
  height = 5
) {
  require_ggplot2()
  normalize <- match.arg(normalize)
  if (length(sample_names) != 2) stop("sample_names must have length 2.")
  if (!is.null(xlim) && (length(xlim) != 2 || !all(is.finite(xlim)))) stop("xlim must be a numeric vector of length 2.")

  if (isTRUE(use_aligned_bins)) {
    dfw <- prep_aligned_spectra_wide(specA, specB, ppm = ppm, normalize = normalize, scale = scale)
    dfA <- tibble::tibble(mz = dfw$mz, intensity = dfw$intA, sample = sample_names[1], shared = dfw$shared, alpha = 1)
    dfB <- tibble::tibble(mz = dfw$mz, intensity = -dfw$intB, sample = sample_names[2], shared = dfw$shared, alpha = 1)
    df <- dplyr::bind_rows(dfA, dfB)
  } else {
    df <- prep_raw_spectra_long(
      specA, specB,
      ppm = ppm,
      normalize = normalize,
      scale = scale,
      sample_names = sample_names,
      alpha_shared = alpha_shared,
      alpha_unique = alpha_unique
    )
  }

  if (nrow(df) == 0) {
    p <- ggplot2::ggplot() + ggplot2::theme_minimal() + ggplot2::labs(title = title, subtitle = subtitle)
    if (!is.null(file)) ggplot2::ggsave(file, p, width = width, height = height)
    return(p)
  }

  if (is.null(label_min)) {
    label_min <- 0.05 * if (normalize == "none") max(abs(df$intensity)) else scale
  }

  label_df <- NULL
  if (is.finite(label_top_n) && label_top_n > 0) {
    label_df <- df %>%
      dplyr::mutate(abs_int = abs(intensity)) %>%
      dplyr::filter(abs_int >= label_min) %>%
      dplyr::group_by(sample) %>%
      dplyr::slice_max(order_by = abs_int, n = label_top_n, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        label = format(round(mz, label_digits), nsmall = label_digits, trim = TRUE),
        vjust = ifelse(intensity >= 0, -0.2, 1.2)
      )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mz)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = mz, y = 0, yend = intensity, color = sample, alpha = alpha),
      linewidth = 0.4
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = "m/z", y = "Relative intensity (mirror)")

  if (!is.null(label_df) && nrow(label_df) > 0) {
    p <- p + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(y = intensity, label = label, color = sample, vjust = vjust),
      size = 3,
      show.legend = FALSE
    )
  }

  if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }

  if (!is.null(file)) {
    ggplot2::ggsave(file, p, width = width, height = height)
  }
  p
}

#' Plot a mirror spectrum comparison for two IDs from a spectra object
#'
#' @param spectra A spectra object returned by `build_spectra()` (or `res$spectra`).
#' @param idA,idB Compound IDs (names of the spectra lists).
#' @param channel Which spectrum list to use.
#' @param params Optional params list; used to provide defaults (e.g., ppm tolerance).
#' @param ppm Optional ppm tolerance; if NULL, uses `params$tol_ppm` or 20.
#' @param ... Passed to `plot_spectrum_mirror()`.
#' @return A ggplot object.
#' @export
plot_pair_spectrum_mirror <- function(
  spectra,
  idA,
  idB,
  channel = c("frag", "loss", "loss_typ", "loss_anchor", "loss_pair", "loss_anchor_typ", "loss_pair_typ"),
  params = NULL,
  ppm = NULL,
  ...
) {
  require_ggplot2()
  channel <- match.arg(channel)
  if (is.null(ppm)) {
    ppm <- if (!is.null(params) && !is.null(params$tol_ppm) && is.finite(params$tol_ppm)) params$tol_ppm else 20
  }
  get_list <- function(ch) {
    if (ch == "frag") return(spectra$frag_list)
    if (ch == "loss") return(spectra$loss_list)
    if (ch == "loss_typ") return(spectra$loss_typ_list)
    if (ch == "loss_anchor") return(spectra$loss_anchor_list)
    if (ch == "loss_pair") return(spectra$loss_pair_list)
    if (ch == "loss_anchor_typ") return(spectra$loss_anchor_typ_list)
    if (ch == "loss_pair_typ") return(spectra$loss_pair_typ_list)
    NULL
  }
  lst <- get_list(channel)
  if (is.null(lst)) stop("Requested channel is not available in 'spectra' (NULL): ", channel)
  if (!idA %in% names(lst)) stop("idA not found in spectra list: ", idA)
  if (!idB %in% names(lst)) stop("idB not found in spectra list: ", idB)
  specA <- lst[[idA]]
  specB <- lst[[idB]]

  subtitle <- NULL
  if (!is.null(spectra$df_spec)) {
    df <- spectra$df_spec
    rowA <- df[df$id == idA, , drop = FALSE]
    rowB <- df[df$id == idB, , drop = FALSE]
    if (nrow(rowA) == 1 && nrow(rowB) == 1) {
      fmt <- function(x) if (is.na(x) || !nzchar(as.character(x))) "NA" else as.character(x)
      subtitle <- paste0(
        idA, ": class=", fmt(rowA$compound_class), ", RI=", fmt(rowA$RI), ", deriv=", fmt(rowA$derivatization_type),
        "  |  ",
        idB, ": class=", fmt(rowB$compound_class), ", RI=", fmt(rowB$RI), ", deriv=", fmt(rowB$derivatization_type)
      )
    }
  }

  plot_spectrum_mirror(
    specA, specB,
    ppm = ppm,
    title = paste0("Mirror spectrum: ", channel),
    subtitle = subtitle,
    sample_names = c(idA, idB),
    ...
  )
}

#' Compare two compounds across multiple channels (panel)
#'
#' Convenience wrapper that generates multiple mirror plots (e.g., fragment and loss)
#' and arranges them into a single figure when `gridExtra` is available.
#'
#' @param res Result object returned by `run_eihrms_similarity()`.
#' @param idA,idB Compound IDs.
#' @param channels Character vector of channels to plot.
#' @param combine If TRUE and `gridExtra` is installed, return a combined grob.
#' @param ncol Number of columns when combining.
#' @param ... Passed to `plot_pair_spectrum_mirror()`.
#' @return A list of ggplot objects, or a combined grob if `combine=TRUE`.
#' @export
plot_pair_spectrum_compare <- function(
  res,
  idA,
  idB,
  channels = c("frag", "loss"),
  combine = TRUE,
  ncol = 1,
  ...
) {
  require_ggplot2()
  if (is.null(res$spectra)) stop("res$spectra is missing. Provide the result from run_eihrms_similarity().")
  params <- res$params
  plots <- lapply(channels, function(ch) {
    plot_pair_spectrum_mirror(res$spectra, idA, idB, channel = ch, params = params, ...)
  })
  names(plots) <- channels
  if (isTRUE(combine) && requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::arrangeGrob(grobs = plots, ncol = ncol))
  }
  plots
}

#' Plot contribution breakdown for a pair
#'
#' Provides a compact, report-friendly visualization of which components drive the
#' (squared) distance for a given pair, using the contribution fields returned by
#' `distance_breakdown_pair()`.
#'
#' By default (`style='auto'`), this function draws a **two-level** view when possible:
#' a total bar (`frag` vs `loss`) plus a second bar showing the *composition of loss*
#' (e.g., anchored vs pairwise, and/or raw vs typical) **normalized within loss**.
#' This tends to remain readable even when `loss` is a small fraction of the total.
#'
#' @param det A list returned by `distance_breakdown_pair()`.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param label_min Minimum fraction (0-1) required to label a segment.
#' @param style Plot style: `auto` (default), `two_level`, or `simple`.
#' @param loss_breakdown_min Minimum `loss` fraction (of total) required to show the loss-composition bar.
#' @param show_legend If TRUE, show a legend.
#' @param base_size Base font size passed to `ggplot2::theme_minimal()`.
#' @return A ggplot object.
#' @export
plot_pair_contribution_bar <- function(det,
                                      title = "Distance contribution (squared)",
                                      subtitle = NULL,
                                      label_min = 0.08,
                                      style = c("auto", "two_level", "simple"),
                                      loss_breakdown_min = 0.05,
                                      show_legend = TRUE,
                                      base_size = 10) {
  require_ggplot2()
  style <- match.arg(style)
  if (is.null(det) || !is.list(det)) stop("det must be a list returned by distance_breakdown_pair().")

  pct <- function(x) paste0(round(x * 100), "%")
  getv <- function(x) {
    if (!is.null(x) && is.finite(x)) return(as.numeric(x))
    0
  }

  # Total contributions (squared-distance fractions that should sum to ~1)
  c_frag <- getv(det$c_frag_total)
  c_loss <- getv(det$c_loss_total)

  # Loss components (fractions of total)
  loss_parts <- character(0)
  loss_vals_total <- numeric(0)

  # Split loss mode (anchored A / pairwise B)
  if (!is.null(det$c_anchor_total) && !is.null(det$c_pair_total)) {
    a_raw <- getv(det$c_anchor_raw_total)
    a_typ <- getv(det$c_anchor_typ_total)
    b_raw <- getv(det$c_pair_raw_total)
    b_typ <- getv(det$c_pair_typ_total)

    loss_parts <- c("loss_anchor_raw", "loss_anchor_typ", "loss_pair_raw", "loss_pair_typ")
    loss_vals_total <- c(a_raw, a_typ, b_raw, b_typ)

    # Drop empty segments
    keep <- (loss_vals_total > 0)
    loss_parts <- loss_parts[keep]
    loss_vals_total <- loss_vals_total[keep]

    # If typical parts are absent, collapse to anchored/pairwise totals
    if (!any(grepl("_typ$", loss_parts))) {
      loss_parts <- c("loss_anchor", "loss_pair")
      loss_vals_total <- c(getv(det$c_anchor_total), getv(det$c_pair_total))
    }

  } else if (!is.null(det$c_loss_raw_total) && !is.null(det$c_loss_typ_total)) {
    # Combined loss mode (raw vs typical)
    loss_parts <- c("loss_raw", "loss_typical")
    loss_vals_total <- c(getv(det$c_loss_raw_total), getv(det$c_loss_typ_total))

    keep <- (loss_vals_total > 0)
    loss_parts <- loss_parts[keep]
    loss_vals_total <- loss_vals_total[keep]

    # If typical part is absent, collapse to overall loss
    if (!any(loss_parts == "loss_typical")) {
      loss_parts <- c("loss")
      loss_vals_total <- c(c_loss)
    }
  } else {
    # Fallback: no deeper loss decomposition available
    loss_parts <- c("loss")
    loss_vals_total <- c(c_loss)
  }

  # Normalize total bar if required
  total_vals <- c(c_frag, c_loss)
  total_vals[!is.finite(total_vals)] <- 0
  s_total <- sum(total_vals)
  if (!is.finite(s_total) || s_total <= 0) {
    total_parts <- c("n/a")
    total_vals <- c(1)
    c_loss <- 0
  } else {
    total_parts <- c("frag", "loss")
    total_vals <- total_vals / s_total
    c_loss <- total_vals[2]
  }

  # Decide whether to show a two-level view
  show_two <- FALSE
  if (style == "two_level") show_two <- TRUE
  if (style == "simple") show_two <- FALSE
  if (style == "auto") {
    show_two <- (c_loss >= loss_breakdown_min) && length(loss_parts) >= 2 && sum(loss_vals_total) > 0
  }

  # Build plotting data
  df_total <- tibble::tibble(
    level = "Total",
    part = factor(total_parts, levels = total_parts),
    value = as.numeric(total_vals)
  )

  dfs <- list(df_total)

  if (isTRUE(show_two) && c_loss > 0) {
    # Loss composition shown as fractions within loss
    loss_vals_total[!is.finite(loss_vals_total)] <- 0
    s_loss_tot <- sum(loss_vals_total)
    if (is.finite(s_loss_tot) && s_loss_tot > 0) {
      loss_within <- loss_vals_total / s_loss_tot
      lvl <- paste0("Loss breakdown (within loss; loss=", pct(c_loss), " of total)")
      df_loss <- tibble::tibble(
        level = lvl,
        part = factor(loss_parts, levels = loss_parts),
        value = as.numeric(loss_within)
      )
      dfs[[length(dfs) + 1]] <- df_loss
    }
  }

  df <- dplyr::bind_rows(dfs) %>%
    dplyr::group_by(level) %>%
    dplyr::mutate(
      cum = cumsum(value),
      mid = cum - value / 2,
      label = ifelse(value >= label_min, pct(value), "")
    ) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(df, ggplot2::aes(x = level, y = value, fill = part)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::coord_flip() +
    ggplot2::geom_text(
      ggplot2::aes(y = mid, label = label),
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(z) paste0(round(z * 100), "%")) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = if (isTRUE(show_legend)) "bottom" else "none"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle)

  # Make legend more compact when present
  if (isTRUE(show_legend)) {
    p <- p + ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))
  }
  p
}

#' Pair summary panel (mirror plots + contribution bar)
#'
#' Generates a compact, report-ready panel for a pair, typically used in the
#' cluster explanation reports.
#'
#' @param spectra A spectra object returned by `build_spectra()`.
#' @param idA,idB Compound IDs.
#' @param det Optional list returned by `distance_breakdown_pair()`. If NULL,
#'   the breakdown is computed.
#' @param channels Channels to include as mirror plots (e.g., c("frag","loss")).
#' @param include_contribution If TRUE, append a contribution bar plot.
#' @param params Optional params list.
#' @param ... Passed to `plot_pair_spectrum_mirror()`.
#' @return A grob (requires gridExtra).
#' @export
plot_pair_summary_panel <- function(spectra,
                                    idA,
                                    idB,
                                    det = NULL,
                                    channels = c("frag", "loss"),
                                    include_contribution = TRUE,
                                    params = NULL,
                                    ...) {
  require_ggplot2()
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required to combine plots into a panel. Install it (install.packages('gridExtra')).")
  }

  if (is.null(det)) {
    det <- distance_breakdown_pair(
      idA, idB, spectra,
      params = if (is.null(params)) eihrms_default_params() else params,
      include_typical = FALSE
    )
  }

  channels <- unique(as.character(channels))
  plots <- lapply(channels, function(ch) {
    plot_pair_spectrum_mirror(
      spectra = spectra,
      idA = idA,
      idB = idB,
      channel = ch,
      params = params,
      ...
    )
  })
  names(plots) <- channels

  if (isTRUE(include_contribution)) {
    subtitle <- paste0(idA, " vs ", idB)
    plots[["contribution"]] <- plot_pair_contribution_bar(det, subtitle = subtitle)
  }

  # Heuristic heights: contribution plot is shorter.
  heights <- rep(1, length(plots))
  if ("contribution" %in% names(plots)) heights[length(plots)] <- 0.85

  gridExtra::arrangeGrob(grobs = plots, ncol = 1, heights = heights)
}
