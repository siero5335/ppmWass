# Derivatization auto-detection evaluation utilities.
# These helpers support threshold sweeps, ROC-style summaries,
# and convenient CSV outputs.

#' Normalize derivatization labels
#'
#' @keywords internal
#' @noRd
normalize_deriv_label <- function(x) {
  if (is.null(x)) return(rep("none", 0))
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "none"
  x <- toupper(x)
  x[x %in% c("NA", "NONE", "NO")] <- "none"
  x[x %in% c("TMS", "TRIMETHYLSILYL", "TRIMETHYLSILYLATION")] <- "TMS"
  x[x %in% c("TBDMS", "T-BDMS", "MTBSTFA")] <- "TBDMS"
  x
}

#' Simulate derivatization detection rule from diagnostic ion intensities
#'
#' Mirrors `detect_derivatization_from_frag()` logic using precomputed i73/i147/i115/i57.
#'
#' @param i73,i147,i115,i57 Numeric vectors of diagnostic-ion relative intensities (0..1).
#' @param min_int Minimum intensity threshold.
#' @param strong_int "Strong" intensity threshold.
#' @param allowed Allowed derivatization types: "TMS", "TBDMS", or "BOTH".
#' @return Character vector of predicted types ("none", "TMS", "TBDMS").
#' @export
predict_derivatization_from_diagnostics <- function(i73, i147, i115, i57,
                                                   min_int = 0.02,
                                                   strong_int = 0.08,
                                                   allowed = "BOTH") {
  allowed <- normalize_deriv_label(allowed)
  allow_tms <- allowed %in% c("TMS", "BOTH")
  allow_tbdms <- allowed %in% c("TBDMS", "BOTH")

  i73 <- as.numeric(i73); i147 <- as.numeric(i147); i115 <- as.numeric(i115); i57 <- as.numeric(i57)
  i73[!is.finite(i73) | i73 < 0] <- 0
  i147[!is.finite(i147) | i147 < 0] <- 0
  i115[!is.finite(i115) | i115 < 0] <- 0
  i57[!is.finite(i57) | i57 < 0] <- 0

  if (!allow_tms) {
    i73 <- 0; i147 <- 0
  }
  if (!allow_tbdms) {
    i115 <- 0; i57 <- 0
  }

  tms_score <- i73 + 0.5 * i147
  tbdms_score <- i115 + 0.2 * i57

  tms_hit <- (i73 >= min_int) & (i147 >= (min_int / 2) | i73 >= strong_int)
  tbdms_hit <- (i115 >= min_int) & (i57 >= (min_int / 2) | i115 >= strong_int)

  out <- rep("none", length(i73))
  both <- which(tms_hit & tbdms_hit)
  if (length(both) > 0) {
    out[both] <- ifelse(tms_score[both] >= tbdms_score[both], "TMS", "TBDMS")
  }
  out[which(tms_hit & !tbdms_hit)] <- "TMS"
  out[which(tbdms_hit & !tms_hit)] <- "TBDMS"
  out
}

#' Confusion-matrix summary (tidy)
#'
#' @param truth True labels.
#' @param pred Predicted labels.
#' @param labels Order of labels.
#' @return A list with `confusion` (table) and `metrics` (overall accuracy + per-class recall).
#' @export
confusion_summary <- function(truth, pred, labels = c("none", "TMS", "TBDMS")) {
  truth <- normalize_deriv_label(truth)
  pred <- normalize_deriv_label(pred)

  truth[!truth %in% labels] <- "none"
  pred[!pred %in% labels] <- "none"

  tab <- table(factor(truth, levels = labels), factor(pred, levels = labels))
  acc <- sum(diag(tab)) / max(1, sum(tab))

  recall <- rep(NA_real_, length(labels))
  names(recall) <- labels
  for (lab in labels) {
    denom <- sum(tab[lab, , drop = TRUE])
    recall[lab] <- ifelse(denom > 0, tab[lab, lab] / denom, NA_real_)
  }

  list(
    confusion = as.data.frame.matrix(tab),
    metrics = data.frame(
      accuracy = acc,
      recall_none = unname(recall["none"]),
      recall_TMS = unname(recall["TMS"]),
      recall_TBDMS = unname(recall["TBDMS"]),
      stringsAsFactors = FALSE
    )
  )
}

#' Binary ROC curve from scores
#'
#' @param scores Numeric score (higher = more likely positive).
#' @param truth_positive Logical vector for true positives.
#' @param thresholds Optional thresholds; default uses unique scores.
#' @return Data frame with threshold, TPR, FPR, precision, recall.
#' @export
roc_curve_binary <- function(scores, truth_positive, thresholds = NULL) {
  s <- as.numeric(scores)
  y <- as.logical(truth_positive)
  ok <- is.finite(s) & !is.na(y)
  s <- s[ok]
  y <- y[ok]

  if (length(s) == 0) {
    return(data.frame(threshold = numeric(0), tpr = numeric(0), fpr = numeric(0),
                      precision = numeric(0), recall = numeric(0),
                      tp = integer(0), fp = integer(0), tn = integer(0), fn = integer(0)))
  }

  if (is.null(thresholds)) {
    thresholds <- sort(unique(s), decreasing = TRUE)
    thresholds <- c(Inf, thresholds, -Inf)
  } else {
    thresholds <- sort(unique(as.numeric(thresholds)), decreasing = TRUE)
    thresholds <- c(Inf, thresholds[is.finite(thresholds)], -Inf)
  }

  out <- vector("list", length(thresholds))
  P <- sum(y)
  N <- sum(!y)

  for (i in seq_along(thresholds)) {
    th <- thresholds[i]
    pred_pos <- s >= th
    tp <- sum(pred_pos & y)
    fp <- sum(pred_pos & !y)
    fn <- sum(!pred_pos & y)
    tn <- sum(!pred_pos & !y)

    tpr <- ifelse(P > 0, tp / P, NA_real_)
    fpr <- ifelse(N > 0, fp / N, NA_real_)
    precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
    recall <- tpr

    out[[i]] <- data.frame(
      threshold = th,
      tpr = tpr,
      fpr = fpr,
      precision = precision,
      recall = recall,
      tp = tp, fp = fp, tn = tn, fn = fn,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

#' AUC (trapezoid) from an ROC curve
#'
#' @param roc_df Output of roc_curve_binary().
#' @return Numeric AUC.
#' @export
roc_auc <- function(roc_df) {
  if (is.null(roc_df) || nrow(roc_df) < 2) return(NA_real_)
  df <- roc_df
  df <- df[is.finite(df$fpr) & is.finite(df$tpr), , drop = FALSE]
  if (nrow(df) < 2) return(NA_real_)
  # Sort by FPR
  df <- df[order(df$fpr, df$tpr), , drop = FALSE]
  x <- df$fpr
  y <- df$tpr
  sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)]) / 2)
}

#' Threshold sweep for derivatization rule (min_int / strong_int)
#'
#' @param df_spec Spectra metadata (needs deriv_i73/deriv_i147/deriv_i115/deriv_i57 columns).
#' @param truth True labels ("none"/"TMS"/"TBDMS").
#' @param min_int_grid Grid of min_int values.
#' @param strong_int_grid Grid of strong_int values.
#' @param allowed Allowed type(s): "TMS"/"TBDMS"/"BOTH".
#' @return Data frame with accuracy and per-class recall over the grid.
#' @export
sweep_derivatization_rule <- function(df_spec, truth,
                                     min_int_grid = seq(0.01, 0.05, by = 0.005),
                                     strong_int_grid = seq(0.05, 0.15, by = 0.01),
                                     allowed = "BOTH") {
  if (is.null(df_spec$deriv_i73) || is.null(df_spec$deriv_i147) ||
      is.null(df_spec$deriv_i115) || is.null(df_spec$deriv_i57)) {
    stop("df_spec must include deriv_i73/deriv_i147/deriv_i115/deriv_i57 (enable derivatization auto mode, or rebuild spectra with updated code).")
  }

  truth <- normalize_deriv_label(truth)
  i73 <- df_spec$deriv_i73
  i147 <- df_spec$deriv_i147
  i115 <- df_spec$deriv_i115
  i57 <- df_spec$deriv_i57

  grid <- expand.grid(min_int = min_int_grid, strong_int = strong_int_grid)

  res <- vector("list", nrow(grid))
  for (k in seq_len(nrow(grid))) {
    mi <- grid$min_int[k]
    si <- grid$strong_int[k]
    pred <- predict_derivatization_from_diagnostics(i73, i147, i115, i57, min_int = mi, strong_int = si, allowed = allowed)
    cs <- confusion_summary(truth, pred)
    res[[k]] <- cbind(grid[k, , drop = FALSE], cs$metrics)
  }
  do.call(rbind, res)
}

#' Evaluate derivatization detection and write CSV outputs
#'
#' @param spectra Output of build_spectra().
#' @param truth Optional named/unnamed vector of true labels.
#' @param output_dir Output folder (created if missing).
#' @param write_outputs If TRUE, write CSV files.
#' @param allowed Allowed type(s) for rule sweep.
#' @return A list with data frames.
#' @export
run_derivatization_evaluation <- function(spectra,
                                         truth = NULL,
                                         output_dir = ".",
                                         write_outputs = TRUE,
                                         allowed = "BOTH") {
  if (is.null(spectra$df_spec)) stop("spectra must be the result of build_spectra().")
  df <- spectra$df_spec

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pred_summary <- df[, c("id", "known", "RI", "derivatization_type", "deriv_score_tms", "deriv_score_tbdms",
                        "deriv_i73", "deriv_i147", "deriv_i115", "deriv_i57"), drop = FALSE]

  out <- list(pred_summary = pred_summary)

  if (isTRUE(write_outputs)) {
    readr::write_csv(pred_summary, file.path(output_dir, "derivatization_scores.csv"))
  }

  if (!is.null(truth)) {
    # Align truth vector to df$id
    if (!is.null(names(truth))) {
      truth <- truth[df$id]
    }
    truth <- normalize_deriv_label(truth)
    out$confusion <- confusion_summary(truth, df$derivatization_type)

    # Rule sweep (min_int/strong_int)
    out$rule_sweep <- sweep_derivatization_rule(df, truth, allowed = allowed)

    # ROC-style curves for score-based detection (one-vs-rest)
    roc_tms <- roc_curve_binary(df$deriv_score_tms, truth_positive = truth %in% c("TMS", "BOTH"))
    roc_tbdms <- roc_curve_binary(df$deriv_score_tbdms, truth_positive = truth %in% c("TBDMS", "BOTH"))
    out$roc_tms <- roc_tms
    out$roc_tbdms <- roc_tbdms
    out$auc <- data.frame(
      auc_tms = roc_auc(roc_tms),
      auc_tbdms = roc_auc(roc_tbdms),
      stringsAsFactors = FALSE
    )

    if (isTRUE(write_outputs)) {
      readr::write_csv(out$confusion$confusion %>% tibble::rownames_to_column("truth"), file.path(output_dir, "derivatization_confusion.csv"))
      readr::write_csv(out$confusion$metrics, file.path(output_dir, "derivatization_confusion_metrics.csv"))
      readr::write_csv(out$rule_sweep, file.path(output_dir, "derivatization_rule_sweep.csv"))
      readr::write_csv(out$roc_tms, file.path(output_dir, "derivatization_roc_tms.csv"))
      readr::write_csv(out$roc_tbdms, file.path(output_dir, "derivatization_roc_tbdms.csv"))
      readr::write_csv(out$auc, file.path(output_dir, "derivatization_auc.csv"))
    }
  }

  out
}
