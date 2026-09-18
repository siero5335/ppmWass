#' @keywords internal
#' @noRd
CLASSIFY_RULES_DEFAULT <- list(
  alkane_min_int = 0.05,
  alkane_count = 4,

  alkene_min_int = 0.05,
  alkene_count = 3,
  alkene_41_min_int = 0.08,

  aromatic_77_min_int = 0.03,
  aromatic_78_min_int = 0.05,

  alkylbenzene_91_min_int = 0.08,

  monoterpene_93_min_int = 0.08,
  monoterpene_136_min_int = 0.02,

  sesquiterpene_161_min_int = 0.03,
  sesquiterpene_204_min_int = 0.02,

  methyl_ketone_43_min_int = 0.20,
  methyl_ketone_58_min_int = 0.05,

  acetyl_43_min_int = 0.10,
  acetyl_ratio = 1.50,

  aldehyde_44_min_int = 0.08,
  aldehyde_29_min_int = 0.05,
  aldehyde_44hydro_min_int = 0.05,
  aldehyde_44hydro_ratio = 0.50,

  alcohol_31_min_int = 0.03,
  alcohol_45_min_int = 0.03,

  acetate_43_min_int = 0.15,
  acetate_61_min_int = 0.03,

  sulfur_47_min_int = 0.05,

  furan_68_min_int = 0.10,
  furan_39_min_int = 0.10,
  furan_68hydro_min_int = 0.10
)

#' Detect Compound Class from EI-HRMS Peaks
#'
#' @param peaks Two-column matrix with m/z and intensity.
#' @param ppm_tol Mass tolerance in ppm.
#' @param threshold Minimum relative intensity threshold.
#' @param rules Optional named list overriding entries in `CLASSIFY_RULES_DEFAULT`.
#' @return A string of inferred class labels.
#' @export
detect_compound_class <- function(peaks, ppm_tol = 15, threshold = 0.03, rules = NULL) {
  if (nrow(peaks) == 0) return("unknown")

  # Classification rule thresholds (defaults can be overridden via `rules`)
  if (is.null(rules)) {
    rules <- CLASSIFY_RULES_DEFAULT
  } else {
    if (!is.list(rules)) stop('rules must be a list when provided.')
    rules <- utils::modifyList(CLASSIFY_RULES_DEFAULT, rules)
  }

  mz <- peaks[, 1]
  it <- peaks[, 2]

  has_peak <- function(target_mz, min_int = threshold, ppm = ppm_tol) {
    tol <- target_mz * ppm * 1e-6
    any(abs(mz - target_mz) < tol & it > min_int)
  }

  get_intensity <- function(target_mz, ppm = ppm_tol) {
    tol <- target_mz * ppm * 1e-6
    idx <- which(abs(mz - target_mz) < tol)
    if (length(idx) == 0) return(0)
    max(it[idx])
  }

  count_matches <- function(target_mzs, min_int = threshold, ppm = ppm_tol) {
    sum(sapply(target_mzs, function(m) has_peak(m, min_int, ppm)))
  }

  classes <- c()
  confidence <- c()

  alkane_ions <- c(43.0548, 57.0704, 71.0861, 85.1017, 99.1174)
  alkane_count <- count_matches(alkane_ions, min_int = rules$alkane_min_int)
  if (alkane_count >= rules$alkane_count) {
    classes <- c(classes, "alkane")
    confidence <- c(confidence, get_intensity(57.0704))
  }

  alkene_ions <- c(41.0391, 55.0548, 69.0704, 83.0861)
  alkene_count <- count_matches(alkene_ions, min_int = rules$alkene_min_int)
  if (alkene_count >= rules$alkene_count && has_peak(41.0391, min_int = rules$alkene_41_min_int)) {
    classes <- c(classes, "alkene")
    confidence <- c(confidence, get_intensity(41.0391))
  }

  if (has_peak(77.0391, min_int = rules$aromatic_77_min_int) || has_peak(78.0470, min_int = rules$aromatic_78_min_int)) {
    classes <- c(classes, "aromatic")
    confidence <- c(confidence, max(get_intensity(77.0391), get_intensity(78.0470)))
  }

  if (has_peak(91.0548, min_int = rules$alkylbenzene_91_min_int)) {
    classes <- c(classes, "alkylbenzene")
    confidence <- c(confidence, get_intensity(91.0548))
  }

  if (has_peak(93.0704, min_int = rules$monoterpene_93_min_int) && has_peak(136.1252, min_int = rules$monoterpene_136_min_int)) {
    classes <- c(classes, "monoterpene_HC")
    confidence <- c(confidence, get_intensity(93.0704))
  }

  if (has_peak(161.1330, min_int = rules$sesquiterpene_161_min_int) && has_peak(204.1878, min_int = rules$sesquiterpene_204_min_int)) {
    classes <- c(classes, "sesquiterpene")
    confidence <- c(confidence, get_intensity(161.1330))
  }

  if (has_peak(43.0184, min_int = rules$methyl_ketone_43_min_int) && has_peak(58.0419, min_int = rules$methyl_ketone_58_min_int)) {
    classes <- c(classes, "methyl_ketone")
    confidence <- c(confidence, get_intensity(43.0184))
  }

  if (has_peak(43.0184, min_int = rules$acetyl_43_min_int) &&
      get_intensity(43.0184) > get_intensity(43.0548) * rules$acetyl_ratio) {
    if (!("methyl_ketone" %in% classes)) {
      classes <- c(classes, "acetyl_compound")
      confidence <- c(confidence, get_intensity(43.0184))
    }
  }

  if (has_peak(44.0262, min_int = rules$aldehyde_44_min_int) && has_peak(29.0027, min_int = rules$aldehyde_29_min_int)) {
    if (get_intensity(44.0262) > get_intensity(44.0626) * rules$aldehyde_44hydro_ratio ||
        !has_peak(44.0626, min_int = rules$aldehyde_44hydro_min_int)) {
      classes <- c(classes, "aldehyde")
      confidence <- c(confidence, get_intensity(44.0262))
    }
  }

  if (has_peak(31.0184, min_int = rules$alcohol_31_min_int) && has_peak(45.0340, min_int = rules$alcohol_45_min_int)) {
    classes <- c(classes, "alcohol")
    confidence <- c(confidence, get_intensity(31.0184) + get_intensity(45.0340))
  }

  if (has_peak(43.0184, min_int = rules$acetate_43_min_int) && has_peak(61.0290, min_int = rules$acetate_61_min_int)) {
    classes <- c(classes, "acetate_ester")
    confidence <- c(confidence, get_intensity(43.0184))
  }

  if (has_peak(46.9955, min_int = rules$sulfur_47_min_int) || has_peak(48.0034, min_int = rules$sulfur_47_min_int)) {
    classes <- c(classes, "sulfur_compound")
    confidence <- c(confidence, max(get_intensity(46.9955), get_intensity(48.0034)))
  }

  if (has_peak(68.0262, min_int = rules$furan_68_min_int) && has_peak(39.0235, min_int = rules$furan_39_min_int)) {
    if (!has_peak(68.0626, min_int = rules$furan_68hydro_min_int) || get_intensity(68.0262) > get_intensity(68.0626)) {
      classes <- c(classes, "furan")
      confidence <- c(confidence, get_intensity(68.0262))
    }
  }

  if (length(classes) == 0) return("unclassified")

  ord <- order(confidence, decreasing = TRUE)
  classes <- classes[ord]
  paste(utils::head(classes, 3), collapse = ";")
}
