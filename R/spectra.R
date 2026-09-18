#' @keywords internal
#' @noRd
parse_ei <- function(x) {
  if (length(x) != 1 || is.na(x) || !nzchar(x)) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity"))))
  }
  pairs <- strsplit(x, " ", fixed = TRUE)[[1]]
  mz <- as.numeric(sub(":.*$", "", pairs))
  int <- as.numeric(sub("^.*:", "", pairs))
  m <- cbind(mz = mz, intensity = int)
  m <- m[is.finite(m[, 1]) & is.finite(m[, 2]) & m[, 1] > 0 & m[, 2] > 0, , drop = FALSE]
  m[order(m[, 1]), , drop = FALSE]
}

#' @keywords internal
#' @noRd
prep_peaks <- function(peaks, params) {
  if (nrow(peaks) == 0) return(peaks)

  mz <- peaks[, 1]
  int <- peaks[, 2]

  keep <- mz >= params$min_mz & mz <= params$max_mz
  mz <- mz[keep]
  int <- int[keep]
  if (length(mz) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity"))))
  }

  noise_threshold <- max(int) * params$noise_thr
  keep <- int >= noise_threshold
  mz <- mz[keep]
  int <- int[keep]
  if (length(mz) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity"))))
  }

  if (length(mz) > 1) {
    ord <- order(mz)
    mz <- mz[ord]
    int <- int[ord]

    n <- length(mz)
    merged_mz <- vector("list", n)
    merged_int <- vector("list", n)
    k <- 0
    i <- 1
    while (i <= n) {
      tol <- mz[i] * params$centroid_ppm * 1e-6
      j <- i
      while (j <= n && abs(mz[j] - mz[i]) <= tol) {
        j <- j + 1
      }
      idx <- i:(j - 1)
      group_int <- int[idx]
      group_mz <- mz[idx]
      k <- k + 1
      merged_mz[[k]] <- sum(group_mz * group_int) / sum(group_int)
      merged_int[[k]] <- sum(group_int)
      i <- j
    }
    mz <- unlist(merged_mz[1:k], use.names = FALSE)
    int <- unlist(merged_int[1:k], use.names = FALSE)
  }

  if (length(mz) > params$topK) {
    ord <- order(int, decreasing = TRUE)[1:params$topK]
    mz <- mz[ord]
    int <- int[ord]
    ord <- order(mz)
    mz <- mz[ord]
    int <- int[ord]
  }

  int <- int / sum(int)

  if (params$use_sqrt) {
    int <- sqrt(int)
    int <- int / sum(int)
  }

  cbind(mz = mz, intensity = int)
}

#' Resolve derivatization handling mode
#'
#' Legacy behavior: if derivatization_mode is NULL, fall back to use_derivatization_losses.
#'
#' @keywords internal
#' @noRd
resolve_derivatization_mode <- function(params) {
  mode <- params$derivatization_mode
  if (!is.null(mode)) {
    mode <- tolower(as.character(mode))
    if (!mode %in% c("off", "manual", "auto")) mode <- "off"
    return(mode)
  }
  if (!is.null(params$use_derivatization_losses) && isTRUE(params$use_derivatization_losses)) {
    return("manual")
  }
  "off"
}

#' @keywords internal
#' @noRd
normalize_derivatization_type <- function(x) {
  if (is.null(x) || !nzchar(as.character(x))) return("TMS")
  y <- toupper(as.character(x))
  if (y %in% c("TMS", "TBDMS", "BOTH")) return(y)
  if (y %in% c("TMS+TBDMS", "TBDMS+TMS")) return("BOTH")
  y
}

#' @keywords internal
#' @noRd
get_derivatization_losses <- function(deriv_type) {
  dt <- normalize_derivatization_type(deriv_type)
  if (dt == "TMS") return(TYPICAL_LOSSES_DERIV_TMS)
  if (dt == "TBDMS") return(TYPICAL_LOSSES_DERIV_TBDMS)
  if (dt == "BOTH") return(c(TYPICAL_LOSSES_DERIV_TMS, TYPICAL_LOSSES_DERIV_TBDMS))
  list()
}

#' Get typical neutral losses for this run / spectrum
#'
#' Returns a combined list of typical neutral losses, optionally including
#' extended losses and derivatization-specific losses (e.g., TMS, TBDMS).
#'
#' In auto mode, derivatization losses can be enabled per spectrum based on the
#' detected derivatization type.
#'
#' @param deriv_type Optional derivatization type for this spectrum ("TMS", "TBDMS", "BOTH", or "none").
#' @param purpose "mref" (for Mref estimation) or "projection" (for typical-loss projection).
#'
#' @keywords internal
#' @noRd
get_typical_losses <- function(params, deriv_type = NULL, purpose = c("mref", "projection")) {
  purpose <- match.arg(purpose)
  losses <- TYPICAL_LOSSES_CORE

  if (!is.null(params$use_extended_losses) && isTRUE(params$use_extended_losses)) {
    losses <- c(losses, TYPICAL_LOSSES_EXTENDED)
  }

  mode <- resolve_derivatization_mode(params)
  allowed <- normalize_derivatization_type(params$derivatization_type)

  if (mode == "manual") {
    # Apply derivatization losses globally (all spectra)
    losses <- c(losses, get_derivatization_losses(allowed))
  } else if (mode == "auto") {
    # Apply derivatization losses per spectrum
    dt <- if (is.null(deriv_type) || !nzchar(as.character(deriv_type))) "none" else toupper(as.character(deriv_type))
    if (dt %in% c("TMS+TBDMS", "TBDMS+TMS")) dt <- "BOTH"

    # Respect allowed type filter
    if (allowed == "TMS" && dt != "TMS") dt <- "none"
    if (allowed == "TBDMS" && dt != "TBDMS") dt <- "none"

    if (dt %in% c("TMS", "TBDMS", "BOTH")) {
      losses <- c(losses, get_derivatization_losses(dt))
    }
  }

  losses
}

#' Universe of typical losses for mapping (class/exact), based on params
#'
#' @keywords internal
#' @noRd
get_typical_losses_universe <- function(params) {
  losses <- TYPICAL_LOSSES_CORE
  if (!is.null(params$use_extended_losses) && isTRUE(params$use_extended_losses)) {
    losses <- c(losses, TYPICAL_LOSSES_EXTENDED)
  }
  mode <- resolve_derivatization_mode(params)
  if (mode %in% c("manual", "auto")) {
    allowed <- normalize_derivatization_type(params$derivatization_type)
    losses <- c(losses, get_derivatization_losses(allowed))
  }
  losses
}

#' Get max intensity near a target m/z within ppm tolerance
#'
#' @keywords internal
#' @noRd
get_intensity_at_mz <- function(peaks, target_mz, ppm = 20) {
  if (is.null(peaks) || nrow(peaks) == 0) return(0)
  mz <- peaks[, 1]
  it <- peaks[, 2]
  tol <- target_mz * ppm * 1e-6
  idx <- which(abs(mz - target_mz) <= tol)
  if (length(idx) == 0) return(0)
  max(it[idx], na.rm = TRUE)
}

#' Auto-detect derivatization (TMS/TBDMS) from fragment ions
#'
#' Uses diagnostic ions with HRMS exact-mass matching:
#' - TMS: m/z 73.0474 and/or 147.0654
#' - TBDMS: m/z 115.0943 (and optional support from 57.0704)
#'
#' @keywords internal
#' @noRd
detect_derivatization_from_frag <- function(frag_peaks, params) {
  ppm <- params$derivatization_auto_ppm
  if (is.null(ppm)) {
    ppm <- if (!is.null(params$class_detection_ppm)) params$class_detection_ppm else params$tol_ppm
  }
  min_int <- if (is.null(params$derivatization_auto_min_int)) 0.02 else params$derivatization_auto_min_int
  strong_int <- if (is.null(params$derivatization_auto_strong_int)) 0.08 else params$derivatization_auto_strong_int

  allowed <- normalize_derivatization_type(params$derivatization_type)
  allow_tms <- allowed %in% c("TMS", "BOTH")
  allow_tbdms <- allowed %in% c("TBDMS", "BOTH")

  i73 <- if (allow_tms) get_intensity_at_mz(frag_peaks, 73.0474, ppm) else 0
  i147 <- if (allow_tms) get_intensity_at_mz(frag_peaks, 147.0654, ppm) else 0
  tms_score <- i73 + 0.5 * i147
  tms_hit <- (i73 >= min_int) && (i147 >= (min_int / 2) || i73 >= strong_int)

  i115 <- if (allow_tbdms) get_intensity_at_mz(frag_peaks, 115.0943, ppm) else 0
  i57 <- if (allow_tbdms) get_intensity_at_mz(frag_peaks, 57.0704, ppm) else 0
  tbdms_score <- i115 + 0.2 * i57
  tbdms_hit <- (i115 >= min_int) && (i57 >= (min_int / 2) || i115 >= strong_int)

  type <- "none"
  if (tms_hit && tbdms_hit) {
    # Choose the stronger signal to keep a single label.
    type <- if (tms_score >= tbdms_score) "TMS" else "TBDMS"
  } else if (tms_hit) {
    type <- "TMS"
  } else if (tbdms_hit) {
    type <- "TBDMS"
  }

  list(
    type = type,
    score_tms = tms_score,
    score_tbdms = tbdms_score,
    i73 = i73,
    i147 = i147,
    i115 = i115,
    i57 = i57
  )
}

#' Downweight intensity vector near target m/z values (no renormalization)
#'
#' This helper is used when downweighting needs to happen *before* peak summarization
#' and truncation (e.g., for neutral-loss spectra). Unlike downweight_mz_matches(),
#' it does not renormalize.
#'
#' @keywords internal
#' @noRd
downweight_intensity_matches <- function(mz_vec, it_vec, target_mz, ppm = 20, weight = 0.5) {
  if (is.null(mz_vec) || length(mz_vec) == 0) return(it_vec)
  if (is.null(it_vec) || length(it_vec) == 0) return(it_vec)
  if (is.null(target_mz) || length(target_mz) == 0) return(it_vec)
  if (!is.finite(weight) || weight < 0 || weight > 1) return(it_vec)

  mz_vec <- as.numeric(mz_vec)
  it_vec <- as.numeric(it_vec)
  it_vec[!is.finite(it_vec) | it_vec < 0] <- 0

  for (m in target_mz) {
    if (!is.finite(m) || m <= 0) next
    tol <- m * ppm * 1e-6
    idx <- which(abs(mz_vec - m) <= tol)
    if (length(idx) > 0) it_vec[idx] <- it_vec[idx] * weight
  }
  it_vec
}

#' Downweight peaks near target m/z values and renormalize
#'
#' @keywords internal
#' @noRd
downweight_mz_matches <- function(spec, target_mz, ppm = 20, weight = 0.5) {
  if (is.null(spec) || nrow(spec) == 0) return(spec)
  if (is.null(target_mz) || length(target_mz) == 0) return(spec)
  if (!is.finite(weight) || weight < 0 || weight > 1) return(spec)

  mz <- spec[, 1]
  it <- spec[, 2]
  it <- downweight_intensity_matches(mz, it, target_mz, ppm = ppm, weight = weight)

  s <- sum(it)
  if (is.finite(s) && s > 0) it <- it / s
  cbind(mz = mz, intensity = it)
}

#' Estimate Mref in HRMS (internal)
#'
#' Chooses the maximum valid candidate, biased toward high m/z peaks.
#'
#' @keywords internal
#' @noRd
estimate_Mref_HRMS <- function(mz, it, ppm = 15, min_rel_int = 0.005, typical_losses = TYPICAL_LOSSES) {
  if (length(mz) == 0) return(list(Mref = NA, detected_losses = list()))

  candidates <- list()
  detected_losses <- list()

  candidates$max_mz <- max(mz)

  hi_region <- mz > max(mz) - 50
  if (any(hi_region & it > min_rel_int)) {
    candidates$hi_strong <- mz[hi_region][which.max(it[hi_region])]
  }

  near_max <- abs(mz - max(mz)) <= 50
  strong_idx <- which(near_max & it > 0.02)

  if (length(strong_idx) > 0) {
    strong_peaks <- mz[strong_idx]
    strong_int <- it[strong_idx]

    for (loss_info in typical_losses) {
      for (i in seq_along(strong_peaks)) {
        potential_M <- strong_peaks[i] + loss_info$exact
        tol <- potential_M * ppm * 1e-6

        if (abs(potential_M - max(mz)) <= max(tol * 3, 1)) {
          candidate_name <- paste0(loss_info$formula, "_from_", round(strong_peaks[i], 4))
          candidates[[candidate_name]] <- potential_M

          detected_losses <- c(detected_losses, list(list(
            formula = loss_info$formula,
            class = loss_info$class,
            exact_loss = loss_info$exact,
            from_peak = strong_peaks[i],
            to_Mref = potential_M,
            intensity = strong_int[i],
            ppm_error = abs(potential_M - max(mz)) / max(mz) * 1e6
          )))
        }
      }
    }
  }

  hi_idx <- which(mz > stats::quantile(mz, 0.85) & it > 0.01)
  if (length(hi_idx) > 0) {
    candidates$quantile_based <- max(mz[hi_idx])
  }

  valid_candidates <- unlist(candidates)
  valid_candidates <- valid_candidates[is.finite(valid_candidates)]

  Mref <- if (length(valid_candidates) == 0) max(mz) else max(valid_candidates)

  list(Mref = Mref, detected_losses = detected_losses)
}

#' @keywords internal
#' @noRd
estimate_Mref <- function(mz, it, min_rel_int = 0.005) {
  result <- estimate_Mref_HRMS(mz, it, ppm = 15, min_rel_int = min_rel_int)
  result$Mref
}

#' @keywords internal
#' @noRd
infer_functional_groups <- function(detected_losses) {
  if (length(detected_losses) == 0) return(character(0))
  classes <- sapply(detected_losses, function(x) x$class)
  unique(classes)
}

#' @keywords internal
#' @noRd
summarize_neutral_losses <- function(detected_losses) {
  if (length(detected_losses) == 0) return("none")
  formulas <- sapply(detected_losses, function(x) x$formula)
  paste(unique(formulas), collapse = ";")
}

#' Compute confidence of Mref estimation
#'
#' Returns a score in the range 0 to 1 based on (i) how many typical-loss-supported candidates were found,
#' (ii) how intense those supporting peaks are, and (iii) whether a strong high-m/z ion exists.
#'
#' @keywords internal
#' @noRd
compute_mref_confidence <- function(mz, it, detected_losses, params) {
  if (length(mz) == 0 || length(it) == 0) return(0)
  mz <- as.numeric(mz)
  it <- as.numeric(it)
  it[!is.finite(it) | it < 0] <- 0

  # Defaults
  nloss_ref <- if (is.null(params$mref_conf_nloss_ref)) 3 else params$mref_conf_nloss_ref
  int_ref <- if (is.null(params$mref_conf_int_ref)) 0.05 else params$mref_conf_int_ref
  hi_win <- if (is.null(params$mref_conf_hi_window)) 50 else params$mref_conf_hi_window
  hi_int_ref <- if (is.null(params$mref_conf_hi_int_ref)) 0.03 else params$mref_conf_hi_int_ref

  w_n <- if (is.null(params$mref_conf_w_nloss)) 0.4 else params$mref_conf_w_nloss
  w_i <- if (is.null(params$mref_conf_w_intensity)) 0.4 else params$mref_conf_w_intensity
  w_h <- if (is.null(params$mref_conf_w_hi)) 0.2 else params$mref_conf_w_hi

  w_sum <- w_n + w_i + w_h
  if (!is.finite(w_sum) || w_sum <= 0) {
    w_n <- 1; w_i <- 0; w_h <- 0; w_sum <- 1
  }
  w_n <- w_n / w_sum
  w_i <- w_i / w_sum
  w_h <- w_h / w_sum

  nloss <- length(detected_losses)
  n_score <- if (is.finite(nloss_ref) && nloss_ref > 0) min(1, nloss / nloss_ref) else 0

  sum_int <- 0
  if (nloss > 0) {
    sum_int <- sum(vapply(detected_losses, function(x) {
      val <- x$intensity
      if (is.null(val) || !is.finite(val)) 0 else val
    }, numeric(1)), na.rm = TRUE)
  }
  i_score <- if (is.finite(int_ref) && int_ref > 0) min(1, sum_int / int_ref) else 0

  hi_int <- 0
  if (is.finite(hi_win) && hi_win > 0) {
    hi_region <- mz > (max(mz, na.rm = TRUE) - hi_win)
    if (any(hi_region)) hi_int <- max(it[hi_region], na.rm = TRUE)
  }
  h_score <- if (is.finite(hi_int_ref) && hi_int_ref > 0) min(1, hi_int / hi_int_ref) else 0

  conf <- w_n * n_score + w_i * i_score + w_h * h_score
  conf <- min(max(conf, 0), 1)
  conf
}

#' @keywords internal
#' @noRd
# Legacy function name retained for API compatibility. The returned secondary
# representation is not a pure neutral-loss spectrum: it pools (A) differences
# from a reference mass estimated heuristically from the observed EI spectrum
# and (B) pairwise absolute differences among high-intensity fragment peaks.
# `reference_mz` is an internal sensitivity-analysis override. Leaving it NULL
# retains the estimated-reference publication method.
build_loss_peaks <- function(peaks, params, return_info = FALSE, deriv_type = NULL,
                             reference_mz = NULL) {
  empty_mat <- matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity")))
  empty_result <- if (return_info) {
    list(
      loss_peaks = empty_mat,
      lossA_peaks = empty_mat,
      lossB_peaks = empty_mat,
      Mref = NA,
      detected_losses = list(),
      functional_groups = character(0),
      mref_confidence = NA_real_
    )
  } else {
    empty_mat
  }

  if (is.null(dim(peaks)) || length(dim(peaks)) != 2L || ncol(peaks) < 2L) {
    stop("peaks must be a matrix-like object with m/z and intensity columns.",
         call. = FALSE)
  }
  if (nrow(peaks) == 0) return(empty_result)
  peak_values <- as.matrix(peaks[, 1L:2L, drop = FALSE])
  if (!is.numeric(peak_values) || any(!is.finite(peak_values))) {
    stop("peaks must contain finite numeric m/z and intensity values; NA removal is not implicit.",
         call. = FALSE)
  }
  if (any(peak_values[, 1L] <= 0) || any(peak_values[, 2L] < 0)) {
    stop("peaks must contain positive m/z and non-negative intensity values.",
         call. = FALSE)
  }
  peaks <- peak_values
  if (!is.null(reference_mz) &&
      (length(reference_mz) != 1L || !is.finite(reference_mz) || reference_mz <= 0)) {
    stop("reference_mz must be NULL or one positive finite value.", call. = FALSE)
  }

  summarize_loss <- function(mz_vec, it_vec) {
    if (length(mz_vec) == 0) return(empty_mat)
    if (sum(it_vec, na.rm = TRUE) <= 0) return(empty_mat)

    loss_df <- tibble::tibble(mz = mz_vec, intensity = it_vec) %>%
      dplyr::mutate(mz_bin = round(mz, 3)) %>%
      dplyr::group_by(mz_bin) %>%
      dplyr::summarise(intensity = sum(intensity), .groups = "drop") %>%
      dplyr::arrange(mz_bin)

    if (nrow(loss_df) > params$loss_max_peaks) {
      loss_df <- loss_df %>%
        dplyr::arrange(dplyr::desc(intensity)) %>%
        utils::head(params$loss_max_peaks) %>%
        dplyr::arrange(mz_bin)
    }

    loss_df <- loss_df %>%
      dplyr::mutate(intensity = intensity / sum(intensity))

    cbind(mz = loss_df$mz_bin, intensity = loss_df$intensity)
  }

  mz <- peaks[, 1]
  it <- peaks[, 2]

  ord <- order(it, decreasing = TRUE)
  k <- min(params$loss_top_peaks, length(ord))
  mz_top <- mz[ord[1:k]]
  it_top <- it[ord[1:k]]

  loss_db <- get_typical_losses(params, deriv_type = deriv_type, purpose = "mref")
  mref_result <- estimate_Mref_HRMS(mz, it, ppm = params$tol_ppm, typical_losses = loss_db)
  Mref <- if (is.null(reference_mz)) mref_result$Mref else as.numeric(reference_mz)
  detected_losses <- mref_result$detected_losses
  functional_groups <- infer_functional_groups(detected_losses)

  mref_conf <- compute_mref_confidence(mz, it, detected_losses, params)

  # A) Reference-anchored mass differences: Mref - fragment m/z
  lossA_mz <- Mref - mz
  lossA_it <- it
  if (isTRUE(params$use_mref_confidence)) {
    # NOTE: In split-loss mode, this uniform scaling is cancelled by within-channel normalization.
    # It still matters for the combined loss spectrum below.
    pwr <- if (is.null(params$mref_conf_power)) 1 else params$mref_conf_power
    if (!is.finite(pwr) || pwr < 0) pwr <- 1
    lossA_it <- lossA_it * (mref_conf^pwr)
  }
  keepA <- is.finite(lossA_mz) & lossA_mz >= params$loss_min & lossA_mz <= params$loss_max
  lossA_mz <- lossA_mz[keepA]
  lossA_it <- lossA_it[keepA]

  # B) Pairwise absolute mass differences among the top fragment peaks. These
  # are invariant to a common shift of the reference anchor because Mref is not
  # used in their construction.
  if (k > 1) {
    dmat <- abs(outer(mz_top, mz_top, "-"))
    wmat <- outer(it_top, it_top, FUN = pmin)
    upper <- upper.tri(dmat, diag = FALSE)
    lossB_mz <- as.numeric(dmat[upper])
    lossB_it <- as.numeric(wmat[upper])
    keepB <- is.finite(lossB_mz) & lossB_mz >= params$loss_min & lossB_mz <= params$loss_max
    lossB_mz <- lossB_mz[keepB]
    lossB_it <- lossB_it[keepB]
  } else {
    lossB_mz <- numeric(0)
    lossB_it <- numeric(0)
  }



  # Optional: downweight derivatization-related neutral losses (e.g., 90 for TMS) before summarization.
  # Doing this before summarize_loss() allows the downweighting to influence peak selection/truncation.
  dw <- params$derivatization_raw_loss_weight
  if (!is.null(dw) && is.finite(dw) && dw < 1) {
    deriv_loss_list <- get_derivatization_losses(deriv_type)
    if (length(deriv_loss_list) > 0) {
      deriv_exacts <- vapply(deriv_loss_list, function(z) z$exact, numeric(1))

      lossA_it <- downweight_intensity_matches(lossA_mz, lossA_it, deriv_exacts, ppm = params$tol_ppm, weight = dw)
      lossB_it <- downweight_intensity_matches(lossB_mz, lossB_it, deriv_exacts, ppm = params$tol_ppm, weight = dw)
    }
  }
  # Build spectra for each channel (normalized separately)
  lossA_peaks <- summarize_loss(lossA_mz, lossA_it)
  lossB_peaks <- summarize_loss(lossB_mz, lossB_it)

  # Combined loss spectrum (backward-compatible, for legacy use)
  all_mz <- c(lossA_mz, lossB_mz)
  all_it <- c(lossA_it, lossB_it)
  loss_peaks <- summarize_loss(all_mz, all_it)

  if (return_info) {
    list(
      loss_peaks = loss_peaks,
      lossA_peaks = lossA_peaks,
      lossB_peaks = lossB_peaks,
      Mref = Mref,
      detected_losses = detected_losses,
      functional_groups = functional_groups,
      mref_confidence = mref_conf
    )
  } else {
    loss_peaks
  }
}



#' Project neutral-loss peaks onto a fixed set of typical losses
#'
#' This creates a compact, interpretable "typical-loss spectrum" that can be blended
#' with the raw neutral-loss spectrum for substructure-oriented similarity.
#'
#' @keywords internal
#' @noRd
project_to_typical_losses <- function(loss_peaks, typical_losses = TYPICAL_LOSSES, ppm = 20, eps = 1e-12, keep_zeros = FALSE, deriv_weight = 1.0) {
  if (is.null(loss_peaks) || nrow(loss_peaks) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity"))))
  }
  if (!is.finite(ppm) || ppm <= 0) stop("ppm must be > 0.")

  exacts <- vapply(typical_losses, function(x) x$exact, numeric(1))
  labels <- vapply(typical_losses, function(x) x$formula, character(1))
  classes <- vapply(typical_losses, function(x) x$class, character(1))
  w_vec <- rep(1, length(exacts))
  if (!is.null(deriv_weight) && is.finite(deriv_weight) && deriv_weight >= 0 && deriv_weight <= 1) {
    w_vec[grepl("^deriv_", classes)] <- deriv_weight
  }
  out <- numeric(length(exacts))

  mz <- loss_peaks[, 1]
  it <- loss_peaks[, 2]
  it[!is.finite(it) | it < 0] <- 0

  for (i in seq_along(mz)) {
    mzi <- mz[i]
    if (!is.finite(mzi) || mzi <= 0) next
    tol <- mzi * ppm * 1e-6
    idx <- which(abs(exacts - mzi) <= tol)
    if (length(idx) == 0) next
    if (length(idx) > 1) {
      idx <- idx[which.min(abs(exacts[idx] - mzi))]
    }
    out[idx] <- out[idx] + it[i] * w_vec[idx]
  }

  s <- sum(out)
  if (!is.finite(s) || s <= eps) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("mz", "intensity"))))
  }
  out <- out / s

  if (!isTRUE(keep_zeros)) {
    keep <- out > 0
    exacts <- exacts[keep]
    out <- out[keep]
    labels <- labels[keep]
  }

  m <- cbind(mz = exacts, intensity = out)
  rownames(m) <- labels
  m[order(m[, 1]), , drop = FALSE]
}


#' @keywords internal
#' @noRd
build_df_spec <- function(df, df_final) {
  as_feature_table(df) %>%
    dplyr::mutate(compound_id = derive_compound_id(.)) %>%
    dplyr::transmute(
      id = compound_id,
      known = dplyr::if_else(`Annotation tag (VS1.0)` == "4", "F", "T"),
      RI = `Average RI`,
      ei = `EI spectrum`
    ) %>%
    dplyr::filter(id %in% df_final$compound_id) %>%
    dplyr::distinct(id, .keep_all = TRUE)
}

#' @keywords internal
#' @noRd
process_single_spectrum <- function(ei_string, params, deriv_mode = resolve_derivatization_mode(params)) {
  peaks <- parse_ei(ei_string)

  deriv_type <- "none"
  deriv_info <- list(type = "none", score_tms = 0, score_tbdms = 0)
  if (deriv_mode == "manual") {
    deriv_type <- normalize_derivatization_type(params$derivatization_type)
    deriv_info$score_tms <- NA_real_
    deriv_info$score_tbdms <- NA_real_
  } else if (deriv_mode == "auto") {
    params_detect <- params
    params_detect$use_sqrt <- FALSE
    frag_detect <- prep_peaks(peaks, params_detect)
    deriv_info <- detect_derivatization_from_frag(frag_detect, params)
    deriv_type <- deriv_info$type
  }

  frag <- prep_peaks(peaks, params)

  fw <- params$derivatization_frag_weight
  if (!is.null(fw) && is.finite(fw) && fw < 1 && deriv_type != "none") {
    diag_mz <- numeric(0)
    if (deriv_type %in% c("TMS", "BOTH")) diag_mz <- c(diag_mz, 73.0474, 147.0654)
    if (deriv_type %in% c("TBDMS", "BOTH")) diag_mz <- c(diag_mz, 57.0704, 115.0943)
    ppm_dw <- if (is.null(params$derivatization_auto_ppm)) params$class_detection_ppm else params$derivatization_auto_ppm
    if (is.null(ppm_dw)) ppm_dw <- params$tol_ppm
    frag <- downweight_mz_matches(frag, unique(diag_mz), ppm = ppm_dw, weight = fw)
  }

  loss_result <- build_loss_peaks(frag, params, return_info = TRUE, deriv_type = deriv_type)

  use_typical_loss <- isTRUE(params$use_typical_loss)
  use_split_loss <- isTRUE(params$use_split_loss)
  loss_typ <- NULL
  loss_anchor_typ <- NULL
  loss_pair_typ <- NULL
  if (use_typical_loss) {
    ppm_typ <- if (is.null(params$loss_typical_ppm)) params$tol_ppm else params$loss_typical_ppm
    typical_losses_this <- get_typical_losses(params, deriv_type = deriv_type, purpose = "projection")
    dw_typ <- if (is.null(params$derivatization_typical_weight)) 1 else params$derivatization_typical_weight
    loss_typ <- project_to_typical_losses(
      loss_result$loss_peaks,
      typical_losses = typical_losses_this,
      ppm = ppm_typ,
      deriv_weight = dw_typ
    )
    if (use_split_loss) {
      loss_anchor_typ <- project_to_typical_losses(
        loss_result$lossA_peaks,
        typical_losses = typical_losses_this,
        ppm = ppm_typ,
        deriv_weight = dw_typ
      )
      loss_pair_typ <- project_to_typical_losses(
        loss_result$lossB_peaks,
        typical_losses = typical_losses_this,
        ppm = ppm_typ,
        deriv_weight = dw_typ
      )
    }
  }

  list(
    frag = frag,
    loss = loss_result$loss_peaks,
    loss_anchor = loss_result$lossA_peaks,
    loss_pair = loss_result$lossB_peaks,
    loss_typ = loss_typ,
    loss_anchor_typ = loss_anchor_typ,
    loss_pair_typ = loss_pair_typ,
    deriv_type = deriv_type,
    deriv_info = deriv_info,
    mref = loss_result$Mref,
    mref_confidence = if (is.null(loss_result$mref_confidence)) NA_real_ else loss_result$mref_confidence,
    neutral_losses = summarize_neutral_losses(loss_result$detected_losses),
    functional_groups = paste(loss_result$functional_groups, collapse = ";"),
    deriv_i73 = if (is.null(deriv_info$i73)) NA_real_ else deriv_info$i73,
    deriv_i147 = if (is.null(deriv_info$i147)) NA_real_ else deriv_info$i147,
    deriv_i115 = if (is.null(deriv_info$i115)) NA_real_ else deriv_info$i115,
    deriv_i57 = if (is.null(deriv_info$i57)) NA_real_ else deriv_info$i57
  )
}

#' Build Spectra Lists and Class Annotations
#'
#' @param df Raw MS-DIAL tibble.
#' @param df_final Cleaned quant table from prepare_quant_data().
#' @param params Parameter list.
#' @param progress Print progress every 100 spectra.
#' @return A list with df_spec, frag_list, loss_list, and summaries.
#' @export
build_spectra <- function(df, df_final, params, progress = TRUE) {
  df_spec <- build_df_spec(df, df_final)

  frag_list <- vector("list", nrow(df_spec))
  loss_list <- vector("list", nrow(df_spec))

  use_typical_loss <- isTRUE(params$use_typical_loss)
  use_split_loss <- isTRUE(params$use_split_loss)

  # Combined loss (legacy) + optional split losses (A = anchored; B = pairwise)
  loss_typ_list <- if (use_typical_loss) vector("list", nrow(df_spec)) else NULL
  loss_anchor_list <- if (use_split_loss) vector("list", nrow(df_spec)) else NULL
  loss_pair_list <- if (use_split_loss) vector("list", nrow(df_spec)) else NULL

  # Typical-loss projections for split channels (optional)
  loss_anchor_typ_list <- if (use_typical_loss && use_split_loss) vector("list", nrow(df_spec)) else NULL
  loss_pair_typ_list <- if (use_typical_loss && use_split_loss) vector("list", nrow(df_spec)) else NULL

  mref_vec <- numeric(nrow(df_spec))
  mref_conf_vec <- numeric(nrow(df_spec))
  neutral_loss_vec <- character(nrow(df_spec))
  functional_group_vec <- character(nrow(df_spec))

  derivatization_type_vec <- character(nrow(df_spec))
  deriv_score_tms_vec <- numeric(nrow(df_spec))
  deriv_score_tbdms_vec <- numeric(nrow(df_spec))

  # Derivatization diagnostic-ion intensities (for evaluation/QC; auto mode)
  deriv_i73_vec <- numeric(nrow(df_spec))
  deriv_i147_vec <- numeric(nrow(df_spec))
  deriv_i115_vec <- numeric(nrow(df_spec))
  deriv_i57_vec <- numeric(nrow(df_spec))

  names(frag_list) <- df_spec$id
  names(loss_list) <- df_spec$id
  if (use_typical_loss) names(loss_typ_list) <- df_spec$id
  if (use_split_loss) {
    names(loss_anchor_list) <- df_spec$id
    names(loss_pair_list) <- df_spec$id
    if (use_typical_loss) {
      names(loss_anchor_typ_list) <- df_spec$id
      names(loss_pair_typ_list) <- df_spec$id
    }
  }
  names(mref_vec) <- df_spec$id
  names(mref_conf_vec) <- df_spec$id
  names(neutral_loss_vec) <- df_spec$id
  names(functional_group_vec) <- df_spec$id
  names(derivatization_type_vec) <- df_spec$id
  names(deriv_score_tms_vec) <- df_spec$id
  names(deriv_score_tbdms_vec) <- df_spec$id
  names(deriv_i73_vec) <- df_spec$id
  names(deriv_i147_vec) <- df_spec$id
  names(deriv_i115_vec) <- df_spec$id
  names(deriv_i57_vec) <- df_spec$id

  deriv_mode <- resolve_derivatization_mode(params)

  for (i in seq_len(nrow(df_spec))) {
    spec_res <- process_single_spectrum(df_spec$ei[i], params, deriv_mode = deriv_mode)

    frag_list[[i]] <- spec_res$frag
    loss_list[[i]] <- spec_res$loss
    if (use_split_loss) {
      loss_anchor_list[[i]] <- spec_res$loss_anchor
      loss_pair_list[[i]] <- spec_res$loss_pair
    }
    if (use_typical_loss) {
      loss_typ_list[[i]] <- spec_res$loss_typ
      if (use_split_loss) {
        loss_anchor_typ_list[[i]] <- spec_res$loss_anchor_typ
        loss_pair_typ_list[[i]] <- spec_res$loss_pair_typ
      }
    }
    mref_vec[i] <- spec_res$mref
    mref_conf_vec[i] <- spec_res$mref_confidence
    neutral_loss_vec[i] <- spec_res$neutral_losses
    functional_group_vec[i] <- spec_res$functional_groups

    derivatization_type_vec[i] <- spec_res$deriv_type
    deriv_score_tms_vec[i] <- if (is.null(spec_res$deriv_info$score_tms)) NA_real_ else spec_res$deriv_info$score_tms
    deriv_score_tbdms_vec[i] <- if (is.null(spec_res$deriv_info$score_tbdms)) NA_real_ else spec_res$deriv_info$score_tbdms
    deriv_i73_vec[i] <- spec_res$deriv_i73
    deriv_i147_vec[i] <- spec_res$deriv_i147
    deriv_i115_vec[i] <- spec_res$deriv_i115
    deriv_i57_vec[i] <- spec_res$deriv_i57

    if (isTRUE(progress) && i %% 100 == 0) {
      message("Spectrum processing: ", i, " / ", nrow(df_spec))
    }
  }

  df_spec$estimated_Mref <- mref_vec
  df_spec$mref_confidence <- mref_conf_vec
  df_spec$detected_neutral_losses <- neutral_loss_vec
  df_spec$inferred_functional_groups <- functional_group_vec
  df_spec$derivatization_type <- derivatization_type_vec
  df_spec$deriv_score_tms <- deriv_score_tms_vec
  df_spec$deriv_score_tbdms <- deriv_score_tbdms_vec
  df_spec$deriv_i73 <- deriv_i73_vec
  df_spec$deriv_i147 <- deriv_i147_vec
  df_spec$deriv_i115 <- deriv_i115_vec
  df_spec$deriv_i57 <- deriv_i57_vec

  df_spec$compound_class <- vapply(seq_len(nrow(df_spec)), function(i) {
    ion_class <- detect_compound_class(frag_list[[i]], ppm_tol = params$class_detection_ppm, rules = params$classify_rules)
    func_groups <- functional_group_vec[i]

    if (ion_class == "unclassified" && nchar(func_groups) > 0) {
      paste0("NL:", func_groups)
    } else if (ion_class != "unclassified" && nchar(func_groups) > 0) {
      ion_class
    } else {
      ion_class
    }
  }, character(1))
  df_spec$compound_class_source <- "inferred"

  class_summary <- df_spec %>%
    dplyr::select(id, known, RI, estimated_Mref, mref_confidence, detected_neutral_losses,
                  inferred_functional_groups, derivatization_type, deriv_score_tms, deriv_score_tbdms,
                  deriv_i73, deriv_i147, deriv_i115, deriv_i57, compound_class) %>%
    dplyr::arrange(compound_class, RI)

  class_summary_simple <- df_spec %>%
    dplyr::select(id, known, RI, derivatization_type, compound_class) %>%
    dplyr::arrange(compound_class, RI)

  func_group_counts <- df_spec %>%
    dplyr::filter(inferred_functional_groups != "") %>%
    tidyr::separate_longer_delim(inferred_functional_groups, delim = ";") %>%
    dplyr::count(inferred_functional_groups, sort = TRUE)

  ri <- df_spec$RI
  names(ri) <- df_spec$id

  list(
    df_spec = df_spec,
    frag_list = frag_list,
    loss_list = loss_list,
    loss_typ_list = loss_typ_list,
    loss_anchor_list = loss_anchor_list,
    loss_pair_list = loss_pair_list,
    loss_anchor_typ_list = loss_anchor_typ_list,
    loss_pair_typ_list = loss_pair_typ_list,
    mref_confidence = mref_conf_vec,
    ri = ri,
    class_summary = class_summary,
    class_summary_simple = class_summary_simple,
    func_group_counts = func_group_counts
  )
}

#' Compute m/z Frequency Map from Spectra
#'
#' @param frag_list Named list of fragment spectra matrices.
#' @param bin_width Bin width in Da for m/z grouping.
#' @return Named numeric vector of m/z-bin frequencies in the range 0 to 1.
#' @export
compute_mz_frequency <- function(frag_list, bin_width = 1.0) {
  n_spectra <- length(frag_list)
  if (n_spectra == 0) return(numeric(0))
  if (!is.finite(bin_width) || bin_width <= 0) stop("bin_width must be > 0.")

  all_bins <- vector("list", n_spectra)
  for (i in seq_along(frag_list)) {
    spec <- frag_list[[i]]
    if (!is.null(spec) && nrow(spec) > 0) {
      bins <- floor(spec[, 1] / bin_width) * bin_width + bin_width / 2
      all_bins[[i]] <- unique(round(bins, 6))
    } else {
      all_bins[[i]] <- numeric(0)
    }
  }

  bin_table <- table(unlist(all_bins, use.names = FALSE))
  freq <- as.numeric(bin_table) / n_spectra
  names(freq) <- names(bin_table)
  freq
}

#' Apply Intensity and m/z Frequency Weighting to a Spectrum
#'
#' @param spec Spectrum matrix with columns (mz, intensity).
#' @param mz_freq Named m/z frequency vector from compute_mz_frequency().
#' @param intensity_power Power for intensity rescaling.
#' @param freq_power Power for m/z frequency weighting.
#' @param bin_width Bin width used to create mz_freq.
#' @return Weighted spectrum matrix with columns (mz, intensity).
#' @export
apply_spectral_weights <- function(spec,
                                   mz_freq = NULL,
                                   intensity_power = 0.5,
                                   freq_power = 0.5,
                                   bin_width = 1.0) {
  if (is.null(spec) || nrow(spec) == 0) return(spec)
  if (!is.finite(bin_width) || bin_width <= 0) stop("bin_width must be > 0.")

  mz <- spec[, 1]
  int <- spec[, 2]
  int[!is.finite(int) | int < 0] <- 0

  if (!is.null(intensity_power) && !isTRUE(all.equal(intensity_power, 1))) {
    if (!is.finite(intensity_power) || intensity_power <= 0) {
      stop("intensity_power must be > 0.")
    }
    int <- int^intensity_power
  }

  if (!is.null(mz_freq) && length(mz_freq) > 0) {
    if (!is.finite(freq_power) || freq_power <= 0) stop("freq_power must be > 0.")
    bins <- as.character(round(floor(mz / bin_width) * bin_width + bin_width / 2, 6))
    freq_values <- mz_freq[bins]
    freq_values[is.na(freq_values)] <- 0
    int <- int * (freq_values^freq_power)
  }

  total <- sum(int)
  if (total > 0) {
    int <- int / total
  } else {
    int <- spec[, 2]
    fallback_total <- sum(int)
    if (fallback_total > 0) {
      int <- int / fallback_total
    }
  }

  cbind(mz = mz, intensity = int)
}

#' Apply Spectral Weighting to Fragment and Loss Lists
#'
#' @param frag_list Named list of fragment spectra matrices.
#' @param loss_list Optional named list of loss spectra matrices.
#' @param mz_freq Named m/z frequency vector.
#' @param intensity_power Power for intensity rescaling.
#' @param freq_power Power for m/z frequency weighting.
#' @param bin_width Bin width used for m/z frequency lookup.
#' @param apply_loss_weights If TRUE, apply intensity weighting to loss spectra.
#' @return List with weighted frag_list and optional loss_list.
#' @export
weight_spectra <- function(frag_list,
                           loss_list = NULL,
                           mz_freq = NULL,
                           intensity_power = 0.5,
                           freq_power = 0.5,
                           bin_width = 1.0,
                           apply_loss_weights = TRUE) {
  weighted_frag <- lapply(
    frag_list,
    apply_spectral_weights,
    mz_freq = mz_freq,
    intensity_power = intensity_power,
    freq_power = freq_power,
    bin_width = bin_width
  )
  names(weighted_frag) <- names(frag_list)

  result <- list(frag_list = weighted_frag)

  if (!is.null(loss_list)) {
    if (isTRUE(apply_loss_weights)) {
      weighted_loss <- lapply(
        loss_list,
        apply_spectral_weights,
        mz_freq = NULL,
        intensity_power = intensity_power,
        freq_power = 1,
        bin_width = bin_width
      )
      names(weighted_loss) <- names(loss_list)
      result$loss_list <- weighted_loss
    } else {
      result$loss_list <- loss_list
    }
  }

  result
}


#' Perturb a Spectrum with Artificial Noise
#'
#' Applies controlled perturbations to a spectrum to simulate measurement
#' variability. Three types of noise can be applied independently or combined:
#' (1) m/z shift, (2) intensity noise, (3) peak dropout.
#'
#' @param spec Two-column data.frame or matrix (mz, intensity).
#' @param mz_noise_ppm Numeric. Standard deviation of Gaussian noise added to
#'   m/z values, in ppm. Default: 0 (no shift).
#' @param intensity_noise_sd Numeric. Standard deviation of multiplicative
#'   log-normal noise on intensity. Default: 0 (no noise).
#' @param dropout_prob Numeric in 0-1. Probability of each peak being
#'   randomly removed. Default: 0 (no dropout).
#' @param seed Optional integer for reproducibility.
#' @return Perturbed spectrum in same format as input.
#' @export
perturb_spectrum <- function(spec, mz_noise_ppm = 0, intensity_noise_sd = 0,
                              dropout_prob = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  is_df <- is.data.frame(spec)
  if (is_df) {
    mz <- spec[, 1]
    int <- spec[, 2]
  } else {
    mz <- spec[, 1]
    int <- spec[, 2]
  }

  n <- length(mz)
  if (n == 0) return(spec)

  # 1. Peak dropout
  if (dropout_prob > 0 && dropout_prob < 1) {
    keep <- stats::runif(n) >= dropout_prob
    if (sum(keep) == 0) {
      # All peaks dropped — return empty spectrum
      if (is_df) {
        return(data.frame(mz = numeric(0), intensity = numeric(0)))
      } else {
        return(matrix(numeric(0), ncol = 2,
                      dimnames = list(NULL, colnames(spec))))
      }
    }
    mz <- mz[keep]
    int <- int[keep]
    n <- length(mz)
  }

  # 2. m/z noise (ppm-proportional Gaussian)
  if (mz_noise_ppm > 0) {
    shift <- stats::rnorm(n, mean = 0, sd = mz_noise_ppm * mz * 1e-6)
    mz <- mz + shift
  }

  # 3. Intensity noise (multiplicative log-normal)
  if (intensity_noise_sd > 0) {
    mult <- exp(stats::rnorm(n, mean = 0, sd = intensity_noise_sd))
    int <- int * mult
    int <- pmax(int, 0)
  }

  if (is_df) {
    result <- data.frame(mz = mz, intensity = int)
    colnames(result) <- colnames(spec)
  } else {
    result <- cbind(mz, int)
    colnames(result) <- colnames(spec)
  }

  result
}


#' Perturb a List of Spectra
#'
#' Applies \code{perturb_spectrum()} to each spectrum in a named list.
#'
#' @param spec_list Named list of spectra (each a 2-column data.frame/matrix).
#' @param mz_noise_ppm Same as \code{perturb_spectrum}.
#' @param intensity_noise_sd Same as \code{perturb_spectrum}.
#' @param dropout_prob Same as \code{perturb_spectrum}.
#' @param seed Base seed. Each spectrum uses seed + index for reproducibility.
#'   NULL for no seed.
#' @return Named list of perturbed spectra (same names as input).
#' @export
perturb_spectra_list <- function(spec_list, mz_noise_ppm = 0,
                                  intensity_noise_sd = 0, dropout_prob = 0,
                                  seed = NULL) {
  n <- length(spec_list)
  result <- vector("list", n)
  names(result) <- names(spec_list)

  for (i in seq_len(n)) {
    s <- if (!is.null(seed)) seed + i else NULL
    result[[i]] <- perturb_spectrum(
      spec_list[[i]],
      mz_noise_ppm = mz_noise_ppm,
      intensity_noise_sd = intensity_noise_sd,
      dropout_prob = dropout_prob,
      seed = s
    )
  }

  result
}
