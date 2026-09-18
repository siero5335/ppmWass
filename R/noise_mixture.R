normalize_tic_spectrum <- function(x) {
  x <- lc_clean_spectrum(x)
  if (!nrow(x)) return(x)
  total <- sum(x[, 2L])
  if (!is.finite(total) || total <= 0) return(lc_empty_spectrum())
  x[, 2L] <- x[, 2L] / total
  x
}

merge_mixture_peaks <- function(x, merge_ppm) {
  x <- lc_clean_spectrum(x)
  if (nrow(x) < 2L || merge_ppm <= 0) return(normalize_tic_spectrum(x))
  out_mz <- numeric(nrow(x))
  out_int <- numeric(nrow(x))
  k <- 0L
  i <- 1L
  while (i <= nrow(x)) {
    anchor <- x[i, 1L]
    tol <- anchor * merge_ppm * 1e-6
    j <- i
    while (j <= nrow(x) && abs(x[j, 1L] - anchor) <= tol) j <- j + 1L
    idx <- i:(j - 1L)
    intensity <- x[idx, 2L]
    k <- k + 1L
    out_mz[[k]] <- sum(x[idx, 1L] * intensity) / sum(intensity)
    out_int[[k]] <- sum(intensity)
    i <- j
  }
  normalize_tic_spectrum(cbind(mz = out_mz[seq_len(k)], intensity = out_int[seq_len(k)]))
}

mix_signal_noise_impl <- function(signal, noise, signal_fraction, dropout_prob,
                                  merge_ppm, seed = NULL) {
  if (length(signal_fraction) != 1L || is.na(signal_fraction) || !is.finite(signal_fraction) ||
      signal_fraction < 0 || signal_fraction > 1) {
    stop("signal_fraction must be between 0 and 1.")
  }
  if (length(dropout_prob) != 1L || is.na(dropout_prob) || !is.finite(dropout_prob) ||
      dropout_prob < 0 || dropout_prob > 1) {
    stop("dropout_prob must be between 0 and 1.")
  }
  if (length(merge_ppm) != 1L || is.na(merge_ppm) || !is.finite(merge_ppm) || merge_ppm < 0) {
    stop("merge_ppm must be a non-negative finite number.")
  }

  run <- function() {
    signal_raw <- normalize_tic_spectrum(signal)
    noise_raw <- normalize_tic_spectrum(noise)
    n_signal_before <- nrow(signal_raw)
    if (nrow(signal_raw) && dropout_prob > 0) {
      signal_raw <- signal_raw[stats::runif(nrow(signal_raw)) >= dropout_prob, , drop = FALSE]
      signal_raw <- normalize_tic_spectrum(signal_raw)
    }
    if (signal_fraction < 1 && !nrow(noise_raw)) {
      stop("noise must contain at least one valid peak when signal_fraction < 1.")
    }
    if (signal_fraction > 0 && !nrow(signal_raw)) {
      signal_weight <- 0
      noise_weight <- 1
    } else if (signal_fraction == 1) {
      signal_weight <- 1
      noise_weight <- 0
    } else if (signal_fraction == 0) {
      signal_weight <- 0
      noise_weight <- 1
    } else {
      signal_weight <- signal_fraction
      noise_weight <- 1 - signal_fraction
    }
    if (nrow(signal_raw)) signal_raw[, 2L] <- signal_raw[, 2L] * signal_weight
    if (nrow(noise_raw)) noise_raw[, 2L] <- noise_raw[, 2L] * noise_weight
    combined <- rbind(signal_raw, noise_raw)
    combined <- combined[combined[, 2L] > 0, , drop = FALSE]
    combined_tic <- sum(combined[, 2L])
    source_labelled_signal_fraction <- if (
      is.finite(combined_tic) && combined_tic > 0
    ) {
      sum(signal_raw[, 2L]) / combined_tic
    } else {
      NA_real_
    }
    spectrum <- merge_mixture_peaks(combined, merge_ppm)
    list(
      spectrum = spectrum,
      diagnostics = c(
        requested_premerge_signal_tic_weight = signal_fraction,
        assigned_premerge_signal_tic_weight = signal_weight,
        # Backward-compatible aliases. These values are assigned before peak
        # coalescing and are not an independently measured post-merge purity.
        requested_signal_fraction = signal_fraction,
        realized_signal_fraction = signal_weight,
        # Source-labelled attribution is computed before source labels are
        # discarded. Coalescing conserves intensity; an overlapping merged
        # peak is allocated in proportion to its signal/noise contributions.
        source_labelled_postmerge_signal_tic_fraction = source_labelled_signal_fraction,
        n_signal_peaks_before = n_signal_before,
        n_signal_peaks_after = nrow(signal_raw),
        n_noise_peaks = nrow(noise_raw),
        n_mixture_peaks = nrow(spectrum)
      )
    )
  }
  with_seed(seed, run())
}

#' Mix a Signal Spectrum with a Noise Spectrum
#'
#' Signal and noise are independently TIC-normalized, weighted, combined, and
#' then normalized again. This preserves the requested signal TIC fraction and
#' avoids the ineffective pattern of scaling a spectrum after normalization.
#'
#' @param signal,noise Two-column spectra matrices.
#' @param signal_fraction Requested signal weight of the pre-merge TIC.
#' @param dropout_prob Independent signal-peak dropout probability.
#' @param merge_ppm Merge tolerance for overlapping signal/noise peaks.
#' @param seed Optional random seed.
#' @return A two-column normalized mixture spectrum.
#' @export
mix_signal_noise_spectrum <- function(signal, noise, signal_fraction = 0.2,
                                      dropout_prob = 0, merge_ppm = 10,
                                      seed = NULL) {
  mix_signal_noise_impl(signal, noise, signal_fraction, dropout_prob,
                        merge_ppm, seed)$spectrum
}

combine_noise_spectra <- function(noise_list, merge_ppm) {
  valid <- lapply(noise_list, normalize_tic_spectrum)
  valid <- valid[vapply(valid, nrow, integer(1L)) > 0L]
  if (!length(valid)) return(lc_empty_spectrum())
  weight <- 1 / length(valid)
  valid <- lapply(valid, function(x) {
    x[, 2L] <- x[, 2L] * weight
    x
  })
  merge_mixture_peaks(do.call(rbind, valid), merge_ppm)
}

#' Create High-noise, Low-signal Query Spectra
#'
#' By default, noise is sampled from other spectra in `signal_list`, providing
#' an interferent/co-isolation stress test. A separate blank or matrix-noise
#' library can be supplied via `noise_pool`.
#'
#' @param signal_list Named list of signal spectra.
#' @param noise_pool Named list of candidate noise spectra. `NULL` reuses
#'   `signal_list`.
#' @param signal_fraction Requested signal weight of the pre-merge TIC.
#' @param dropout_prob Signal-peak dropout probability.
#' @param noise_components Number of independently sampled noise spectra pooled
#'   for each query.
#' @param merge_ppm Peak merge tolerance.
#' @param avoid_self Avoid noise entries whose ID equals the signal ID.
#' @param seed Optional random seed.
#' @return A list with `spectra` and a per-query `diagnostics` data frame.
#' @export
mix_signal_noise_spectra <- function(signal_list, noise_pool = NULL,
                                     signal_fraction = 0.2,
                                     dropout_prob = 0,
                                     noise_components = 1L,
                                     merge_ppm = 10,
                                     avoid_self = TRUE,
                                     seed = NULL) {
  if (!is.list(signal_list) || length(signal_list) < 1L) {
    stop("signal_list must contain at least one spectrum.")
  }
  if (is.null(names(signal_list))) names(signal_list) <- as.character(seq_along(signal_list))
  if (is.null(noise_pool)) noise_pool <- signal_list
  if (!is.list(noise_pool) || length(noise_pool) < 1L) {
    stop("noise_pool must contain at least one spectrum.")
  }
  if (is.null(names(noise_pool))) names(noise_pool) <- as.character(seq_along(noise_pool))
  noise_components <- as.integer(noise_components)
  if (length(noise_components) != 1L || is.na(noise_components) || noise_components < 1L) {
    stop("noise_components must be >= 1.")
  }

  run <- function() {
    out <- vector("list", length(signal_list))
    names(out) <- names(signal_list)
    diagnostics <- vector("list", length(signal_list))
    for (i in seq_along(signal_list)) {
      eligible <- seq_along(noise_pool)
      if (isTRUE(avoid_self)) {
        eligible <- eligible[names(noise_pool)[eligible] != names(signal_list)[i]]
      }
      if (signal_fraction < 1 && !length(eligible)) {
        stop("No eligible noise spectra remain for signal ID: ", names(signal_list)[i])
      }
      if (signal_fraction == 1) {
        chosen <- integer()
        noise <- lc_empty_spectrum()
      } else {
        chosen <- sample(eligible, noise_components, replace = length(eligible) < noise_components)
        noise <- combine_noise_spectra(noise_pool[chosen], merge_ppm)
      }
      mixed <- mix_signal_noise_impl(
        signal = signal_list[[i]],
        noise = noise,
        signal_fraction = signal_fraction,
        dropout_prob = dropout_prob,
        merge_ppm = merge_ppm,
        seed = NULL
      )
      out[[i]] <- mixed$spectrum
      d <- as.list(mixed$diagnostics)
      d$signal_id <- names(signal_list)[i]
      d$noise_ids <- paste(names(noise_pool)[chosen], collapse = ";")
      d$noise_components <- length(chosen)
      diagnostics[[i]] <- as.data.frame(d, stringsAsFactors = FALSE)
    }
    diagnostics <- do.call(rbind, diagnostics)
    diagnostics <- diagnostics[, c(
      "signal_id", "noise_ids", "requested_premerge_signal_tic_weight",
      "assigned_premerge_signal_tic_weight", "requested_signal_fraction",
      "realized_signal_fraction",
      "source_labelled_postmerge_signal_tic_fraction", "noise_components",
      "n_signal_peaks_before", "n_signal_peaks_after", "n_noise_peaks",
      "n_mixture_peaks"
    ), drop = FALSE]
    list(spectra = out, diagnostics = diagnostics)
  }
  with_seed(seed, run())
}
