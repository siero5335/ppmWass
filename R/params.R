#' @keywords internal
#' @noRd
validate_scalar_param <- function(params, name,
                                  lower = NULL, upper = NULL,
                                  lower_inclusive = TRUE, upper_inclusive = TRUE,
                                  integerish = FALSE, allow_null = TRUE) {
  value <- params[[name]]

  if (is.null(value)) {
    if (allow_null) return(invisible(NULL))
    stop(name, " must not be NULL.")
  }

  if (length(value) != 1 || is.na(value) || !is.finite(value)) {
    stop(name, " must be a single finite numeric value.")
  }

  if (integerish && abs(value - round(value)) > sqrt(.Machine$double.eps)) {
    stop(name, " must be an integer-like value.")
  }

  if (!is.null(lower)) {
    ok_lower <- if (lower_inclusive) value >= lower else value > lower
    if (!ok_lower) {
      op <- if (lower_inclusive) ">=" else ">"
      stop(name, " must be ", op, " ", lower, ".")
    }
  }

  if (!is.null(upper)) {
    ok_upper <- if (upper_inclusive) value <= upper else value < upper
    if (!ok_upper) {
      op <- if (upper_inclusive) "<=" else "<"
      stop(name, " must be ", op, " ", upper, ".")
    }
  }

  invisible(NULL)
}

#' @keywords internal
#' @noRd
validate_flag_param <- function(params, name, allow_null = TRUE) {
  value <- params[[name]]
  if (is.null(value)) {
    if (allow_null) return(invisible(NULL))
    stop(name, " must not be NULL.")
  }
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    stop(name, " must be TRUE/FALSE.")
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
validate_choice_param <- function(params, name, choices, allow_null = TRUE) {
  value <- params[[name]]
  if (is.null(value)) {
    if (allow_null) return(invisible(NULL))
    stop(name, " must not be NULL.")
  }
  if (!is.character(value) || length(value) != 1 || is.na(value) || !value %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "))
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
check_common_param_typos <- function(params) {
  typo_map <- c(
    backedn = "backend",
    binppm = "bin_ppm",
    smearppm = "smear_ppm",
    smear_knd = "smear_kind",
    max_dense_cell = "max_dense_cells",
    prefilter_topk = "prefilter_top_k"
  )
  hit <- intersect(names(typo_map), names(params))
  if (length(hit) > 0L) {
    stop("Unknown parameter '", hit[1], "'; did you mean '", typo_map[[hit[1]]], "'?")
  }
  invisible(NULL)
}

#' Default Parameters for EI-HRMS Similarity
#'
#' @return A list of default parameters.
#' @export
#'
#' @examples
#' params <- eihrms_default_params()
#' params$distance_method
#'
#' @details
#' The returned list can be modified and passed to other functions in this package.
#' Optional values:
#' - normalize_ref: NULL (no normalization) or an Alignment ID / Metabolite name / row index
#' - normalize_cols: NULL (all sample columns) or a character vector of column names
#'
#' Optimal-transport backend (applies when \code{distance_method = "ppm_wasserstein"}):
#' - ot_method: \code{"exact"} (network simplex via
#'   \code{transport::transport()}), \code{"sinkhorn"} (finite-iteration
#'   approximation, requires \pkg{approxOT}), or \code{"greenkhorn"}
#'   (approximate variant, requires \pkg{approxOT}), or \code{"exact_sparse"}
#'   (unregularized OT using sparse connected components). Default: \code{"exact"}.
#'   The exact backend is the publication-grade default because it directly
#'   evaluates the unregularized transport objective and is invariant to
#'   exchanging the two spectra. Finite-iteration approximate backends are
#'   retained for explicitly requested sensitivity and performance work.
#' - sinkhorn_epsilon: regularisation strength for Sinkhorn/Greenkhorn (default 0.05).
#'   Smaller values improve accuracy but slow convergence.
#' - sinkhorn_niter: maximum iterations for Sinkhorn/Greenkhorn (default 100).
#'
#' The legacy parameter \code{tol_ppm} has method-specific semantics. For
#' hard-alignment comparators it is a peak-matching tolerance. For
#' \code{ppm_wasserstein} it is the base scale of the ground cost, not a
#' zero-cost window: the cost saturates at
#' \code{tol_ppm * wasserstein_transition_mult}. The names \code{w_loss} and
#' \code{use_split_loss} are also retained for compatibility; they refer to the
#' pooled or split derived mass-difference representation.
eihrms_default_params <- function() {
  list(
    # Internal standard normalization
    normalize_ref = NULL,
    normalize_cols = NULL,

    # EI peak preprocessing
    min_mz = 35,
    max_mz = 650,
    centroid_ppm = 10,
    noise_thr = 0.01,
    topK = 200,
    use_sqrt = FALSE,

    # Neutral loss generation
    loss_top_peaks = 40,
    loss_max_peaks = 250,
    loss_min = 5,
    loss_max = 350,

    # Typical loss projection (optional, for substructure-oriented similarity)
    use_typical_loss = FALSE,
    # Blend within loss channel: beta=1 -> use raw loss only; beta=0 -> use typical-loss only
    loss_raw_beta = 0.80,
    # ppm tolerance for matching observed losses to TYPICAL_LOSSES (defaults to tol_ppm if NULL)
    loss_typical_ppm = NULL,

    # Split neutral-loss channels (optional): anchored (A) vs pairwise (B)
    use_split_loss = FALSE,
    # Within the loss channel, weight for anchored losses (0..1); pairwise weight = 1 - this.
    loss_anchor_weight = 0.6,
    # If Mref-confidence is enabled, modulate anchored-loss weight by pair-confidence^gamma
    loss_anchor_conf_gamma = 1.0,

    # Mref confidence (optional): down-weight anchored neutral losses and typical-loss signal when Mref is unreliable
    use_mref_confidence = FALSE,
    # Reference number of detected typical-loss-supported Mref candidates corresponding to high confidence
    mref_conf_nloss_ref = 3,
    # Reference summed intensity (in normalized spectrum units) of peaks supporting detected_losses
    mref_conf_int_ref = 0.05,
    # High-m/z window (Da) near max(m/z) for checking whether a strong high-mass ion exists
    mref_conf_hi_window = 50,
    # Reference high-m/z max intensity corresponding to high confidence
    mref_conf_hi_int_ref = 0.03,
    # Component weights (will be normalized internally)
    mref_conf_w_nloss = 0.4,
    mref_conf_w_intensity = 0.4,
    mref_conf_w_hi = 0.2,
    # Power for scaling anchored-loss intensities by confidence (conf^power)
    mref_conf_power = 1.0,
    # If typical-loss is enabled, scale its contribution by pair-confidence^gamma
    loss_typical_conf_gamma = 1.0,

    # Typical-loss library options (for Mref estimation and optional substructure-oriented features)
    use_extended_losses = FALSE,
    # Derivatization-related losses (e.g., TMS, TBDMS) are OFF by default; enable explicitly when appropriate.
    use_derivatization_losses = FALSE,
    # "TMS", "TBDMS", or "BOTH"
    derivatization_type = "TMS",

    # Derivatization handling (optional).
    # derivatization_mode: NULL (legacy behavior using use_derivatization_losses),
    #   or one of "off", "manual", "auto".
    derivatization_mode = NULL,
    # Auto-detection settings (used when derivatization_mode == "auto")
    derivatization_auto_ppm = NULL,
    derivatization_auto_min_int = 0.02,
    derivatization_auto_strong_int = 0.08,
    # Downweight derivatization-related signals to avoid ubiquitous TMS/TBDMS patterns dominating similarity
    derivatization_typical_weight = 0.5,
    derivatization_raw_loss_weight = 0.5,
    derivatization_frag_weight = 1.0,

    # Distance method
    # "hellinger", "wasserstein", "ppm_wasserstein", "cosine", "entropy" (= weighted),
    # "entropy_weighted" (alias of "entropy"),
    # "entropy_unweighted", "weighted_cosine", "composite"
    distance_method = "wasserstein",
    wasserstein_align = FALSE,
    wasserstein_transition_mult = 3,
    # Optimal transport solver for ppm_wasserstein:
    #   "sinkhorn"   - finite-iteration approximate (requires approxOT)
    #   "greenkhorn" - approximate variant (requires approxOT)
    #   "exact"      - precise network simplex (uses transport::transport())
    # Approximate solvers use approxOT's C++ backend; exact always uses
    # the transport package for numerical stability and symmetry.
    ot_method = "exact",
    sinkhorn_epsilon = 0.05,
    sinkhorn_niter = 100L,
    # Method-specific ppm value. For ppm_wasserstein this is a base ground-cost
    # scale; for alignment-based comparators it is a hard matching tolerance.
    tol_ppm = 20,
    mass_power = 3,
    intensity_power = 0.5,

    # Blend weights. w_loss is a legacy name for the pooled derived
    # mass-difference branch (anchored plus pairwise differences).
    w_frag = 0.70,
    w_loss = 0.30,

    # RI constraints
    ri_mode_cluster = "none",
    ri_mode_analog = "none",
    ri_sigma = 30,
    ri_window_hard = 100,

    # Analog search
    analog_k = 30,
    analog_sim_threshold = 0.5,
    analog_ri_gap = 100,

    # Parallelization (compute_distance_matrix() uses mclapply on Unix and PSOCK clusters on Windows)
    use_parallel = FALSE,
    n_cores = 4,

    # Sparse-bin backend (experimental, opt-in)
    backend = "pair_loop",
    bin_ppm = 2,
    smear_ppm = NULL,
    smear_kind = "tri",
    max_dense_cells = 1e7,
    prefilter_method = NULL,
    prefilter_top_k = NULL,

    # HRMS class detection
    class_detection_ppm = 15,
    classify_rules = NULL,



    # Distance component diagnostics
    # If TRUE, compute_similarity_matrices() will also return per-component distance matrices.
    return_distance_components = FALSE
  )
}

#' @keywords internal
#' @noRd
validate_params <- function(params) {
  check_common_param_typos(params)

  valid_methods <- c(
    "hellinger", "wasserstein", "ppm_wasserstein", "cosine", "entropy",
    "entropy_weighted",
    "entropy_unweighted", "weighted_cosine", "composite"
  )
  if (!params$distance_method %in% valid_methods) {
    stop("distance_method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  valid_ri_modes <- c("none", "soft", "hard")
  if (!params$ri_mode_cluster %in% valid_ri_modes) {
    stop("ri_mode_cluster must be one of: ", paste(valid_ri_modes, collapse = ", "))
  }
  if (!params$ri_mode_analog %in% valid_ri_modes) {
    stop("ri_mode_analog must be one of: ", paste(valid_ri_modes, collapse = ", "))
  }

  validate_scalar_param(params, "min_mz", lower = 0, allow_null = FALSE)
  validate_scalar_param(params, "max_mz", lower = 0, allow_null = FALSE)
  if (params$max_mz <= params$min_mz) {
    stop("max_mz must be greater than min_mz.")
  }
  validate_scalar_param(params, "centroid_ppm", lower = 0, lower_inclusive = FALSE, allow_null = FALSE)
  validate_scalar_param(params, "noise_thr", lower = 0, upper = 1, allow_null = FALSE)
  validate_scalar_param(params, "topK", lower = 1, integerish = TRUE, allow_null = FALSE)
  validate_flag_param(params, "use_sqrt", allow_null = FALSE)
  validate_scalar_param(params, "loss_top_peaks", lower = 1, integerish = TRUE, allow_null = FALSE)
  validate_scalar_param(params, "loss_max_peaks", lower = 1, integerish = TRUE, allow_null = FALSE)
  validate_scalar_param(params, "loss_min", lower = 0, allow_null = FALSE)
  validate_scalar_param(params, "loss_max", lower = 0, allow_null = FALSE)
  if (params$loss_max <= params$loss_min) {
    stop("loss_max must be greater than loss_min.")
  }
  if (params$loss_max_peaks < params$loss_top_peaks) {
    warning("loss_max_peaks is smaller than loss_top_peaks; pairwise-loss generation may truncate early.")
  }

  validate_scalar_param(params, "tol_ppm", lower = 0, lower_inclusive = FALSE, allow_null = FALSE)
  validate_scalar_param(params, "wasserstein_transition_mult", lower = 0, lower_inclusive = FALSE, allow_null = FALSE)

  # Optimal transport backend params
  if (!is.null(params$ot_method)) {
    valid_ot <- c("sinkhorn", "exact", "greenkhorn", "exact_sparse")
    if (!params$ot_method %in% valid_ot) {
      stop("ot_method must be one of: ", paste(valid_ot, collapse = ", "))
    }
  }
  validate_scalar_param(params, "sinkhorn_epsilon", lower = 0, lower_inclusive = FALSE, allow_null = TRUE)
  validate_scalar_param(params, "sinkhorn_niter", lower = 1, integerish = TRUE, allow_null = TRUE)

  validate_scalar_param(params, "mass_power", lower = 0, allow_null = FALSE)
  validate_scalar_param(params, "intensity_power", lower = 0, allow_null = FALSE)
  validate_scalar_param(params, "w_frag", lower = 0, allow_null = FALSE)
  validate_scalar_param(params, "w_loss", lower = 0, allow_null = FALSE)
  if ((params$w_frag + params$w_loss) <= 0) {
    stop("w_frag + w_loss must be greater than 0.")
  }
  if (abs(params$w_frag + params$w_loss - 1) > 1e-6) {
    warning("w_frag + w_loss is not 1.0; combined distances will be rescaled by the supplied weights.")
  }

  validate_scalar_param(params, "ri_sigma", lower = 0, lower_inclusive = FALSE, allow_null = FALSE)
  validate_scalar_param(params, "ri_window_hard", lower = 0, allow_null = FALSE)
  validate_scalar_param(params, "analog_k", lower = 1, integerish = TRUE, allow_null = FALSE)
  validate_scalar_param(params, "analog_sim_threshold", lower = 0, upper = 1, allow_null = FALSE)
  validate_scalar_param(params, "analog_ri_gap", lower = 0, allow_null = FALSE)
  validate_flag_param(params, "use_parallel", allow_null = FALSE)
  validate_scalar_param(params, "n_cores", lower = 1, integerish = TRUE, allow_null = FALSE)
  validate_scalar_param(params, "class_detection_ppm", lower = 0, lower_inclusive = FALSE, allow_null = FALSE)

  # Sparse-bin backend (experimental, opt-in)
  validate_choice_param(params, "backend", choices = c("pair_loop", "sparse_bins"), allow_null = TRUE)
  validate_scalar_param(params, "bin_ppm", lower = 0, lower_inclusive = FALSE, allow_null = TRUE)
  validate_scalar_param(params, "smear_ppm", lower = 0, lower_inclusive = FALSE, allow_null = TRUE)
  validate_choice_param(params, "smear_kind", choices = c("tri", "rect"), allow_null = TRUE)
  validate_scalar_param(params, "max_dense_cells", lower = 1, integerish = TRUE, allow_null = TRUE)
  validate_choice_param(params, "prefilter_method", choices = valid_methods, allow_null = TRUE)
  validate_scalar_param(params, "prefilter_top_k", lower = 1, integerish = TRUE, allow_null = TRUE)

  if (!is.null(params$normalize_cols)) {
    if (!is.character(params$normalize_cols) || anyNA(params$normalize_cols) || any(params$normalize_cols == "")) {
      stop("normalize_cols must be a character vector of non-empty column names or NULL.")
    }
  }

  # Typical loss projection
  validate_flag_param(params, "use_typical_loss", allow_null = TRUE)

  # Split neutral-loss channels
  validate_flag_param(params, "use_split_loss", allow_null = TRUE)
  if (!is.null(params$loss_anchor_weight)) {
    validate_scalar_param(params, "loss_anchor_weight", lower = 0, upper = 1, allow_null = TRUE)
  }
  if (!is.null(params$loss_anchor_conf_gamma)) {
    validate_scalar_param(params, "loss_anchor_conf_gamma", lower = 0, allow_null = TRUE)
  }

  # Mref confidence
  validate_flag_param(params, "use_mref_confidence", allow_null = TRUE)
  for (nm in c("mref_conf_nloss_ref", "mref_conf_int_ref", "mref_conf_hi_window", "mref_conf_hi_int_ref",
               "mref_conf_w_nloss", "mref_conf_w_intensity", "mref_conf_w_hi", "mref_conf_power",
               "loss_typical_conf_gamma")) {
    if (!is.null(params[[nm]])) {
      validate_scalar_param(params, nm, lower = 0, allow_null = TRUE)
    }
  }
  mref_weight_sum <- sum(c(params$mref_conf_w_nloss, params$mref_conf_w_intensity, params$mref_conf_w_hi), na.rm = TRUE)
  if (is.finite(mref_weight_sum) && mref_weight_sum <= 0) {
    stop("mref_conf_w_nloss + mref_conf_w_intensity + mref_conf_w_hi must be greater than 0.")
  }
  if (is.finite(mref_weight_sum) && abs(mref_weight_sum - 1) > 1e-6) {
    warning("Mref confidence component weights do not sum to 1; they will be normalized internally.")
  }
  if (!is.null(params$loss_typical_conf_gamma) && params$loss_typical_conf_gamma == 0) {
    warning("loss_typical_conf_gamma is 0; typical-loss will not contribute even at high confidence.")
  }
  if (!is.null(params$loss_raw_beta)) {
    validate_scalar_param(params, "loss_raw_beta", lower = 0, upper = 1, allow_null = TRUE)
  }
  if (!is.null(params$loss_typical_ppm)) {
    validate_scalar_param(params, "loss_typical_ppm", lower = 0, lower_inclusive = FALSE, allow_null = TRUE)
  }


  # Typical-loss library flags
  validate_flag_param(params, "use_extended_losses", allow_null = TRUE)
  validate_flag_param(params, "use_derivatization_losses", allow_null = TRUE)
  if (!is.null(params$derivatization_type)) {
    deriv <- toupper(as.character(params$derivatization_type))
    valid_deriv <- c("TMS", "TBDMS", "BOTH", "TMS+TBDMS", "TBDMS+TMS")
    if (!deriv %in% valid_deriv) {
      stop("derivatization_type must be one of: ", paste(valid_deriv, collapse = ", "))
    }
  }

  # Derivatization mode and weights
  if (!is.null(params$derivatization_mode)) {
    mode <- tolower(as.character(params$derivatization_mode))
    valid_modes <- c("off", "manual", "auto")
    if (!mode %in% valid_modes) {
      stop("derivatization_mode must be one of: ", paste(valid_modes, collapse = ", "), " (or NULL for legacy behavior).")
    }
  }
  if (!is.null(params$derivatization_auto_ppm)) {
    validate_scalar_param(params, "derivatization_auto_ppm", lower = 0, lower_inclusive = FALSE, allow_null = TRUE)
  }
  for (nm in c("derivatization_auto_min_int", "derivatization_auto_strong_int",
               "derivatization_typical_weight", "derivatization_raw_loss_weight", "derivatization_frag_weight")) {
    if (!is.null(params[[nm]])) {
      validate_scalar_param(params, nm, lower = 0, upper = 1, allow_null = TRUE)
    }
  }
  if (!is.null(params$derivatization_auto_strong_int) &&
      !is.null(params$derivatization_auto_min_int) &&
      params$derivatization_auto_strong_int < params$derivatization_auto_min_int) {
    warning("derivatization_auto_strong_int is smaller than derivatization_auto_min_int; auto-detection may become overly permissive.")
  }

  # Compound-class detection rules
  if (!is.null(params$classify_rules) && !is.list(params$classify_rules)) {
    stop("classify_rules must be a list when provided (or NULL).")
  }



  # Distance component diagnostics
  validate_flag_param(params, "return_distance_components", allow_null = TRUE)

  params

}
