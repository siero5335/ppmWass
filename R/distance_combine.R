#' @keywords internal
#' @noRd
combined_distance <- function(fragA, fragB, lossA, lossB, params,
                              lossA_typ = NULL, lossB_typ = NULL,
                              confA = NULL, confB = NULL) {
  mass_power <- if (is.null(params$mass_power)) 3 else params$mass_power
  intensity_power <- if (is.null(params$intensity_power)) 0.5 else params$intensity_power
  trans_mult <- if (is.null(params$wasserstein_transition_mult)) 3 else params$wasserstein_transition_mult
  ot_method_val <- if (is.null(params$ot_method)) "exact" else params$ot_method
  sink_eps <- if (is.null(params$sinkhorn_epsilon)) 0.05 else params$sinkhorn_epsilon
  sink_niter <- if (is.null(params$sinkhorn_niter)) 100L else params$sinkhorn_niter

  d_frag <- compute_distance(
    fragA, fragB, params$distance_method, params$tol_ppm, params$wasserstein_align,
    mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
    ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
  )

  d_loss_raw <- compute_distance(
    lossA, lossB, params$distance_method, params$tol_ppm, params$wasserstein_align,
    mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
    ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
  )

  d_loss <- d_loss_raw

  if (isTRUE(params$use_typical_loss) &&
      !is.null(lossA_typ) && !is.null(lossB_typ) &&
      nrow(lossA_typ) > 0 && nrow(lossB_typ) > 0) {

    beta_base <- if (is.null(params$loss_raw_beta)) 0.80 else params$loss_raw_beta
    beta_base <- min(max(beta_base, 0), 1)

    d_loss_typ <- compute_distance(
      lossA_typ, lossB_typ, params$distance_method, params$tol_ppm, params$wasserstein_align,
      mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
      ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
    )

    typ_w <- 1 - beta_base

    if (isTRUE(params$use_mref_confidence) && !is.null(confA) && !is.null(confB)) {
      a <- if (is.finite(confA)) confA else 0
      b <- if (is.finite(confB)) confB else 0
      a <- min(max(a, 0), 1)
      b <- min(max(b, 0), 1)
      conf_pair <- sqrt(a * b)
      gamma <- if (is.null(params$loss_typical_conf_gamma)) 1 else params$loss_typical_conf_gamma
      if (!is.finite(gamma) || gamma < 0) gamma <- 1
      typ_w <- typ_w * (conf_pair^gamma)
    }

    raw_w <- 1 - typ_w
    d_loss <- sqrt(raw_w * d_loss_raw^2 + typ_w * d_loss_typ^2)
  }

  sqrt(params$w_frag * d_frag^2 + params$w_loss * d_loss^2)
}

#' Combine fragment and derived-spectrum distances using one raw requested
#' approximate-OT call per nonempty channel
#'
#' This helper is intentionally narrow: it supports the publication Sinkhorn
#' sensitivity diagnostic, where retries and exact fallback would change the
#' estimand. The ordinary exact-primary and approximate fallback paths continue
#' to use `combined_distance()`.
#' @keywords internal
#' @noRd
combined_distance_requested_approx <- function(
    fragA, fragB, lossA, lossB, params, diagnostics = NULL,
    context = NULL) {
  if (!identical(params$distance_method, "ppm_wasserstein")) {
    stop(
      "combined_distance_requested_approx requires distance_method='ppm_wasserstein'.",
      call. = FALSE
    )
  }
  if (!params$ot_method %in% c("sinkhorn", "greenkhorn")) {
    stop(
      "combined_distance_requested_approx requires sinkhorn or greenkhorn.",
      call. = FALSE
    )
  }
  if (isTRUE(params$use_typical_loss) || isTRUE(params$use_split_loss)) {
    stop(
      "Requested-only publication sensitivity does not support typical/split loss channels.",
      call. = FALSE
    )
  }

  context_for <- function(channel) {
    base <- if (is.null(context)) list() else context
    c(base, list(channel = channel))
  }
  raw_distance <- function(a, b, channel) {
    ppm_wasserstein_requested_approx(
      a, b, ppm = params$tol_ppm,
      transition_mult = params$wasserstein_transition_mult,
      align = params$wasserstein_align,
      ot_method = params$ot_method,
      sinkhorn_epsilon = params$sinkhorn_epsilon,
      sinkhorn_niter = params$sinkhorn_niter,
      .diagnostics = diagnostics,
      .context = context_for(channel)
    )
  }

  d_frag <- raw_distance(fragA, fragB, "fragment")
  d_loss <- raw_distance(lossA, lossB, "derived")
  if (!is.finite(d_frag) || !is.finite(d_loss)) return(NA_real_)
  sqrt(params$w_frag * d_frag^2 + params$w_loss * d_loss^2)
}

#' @keywords internal
#' @noRd
combined_distance_split <- function(fragA, fragB,
                                    lossA_anchor, lossB_anchor,
                                    lossA_pair, lossB_pair,
                                    params,
                                    lossA_anchor_typ = NULL, lossB_anchor_typ = NULL,
                                    lossA_pair_typ = NULL, lossB_pair_typ = NULL,
                                    confA = NULL, confB = NULL) {
  mass_power <- if (is.null(params$mass_power)) 3 else params$mass_power
  intensity_power <- if (is.null(params$intensity_power)) 0.5 else params$intensity_power
  trans_mult <- if (is.null(params$wasserstein_transition_mult)) 3 else params$wasserstein_transition_mult
  ot_method_val <- if (is.null(params$ot_method)) "exact" else params$ot_method
  sink_eps <- if (is.null(params$sinkhorn_epsilon)) 0.05 else params$sinkhorn_epsilon
  sink_niter <- if (is.null(params$sinkhorn_niter)) 100L else params$sinkhorn_niter

  d_frag <- compute_distance(
    fragA, fragB, params$distance_method, params$tol_ppm, params$wasserstein_align,
    mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
    ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
  )

  channel_distance <- function(a, b) {
    if (is.null(a) || is.null(b)) return(1)
    if (nrow(a) == 0 && nrow(b) == 0) return(0)
    compute_distance(
      a, b, params$distance_method, params$tol_ppm, params$wasserstein_align,
      mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
      ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
    )
  }

  d_anchor_raw <- channel_distance(lossA_anchor, lossB_anchor)
  d_anchor <- d_anchor_raw

  if (isTRUE(params$use_typical_loss) &&
      !is.null(lossA_anchor_typ) && !is.null(lossB_anchor_typ) &&
      nrow(lossA_anchor_typ) > 0 && nrow(lossB_anchor_typ) > 0) {

    beta_base <- if (is.null(params$loss_raw_beta)) 0.80 else params$loss_raw_beta
    beta_base <- min(max(beta_base, 0), 1)

    d_anchor_typ <- channel_distance(lossA_anchor_typ, lossB_anchor_typ)

    typ_w <- 1 - beta_base

    if (isTRUE(params$use_mref_confidence) && !is.null(confA) && !is.null(confB)) {
      a <- if (is.finite(confA)) confA else 0
      b <- if (is.finite(confB)) confB else 0
      a <- min(max(a, 0), 1)
      b <- min(max(b, 0), 1)
      conf_pair <- sqrt(a * b)
      gamma <- if (is.null(params$loss_typical_conf_gamma)) 1 else params$loss_typical_conf_gamma
      if (!is.finite(gamma) || gamma < 0) gamma <- 1
      typ_w <- typ_w * (conf_pair^gamma)
    }

    raw_w <- 1 - typ_w
    d_anchor <- sqrt(raw_w * d_anchor_raw^2 + typ_w * d_anchor_typ^2)
  }

  d_pair_raw <- channel_distance(lossA_pair, lossB_pair)
  d_pair <- d_pair_raw

  if (isTRUE(params$use_typical_loss) &&
      !is.null(lossA_pair_typ) && !is.null(lossB_pair_typ) &&
      nrow(lossA_pair_typ) > 0 && nrow(lossB_pair_typ) > 0) {

    beta_base <- if (is.null(params$loss_raw_beta)) 0.80 else params$loss_raw_beta
    beta_base <- min(max(beta_base, 0), 1)

    d_pair_typ <- channel_distance(lossA_pair_typ, lossB_pair_typ)
    typ_w <- 1 - beta_base
    raw_w <- 1 - typ_w
    d_pair <- sqrt(raw_w * d_pair_raw^2 + typ_w * d_pair_typ^2)
  }

  w_anchor_base <- if (is.null(params$loss_anchor_weight)) 0.6 else params$loss_anchor_weight
  w_anchor_base <- min(max(w_anchor_base, 0), 1)
  w_pair_base <- 1 - w_anchor_base

  w_anchor <- w_anchor_base
  w_pair <- w_pair_base

  if (isTRUE(params$use_mref_confidence) && !is.null(confA) && !is.null(confB)) {
    a <- if (is.finite(confA)) confA else 0
    b <- if (is.finite(confB)) confB else 0
    a <- min(max(a, 0), 1)
    b <- min(max(b, 0), 1)
    conf_pair <- sqrt(a * b)
    gamma2 <- if (is.null(params$loss_anchor_conf_gamma)) 1 else params$loss_anchor_conf_gamma
    if (!is.finite(gamma2) || gamma2 < 0) gamma2 <- 1
    w_anchor <- w_anchor * (conf_pair^gamma2)
  }

  if (!is.null(lossA_anchor) && !is.null(lossB_anchor) &&
      nrow(lossA_anchor) == 0 && nrow(lossB_anchor) == 0) {
    w_anchor <- 0
  }
  if (!is.null(lossA_pair) && !is.null(lossB_pair) &&
      nrow(lossA_pair) == 0 && nrow(lossB_pair) == 0) {
    w_pair <- 0
  }

  s <- w_anchor + w_pair
  if (!is.finite(s) || s <= 0) {
    d_loss <- 1
  } else {
    w_anchor <- w_anchor / s
    w_pair <- w_pair / s
    d_loss <- sqrt(w_anchor * d_anchor^2 + w_pair * d_pair^2)
  }

  sqrt(params$w_frag * d_frag^2 + params$w_loss * d_loss^2)
}

#' @keywords internal
#' @noRd
combined_distance_details <- function(fragA, fragB, lossA, lossB, params,
                                      lossA_typ = NULL, lossB_typ = NULL,
                                      confA = NULL, confB = NULL) {
  mass_power <- if (is.null(params$mass_power)) 3 else params$mass_power
  intensity_power <- if (is.null(params$intensity_power)) 0.5 else params$intensity_power
  trans_mult <- if (is.null(params$wasserstein_transition_mult)) 3 else params$wasserstein_transition_mult
  ot_method_val <- if (is.null(params$ot_method)) "exact" else params$ot_method
  sink_eps <- if (is.null(params$sinkhorn_epsilon)) 0.05 else params$sinkhorn_epsilon
  sink_niter <- if (is.null(params$sinkhorn_niter)) 100L else params$sinkhorn_niter

  d_frag <- compute_distance(
    fragA, fragB, params$distance_method, params$tol_ppm, params$wasserstein_align,
    mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
    ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
  )

  d_loss_raw <- compute_distance(
    lossA, lossB, params$distance_method, params$tol_ppm, params$wasserstein_align,
    mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
    ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
  )

  d_loss <- d_loss_raw
  d_loss_typ <- NA_real_
  w_loss_typical <- 0
  w_loss_raw <- 1

  if (isTRUE(params$use_typical_loss) &&
      !is.null(lossA_typ) && !is.null(lossB_typ) &&
      nrow(lossA_typ) > 0 && nrow(lossB_typ) > 0) {

    beta_base <- if (is.null(params$loss_raw_beta)) 0.80 else params$loss_raw_beta
    beta_base <- min(max(beta_base, 0), 1)

    d_loss_typ <- compute_distance(
      lossA_typ, lossB_typ, params$distance_method, params$tol_ppm, params$wasserstein_align,
      mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
      ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
    )

    w_loss_typical <- 1 - beta_base

    if (isTRUE(params$use_mref_confidence) && !is.null(confA) && !is.null(confB)) {
      a <- if (is.finite(confA)) confA else 0
      b <- if (is.finite(confB)) confB else 0
      a <- min(max(a, 0), 1)
      b <- min(max(b, 0), 1)
      conf_pair <- sqrt(a * b)
      gamma <- if (is.null(params$loss_typical_conf_gamma)) 1 else params$loss_typical_conf_gamma
      if (!is.finite(gamma) || gamma < 0) gamma <- 1
      w_loss_typical <- w_loss_typical * (conf_pair^gamma)
    }

    w_loss_typical <- min(max(w_loss_typical, 0), 1)
    w_loss_raw <- 1 - w_loss_typical

    d_loss <- sqrt(w_loss_raw * d_loss_raw^2 + w_loss_typical * d_loss_typ^2)
  }

  d_total <- sqrt(params$w_frag * d_frag^2 + params$w_loss * d_loss^2)

  list(
    d_total = d_total,
    d_frag = d_frag,
    d_loss = d_loss,
    d_loss_raw = d_loss_raw,
    d_loss_typ = d_loss_typ,
    w_loss_raw = w_loss_raw,
    w_loss_typical = w_loss_typical
  )
}

#' @keywords internal
#' @noRd
combined_distance_split_details <- function(fragA, fragB,
                                            lossA_anchor, lossB_anchor,
                                            lossA_pair, lossB_pair,
                                            params,
                                            lossA_anchor_typ = NULL, lossB_anchor_typ = NULL,
                                            lossA_pair_typ = NULL, lossB_pair_typ = NULL,
                                            confA = NULL, confB = NULL) {
  mass_power <- if (is.null(params$mass_power)) 3 else params$mass_power
  intensity_power <- if (is.null(params$intensity_power)) 0.5 else params$intensity_power
  trans_mult <- if (is.null(params$wasserstein_transition_mult)) 3 else params$wasserstein_transition_mult
  ot_method_val <- if (is.null(params$ot_method)) "exact" else params$ot_method
  sink_eps <- if (is.null(params$sinkhorn_epsilon)) 0.05 else params$sinkhorn_epsilon
  sink_niter <- if (is.null(params$sinkhorn_niter)) 100L else params$sinkhorn_niter

  d_frag <- compute_distance(
    fragA, fragB, params$distance_method, params$tol_ppm, params$wasserstein_align,
    mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
    ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
  )

  channel_distance <- function(a, b) {
    if (is.null(a) || is.null(b)) return(1)
    if (nrow(a) == 0 && nrow(b) == 0) return(0)
    compute_distance(
      a, b, params$distance_method, params$tol_ppm, params$wasserstein_align,
      mass_power, intensity_power, wasserstein_transition_mult = trans_mult,
      ot_method = ot_method_val, sinkhorn_epsilon = sink_eps, sinkhorn_niter = sink_niter
    )
  }

  conf_pair <- NA_real_
  if (isTRUE(params$use_mref_confidence) && !is.null(confA) && !is.null(confB)) {
    a <- if (is.finite(confA)) confA else 0
    b <- if (is.finite(confB)) confB else 0
    a <- min(max(a, 0), 1)
    b <- min(max(b, 0), 1)
    conf_pair <- sqrt(a * b)
  }

  d_anchor_raw <- channel_distance(lossA_anchor, lossB_anchor)
  d_anchor <- d_anchor_raw
  d_anchor_typ <- NA_real_
  w_anchor_typical <- 0
  w_anchor_raw <- 1

  if (isTRUE(params$use_typical_loss) &&
      !is.null(lossA_anchor_typ) && !is.null(lossB_anchor_typ) &&
      nrow(lossA_anchor_typ) > 0 && nrow(lossB_anchor_typ) > 0) {

    beta_base <- if (is.null(params$loss_raw_beta)) 0.80 else params$loss_raw_beta
    beta_base <- min(max(beta_base, 0), 1)

    d_anchor_typ <- channel_distance(lossA_anchor_typ, lossB_anchor_typ)

    w_anchor_typical <- 1 - beta_base

    if (isTRUE(params$use_mref_confidence) && is.finite(conf_pair)) {
      gamma <- if (is.null(params$loss_typical_conf_gamma)) 1 else params$loss_typical_conf_gamma
      if (!is.finite(gamma) || gamma < 0) gamma <- 1
      w_anchor_typical <- w_anchor_typical * (conf_pair^gamma)
    }

    w_anchor_typical <- min(max(w_anchor_typical, 0), 1)
    w_anchor_raw <- 1 - w_anchor_typical

    d_anchor <- sqrt(w_anchor_raw * d_anchor_raw^2 + w_anchor_typical * d_anchor_typ^2)
  }

  d_pair_raw <- channel_distance(lossA_pair, lossB_pair)
  d_pair <- d_pair_raw
  d_pair_typ <- NA_real_
  w_pair_typical <- 0
  w_pair_raw <- 1

  if (isTRUE(params$use_typical_loss) &&
      !is.null(lossA_pair_typ) && !is.null(lossB_pair_typ) &&
      nrow(lossA_pair_typ) > 0 && nrow(lossB_pair_typ) > 0) {

    beta_base <- if (is.null(params$loss_raw_beta)) 0.80 else params$loss_raw_beta
    beta_base <- min(max(beta_base, 0), 1)

    d_pair_typ <- channel_distance(lossA_pair_typ, lossB_pair_typ)

    w_pair_typical <- 1 - beta_base
    w_pair_typical <- min(max(w_pair_typical, 0), 1)
    w_pair_raw <- 1 - w_pair_typical

    d_pair <- sqrt(w_pair_raw * d_pair_raw^2 + w_pair_typical * d_pair_typ^2)
  }

  w_anchor_base <- if (is.null(params$loss_anchor_weight)) 0.6 else params$loss_anchor_weight
  w_anchor_base <- min(max(w_anchor_base, 0), 1)
  w_pair_base <- 1 - w_anchor_base

  w_anchor_final <- w_anchor_base
  w_pair_final <- w_pair_base

  if (isTRUE(params$use_mref_confidence) && is.finite(conf_pair)) {
    gamma2 <- if (is.null(params$loss_anchor_conf_gamma)) 1 else params$loss_anchor_conf_gamma
    if (!is.finite(gamma2) || gamma2 < 0) gamma2 <- 1
    w_anchor_final <- w_anchor_final * (conf_pair^gamma2)
  }

  if (!is.null(lossA_anchor) && !is.null(lossB_anchor) &&
      nrow(lossA_anchor) == 0 && nrow(lossB_anchor) == 0) {
    w_anchor_final <- 0
  }
  if (!is.null(lossA_pair) && !is.null(lossB_pair) &&
      nrow(lossA_pair) == 0 && nrow(lossB_pair) == 0) {
    w_pair_final <- 0
  }

  s <- w_anchor_final + w_pair_final
  if (!is.finite(s) || s <= 0) {
    d_loss <- 1
    w_anchor_final <- 0
    w_pair_final <- 0
  } else {
    w_anchor_final <- w_anchor_final / s
    w_pair_final <- w_pair_final / s
    d_loss <- sqrt(w_anchor_final * d_anchor^2 + w_pair_final * d_pair^2)
  }

  d_total <- sqrt(params$w_frag * d_frag^2 + params$w_loss * d_loss^2)

  list(
    d_total = d_total,
    d_frag = d_frag,
    d_loss = d_loss,
    d_anchor = d_anchor,
    d_anchor_raw = d_anchor_raw,
    d_anchor_typ = d_anchor_typ,
    w_anchor_raw = w_anchor_raw,
    w_anchor_typical = w_anchor_typical,
    d_pair = d_pair,
    d_pair_raw = d_pair_raw,
    d_pair_typ = d_pair_typ,
    w_pair_raw = w_pair_raw,
    w_pair_typical = w_pair_typical,
    w_anchor_final = w_anchor_final,
    w_pair_final = w_pair_final,
    conf_pair = conf_pair
  )
}
