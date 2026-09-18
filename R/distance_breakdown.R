#' Typical-loss matches between two spectra
#'
#' Creates a ranked table of typical neutral losses that are present in **both** spectra.
#' This is intended for interpretability / substructure hints (e.g., "H2O loss overlaps strongly").
#'
#' @param typA,typB Typical-loss spectra matrices (mz,intensity) with rownames = formula
#'   as produced by project_to_typical_losses().
#' @param top_n Maximum number of matches to return.
#' @param min_overlap Minimum overlap (min(intensity_A, intensity_B)) required.
#' @param sort_by Sort key: "overlap" (default) or "score" (geometric mean).
#' @param typical_losses Optional typical-loss definition list used for class and exact-mass annotations.
#' @return A data.frame with columns: formula, class, exact, intensity_A, intensity_B,
#'   overlap, score, delta.
#' @export
typical_loss_match_table <- function(typA, typB, top_n = 10, min_overlap = 0,
                                     sort_by = c("overlap", "score"), typical_losses = NULL) {
  sort_by <- match.arg(sort_by)
  if (is.null(typical_losses)) typical_losses <- TYPICAL_LOSSES
  if (is.null(typA) || is.null(typB) || nrow(typA) == 0 || nrow(typB) == 0) {
    return(data.frame(
      formula = character(0),
      class = character(0),
      exact = numeric(0),
      intensity_A = numeric(0),
      intensity_B = numeric(0),
      overlap = numeric(0),
      score = numeric(0),
      delta = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  fA <- rownames(typA)
  fB <- rownames(typB)
  if (is.null(fA) || length(fA) != nrow(typA)) fA <- rep("", nrow(typA))
  if (is.null(fB) || length(fB) != nrow(typB)) fB <- rep("", nrow(typB))

  a <- as.numeric(typA[, 2])
  b <- as.numeric(typB[, 2])
  a[!is.finite(a) | a < 0] <- 0
  b[!is.finite(b) | b < 0] <- 0

  vecA <- stats::setNames(a, fA)
  vecB <- stats::setNames(b, fB)

  formulas <- sort(unique(c(names(vecA), names(vecB))))
  if (length(formulas) == 0) {
    return(data.frame(
      formula = character(0),
      class = character(0),
      exact = numeric(0),
      intensity_A = numeric(0),
      intensity_B = numeric(0),
      overlap = numeric(0),
      score = numeric(0),
      delta = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  ia <- vecA[formulas]
  ib <- vecB[formulas]
  ia[is.na(ia)] <- 0
  ib[is.na(ib)] <- 0

  overlap <- pmin(ia, ib)
  score <- sqrt(ia * ib)
  delta <- abs(ia - ib)

  keep <- is.finite(overlap) & overlap > min_overlap
  formulas <- formulas[keep]
  ia <- ia[keep]
  ib <- ib[keep]
  overlap <- overlap[keep]
  score <- score[keep]
  delta <- delta[keep]

  if (length(formulas) == 0) {
    return(data.frame(
      formula = character(0),
      class = character(0),
      exact = numeric(0),
      intensity_A = numeric(0),
      intensity_B = numeric(0),
      overlap = numeric(0),
      score = numeric(0),
      delta = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  class_map <- stats::setNames(vapply(typical_losses, function(z) z$class, character(1)),
                               vapply(typical_losses, function(z) z$formula, character(1)))
  exact_map <- stats::setNames(vapply(typical_losses, function(z) z$exact, numeric(1)),
                               vapply(typical_losses, function(z) z$formula, character(1)))

  cls <- unname(class_map[formulas])
  ex <- unname(exact_map[formulas])
  cls[is.na(cls)] <- ""
  ex[is.na(ex)] <- NA_real_

  df <- data.frame(
    formula = formulas,
    class = cls,
    exact = ex,
    intensity_A = as.numeric(ia),
    intensity_B = as.numeric(ib),
    overlap = as.numeric(overlap),
    score = as.numeric(score),
    delta = as.numeric(delta),
    stringsAsFactors = FALSE
  )

  ord <- if (sort_by == "score") order(df$score, decreasing = TRUE) else order(df$overlap, decreasing = TRUE)
  df <- df[ord, , drop = FALSE]
  if (nrow(df) > top_n) df <- df[seq_len(top_n), , drop = FALSE]

  df
}

#' Distance breakdown for a specific pair (helper)
#'
#' @param idA,idB Compound IDs.
#' @param spectra List returned by build_spectra().
#' @param params Parameter list.
#' @param include_typical If TRUE, also include typical-loss tables (top_n per side).
#' @param top_n Maximum rows per typical-loss table.
#' @return A named list with distances/weights (and optionally typical-loss tables).
#' @export
distance_breakdown_pair <- function(idA, idB, spectra, params = eihrms_default_params(),
                                    include_typical = TRUE, top_n = 10) {
  params <- validate_params(params)
  if (is.null(spectra$frag_list) || is.null(spectra$loss_list)) {
    stop("spectra must be the result of build_spectra() (needs frag_list and loss_list).")
  }
  if (!idA %in% names(spectra$frag_list)) stop("idA not found in spectra$frag_list: ", idA)
  if (!idB %in% names(spectra$frag_list)) stop("idB not found in spectra$frag_list: ", idB)

  confA <- confB <- NULL
  if (isTRUE(params$use_mref_confidence) && !is.null(spectra$mref_confidence)) {
    confA <- unname(spectra$mref_confidence[idA])
    confB <- unname(spectra$mref_confidence[idB])
  }

  typ_univ <- get_typical_losses_universe(params)
  typ_class_map <- stats::setNames(vapply(typ_univ, function(z) z$class, character(1)),
                                   vapply(typ_univ, function(z) z$formula, character(1)))

  typical_table <- function(x) {
    if (is.null(x) || nrow(x) == 0) {
      return(data.frame(formula = character(0), class = character(0),
                        exact = numeric(0), intensity = numeric(0),
                        stringsAsFactors = FALSE))
    }
    formula <- rownames(x)
    if (is.null(formula) || any(!nzchar(formula))) {
      formula <- rep("", nrow(x))
    }
    cls <- unname(typ_class_map[formula])
    cls[is.na(cls)] <- ""
    df <- data.frame(
      formula = formula,
      class = cls,
      exact = as.numeric(x[, 1]),
      intensity = as.numeric(x[, 2]),
      stringsAsFactors = FALSE
    )
    df <- df[order(df$intensity, decreasing = TRUE), , drop = FALSE]
    if (nrow(df) > top_n) df <- df[seq_len(top_n), , drop = FALSE]
    df
  }

  add_contributions <- function(det, split = FALSE) {
    D2 <- params$w_frag * (det$d_frag^2) + params$w_loss * (det$d_loss^2)
    if (!is.finite(D2) || D2 <= 0) {
      det$c_frag_total <- 0
      det$c_loss_total <- 0
    } else {
      det$c_frag_total <- params$w_frag * (det$d_frag^2) / D2
      det$c_loss_total <- params$w_loss * (det$d_loss^2) / D2
    }

    if (!isTRUE(split)) {
      w_raw <- if (!is.null(det$w_loss_raw) && is.finite(det$w_loss_raw)) det$w_loss_raw else 1
      w_typ <- if (!is.null(det$w_loss_typical) && is.finite(det$w_loss_typical)) det$w_loss_typical else 0

      d_raw2 <- if (!is.null(det$d_loss_raw) && is.finite(det$d_loss_raw)) det$d_loss_raw^2 else det$d_loss^2
      d_typ2 <- if (!is.null(det$d_loss_typ) && is.finite(det$d_loss_typ)) det$d_loss_typ^2 else 0

      loss2 <- det$d_loss^2
      det$c_loss_raw_in_loss <- if (is.finite(loss2) && loss2 > 0) w_raw * d_raw2 / loss2 else 0
      det$c_loss_typ_in_loss <- if (is.finite(loss2) && loss2 > 0) w_typ * d_typ2 / loss2 else 0

      det$c_loss_raw_total <- if (is.finite(D2) && D2 > 0) params$w_loss * w_raw * d_raw2 / D2 else 0
      det$c_loss_typ_total <- if (is.finite(D2) && D2 > 0) params$w_loss * w_typ * d_typ2 / D2 else 0

      return(det)
    }

    wA <- if (!is.null(det$w_anchor_final) && is.finite(det$w_anchor_final)) det$w_anchor_final else 0
    wB <- if (!is.null(det$w_pair_final) && is.finite(det$w_pair_final)) det$w_pair_final else 0

    dA2 <- if (!is.null(det$d_anchor) && is.finite(det$d_anchor)) det$d_anchor^2 else 0
    dB2 <- if (!is.null(det$d_pair) && is.finite(det$d_pair)) det$d_pair^2 else 0

    loss2 <- det$d_loss^2
    det$c_anchor_in_loss <- if (is.finite(loss2) && loss2 > 0) wA * dA2 / loss2 else 0
    det$c_pair_in_loss <- if (is.finite(loss2) && loss2 > 0) wB * dB2 / loss2 else 0

    det$c_anchor_total <- if (is.finite(D2) && D2 > 0) params$w_loss * wA * dA2 / D2 else 0
    det$c_pair_total <- if (is.finite(D2) && D2 > 0) params$w_loss * wB * dB2 / D2 else 0

    wA_raw <- if (!is.null(det$w_anchor_raw) && is.finite(det$w_anchor_raw)) det$w_anchor_raw else 1
    wA_typ <- if (!is.null(det$w_anchor_typical) && is.finite(det$w_anchor_typical)) det$w_anchor_typical else 0

    dA_raw2 <- if (!is.null(det$d_anchor_raw) && is.finite(det$d_anchor_raw)) det$d_anchor_raw^2 else dA2
    dA_typ2 <- if (!is.null(det$d_anchor_typ) && is.finite(det$d_anchor_typ)) det$d_anchor_typ^2 else 0

    anchor2 <- if (!is.null(det$d_anchor) && is.finite(det$d_anchor)) det$d_anchor^2 else 0
    det$c_anchor_raw_in_anchor <- if (is.finite(anchor2) && anchor2 > 0) wA_raw * dA_raw2 / anchor2 else 0
    det$c_anchor_typ_in_anchor <- if (is.finite(anchor2) && anchor2 > 0) wA_typ * dA_typ2 / anchor2 else 0

    det$c_anchor_raw_total <- if (is.finite(D2) && D2 > 0) params$w_loss * wA * wA_raw * dA_raw2 / D2 else 0
    det$c_anchor_typ_total <- if (is.finite(D2) && D2 > 0) params$w_loss * wA * wA_typ * dA_typ2 / D2 else 0

    wB_raw <- if (!is.null(det$w_pair_raw) && is.finite(det$w_pair_raw)) det$w_pair_raw else 1
    wB_typ <- if (!is.null(det$w_pair_typical) && is.finite(det$w_pair_typical)) det$w_pair_typical else 0

    dB_raw2 <- if (!is.null(det$d_pair_raw) && is.finite(det$d_pair_raw)) det$d_pair_raw^2 else dB2
    dB_typ2 <- if (!is.null(det$d_pair_typ) && is.finite(det$d_pair_typ)) det$d_pair_typ^2 else 0

    pair2 <- if (!is.null(det$d_pair) && is.finite(det$d_pair)) det$d_pair^2 else 0
    det$c_pair_raw_in_pair <- if (is.finite(pair2) && pair2 > 0) wB_raw * dB_raw2 / pair2 else 0
    det$c_pair_typ_in_pair <- if (is.finite(pair2) && pair2 > 0) wB_typ * dB_typ2 / pair2 else 0

    det$c_pair_raw_total <- if (is.finite(D2) && D2 > 0) params$w_loss * wB * wB_raw * dB_raw2 / D2 else 0
    det$c_pair_typ_total <- if (is.finite(D2) && D2 > 0) params$w_loss * wB * wB_typ * dB_typ2 / D2 else 0

    det
  }

  if (isTRUE(params$use_split_loss) &&
      !is.null(spectra$loss_anchor_list) && !is.null(spectra$loss_pair_list)) {

    det <- combined_distance_split_details(
      spectra$frag_list[[idA]], spectra$frag_list[[idB]],
      spectra$loss_anchor_list[[idA]], spectra$loss_anchor_list[[idB]],
      spectra$loss_pair_list[[idA]], spectra$loss_pair_list[[idB]],
      params,
      lossA_anchor_typ = if (!is.null(spectra$loss_anchor_typ_list)) spectra$loss_anchor_typ_list[[idA]] else NULL,
      lossB_anchor_typ = if (!is.null(spectra$loss_anchor_typ_list)) spectra$loss_anchor_typ_list[[idB]] else NULL,
      lossA_pair_typ = if (!is.null(spectra$loss_pair_typ_list)) spectra$loss_pair_typ_list[[idA]] else NULL,
      lossB_pair_typ = if (!is.null(spectra$loss_pair_typ_list)) spectra$loss_pair_typ_list[[idB]] else NULL,
      confA = confA, confB = confB
    )

    det$idA <- idA
    det$idB <- idB
    det$confA <- confA
    det$confB <- confB

    if (isTRUE(include_typical) && isTRUE(params$use_typical_loss)) {
      det$typical_anchor_A <- typical_table(if (!is.null(spectra$loss_anchor_typ_list)) spectra$loss_anchor_typ_list[[idA]] else NULL)
      det$typical_anchor_B <- typical_table(if (!is.null(spectra$loss_anchor_typ_list)) spectra$loss_anchor_typ_list[[idB]] else NULL)
      det$typical_pair_A <- typical_table(if (!is.null(spectra$loss_pair_typ_list)) spectra$loss_pair_typ_list[[idA]] else NULL)
      det$typical_pair_B <- typical_table(if (!is.null(spectra$loss_pair_typ_list)) spectra$loss_pair_typ_list[[idB]] else NULL)

      det$typical_anchor_matches <- typical_loss_match_table(
        if (!is.null(spectra$loss_anchor_typ_list)) spectra$loss_anchor_typ_list[[idA]] else NULL,
        if (!is.null(spectra$loss_anchor_typ_list)) spectra$loss_anchor_typ_list[[idB]] else NULL,
        top_n = top_n,
        typical_losses = typ_univ
      )
      det$typical_pair_matches <- typical_loss_match_table(
        if (!is.null(spectra$loss_pair_typ_list)) spectra$loss_pair_typ_list[[idA]] else NULL,
        if (!is.null(spectra$loss_pair_typ_list)) spectra$loss_pair_typ_list[[idB]] else NULL,
        top_n = top_n,
        typical_losses = typ_univ
      )
    }

    det <- add_contributions(det, split = TRUE)
    return(det)
  }

  det <- combined_distance_details(
    spectra$frag_list[[idA]], spectra$frag_list[[idB]],
    spectra$loss_list[[idA]], spectra$loss_list[[idB]], params,
    lossA_typ = if (!is.null(spectra$loss_typ_list)) spectra$loss_typ_list[[idA]] else NULL,
    lossB_typ = if (!is.null(spectra$loss_typ_list)) spectra$loss_typ_list[[idB]] else NULL,
    confA = confA, confB = confB
  )

  det$idA <- idA
  det$idB <- idB
  det$confA <- confA
  det$confB <- confB

  if (isTRUE(include_typical) && isTRUE(params$use_typical_loss)) {
    det$typical_loss_A <- typical_table(if (!is.null(spectra$loss_typ_list)) spectra$loss_typ_list[[idA]] else NULL)
    det$typical_loss_B <- typical_table(if (!is.null(spectra$loss_typ_list)) spectra$loss_typ_list[[idB]] else NULL)

    det$typical_loss_matches <- typical_loss_match_table(
      if (!is.null(spectra$loss_typ_list)) spectra$loss_typ_list[[idA]] else NULL,
      if (!is.null(spectra$loss_typ_list)) spectra$loss_typ_list[[idB]] else NULL,
      top_n = top_n,
      typical_losses = typ_univ
    )
  }

  det <- add_contributions(det, split = FALSE)

  det
}

#' Convert component distance matrices into a long table (for plotting)
#'
#' @param dist_res Result of compute_distance_matrices().
#' @param upper_only If TRUE, return only i<j pairs (recommended).
#' @param include_diag If TRUE, include diagonal pairs.
#' @return A data.frame with columns idA, idB, dist, and component columns.
#' @export
distance_components_to_long <- function(dist_res, upper_only = TRUE, include_diag = FALSE) {
  if (is.null(dist_res$dist) || is.null(dist_res$components)) {
    stop("dist_res must be the output of compute_distance_matrices().")
  }
  dist <- dist_res$dist
  comps <- dist_res$components
  ids <- rownames(dist)
  n <- nrow(dist)

  grid <- expand.grid(i = seq_len(n), j = seq_len(n))
  if (isTRUE(upper_only)) {
    if (isTRUE(include_diag)) {
      grid <- grid[grid$i <= grid$j, , drop = FALSE]
    } else {
      grid <- grid[grid$i < grid$j, , drop = FALSE]
    }
  } else if (!isTRUE(include_diag)) {
    grid <- grid[grid$i != grid$j, , drop = FALSE]
  }

  out <- data.frame(
    idA = ids[grid$i],
    idB = ids[grid$j],
    dist = dist[cbind(grid$i, grid$j)],
    stringsAsFactors = FALSE
  )

  for (nm in names(comps)) {
    out[[nm]] <- comps[[nm]][cbind(grid$i, grid$j)]
  }
  out
}
