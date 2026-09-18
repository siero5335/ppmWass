#' Default Parameters for LC-HRMS/MS Fragment Search
#'
#' Returns conservative fragment-only parameters for the experimental LC
#' workflow. Neutral-loss scoring remains disabled until precursor/adduct
#' handling has been validated for the target acquisition method.
#'
#' @return A parameter list compatible with [search_topk_rerank()].
#' @export
lc_default_params <- function() {
  params <- eihrms_default_params()
  params$min_mz <- 20
  params$max_mz <- 2000
  params$topK <- 300L
  params$w_frag <- 1
  params$w_loss <- 0
  params$use_typical_loss <- FALSE
  params$use_split_loss <- FALSE
  params$use_mref_confidence <- FALSE
  params
}

#' Define an LC Library-search Plan
#'
#' The plan separates metadata-based candidate gates from spectral scoring so
#' that every reduction in the candidate set can be audited. A `NULL` numeric
#' gate disables that gate.
#'
#' @param precursor_ppm Maximum precursor-mass difference in ppm, or `NULL`.
#' @param rt_window Maximum absolute retention-time difference, in the same
#'   unit used by query and library metadata, or `NULL`.
#' @param require_ion_mode Require matching positive/negative ion mode.
#' @param require_adduct Require matching normalized adduct labels.
#' @param collision_energy_tolerance Maximum absolute collision-energy
#'   difference, or `NULL`.
#' @param missing_metadata Policy for values missing from an active gate:
#'   `"allow"`, `"exclude"`, or `"error"`.
#' @param exclude_same_id Remove library entries whose ID equals the query ID.
#' @param first_pass_method,rerank_method Distance methods passed to
#'   [search_topk_rerank()].
#' @param first_pass_top_k,final_top_k Candidate and returned hit counts.
#' @return An object of class `LCSearchPlan`.
#' @export
lc_search_plan <- function(
    precursor_ppm = 10,
    rt_window = NULL,
    require_ion_mode = TRUE,
    require_adduct = FALSE,
    collision_energy_tolerance = NULL,
    missing_metadata = c("allow", "exclude", "error"),
    exclude_same_id = FALSE,
    first_pass_method = "entropy_weighted",
    rerank_method = "ppm_wasserstein",
    first_pass_top_k = 50L,
    final_top_k = 10L) {
  missing_metadata <- match.arg(missing_metadata)
  validate_optional_nonnegative <- function(x, name, strictly_positive = FALSE) {
    if (is.null(x)) return(invisible(NULL))
    if (length(x) != 1L || is.na(x) || !is.finite(x) ||
        (strictly_positive && x <= 0) || (!strictly_positive && x < 0)) {
      stop(name, if (strictly_positive) " must be a positive finite number or NULL." else
        " must be a non-negative finite number or NULL.")
    }
  }
  validate_optional_nonnegative(precursor_ppm, "precursor_ppm", strictly_positive = TRUE)
  validate_optional_nonnegative(rt_window, "rt_window")
  validate_optional_nonnegative(collision_energy_tolerance, "collision_energy_tolerance")
  if (!is.logical(require_ion_mode) || length(require_ion_mode) != 1L || is.na(require_ion_mode)) {
    stop("require_ion_mode must be TRUE/FALSE.")
  }
  if (!is.logical(require_adduct) || length(require_adduct) != 1L || is.na(require_adduct)) {
    stop("require_adduct must be TRUE/FALSE.")
  }
  if (!is.logical(exclude_same_id) || length(exclude_same_id) != 1L || is.na(exclude_same_id)) {
    stop("exclude_same_id must be TRUE/FALSE.")
  }
  validate_distance_method_name(first_pass_method, "first_pass_method")
  validate_distance_method_name(rerank_method, "rerank_method")
  first_pass_top_k <- as.integer(first_pass_top_k)
  final_top_k <- as.integer(final_top_k)
  if (length(first_pass_top_k) != 1L || is.na(first_pass_top_k) || first_pass_top_k < 1L) {
    stop("first_pass_top_k must be >= 1.")
  }
  if (length(final_top_k) != 1L || is.na(final_top_k) || final_top_k < 1L) {
    stop("final_top_k must be >= 1.")
  }
  if (final_top_k > first_pass_top_k) {
    stop("final_top_k must be <= first_pass_top_k.")
  }

  structure(
    list(
      precursor_ppm = precursor_ppm,
      rt_window = rt_window,
      require_ion_mode = require_ion_mode,
      require_adduct = require_adduct,
      collision_energy_tolerance = collision_energy_tolerance,
      missing_metadata = missing_metadata,
      exclude_same_id = exclude_same_id,
      first_pass_method = first_pass_method,
      rerank_method = rerank_method,
      first_pass_top_k = first_pass_top_k,
      final_top_k = final_top_k
    ),
    class = "LCSearchPlan"
  )
}

lc_empty_spectrum <- function() {
  matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("mz", "intensity")))
}

lc_clean_spectrum <- function(x, id = NULL) {
  if (is.null(x) || length(x) == 0L) return(lc_empty_spectrum())
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || ncol(x) < 2L) {
    stop("Spectrum", if (!is.null(id)) paste0(" '", id, "'") else "",
         " must be a matrix or data frame with at least two columns.")
  }
  out <- cbind(mz = suppressWarnings(as.numeric(x[, 1L])),
               intensity = suppressWarnings(as.numeric(x[, 2L])))
  keep <- is.finite(out[, 1L]) & is.finite(out[, 2L]) &
    out[, 1L] > 0 & out[, 2L] > 0
  out <- out[keep, , drop = FALSE]
  out[order(out[, 1L]), , drop = FALSE]
}

lc_normalize_ion_mode <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[y %in% c("+", "pos", "positive", "positive ion", "esi+")] <- "positive"
  y[y %in% c("-", "neg", "negative", "negative ion", "esi-")] <- "negative"
  y[is.na(x) | !nzchar(y)] <- NA_character_
  y
}

lc_normalize_adduct <- function(x) {
  y <- toupper(gsub("\\s+", "", trimws(as.character(x))))
  y[is.na(x) | !nzchar(y)] <- NA_character_
  y
}

lc_parse_numeric <- function(x) {
  x <- as.character(x)
  hit <- regexpr("[-+]?[0-9]*\\.?[0-9]+", x, perl = TRUE)
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x) & hit > 0L
  out[ok] <- suppressWarnings(as.numeric(regmatches(x, hit)[ok]))
  out
}

lc_canonical_metadata <- function(metadata, ids) {
  n <- length(ids)
  if (is.null(metadata)) {
    metadata <- data.frame(id = ids, stringsAsFactors = FALSE)
  }
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (nrow(metadata) != n) {
    stop("metadata must have one row per spectrum.")
  }
  if (!"id" %in% names(metadata)) metadata$id <- ids
  metadata$id <- as.character(metadata$id)
  if (anyNA(metadata$id) || any(!nzchar(metadata$id)) || anyDuplicated(metadata$id)) {
    stop("metadata$id must contain unique, non-missing IDs.")
  }
  if (!setequal(metadata$id, ids)) {
    stop("metadata$id must match the spectrum IDs.")
  }
  metadata <- metadata[match(ids, metadata$id), , drop = FALSE]

  first_column <- function(candidates, default = NULL) {
    hit <- candidates[candidates %in% names(metadata)]
    if (length(hit)) metadata[[hit[[1L]]]] else default
  }
  metadata$precursor_mz <- suppressWarnings(as.numeric(first_column(
    c("precursor_mz", "precursormz", "mw"), rep(NA_real_, n)
  )))
  metadata$rt <- suppressWarnings(as.numeric(first_column(
    c("rt", "retention_time"), rep(NA_real_, n)
  )))
  metadata$ion_mode <- lc_normalize_ion_mode(first_column(
    c("ion_mode", "ionmode"), rep(NA_character_, n)
  ))
  metadata$adduct <- lc_normalize_adduct(first_column(
    c("adduct", "precursor_type", "precursortype"), rep(NA_character_, n)
  ))
  metadata$collision_energy <- lc_parse_numeric(first_column(
    c("collision_energy", "collisionenergy", "ce"), rep(NA_character_, n)
  ))

  canonical <- c("id", "precursor_mz", "rt", "ion_mode", "adduct", "collision_energy")
  metadata[, c(canonical, setdiff(names(metadata), canonical)), drop = FALSE]
}

#' Create a Validated LC Spectra Container
#'
#' @param frag_list Named list of two-column fragment spectra.
#' @param metadata Data frame with one row per spectrum. Canonical columns are
#'   `id`, `precursor_mz`, `rt`, `ion_mode`, `adduct`, and `collision_energy`.
#'   Common MSP aliases are accepted.
#' @param ids IDs aligned with `frag_list`.
#' @param preprocess If `TRUE`, run the package peak preprocessing using
#'   `params`; otherwise retain all valid peaks and normalize each spectrum.
#' @param params Parameters from [lc_default_params()].
#' @param strict If `TRUE`, reject empty spectra after processing.
#' @return An object of class `LCSpectra` containing spectra, metadata, and QC.
#' @export
as_lc_spectra <- function(frag_list, metadata = NULL, ids = names(frag_list),
                          preprocess = FALSE, params = lc_default_params(),
                          strict = TRUE) {
  if (!is.list(frag_list) || length(frag_list) < 1L) {
    stop("frag_list must contain at least one spectrum.")
  }
  if (is.null(ids)) ids <- as.character(seq_along(frag_list))
  ids <- as.character(ids)
  if (length(ids) != length(frag_list) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("ids must be unique, non-missing, and aligned with frag_list.")
  }

  cleaned <- lapply(seq_along(frag_list), function(i) {
    x <- lc_clean_spectrum(frag_list[[i]], ids[[i]])
    if (isTRUE(preprocess)) {
      prep_peaks(x, params)
    } else if (nrow(x) && sum(x[, 2L]) > 0) {
      x[, 2L] <- x[, 2L] / sum(x[, 2L])
      x
    } else {
      lc_empty_spectrum()
    }
  })
  names(cleaned) <- ids
  metadata <- lc_canonical_metadata(metadata, ids)
  qc <- data.frame(
    id = ids,
    n_peaks = vapply(cleaned, nrow, integer(1L)),
    min_mz = vapply(cleaned, function(x) if (nrow(x)) min(x[, 1L]) else NA_real_, numeric(1L)),
    max_mz = vapply(cleaned, function(x) if (nrow(x)) max(x[, 1L]) else NA_real_, numeric(1L)),
    tic = vapply(cleaned, function(x) if (nrow(x)) sum(x[, 2L]) else 0, numeric(1L)),
    valid = vapply(cleaned, function(x) nrow(x) > 0L && all(is.finite(x)), logical(1L)),
    stringsAsFactors = FALSE
  )
  if (isTRUE(strict) && any(!qc$valid)) {
    stop("Invalid or empty spectra: ", paste(qc$id[!qc$valid], collapse = ", "))
  }

  structure(list(frag_list = cleaned, metadata = metadata, qc = qc), class = "LCSpectra")
}

#' Read an LC-MS/MS MSP Library
#'
#' Reads MSP records without invoking the GC-specific Mref, derivatization, or
#' neutral-loss pipeline. Duplicate entry names receive stable unique suffixes.
#'
#' @param file Path to an MSP file.
#' @param params Parameters from [lc_default_params()].
#' @param id_column MSP metadata column used as the entry ID.
#' @param progress Progress interval passed to [read_msp()].
#' @param drop_invalid Drop empty spectra after preprocessing.
#' @return An `LCSpectra` object.
#' @export
read_lc_msp <- function(file, params = lc_default_params(), id_column = "name",
                        progress = 1000L, drop_invalid = TRUE) {
  raw <- read_msp(file, progress = progress)
  if (!id_column %in% names(raw)) stop("id_column not found in MSP data: ", id_column)
  ids <- as.character(raw[[id_column]])
  missing_id <- is.na(ids) | !nzchar(ids)
  ids[missing_id] <- paste0("entry_", which(missing_id))
  ids <- make.unique(ids, sep = "__")
  frag <- lapply(raw$spectrum, function(x) prep_peaks(parse_ei(x), params))
  names(frag) <- ids
  raw$id <- ids
  metadata <- lc_canonical_metadata(raw, ids)
  valid <- vapply(frag, nrow, integer(1L)) > 0L
  if (any(!valid)) {
    if (!isTRUE(drop_invalid)) {
      stop("MSP file contains ", sum(!valid), " empty spectra after preprocessing.")
    }
    warning("Dropped ", sum(!valid), " empty spectra after preprocessing.", call. = FALSE)
    frag <- frag[valid]
    metadata <- metadata[valid, , drop = FALSE]
  }
  as_lc_spectra(frag, metadata = metadata, preprocess = FALSE, strict = TRUE)
}

lc_gate_missing <- function(query_value, library_value, comparison, policy, field) {
  missing <- is.na(query_value) | is.na(library_value)
  if (policy == "error" && any(missing)) {
    stop("Missing metadata in active LC candidate gate: ", field)
  }
  matched <- comparison(query_value, library_value)
  matched[is.na(matched)] <- FALSE
  if (policy == "allow") matched[missing] <- TRUE
  matched
}

#' Build Auditable LC Candidate Sets
#'
#' @param query,library `LCSpectra` objects.
#' @param plan An object from [lc_search_plan()].
#' @return A list with library row indices per query and an audit table.
#' @export
lc_candidate_sets <- function(query, library, plan = lc_search_plan()) {
  if (!inherits(query, "LCSpectra") || !inherits(library, "LCSpectra")) {
    stop("query and library must be LCSpectra objects.")
  }
  if (!inherits(plan, "LCSearchPlan")) stop("plan must be an LCSearchPlan object.")

  qmeta <- query$metadata
  lmeta <- library$metadata
  indices <- vector("list", nrow(qmeta))
  names(indices) <- qmeta$id
  audit_rows <- vector("list", nrow(qmeta))

  for (i in seq_len(nrow(qmeta))) {
    keep <- rep(TRUE, nrow(lmeta))
    counts <- c(initial = sum(keep))
    apply_gate <- function(gate, name) {
      keep <<- keep & gate
      counts[[name]] <<- sum(keep)
    }

    if (!is.null(plan$precursor_ppm)) {
      q <- qmeta$precursor_mz[[i]]
      gate <- lc_gate_missing(q, lmeta$precursor_mz, function(a, b) {
        is.finite(a) & a > 0 & is.finite(b) & abs(b - a) / a * 1e6 <= plan$precursor_ppm
      }, plan$missing_metadata, "precursor_mz")
      apply_gate(gate, "after_precursor")
    }
    if (isTRUE(plan$require_ion_mode)) {
      gate <- lc_gate_missing(qmeta$ion_mode[[i]], lmeta$ion_mode, `==`,
                              plan$missing_metadata, "ion_mode")
      apply_gate(gate, "after_ion_mode")
    }
    if (isTRUE(plan$require_adduct)) {
      gate <- lc_gate_missing(qmeta$adduct[[i]], lmeta$adduct, `==`,
                              plan$missing_metadata, "adduct")
      apply_gate(gate, "after_adduct")
    }
    if (!is.null(plan$rt_window)) {
      q <- qmeta$rt[[i]]
      gate <- lc_gate_missing(q, lmeta$rt, function(a, b) abs(b - a) <= plan$rt_window,
                              plan$missing_metadata, "rt")
      apply_gate(gate, "after_rt")
    }
    if (!is.null(plan$collision_energy_tolerance)) {
      q <- qmeta$collision_energy[[i]]
      gate <- lc_gate_missing(q, lmeta$collision_energy,
                              function(a, b) abs(b - a) <= plan$collision_energy_tolerance,
                              plan$missing_metadata, "collision_energy")
      apply_gate(gate, "after_collision_energy")
    }
    if (isTRUE(plan$exclude_same_id)) {
      apply_gate(lmeta$id != qmeta$id[[i]], "after_same_id")
    }
    indices[[i]] <- which(keep)
    counts[["final"]] <- sum(keep)
    row <- as.data.frame(as.list(counts), stringsAsFactors = FALSE)
    row$query_id <- qmeta$id[[i]]
    audit_rows[[i]] <- row[, c("query_id", setdiff(names(row), "query_id")), drop = FALSE]
  }

  all_names <- unique(unlist(lapply(audit_rows, names), use.names = FALSE))
  audit_rows <- lapply(audit_rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA_integer_
    x[, all_names, drop = FALSE]
  })
  list(indices = indices, audit = do.call(rbind, audit_rows), plan = plan)
}

#' Search an LC Library with Metadata Candidate Gates
#'
#' This experimental outer workflow performs metadata gating first and then
#' calls [search_topk_rerank()] independently for each query's candidate set.
#' Scoring is fragment-only by design.
#'
#' @param query,library `LCSpectra` objects.
#' @param plan An object from [lc_search_plan()].
#' @param params Parameters from [lc_default_params()].
#' @param progress Emit per-query progress.
#' @param use_parallel,n_cores Parallel settings for spectral scoring.
#' @return A list containing long-form results, candidate audit, plan, and time.
#' @export
search_lc_library <- function(query, library, plan = lc_search_plan(),
                              params = lc_default_params(), progress = TRUE,
                              use_parallel = params$use_parallel %||% FALSE,
                              n_cores = params$n_cores %||% 1L) {
  started <- proc.time()[["elapsed"]]
  candidate_info <- lc_candidate_sets(query, library, plan)
  result_rows <- vector("list", length(query$frag_list))

  for (i in seq_along(query$frag_list)) {
    idx <- candidate_info$indices[[i]]
    if (!length(idx)) {
      result_rows[[i]] <- NULL
      next
    }
    first_k <- min(plan$first_pass_top_k, length(idx))
    final_k <- min(plan$final_top_k, first_k)
    if (isTRUE(progress)) {
      message("[LC search] query ", i, "/", length(query$frag_list),
              ": ", length(idx), " candidates")
    }
    one <- search_topk_rerank(
      query_frag_list = query$frag_list[i],
      query_loss_list = NULL,
      lib_frag_list = library$frag_list[idx],
      lib_loss_list = NULL,
      params = params,
      first_pass_method = plan$first_pass_method,
      rerank_method = plan$rerank_method,
      first_pass_top_k = first_k,
      final_top_k = final_k,
      query_ids = query$metadata$id[i],
      lib_ids = library$metadata$id[idx],
      exclude_self = FALSE,
      use_parallel = use_parallel,
      n_cores = n_cores,
      progress = FALSE
    )
    one$results$n_gated_candidates <- length(idx)
    result_rows[[i]] <- one$results
  }

  nonempty <- result_rows[lengths(result_rows) > 0L]
  results <- if (length(nonempty)) do.call(rbind, nonempty) else data.frame()
  rownames(results) <- NULL
  summary <- candidate_info$audit
  returned <- if (nrow(results)) table(results$query_id) else integer()
  summary$n_returned <- as.integer(returned[match(summary$query_id, names(returned))])
  summary$n_returned[is.na(summary$n_returned)] <- 0L

  list(
    results = results,
    candidate_summary = summary,
    plan = plan,
    params = params,
    timing = list(total_sec = proc.time()[["elapsed"]] - started)
  )
}
