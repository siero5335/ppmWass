#' @keywords internal
#' @noRd
normalize_benchmark_grid <- function(grid) {
  if (is.data.frame(grid)) {
    if (nrow(grid) == 0 || ncol(grid) == 0) {
      stop("grid must have at least one row and one column.")
    }
    if (is.null(names(grid)) || any(names(grid) == "")) {
      stop("grid must have non-empty column names.")
    }
    return(grid)
  }

  if (is.list(grid)) {
    if (length(grid) == 0 || is.null(names(grid)) || any(names(grid) == "")) {
      stop("grid must be a non-empty named list or data frame.")
    }
    empty <- vapply(grid, length, integer(1)) == 0
    if (any(empty)) {
      stop("grid entries must contain at least one value: ",
           paste(names(grid)[empty], collapse = ", "))
    }
    return(expand.grid(grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  }

  stop("grid must be a data.frame or named list.")
}

#' @keywords internal
#' @noRd
coerce_grid_value <- function(x) {
  if (is.factor(x)) {
    return(as.character(x))
  }
  x
}

#' @keywords internal
#' @noRd
apply_param_overrides <- function(params, overrides) {
  unknown <- setdiff(names(overrides), names(params))
  if (length(unknown) > 0) {
    stop("Unknown parameter name(s) in grid: ", paste(unknown, collapse = ", "))
  }

  for (nm in names(overrides)) {
    params[[nm]] <- coerce_grid_value(overrides[[nm]])
  }

  validate_params(params)
}

#' @keywords internal
#' @noRd
subset_named_list_by_ids <- function(x, ids) {
  if (is.null(x)) return(NULL)
  if (is.null(names(x))) {
    stop("Expected a named list or vector when subsetting spectra_result.")
  }
  x[ids]
}

#' @keywords internal
#' @noRd
subset_ground_truth_by_ids <- function(ground_truth, ids) {
  if (is.null(ground_truth)) {
    return(NULL)
  }
  if (!is.data.frame(ground_truth) || !"id" %in% names(ground_truth)) {
    stop("ground_truth must be a data frame containing an 'id' column.")
  }

  ground_truth[ground_truth$id %in% ids, , drop = FALSE]
}

#' @keywords internal
#' @noRd
subset_spectra_result_by_ids <- function(spectra_result, ids) {
  out <- spectra_result

  list_keys <- c(
    "frag_list", "loss_list", "loss_typ_list",
    "loss_anchor_list", "loss_pair_list",
    "loss_anchor_typ_list", "loss_pair_typ_list"
  )
  for (nm in list_keys) {
    if (!is.null(out[[nm]])) {
      out[[nm]] <- subset_named_list_by_ids(out[[nm]], ids)
    }
  }

  vector_keys <- c("ri", "mref_confidence")
  for (nm in vector_keys) {
    if (!is.null(out[[nm]])) {
      out[[nm]] <- subset_named_list_by_ids(out[[nm]], ids)
    }
  }

  if (!is.null(out$df_spec)) {
    idx <- match(ids, out$df_spec$id)
    out$df_spec <- out$df_spec[idx, , drop = FALSE]
  }

  out
}

#' Sweep Distance Parameters Across a Benchmark Grid
#'
#' Runs `run_benchmark()` repeatedly over a grid of distance-parameter settings
#' and collects the summary metrics into one tidy table.
#'
#' @param spectra_result Output from `build_spectra()` or `build_spectra_from_msp()`.
#' @param grid Parameter combinations as either a data frame or a named list
#'   that will be expanded with `expand.grid()`.
#' @param methods Character vector of distance methods to compare.
#' @param ground_truth Optional data frame with true groupings.
#' @param params Base parameter list.
#' @param metrics Optional subset of summary metric column names to keep.
#'   `NULL` keeps all metrics returned by `run_benchmark()`.
#' @param n_replicates Number of repeated benchmark runs per grid row.
#' @param sample_frac Fraction of spectra to sample without replacement for each
#'   replicate. Use values below 1 to estimate stability by repeated
#'   subsampling; `1` reuses the full dataset.
#' @param seed Random seed used for replicate subsampling.
#' @param keep_results If `TRUE`, also return the full benchmark object for each
#'   grid row.
#' @param progress If `TRUE`, print sweep progress.
#' @param ... Additional arguments passed to `run_benchmark()`.
#' @return A list with `summary`, `grid`, and optionally `results`.
#'
#' @details
#' By default each grid row is evaluated once on the full dataset. Set
#' `n_replicates > 1` together with `sample_frac < 1` to obtain repeated
#' subsampled estimates that can be summarized into confidence intervals or
#' stability plots for manuscript reporting.
#' @export
sweep_distance_params <- function(spectra_result,
                                  grid,
                                  methods = c("entropy", "hellinger", "cosine"),
                                  ground_truth = NULL,
                                  params = eihrms_default_params(),
                                  metrics = NULL,
                                  n_replicates = 1,
                                  sample_frac = 1,
                                  seed = 42,
                                  keep_results = FALSE,
                                  progress = TRUE,
                                  ...) {
  grid_df <- normalize_benchmark_grid(grid)
  if (!is.numeric(n_replicates) || length(n_replicates) != 1 || n_replicates < 1) {
    stop("n_replicates must be >= 1.")
  }
  if (!is.numeric(sample_frac) || length(sample_frac) != 1 || sample_frac <= 0 || sample_frac > 1) {
    stop("sample_frac must be in (0, 1].")
  }

  ids <- names(spectra_result$frag_list)
  if (length(ids) < 2) {
    stop("At least two spectra are required for benchmark sweeping.")
  }

  total_runs <- nrow(grid_df) * n_replicates
  out_rows <- vector("list", total_runs)
  full_results <- if (isTRUE(keep_results)) vector("list", total_runs) else NULL
  run_idx <- 0L

  for (i in seq_len(nrow(grid_df))) {
    row_overrides <- as.list(grid_df[i, , drop = FALSE])
    params_i <- apply_param_overrides(params, row_overrides)

    for (rep_idx in seq_len(n_replicates)) {
      run_idx <- run_idx + 1L
      ids_sub <- ids
      if (sample_frac < 1) {
        set.seed(seed + i * 10000 + rep_idx)
        n_sub <- max(2L, floor(length(ids) * sample_frac))
        ids_sub <- sample(ids, size = n_sub, replace = FALSE)
      }

      spectra_i <- subset_spectra_result_by_ids(spectra_result, ids_sub)
      ground_truth_i <- subset_ground_truth_by_ids(ground_truth, ids_sub)

      if (isTRUE(progress)) {
        message("Distance-parameter sweep: ", run_idx, " / ", total_runs)
      }

      bench_i <- run_benchmark(
        spectra_result = spectra_i,
        methods = methods,
        ground_truth = ground_truth_i,
        params = params_i,
        ...
      )

      summary_i <- bench_i$summary
      summary_i$sweep_id <- i
      summary_i$replicate <- rep_idx
      summary_i$n_spectra <- length(ids_sub)
      for (nm in names(row_overrides)) {
        summary_i[[nm]] <- coerce_grid_value(row_overrides[[nm]])
      }

      base_cols <- c("sweep_id", "replicate", "n_spectra", names(row_overrides), "method")
      metric_cols <- setdiff(names(summary_i), base_cols)
      if (!is.null(metrics)) {
        missing_metrics <- setdiff(metrics, metric_cols)
        if (length(missing_metrics) > 0) {
          stop("Unknown metric column(s) requested: ", paste(missing_metrics, collapse = ", "))
        }
        keep_cols <- c(base_cols, metrics)
        summary_i <- summary_i[, keep_cols, drop = FALSE]
      }

      out_rows[[run_idx]] <- summary_i
      if (isTRUE(keep_results)) {
        full_results[[run_idx]] <- bench_i
      }
    }
  }

  summary_df <- do.call(rbind, out_rows)
  rownames(summary_df) <- NULL

  list(
    summary = summary_df,
    grid = grid_df,
    results = full_results
  )
}

#' Benchmark Distance Runtime Across Dataset Sizes
#'
#' Measures wall-clock time for `compute_distance_matrix()` over repeated random
#' subsets of increasing size.
#'
#' @param spectra_result Output from `build_spectra()` or `build_spectra_from_msp()`.
#' @param methods Character vector of distance methods to benchmark.
#' @param sizes Integer vector of subset sizes to evaluate.
#' @param n_replicates Number of random subsets per size.
#' @param params Base parameter list.
#' @param seed Random seed for reproducible subset selection.
#' @param progress If `TRUE`, print progress messages.
#' @return A list with `results` and `summary` data frames.
#' @export
benchmark_distance_runtime <- function(spectra_result,
                                       methods = c("entropy", "hellinger", "cosine"),
                                       sizes = NULL,
                                       n_replicates = 3,
                                       params = eihrms_default_params(),
                                       seed = 42,
                                       progress = TRUE) {
  ids <- names(spectra_result$frag_list)
  n_total <- length(ids)
  if (n_total < 2) {
    stop("At least two spectra are required for runtime benchmarking.")
  }

  if (is.null(sizes)) {
    sizes <- unique(pmin(c(100, 500, 1000, n_total), n_total))
  }
  sizes <- sort(unique(as.integer(sizes)))
  if (any(!is.finite(sizes)) || any(sizes < 2)) {
    stop("sizes must contain integers >= 2.")
  }
  if (any(sizes > n_total)) {
    stop("sizes cannot exceed the number of available spectra (", n_total, ").")
  }
  if (!is.numeric(n_replicates) || length(n_replicates) != 1 || n_replicates < 1) {
    stop("n_replicates must be >= 1.")
  }

  row_idx <- 0
  res_rows <- vector("list", length(methods) * length(sizes) * n_replicates)

  for (method in methods) {
    params_i <- params
    params_i$distance_method <- method
    params_i <- validate_params(params_i)

    for (size in sizes) {
      for (rep_idx in seq_len(n_replicates)) {
        row_idx <- row_idx + 1
        set.seed(seed + match(method, methods) * 10000 + size * 100 + rep_idx)
        ids_sub <- if (size == n_total) ids else sample(ids, size, replace = FALSE)
        spectra_sub <- subset_spectra_result_by_ids(spectra_result, ids_sub)

        if (isTRUE(progress)) {
          message("Runtime benchmark: method=", method,
                  " size=", size,
                  " rep=", rep_idx, "/", n_replicates)
        }

        elapsed <- system.time(
          compute_distance_matrix(
            spectra_sub$frag_list,
            spectra_sub$loss_list,
            params_i,
            progress = FALSE,
            loss_typ_list = spectra_sub$loss_typ_list,
            mref_conf = spectra_sub$mref_confidence,
            loss_anchor_list = spectra_sub$loss_anchor_list,
            loss_pair_list = spectra_sub$loss_pair_list,
            loss_anchor_typ_list = spectra_sub$loss_anchor_typ_list,
            loss_pair_typ_list = spectra_sub$loss_pair_typ_list
          )
        )[["elapsed"]]

        res_rows[[row_idx]] <- data.frame(
          method = method,
          size = size,
          replicate = rep_idx,
          n_pairs = size * (size - 1) / 2,
          elapsed_sec = unname(elapsed),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  results_df <- do.call(rbind, res_rows)
  rownames(results_df) <- NULL

  split_key <- interaction(results_df$method, results_df$size, results_df$n_pairs, drop = TRUE)
  summary_rows <- lapply(split(results_df, split_key), function(df) {
    data.frame(
      method = df$method[1],
      size = df$size[1],
      n_pairs = df$n_pairs[1],
      mean = mean(df$elapsed_sec),
      sd = stats::sd(df$elapsed_sec),
      median = stats::median(df$elapsed_sec),
      min = min(df$elapsed_sec),
      max = max(df$elapsed_sec),
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  list(
    results = results_df,
    summary = summary_df
  )
}
