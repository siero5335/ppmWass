#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]])
}

repo_dir <- normalizePath(get_arg("repo", getwd()), mustWork = TRUE)
input_dir <- normalizePath(get_arg("input-dir", "."), mustWork = TRUE)
output_dir <- get_arg("output-dir", file.path(dirname(repo_dir), "audit"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
pkgload::load_all(repo_dir, quiet = TRUE)

inputs <- c(
  RECETOX = file.path(input_dir, "RECETOX_merged.msp"),
  `HREI-MSDB` = file.path(input_dir, "HREI-MSDB.msp"),
  NOMINAL = file.path(input_dir, "combine_ei_matched_v3.msp")
)
stopifnot(all(file.exists(inputs)))

split_spectrum <- function(x) {
  if (is.na(x) || !nchar(x)) {
    return(data.frame(mz = numeric(), intensity = numeric()))
  }
  token <- strsplit(x, "[[:space:]]+")[[1]]
  pieces <- strsplit(token, ":", fixed = TRUE)
  data.frame(
    mz = vapply(pieces, function(z) as.numeric(z[[1]]), numeric(1)),
    intensity = vapply(pieces, function(z) as.numeric(z[[2]]), numeric(1))
  )
}

count_peak_encodings <- function(path) {
  lines <- readLines(path, warn = FALSE)
  in_peaks <- FALSE
  colon <- 0L
  whitespace <- 0L
  annotated <- 0L
  for (line in lines) {
    txt <- trimws(line)
    if (!nchar(txt)) {
      in_peaks <- FALSE
      next
    }
    if (grepl("^Num[ _]*Peaks[[:space:]]*:", txt, ignore.case = TRUE)) {
      in_peaks <- TRUE
      next
    }
    if (!in_peaks) next
    if (grepl("^(Name|Title)[[:space:]]*:", txt, ignore.case = TRUE)) {
      in_peaks <- FALSE
      next
    }
    if (grepl("^[+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?[[:space:]]*:", txt,
              perl = TRUE)) {
      colon <- colon + 1L
    } else if (grepl("^[+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?[[:space:]]+", txt,
                     perl = TRUE)) {
      whitespace <- whitespace + 1L
    }
    if (grepl("[\"']", txt)) annotated <- annotated + 1L
  }
  c(colon_peak_lines = colon, whitespace_peak_lines = whitespace,
    annotated_peak_lines = annotated)
}

scan_peak_records <- function(path) {
  lines <- readLines(path, warn = FALSE)
  repair <- getFromNamespace("repair_embedded_msp_record_starts", "ppmWass")
  lines <- repair(lines)$lines
  record_index <- 0L
  in_peaks <- FALSE
  rows <- list()
  ensure_record <- function() {
    if (record_index < 1L || length(rows) >= record_index) return()
    rows[[record_index]] <<- list(
      colon_lines = 0L, whitespace_lines = 0L,
      odd_numeric_token_lines = 0L, malformed_peak_lines = 0L
    )
  }
  number <- "[+-]?(?:(?:[0-9]+(?:\\.[0-9]*)?)|(?:\\.[0-9]+))(?:[eE][+-]?[0-9]+)?"
  for (line in lines) {
    txt <- trimws(line)
    if (!nchar(txt)) {
      in_peaks <- FALSE
      next
    }
    if (grepl("^(Name|Title)[[:space:]]*:", txt, ignore.case = TRUE)) {
      record_index <- record_index + 1L
      ensure_record()
      in_peaks <- FALSE
      next
    }
    if (grepl("^Num[ _]*Peaks[[:space:]]*:", txt, ignore.case = TRUE)) {
      if (record_index < 1L) {
        record_index <- 1L
        ensure_record()
      }
      in_peaks <- TRUE
      next
    }
    if (!in_peaks || record_index < 1L) next
    ensure_record()
    numeric_part <- sub("[\"'].*$", "", txt)
    numeric_tokens <- regmatches(
      numeric_part, gregexpr(number, numeric_part, perl = TRUE)
    )[[1L]]
    numeric_tokens <- numeric_tokens[nzchar(numeric_tokens)]
    starts_numeric <- grepl(paste0("^[[:space:]]*\\(?[[:space:]]*", number),
                            txt, perl = TRUE)
    colon_line <- grepl(
      paste0("^[[:space:]]*\\(?[[:space:]]*", number,
             "[[:space:]]*:[[:space:]]*", number),
      txt, perl = TRUE
    )
    whitespace_line <- !colon_line && grepl(
      paste0("^[[:space:]]*\\(?[[:space:]]*", number,
             "[[:space:]]+", number),
      txt, perl = TRUE
    )
    if (colon_line) rows[[record_index]]$colon_lines <-
      rows[[record_index]]$colon_lines + 1L
    if (whitespace_line) rows[[record_index]]$whitespace_lines <-
      rows[[record_index]]$whitespace_lines + 1L
    if (starts_numeric && length(numeric_tokens) %% 2L == 1L) {
      rows[[record_index]]$odd_numeric_token_lines <-
        rows[[record_index]]$odd_numeric_token_lines + 1L
    }
    if (starts_numeric && (!colon_line && !whitespace_line ||
                           length(numeric_tokens) < 2L)) {
      rows[[record_index]]$malformed_peak_lines <-
        rows[[record_index]]$malformed_peak_lines + 1L
    }
  }
  if (!length(rows)) {
    return(data.frame(
      record_index = integer(), parser_branch = character(),
      colon_lines = integer(), whitespace_lines = integer(),
      odd_numeric_token_lines = integer(), malformed_peak_lines = integer()
    ))
  }
  do.call(rbind, lapply(seq_along(rows), function(i) {
    x <- rows[[i]]
    branch <- if (x$colon_lines > 0L && x$whitespace_lines > 0L) {
      "mixed_colon_and_whitespace"
    } else if (x$colon_lines > 0L) {
      "colon_pairs"
    } else if (x$whitespace_lines > 0L) {
      "whitespace_pairs"
    } else {
      "no_recognized_peak_line"
    }
    data.frame(
      record_index = i, parser_branch = branch,
      colon_lines = x$colon_lines, whitespace_lines = x$whitespace_lines,
      odd_numeric_token_lines = x$odd_numeric_token_lines,
      malformed_peak_lines = x$malformed_peak_lines,
      stringsAsFactors = FALSE
    )
  }))
}

precision_rows <- list()
record_rows <- list()
problem_text <- character()

for (dataset in names(inputs)) {
  path <- inputs[[dataset]]
  message("Auditing ", dataset, ": ", path)
  lib <- read_msp(path, progress = 0)
  parser_diag <- attr(lib, "msp_parser_diagnostics")
  enc <- count_peak_encodings(path)
  raw_record_scan <- scan_peak_records(path)

  peak_list <- lapply(lib$spectrum, split_spectrum)
  parsed_counts <- vapply(peak_list, nrow, integer(1))
  declared <- if ("num_peaks" %in% names(lib)) lib$num_peaks else rep(NA_integer_, nrow(lib))
  mismatch <- is.finite(declared) & declared != parsed_counts

  issue_by_record <- if (nrow(parser_diag$record_issues)) {
    split(parser_diag$record_issues$issue_type,
          parser_diag$record_issues$record_index)
  } else list()
  issue_text <- vapply(seq_len(nrow(lib)), function(i) {
    value <- unique(issue_by_record[[as.character(i)]])
    if (is.null(value) || !length(value)) "" else paste(value, collapse = ";")
  }, character(1L))
  raw_record_scan <- raw_record_scan[match(seq_len(nrow(lib)),
                                           raw_record_scan$record_index), , drop = FALSE]
  record_rows[[dataset]] <- data.frame(
    dataset = dataset,
    record_index = seq_len(nrow(lib)),
    name = lib$name,
    declared_num_peaks = declared,
    parsed_num_peaks = parsed_counts,
    declared_parsed_mismatch = mismatch,
    parser_branch = raw_record_scan$parser_branch,
    colon_peak_lines = raw_record_scan$colon_lines,
    whitespace_peak_lines = raw_record_scan$whitespace_lines,
    odd_numeric_token_lines = raw_record_scan$odd_numeric_token_lines,
    malformed_peak_lines = raw_record_scan$malformed_peak_lines,
    parser_issue_types = issue_text,
    stringsAsFactors = FALSE
  )

  all_peaks <- do.call(rbind, peak_list)
  old_mz <- as.numeric(sprintf("%.4f", all_peaks$mz))
  old_intensity <- as.numeric(sprintf("%.0f", all_peaks$intensity))
  mz_abs <- abs(old_mz - all_peaks$mz)
  mz_ppm <- mz_abs / all_peaks$mz * 1e6
  int_abs <- abs(old_intensity - all_peaks$intensity)

  q <- function(x, p) {
    if (!length(x)) return(NA_real_)
    unname(stats::quantile(x, p, na.rm = TRUE, type = 8))
  }
  converted <- convert_msp_to_internal(lib, require_ri = FALSE)
  mz_ranges <- cut(
    all_peaks$mz,
    breaks = c(-Inf, 100, 200, 300, 400, 500, Inf),
    right = FALSE,
    labels = c("<100", "100-<200", "200-<300", "300-<400",
               "400-<500", ">=500")
  )
  make_precision_row <- function(index, mz_range) data.frame(
    dataset = dataset, mz_range = mz_range,
    input_file = normalizePath(path),
    legacy_name_or_title_lines = sum(grepl(
      "^[[:space:]]*(Name|Title)[[:space:]]*:",
      readLines(path, warn = FALSE), ignore.case = TRUE
    )),
    fixed_records = nrow(lib),
    converted_spectra = nrow(converted),
    excluded_during_conversion = nrow(lib) - nrow(converted),
    num_peaks_fields = sum(grepl(
      "^[[:space:]]*Num[ _]*Peaks[[:space:]]*:",
      readLines(path, warn = FALSE), ignore.case = TRUE
    )),
    embedded_record_repairs = parser_diag$n_embedded_record_repairs,
    colon_peak_lines = enc[["colon_peak_lines"]],
    whitespace_peak_lines = enc[["whitespace_peak_lines"]],
    annotated_peak_lines = enc[["annotated_peak_lines"]],
    records_colon_only = sum(raw_record_scan$parser_branch == "colon_pairs", na.rm = TRUE),
    records_whitespace_only = sum(raw_record_scan$parser_branch == "whitespace_pairs", na.rm = TRUE),
    records_mixed_encoding = sum(raw_record_scan$parser_branch ==
                                   "mixed_colon_and_whitespace", na.rm = TRUE),
    odd_numeric_token_lines = sum(raw_record_scan$odd_numeric_token_lines, na.rm = TRUE),
    malformed_peak_lines = sum(raw_record_scan$malformed_peak_lines, na.rm = TRUE),
    records_with_duplicate_metadata = length(unique(
      parser_diag$record_issues$record_index[
        parser_diag$record_issues$issue_type == "duplicate_metadata_field"
      ]
    )),
    parsed_peaks = sum(index),
    mz_changed_by_legacy_format = sum(mz_abs[index] > 0),
    mz_max_abs_da = if (any(index)) max(mz_abs[index], na.rm = TRUE) else NA_real_,
    mz_median_abs_da = if (any(index)) stats::median(mz_abs[index], na.rm = TRUE) else NA_real_,
    mz_p95_abs_da = q(mz_abs[index], 0.95),
    mz_p99_abs_da = q(mz_abs[index], 0.99),
    ppm_median = if (any(index)) stats::median(mz_ppm[index], na.rm = TRUE) else NA_real_,
    ppm_p95 = q(mz_ppm[index], 0.95),
    ppm_p99 = q(mz_ppm[index], 0.99),
    ppm_max = if (any(index)) max(mz_ppm[index], na.rm = TRUE) else NA_real_,
    intensity_changed_by_legacy_format = sum(int_abs[index] > 0),
    intensity_max_abs = if (any(index)) max(int_abs[index], na.rm = TRUE) else NA_real_,
    positive_intensity_rounded_to_zero = sum(
      index & all_peaks$intensity > 0 & old_intensity == 0
    ),
    records_with_declared_count_mismatch = sum(mismatch, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  precision_rows[[dataset]] <- rbind(
    make_precision_row(rep(TRUE, nrow(all_peaks)), "all_mz"),
    do.call(rbind, lapply(levels(mz_ranges), function(label) {
      make_precision_row(!is.na(mz_ranges) & mz_ranges == label, label)
    }))
  )

  if (parser_diag$n_embedded_record_repairs > 0) {
    repair_lines <- apply(parser_diag$embedded_record_repairs, 1, function(z) {
      paste0(
        dataset, " embedded record repair at source line ", z[["source_line"]],
        ": ", z[["peak_prefix"]], " | ", z[["record_start"]]
      )
    })
    problem_text <- c(problem_text, repair_lines)
  }
  if (any(mismatch, na.rm = TRUE)) {
    idx <- which(mismatch)
    problem_text <- c(
      problem_text,
      paste0(
        dataset, " declared/parsed mismatch record ", idx,
        " name=", lib$name[idx],
        " declared=", declared[idx],
        " parsed=", parsed_counts[idx]
      )
    )
  }
}

precision <- do.call(rbind, precision_rows)
records <- do.call(rbind, record_rows)
utils::write.csv(precision, file.path(output_dir, "msp_precision_qc.csv"), row.names = FALSE)
utils::write.csv(records, file.path(output_dir, "msp_record_qc.csv"), row.names = FALSE)
writeLines(problem_text, file.path(output_dir, "msp_problem_records.txt"), useBytes = TRUE)
saveRDS(list(precision = precision, records = records),
        file.path(output_dir, "msp_parser_audit.rds"))
git_commit <- tryCatch(
  trimws(system2("git", c("-C", repo_dir, "rev-parse", "HEAD"),
                  stdout = TRUE, stderr = FALSE)),
  error = function(e) NA_character_
)
git_status <- tryCatch(
  system2("git", c("-C", repo_dir, "status", "--porcelain"),
          stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
utils::write.csv(
  data.frame(
    parameter = c("timestamp_utc", "package_commit", "package_tree_dirty",
                  "command", "legacy_rounding_model"),
    value = c(format(Sys.time(), tz = "UTC"), git_commit,
              length(git_status) > 0L, paste(commandArgs(), collapse = " "),
              "whitespace branch sprintf(%.4f:%.0f)"),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "msp_parser_audit_parameters.csv"), row.names = FALSE
)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "msp_parser_audit_sessionInfo.txt"))
message("Parser audit complete: ", output_dir)
