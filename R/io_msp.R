#' Read MSP (NIST format) Spectral Library
#'
#' Reads .msp files from NIST, MoNA, MassBank, or other sources.
#' Supports both GC-EI and LC-MS/MS spectra.
#'
#' @param file Path to .msp file.
#' @param progress Show progress every N spectra (0 = silent).
#' @return A tibble with columns: name, cas, inchikey, formula, mw, ri, 
#'         num_peaks, spectrum (as "mz:intensity" string), and other metadata.
#' @export
#' @examples
#' \dontrun{
#' lib <- read_msp("nist_gcms.msp")
#' head(lib)
#' }
read_msp <- function(file, progress = 1000) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  # Read all lines. Some public MSP exports contain a missing newline between
  # the final peak of one record and the NAME/TITLE field of the next record.
  # Repair only markers embedded in a declared peak block and retain a
  # diagnostic attribute so the repair is auditable.
  raw_lines <- readLines(file, warn = FALSE)
  repaired <- repair_embedded_msp_record_starts(raw_lines)
  lines <- repaired$lines
  n_lines <- length(lines)
  
  # Initialize storage
  records <- list()
  record_idx <- 0
  
  # Current record being parsed
  current <- list()
  in_peaks <- FALSE
  peak_lines <- list()
  peak_idx <- 0
  current_field_counts <- integer(0)
  record_issues <- list()

  add_record_issue <- function(issue_type, field = NA_character_,
                               declared = NA_integer_, parsed = NA_integer_,
                               value = NA_character_) {
    record_issues[[length(record_issues) + 1L]] <<- data.frame(
      record_index = record_idx + 1L,
      record_name = if (is.null(current$name)) NA_character_ else current$name,
      issue_type = issue_type,
      field = field,
      declared_peak_count = declared,
      parsed_peak_count = parsed,
      value = value,
      stringsAsFactors = FALSE
    )
  }

  finalize_current_record <- function() {
    if (length(current) == 0) return(FALSE)

    if (length(peak_lines) > 0) {
      current$spectrum <<- parse_msp_peaks(unlist(peak_lines, use.names = FALSE))
    }
    parsed_peak_count <- if (
      is.null(current$spectrum) || is.na(current$spectrum) ||
      !nzchar(current$spectrum)
    ) 0L else length(strsplit(current$spectrum, " ", fixed = TRUE)[[1L]])
    if (!is.null(current$num_peaks) && is.finite(current$num_peaks) &&
        current$num_peaks != parsed_peak_count) {
      add_record_issue(
        "declared_peak_count_mismatch", field = "num_peaks",
        declared = as.integer(current$num_peaks), parsed = parsed_peak_count
      )
    }
    record_idx <<- record_idx + 1
    records[[record_idx]] <<- current

    if (progress > 0 && record_idx %% progress == 0) {
      message("Read ", record_idx, " spectra...")
    }

    current <<- list()
    in_peaks <<- FALSE
    peak_lines <<- list()
    peak_idx <<- 0
    current_field_counts <<- integer(0)
    TRUE
  }

  looks_like_peak_line <- function(x) {
    if (nchar(x) == 0) return(FALSE)
    txt <- gsub("\\(|\\)", "", x)
    if (grepl("[A-Za-z]", txt)) return(FALSE)
    txt <- gsub(":", " ", txt)
    txt <- gsub("[,;]", " ", txt)
    tokens <- unlist(strsplit(txt, "\\s+"))
    tokens <- tokens[nchar(tokens) > 0]
    if (length(tokens) < 2) return(FALSE)
    nums <- suppressWarnings(as.numeric(tokens))
    sum(is.finite(nums)) == length(tokens)
  }
  
  for (i in seq_len(n_lines)) {
    line <- trimws(lines[i])
    
    # Empty line or end of record
    if (nchar(line) == 0) {
      finalize_current_record()
      next
    }

    line_key_lower <- NA_character_
    if (grepl(":", line, fixed = TRUE)) {
      colon_pos <- regexpr(":", line, fixed = TRUE)
      key <- trimws(substr(line, 1, colon_pos - 1))
      line_key_lower <- tolower(key)
      if (line_key_lower %in% c("name", "title") &&
          (in_peaks || length(peak_lines) > 0)) {
        finalize_current_record()
      }
    }

    # A second Num Peaks field inside the peak block is a malformed duplicate,
    # not a peak line. Preserve the first declaration and report the record.
    if (in_peaks && !is.na(line_key_lower) &&
        line_key_lower %in% c("num peaks", "numpeaks", "num_peaks")) {
      add_record_issue(
        "duplicate_metadata_field", field = "num_peaks",
        value = trimws(substr(line, regexpr(":", line, fixed = TRUE) + 1L,
                              nchar(line)))
      )
      next
    }

    # Defensive: treat numeric-only lines as peaks if Num Peaks is missing
    if (!in_peaks && !is.null(current$name) && looks_like_peak_line(line)) {
      in_peaks <- TRUE
      peak_idx <- peak_idx + 1
      peak_lines[[peak_idx]] <- line
      next
    }
    
    # Check for key-value pairs
    if (grepl(":", line, fixed = TRUE) && !in_peaks) {
      # Split on first colon only
      colon_pos <- regexpr(":", line, fixed = TRUE)
      key <- trimws(substr(line, 1, colon_pos - 1))
      value <- trimws(substr(line, colon_pos + 1, nchar(line)))
      
      # Normalize key names
      key_lower <- tolower(key)

      # Repeated synonyms are an intentional multi-valued MSP convention.
      # Other repeated fields retain the historical last-value behavior but
      # are now explicitly auditable.
      count <- if (key_lower %in% names(current_field_counts)) {
        current_field_counts[[key_lower]]
      } else 0L
      if (count > 0L && key_lower != "synon") {
        add_record_issue(
          "duplicate_metadata_field", field = key_lower, value = value
        )
      }
      current_field_counts[[key_lower]] <- count + 1L
      
      if (key_lower %in% c("name", "title")) {
        current$name <- value
      } else if (key_lower %in% c("cas#", "casno", "cas")) {
        current$cas <- value
      } else if (key_lower == "inchikey") {
        current$inchikey <- value
      } else if (key_lower == "inchi") {
        current$inchi <- value
      } else if (key_lower == "smiles") {
        current$smiles <- value
      } else if (key_lower %in% c("formula", "molecular formula")) {
        current$formula <- value
      } else if (key_lower %in% c("mw", "molweight", "exactmass", "precursormz")) {
        current$mw <- suppressWarnings(as.numeric(value))
      } else if (key_lower %in% c("ri", "retentionindex", "retention_index", "kovats")) {
        current$ri <- suppressWarnings(as.numeric(value))
      } else if (key_lower %in% c("rt", "retentiontime", "retention_time")) {
        current$rt <- suppressWarnings(as.numeric(value))
      } else if (key_lower %in% c("num peaks", "numpeaks", "num_peaks")) {
        current$num_peaks <- suppressWarnings(as.integer(value))
        in_peaks <- TRUE  # Next lines are peak data
      } else if (key_lower == "synon") {
        # Synonyms - append to existing
        if (is.null(current$synonyms)) {
          current$synonyms <- value
        } else {
          current$synonyms <- paste(current$synonyms, value, sep = "; ")
        }
      } else if (key_lower %in% c("comment", "comments")) {
        current$comment <- value
      } else if (key_lower %in% c("instrument", "instrumenttype", "instrument_type")) {
        current$instrument <- value
      } else if (key_lower %in% c("ionmode", "ion_mode", "ionization")) {
        current$ion_mode <- value
      } else if (key_lower == "precursortype") {
        current$precursor_type <- value
      } else if (key_lower %in% c("collisionenergy", "collision_energy", "ce")) {
        current$collision_energy <- value
      } else if (key_lower %in% c("source", "db_source")) {
        current$source <- value
      } else if (key_lower == "splash") {
        current$splash <- value
      } else {
        # Store unknown keys in a generic metadata field
        if (is.null(current$other_metadata)) {
          current$other_metadata <- paste0(key, "=", value)
        } else {
          current$other_metadata <- paste(current$other_metadata, 
                                          paste0(key, "=", value), sep = "; ")
        }
      }
    } else if (in_peaks) {
      # Peak data line
      peak_idx <- peak_idx + 1
      peak_lines[[peak_idx]] <- line
    }
  }
  
  # Don't forget the last record
  finalize_current_record()
  
  message("Total spectra read: ", record_idx)
  
  # Convert to tibble
  out <- msp_list_to_tibble(records)
  attr(out, "msp_parser_diagnostics") <- list(
    input_file = normalizePath(file, mustWork = TRUE),
    input_lines = length(raw_lines),
    parsed_lines = length(lines),
    embedded_record_repairs = repaired$diagnostics,
    n_embedded_record_repairs = nrow(repaired$diagnostics),
    record_issues = if (length(record_issues)) do.call(rbind, record_issues) else
      data.frame(
        record_index = integer(), record_name = character(),
        issue_type = character(), field = character(),
        declared_peak_count = integer(), parsed_peak_count = integer(),
        value = character(), stringsAsFactors = FALSE
      ),
    n_records = nrow(out)
  )
  out
}

#' Repair embedded MSP record starts
#'
#' A narrow repair for malformed lines such as
#' `429.08865 1533679NAME: next record`, where the newline before `NAME:` was
#' lost. Markers are split only while inside a declared peak block and only
#' when immediately preceded by a digit. This avoids splitting ordinary
#' metadata values or comments containing the word "name".
#'
#' @param lines Character vector returned by [readLines()].
#' @return A list with repaired `lines` and a diagnostics data frame.
#' @keywords internal
#' @noRd
repair_embedded_msp_record_starts <- function(lines) {
  # At most one extra output line is created for each input line. A preallocated
  # list avoids quadratic copying for large commercial-library MSP files.
  out <- vector("list", length(lines) * 2L)
  out_idx <- 0L
  diagnostics <- list()
  in_peaks <- FALSE

  append_line <- function(value) {
    out_idx <<- out_idx + 1L
    out[[out_idx]] <<- value
  }

  for (i in seq_along(lines)) {
    line <- lines[[i]]
    trimmed <- trimws(line)

    if (!nchar(trimmed)) {
      append_line(line)
      in_peaks <- FALSE
      next
    }

    if (grepl("^[[:space:]]*(Num[[:space:]_]*Peaks)[[:space:]]*:", line,
              ignore.case = TRUE, perl = TRUE)) {
      in_peaks <- TRUE
      append_line(line)
      next
    }

    if (in_peaks) {
      pos <- regexpr("(?<=[0-9])(?:(?:NAME)|(?:TITLE)):", line,
                     ignore.case = TRUE, perl = TRUE)
      if (pos[[1]] > 1L) {
        prefix <- substr(line, 1L, pos[[1]] - 1L)
        suffix <- substr(line, pos[[1]], nchar(line))
        append_line(prefix)
        append_line(suffix)
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
          source_line = i,
          marker = sub(":.*$", "", suffix),
          peak_prefix = prefix,
          record_start = suffix,
          stringsAsFactors = FALSE
        )
        in_peaks <- FALSE
        next
      }
    }

    # A normal NAME/TITLE line also closes any prior peak block.
    if (grepl("^[[:space:]]*(NAME|TITLE)[[:space:]]*:", line,
              ignore.case = TRUE, perl = TRUE)) {
      in_peaks <- FALSE
    }
    append_line(line)
  }

  diag_df <- if (length(diagnostics)) {
    do.call(rbind, diagnostics)
  } else {
    data.frame(
      source_line = integer(), marker = character(), peak_prefix = character(),
      record_start = character(), stringsAsFactors = FALSE
    )
  }
  list(lines = unlist(out[seq_len(out_idx)], use.names = FALSE), diagnostics = diag_df)
}

#' Parse MSP Peak Lines to MS-DIAL Format
#'
#' @param peak_lines Character vector of peak data lines.
#' @return String in "mz:intensity mz:intensity ..." format.
#' @keywords internal
#' @noRd
parse_msp_peaks <- function(peak_lines) {
  if (length(peak_lines) == 0) return("")

  # Parse each physical line independently. This prevents numeric values in a
  # quoted annotation from being paired with tokens on a later line. Numeric
  # token strings are retained verbatim after validation, so no precision is
  # lost through fixed-width formatting.
  number <- "[+-]?(?:(?:[0-9]+(?:\\.[0-9]*)?)|(?:\\.[0-9]+))(?:[eE][+-]?[0-9]+)?"
  colon_pair <- paste0("^[[:space:]]*\\(?[[:space:]]*(", number,
                       ")[[:space:]]*:[[:space:]]*(", number, ")")
  space_pair <- paste0("^[[:space:]]*\\(?[[:space:]]*(", number,
                       ")[[:space:]]+(", number, ")")
  starts_number <- paste0("^[[:space:]]*\\(?[[:space:]]*", number)
  contains_missing_value <- paste0(
    "(?i)(?:^|[[:space:]:;,()])(?:NA|N/A|NaN|NULL|\\.)",
    "(?=$|[[:space:]:;,()])"
  )
  peaks <- character(0)
  problem_lines <- integer(0)
  invalid_value_lines <- integer(0)

  for (line_idx in seq_along(peak_lines)) {
    rest <- peak_lines[[line_idx]]
    parsed_on_line <- 0L

    repeat {
      rest <- sub("^[[:space:];,)]*", "", rest, perl = TRUE)
      if (!nchar(rest)) break

      match <- regexec(colon_pair, rest, perl = TRUE)
      hit <- regmatches(rest, match)[[1]]
      if (!length(hit)) {
        match <- regexec(space_pair, rest, perl = TRUE)
        hit <- regmatches(rest, match)[[1]]
      }

      if (!length(hit)) {
        # A numeric prefix without a complete pair is malformed even when it is
        # the first token on the physical line (for example, "200" or
        # "100 NA"). Explicit missing-value markers are reported as well.
        if (grepl(starts_number, rest, perl = TRUE) ||
            grepl(contains_missing_value, rest, perl = TRUE)) {
          problem_lines <- c(problem_lines, line_idx)
        }
        break
      }

      mz_token <- hit[[2]]
      intensity_token <- hit[[3]]
      mz <- suppressWarnings(as.numeric(mz_token))
      intensity <- suppressWarnings(as.numeric(intensity_token))
      if (is.finite(mz) && is.finite(intensity) && mz > 0 && intensity > 0) {
        peaks <- c(peaks, paste0(mz_token, ":", intensity_token))
      } else {
        invalid_value_lines <- c(invalid_value_lines, line_idx)
      }
      parsed_on_line <- parsed_on_line + 1L
      rest <- substr(rest, nchar(hit[[1]]) + 1L, nchar(rest))

      # Quoted or free-text annotations terminate peak parsing for this line.
      annotation_probe <- sub("^[[:space:];,)]*", "", rest, perl = TRUE)
      if (nchar(annotation_probe) &&
          !grepl(starts_number, annotation_probe, perl = TRUE)) break
    }
  }

  if (length(problem_lines)) {
    warning(
      "Unpaired numeric token(s) or missing value(s) in MSP peak line(s): ",
      paste(unique(problem_lines), collapse = ", "),
      call. = FALSE
    )
  }
  if (length(invalid_value_lines)) {
    warning(
      "Invalid non-positive or nonfinite value(s) in MSP peak line(s): ",
      paste(unique(invalid_value_lines), collapse = ", "),
      call. = FALSE
    )
  }
  paste(peaks, collapse = " ")
}

#' Convert List of Records to Tibble
#' @keywords internal
#' @noRd
msp_list_to_tibble <- function(records) {
  if (length(records) == 0) {
    return(tibble::tibble(
      name = character(),
      cas = character(),
      inchikey = character(),
      inchi = character(),
      smiles = character(),
      formula = character(),
      mw = numeric(),
      ri = numeric(),
      rt = numeric(),
      num_peaks = integer(),
      spectrum = character()
    ))
  }
  
  # Define all possible columns and types
  all_cols <- c(
    "name", "cas", "inchikey", "inchi", "smiles", "formula",
    "mw", "ri", "rt", "num_peaks", "spectrum", "synonyms",
    "comment", "instrument", "ion_mode", "precursor_type",
    "collision_energy", "source", "splash", "other_metadata"
  )
  col_types <- c(
    name = "character",
    cas = "character",
    inchikey = "character",
    inchi = "character",
    smiles = "character",
    formula = "character",
    mw = "numeric",
    ri = "numeric",
    rt = "numeric",
    num_peaks = "integer",
    spectrum = "character",
    synonyms = "character",
    comment = "character",
    instrument = "character",
    ion_mode = "character",
    precursor_type = "character",
    collision_energy = "character",
    source = "character",
    splash = "character",
    other_metadata = "character"
  )

  n <- length(records)
  result <- vector("list", length(all_cols))
  names(result) <- all_cols

  for (col in all_cols) {
    type <- col_types[[col]]
    if (type == "numeric") {
      result[[col]] <- rep(NA_real_, n)
    } else if (type == "integer") {
      result[[col]] <- rep(NA_integer_, n)
    } else {
      result[[col]] <- rep(NA_character_, n)
    }
  }

  for (i in seq_len(n)) {
    rec <- records[[i]]
    if (length(rec) == 0) next
    rec_names <- intersect(names(rec), all_cols)
    if (length(rec_names) == 0) next
    for (nm in rec_names) {
      val <- rec[[nm]]
      if (is.null(val)) next
      type <- col_types[[nm]]
      if (type == "numeric") {
        result[[nm]][i] <- suppressWarnings(as.numeric(val))
      } else if (type == "integer") {
        result[[nm]][i] <- suppressWarnings(as.integer(val))
      } else {
        result[[nm]][i] <- as.character(val)
      }
    }
  }

  df <- tibble::as_tibble(result)
  
  # Remove columns that are all NA
  all_na <- sapply(df, function(x) all(is.na(x)))
  df <- df[, !all_na, drop = FALSE]
  
  df
}

#' Convert MSP Library to Package Format
#'
#' Converts a tibble from read_msp() to the format used by build_spectra().
#'
#' @param msp_lib Tibble from read_msp().
#' @param id_column Column to use as compound ID (default: "name").
#' @param require_ri If TRUE, filter to entries with RI values.
#' @param require_inchikey If TRUE, filter to entries with InChIKey.
#' @return A tibble compatible with the package workflow.
#' @export
convert_msp_to_internal <- function(msp_lib, 
                                    id_column = "name",
                                    require_ri = FALSE,
                                    require_inchikey = FALSE) {
  
  # Validate

  if (!id_column %in% names(msp_lib)) {
    stop("id_column '", id_column, "' not found in MSP library")
  }
  
  if (!"spectrum" %in% names(msp_lib)) {
    stop("No 'spectrum' column found in MSP library")
  }
  
  # Filter
  df <- msp_lib
  
  if (require_ri) {
    if (!"ri" %in% names(df)) {
      stop("require_ri=TRUE but no 'ri' column in library")
    }
    df <- df[!is.na(df$ri), ]
  }
  
  if (require_inchikey) {
    if (!"inchikey" %in% names(df)) {
      stop("require_inchikey=TRUE but no 'inchikey' column in library")
    }
    df <- df[!is.na(df$inchikey) & nchar(df$inchikey) > 0, ]
  }
  
  # Remove duplicates (keep first occurrence)
  df <- df[!duplicated(df[[id_column]]), ]
  
  # Create output format matching build_df_spec output
  result <- tibble::tibble(
    id = df[[id_column]],
    known = "T",  # Library entries are "known"
    RI = if ("ri" %in% names(df)) df$ri else NA_real_,
    ei = df$spectrum
  )
  
  # Add optional metadata
  if ("name" %in% names(df)) result$name <- df$name
  if ("inchikey" %in% names(df)) result$inchikey <- df$inchikey
  if ("formula" %in% names(df)) result$formula <- df$formula
  if ("mw" %in% names(df)) result$mw <- df$mw
  if ("cas" %in% names(df)) result$cas <- df$cas
  if ("smiles" %in% names(df)) result$smiles <- df$smiles
  
  # Remove entries with empty spectra
  result <- result[nchar(result$ei) > 0 & !is.na(result$ei), ]
  
  message("Converted ", nrow(result), " spectra to internal format")
  
  result
}

#' Build Spectra from MSP Library
#'
#' Convenience function to process an MSP library through the full pipeline.
#'
#' @param msp_file Path to .msp file.
#' @param params Parameter list from eihrms_default_params().
#' @param require_ri If TRUE, filter to entries with RI values.
#' @param progress Show progress.
#' @return A list similar to build_spectra() output.
#' @export
build_spectra_from_msp <- function(msp_file, params, require_ri = FALSE, progress = TRUE) {
  
  # Read MSP
  message("Reading MSP file...")
  msp_lib <- read_msp(msp_file, progress = if (progress) 1000 else 0)
  
  # Convert to internal format
  message("Converting to internal format...")
  df_spec <- convert_msp_to_internal(msp_lib, require_ri = require_ri)
  
  # Process spectra
  message("Building fragment and loss spectra...")
  
  frag_list <- vector("list", nrow(df_spec))
  loss_list <- vector("list", nrow(df_spec))

  use_typical_loss <- isTRUE(params$use_typical_loss)
  use_split_loss <- isTRUE(params$use_split_loss)

  loss_typ_list <- if (use_typical_loss) vector("list", nrow(df_spec)) else NULL
  loss_anchor_list <- if (use_split_loss) vector("list", nrow(df_spec)) else NULL
  loss_pair_list <- if (use_split_loss) vector("list", nrow(df_spec)) else NULL

  loss_anchor_typ_list <- if (use_typical_loss && use_split_loss) vector("list", nrow(df_spec)) else NULL
  loss_pair_typ_list <- if (use_typical_loss && use_split_loss) vector("list", nrow(df_spec)) else NULL

  mref_vec <- numeric(nrow(df_spec))
  mref_conf_vec <- numeric(nrow(df_spec))
  neutral_loss_vec <- character(nrow(df_spec))
  functional_group_vec <- character(nrow(df_spec))

  derivatization_type_vec <- character(nrow(df_spec))
  deriv_score_tms_vec <- numeric(nrow(df_spec))
  deriv_score_tbdms_vec <- numeric(nrow(df_spec))
  
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
  
  # Compound class detection
  df_spec$compound_class <- vapply(seq_len(nrow(df_spec)), function(i) {
    ion_class <- detect_compound_class(frag_list[[i]], ppm_tol = params$class_detection_ppm, rules = params$classify_rules)
    func_groups <- functional_group_vec[i]
    
    if (ion_class == "unclassified" && nchar(func_groups) > 0) {
      paste0("NL:", func_groups)
    } else {
      ion_class
    }
  }, character(1))
  df_spec$compound_class_source <- "inferred"
  
  ri <- df_spec$RI
  names(ri) <- df_spec$id
  
  list(
    df_spec = df_spec,
    frag_list = frag_list,
    # `loss_*` names are retained as backward-compatible aliases. The accurate
    # terminology is pooled/anchored/pairwise derived mass differences.
    loss_list = loss_list,
    derived_list = loss_list,
    loss_typ_list = loss_typ_list,
    loss_anchor_list = loss_anchor_list,
    loss_pair_list = loss_pair_list,
    derived_anchor_list = loss_anchor_list,
    derived_pair_list = loss_pair_list,
    loss_anchor_typ_list = loss_anchor_typ_list,
    loss_pair_typ_list = loss_pair_typ_list,
    mref_confidence = mref_conf_vec,
    ri = ri,
    source = "msp",
    msp_metadata = msp_lib
  )
}

#' Write Spectra to MSP Format
#'
#' Export spectra to MSP format for sharing or use in other software.
#'
#' @param df_spec Data frame with spectral information.
#' @param frag_list Named list of fragment spectra matrices.
#' @param file Output file path.
#' @param include_losses If TRUE, also write loss spectra.
#' @export
write_msp <- function(df_spec, frag_list, file, include_losses = FALSE) {
  
  con <- file(file, "w")
  on.exit(close(con))

  for (i in seq_len(nrow(df_spec))) {
    id <- df_spec$id[i]
    spec <- frag_list[[id]]
    
    if (is.null(spec) || nrow(spec) == 0) next
    
    # Write metadata
    writeLines(paste0("NAME: ", id), con)
    
    if ("ri" %in% names(df_spec) || "RI" %in% names(df_spec)) {
      ri_val <- if ("RI" %in% names(df_spec)) df_spec$RI[i] else df_spec$ri[i]
      if (!is.na(ri_val)) {
        writeLines(paste0("RI: ", round(ri_val)), con)
      }
    }
    
    if ("inchikey" %in% names(df_spec) && !is.na(df_spec$inchikey[i])) {
      writeLines(paste0("INCHIKEY: ", df_spec$inchikey[i]), con)
    }
    
    if ("formula" %in% names(df_spec) && !is.na(df_spec$formula[i])) {
      writeLines(paste0("FORMULA: ", df_spec$formula[i]), con)
    }
    
    if ("compound_class" %in% names(df_spec)) {
      writeLines(paste0("COMMENT: CompoundClass=", df_spec$compound_class[i]), con)
    }
    
    writeLines(paste0("Num Peaks: ", nrow(spec)), con)
    
    # Write peaks (mz intensity format), scale to max = 999
    max_int <- max(spec[, 2])
    if (max_int > 0) {
      scaled_int <- spec[, 2] / max_int * 999
    } else {
      scaled_int <- spec[, 2]
    }
    for (j in seq_len(nrow(spec))) {
      writeLines(sprintf("%.4f %.0f", spec[j, 1], scaled_int[j]), con)
    }
    
    writeLines("", con)  # Empty line between records
  }
  
  message("Wrote ", nrow(df_spec), " spectra to ", file)
}
