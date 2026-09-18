#' @keywords internal
#' @noRd
ms_dial_filter <- function(df, column_source = df, deduplicate = TRUE) {
  df_tbl <- as_feature_table(df)
  quant_cols <- intersect(get_quant_cols(column_source), names(df_tbl))
  signal_cols <- intersect(get_signal_cols(column_source, include_qc = FALSE), names(df_tbl))
  if (length(signal_cols) == 0) signal_cols <- quant_cols

  df_pre <- df_tbl %>% dplyr::select(`Alignment ID`, `Metabolite name`, `Annotation tag (VS1.0)`, `Total score`,
                       `RT similarity`, `Average Rt(min)`, INCHIKEY, dplyr::all_of(quant_cols)) %>%
    dplyr::mutate(mean = rowMeans(dplyr::select(., dplyr::all_of(signal_cols)), na.rm = TRUE)) %>%
    dplyr::filter(mean > 0) %>%
    dplyr::arrange(dplyr::desc(mean))

  if (!isTRUE(deduplicate)) {
    return(df_pre)
  }

  df_pre <- df_pre %>%
    dplyr::distinct(`Metabolite name`, .keep_all = TRUE)

  df_non_na <- df_pre %>%
    dplyr::filter(!is.na(INCHIKEY) & INCHIKEY != "") %>%
    dplyr::distinct(INCHIKEY, .keep_all = TRUE)

  df_na <- df_pre %>%
    dplyr::filter(is.na(INCHIKEY) | INCHIKEY == "")

  dplyr::bind_rows(df_non_na, df_na)
}

#' @keywords internal
#' @noRd
ms_dial_filter_unknown <- function(df, column_source = df, deduplicate = FALSE) {
  df_tbl <- as_feature_table(df)
  quant_cols <- intersect(get_quant_cols(column_source), names(df_tbl))
  signal_cols <- intersect(get_signal_cols(column_source, include_qc = FALSE), names(df_tbl))
  if (length(signal_cols) == 0) signal_cols <- quant_cols

  df_out <- df_tbl %>% dplyr::select(`Alignment ID`, `Metabolite name`, `Annotation tag (VS1.0)`, `Total score`,
                       `RT similarity`, `Average Rt(min)`, INCHIKEY, dplyr::all_of(quant_cols)) %>%
    dplyr::mutate(mean = rowMeans(dplyr::select(., dplyr::all_of(signal_cols)), na.rm = TRUE)) %>%
    dplyr::filter(mean > 0)

  if (isTRUE(deduplicate)) {
    df_out <- df_out %>%
      dplyr::arrange(dplyr::desc(mean)) %>%
      dplyr::distinct(`Metabolite name`, `Alignment ID`, .keep_all = TRUE)
  }

  df_out
}

#' Prepare Quantification Data
#'
#' @param df MS-DIAL tibble.
#' @param params Parameter list from eihrms_default_params().
#' @param blank_factor Multiplier for blank threshold.
#' @param nonzero_ratio Minimum nonzero ratio for filtering.
#' @param write_cleaned_csv If TRUE, write cleaned CSV.
#' @param cleaned_csv_path Output path for cleaned CSV.
#' @param export_unknown_as_id If TRUE, write unknowns using compound_id in cleaned CSV.
#' @param deduplicate If TRUE, keep the current prepare-stage deduplication behavior for known compounds.
#' @return A list with normalized, cleaned, final tables and column metadata.
#' @export
prepare_quant_data <- function(
  df,
  params,
  blank_factor = 3,
  nonzero_ratio = 0.2,
  write_cleaned_csv = FALSE,
  cleaned_csv_path = "cleaned_voc.csv",
  export_unknown_as_id = TRUE,
  deduplicate = TRUE
) {
  df_tbl <- as_feature_table(df)
  quant_cols <- intersect(get_quant_cols(df), names(df_tbl))
  sample_cols <- intersect(get_sample_cols(df), names(df_tbl))
  blank_cols <- intersect(get_blank_cols(df), names(df_tbl))
  qc_cols <- intersect(get_qc_cols(df), names(df_tbl))
  std_cols <- intersect(get_std_cols(df), names(df_tbl))
  signal_cols <- intersect(get_signal_cols(df, include_qc = FALSE), names(df_tbl))
  if (length(signal_cols) == 0) {
    signal_cols <- sample_cols
  }
  if (length(signal_cols) == 0) {
    signal_cols <- intersect(get_signal_cols(df, include_qc = TRUE), names(df_tbl))
  }
  if (length(signal_cols) == 0) {
    signal_cols <- quant_cols
  }
  if (length(sample_cols) == 0) stop("No sample columns detected. Check input format.")

  # MS-DIAL convention: Annotation tag "4" is treated as unknown.
  known <- df_tbl %>%
    dplyr::filter(`Annotation tag (VS1.0)` != "4") %>%
    ms_dial_filter(column_source = df, deduplicate = deduplicate) %>%
    dplyr::mutate(known = "T")

  unknown <- df_tbl %>%
    dplyr::filter(`Annotation tag (VS1.0)` == "4") %>%
    ms_dial_filter_unknown(column_source = df, deduplicate = FALSE) %>%
    dplyr::mutate(known = "F")

  df2 <- dplyr::bind_rows(known, unknown) %>%
    dplyr::mutate(compound_id = derive_compound_id(.))

  df2 <- df2 %>% dplyr::mutate(dplyr::across(dplyr::all_of(quant_cols), as.numeric))

  normalize_cols <- quant_cols
  if (!is.null(params$normalize_cols)) {
    missing_cols <- setdiff(params$normalize_cols, quant_cols)
    if (length(missing_cols) > 0) {
      stop("normalize_cols not found: ", paste(missing_cols, collapse = ", "))
    }
    normalize_cols <- params$normalize_cols
  }

  if (!is.null(params$normalize_ref)) {
    ref_row <- resolve_reference_row(df2, params$normalize_ref)
    ref_vals <- as.numeric(unlist(df2[ref_row, normalize_cols], use.names = FALSE))
    ok <- is.finite(ref_vals) & ref_vals != 0
    if (!any(ok)) {
      stop("normalize_ref has only NA/0 values in selected columns.")
    }
    if (!all(ok)) {
      warning("normalize_ref has NA/0 values in some columns; those columns are skipped.")
    }
    df2[, normalize_cols[ok]] <- sweep(
      as.matrix(df2[, normalize_cols[ok]]),
      2,
      ref_vals[ok],
      "/"
    )
  }

  df3 <- df2 %>%
    dplyr::mutate(blank_mean = if (length(blank_cols) > 0) {
      blank_mean <- rowMeans(dplyr::select(., dplyr::all_of(blank_cols)), na.rm = TRUE)
      blank_mean[is.nan(blank_mean)] <- 0
      blank_mean
    } else 0) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(quant_cols),
                                ~ dplyr::case_when(. > blank_mean * blank_factor ~ .,
                                                   . <= blank_mean * blank_factor ~ 0,
                                                   TRUE ~ .))) %>%
    dplyr::mutate(mean = rowMeans(dplyr::select(., dplyr::all_of(signal_cols)), na.rm = TRUE)) %>%
    dplyr::filter(mean > 0) %>%
    dplyr::select(!dplyr::any_of(c(blank_cols, std_cols))) %>%
    dplyr::select(!`RT similarity`) %>%
    dplyr::mutate(nonzero = rowSums(dplyr::select(., dplyr::all_of(signal_cols)) > 0)) %>%
    dplyr::filter(nonzero > length(signal_cols) * nonzero_ratio)

  df4 <- df3 %>%
    dplyr::select(!c(`Annotation tag (VS1.0)`, `Total score`, `Average Rt(min)`, INCHIKEY, nonzero))

  if (isTRUE(write_cleaned_csv)) {
    df3_export <- df3
    if (isTRUE(export_unknown_as_id) && "compound_id" %in% names(df3_export)) {
      df3_export <- df3_export %>%
        dplyr::mutate(`Metabolite name` = dplyr::if_else(
          `Annotation tag (VS1.0)` == "4",
          compound_id,
          `Metabolite name`
        ))
    }
    readr::write_csv(df3_export, cleaned_csv_path)
  }

  list(
    normalized = df2,
    cleaned = df3,
    final = df4,
    sample_cols = sample_cols,
    quant_cols = quant_cols,
    signal_cols = signal_cols,
    blank_cols = blank_cols,
    qc_cols = qc_cols,
    std_cols = std_cols
  )
}
