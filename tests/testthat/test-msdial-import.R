write_msdial_fixture <- function(path) {
  cols <- c(
    "Alignment ID",
    "Average Rt(min)",
    "Average RI",
    "Metabolite name",
    "Annotation tag (VS1.0)",
    "Total score",
    "RT similarity",
    "INCHIKEY",
    "EI spectrum",
    "BLK_01",
    "QC_01",
    "STD_01",
    "S01",
    "S02",
    "Average",
    "Stdev",
    "1",
    "1"
  )

  make_meta_row <- function(label, values) {
    out <- rep("", length(cols))
    names(out) <- cols
    out[names(values)] <- values
    paste(c(label, out), collapse = "\t")
  }

  lines <- c(
    make_meta_row("Class", c(BLK_01 = "Blank", QC_01 = "QC", STD_01 = "STD", S01 = "Sample", S02 = "Sample")),
    make_meta_row("File type", c(BLK_01 = "Blank", QC_01 = "QC", STD_01 = "Standard", S01 = "Sample", S02 = "Sample")),
    make_meta_row("Injection order", c(BLK_01 = "1", QC_01 = "2", STD_01 = "3", S01 = "4", S02 = "5")),
    make_meta_row("Batch ID", c(BLK_01 = "B1", QC_01 = "B1", STD_01 = "B1", S01 = "B1", S02 = "B1")),
    paste(cols, collapse = "\t"),
    paste(c("1", "5.00", "600", "Compound_A", "1", "95", "99", "KEY1", "50:100 73.0474:300 147.0654:120", "1", "25", "15", "20", "22", "21", "1.4", "21", "1.4"), collapse = "\t"),
    paste(c("2", "6.00", "650", "Unknown", "4", "80", "90", "", "60:120 73.0474:150 115.0943:90", "0.5", "8", "4", "5", "6", "5.5", "0.7", "5.5", "0.7"), collapse = "\t"),
    paste(c("3", "7.20", "710", "Compound_B", "1", "88", "92", "KEY2", "57.0704:240 73.0474:180 115.0943:120", "0.7", "19", "12", "16", "18", "17", "1.1", "17", "1.1"), collapse = "\t"),
    paste(c("4", "9.10", "780", "Unknown", "4", "75", "85", "", "43.0180:210 55.0548:160 91.0543:90", "0.2", "7", "3", "9", "11", "10", "1.4", "10", "1.4"), collapse = "\t"),
    paste(c("5", "10.40", "820", "Compound_C", "1", "90", "94", "KEY3", "59.0491:220 73.0474:140 147.0654:100", "0.4", "21", "11", "14", "15", "14.5", "0.7", "14.5", "0.7"), collapse = "\t"),
    paste(c("6", "11.60", "860", "Unknown", "4", "70", "82", "", "45.0336:190 77.0387:130 91.0543:110", "0.3", "9", "2", "7", "8", "7.5", "0.7", "7.5", "0.7"), collapse = "\t")
  )

  writeLines(lines, path)
}

test_that("read_ms_dial preserves raw header and sample metadata", {
  tmp <- tempfile(fileext = ".txt")
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  write_msdial_fixture(tmp)

  imported <- ppmWass::read_ms_dial(tmp, show_col_types = FALSE)

  expect_s3_class(imported, "msdial_import")
  expect_equal(imported$raw_header, readr::read_lines(tmp, n_max = 4))
  expect_true(all(c("data", "sample_meta", "raw_header", "source") %in% names(imported)))

  sample_meta <- imported$sample_meta
  expect_true(all(c("column", "column_raw", "index", "class", "file_type", "injection_order", "batch_id", "role") %in% names(sample_meta)))
  expect_equal(sample_meta$role[sample_meta$column == "BLK_01"], "blank")
  expect_equal(sample_meta$role[sample_meta$column == "QC_01"], "qc")
  expect_equal(sample_meta$role[sample_meta$column == "STD_01"], "standard")
  expect_equal(sample_meta$role[sample_meta$column == "S01"], "sample")
  expect_equal(sample_meta$role[sample_meta$column == "Average"], "summary")

  repaired_numeric <- sample_meta[sample_meta$column_raw == "1", ]
  expect_equal(nrow(repaired_numeric), 2)
  expect_true(all(repaired_numeric$role == "summary"))
})

test_that("metadata-driven column helpers classify roles and keep QC explicit", {
  tmp <- tempfile(fileext = ".txt")
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  write_msdial_fixture(tmp)

  imported <- ppmWass::read_ms_dial(tmp, show_col_types = FALSE)

  expect_equal(ppmWass:::get_sample_cols(imported), c("S01", "S02"))
  expect_equal(ppmWass:::get_blank_cols(imported), "BLK_01")
  expect_equal(ppmWass:::get_qc_cols(imported), "QC_01")
  expect_equal(ppmWass:::get_std_cols(imported), "STD_01")
  expect_equal(ppmWass:::get_signal_cols(imported), c("S01", "S02"))
  expect_equal(ppmWass:::get_signal_cols(imported, include_qc = TRUE), c("QC_01", "S01", "S02"))
})

test_that("plain fallback detects suffix-style STD names without matching Stdev", {
  cols <- c("20230208STD", "STD_01", "Stdev", "S01")
  expect_equal(ppmWass:::get_std_cols(cols), c("20230208STD", "STD_01"))
})

test_that("process_single_spectrum returns projected channels and derivatization summary", {
  params <- ppmWass::eihrms_default_params()
  params$use_typical_loss <- TRUE
  params$use_split_loss <- TRUE

  res <- ppmWass:::process_single_spectrum(
    "50:100 73.0474:300 147.0654:120 200:80",
    params
  )

  expect_true(is.matrix(res$frag))
  expect_true(is.matrix(res$loss))
  expect_true(is.matrix(res$loss_anchor))
  expect_true(is.matrix(res$loss_pair))
  expect_true(is.matrix(res$loss_typ))
  expect_true(is.matrix(res$loss_anchor_typ))
  expect_true(is.matrix(res$loss_pair_typ))
  expect_true(res$deriv_type %in% c("none", "TMS", "TBDMS", "BOTH"))
  expect_true(is.character(res$neutral_losses))
  expect_true(is.character(res$functional_groups))
})

test_that("resolve_named_options merges valid overrides and rejects unknown names", {
  defaults <- list(alpha = 1, beta = 2)

  merged <- ppmWass:::resolve_named_options(
    defaults,
    list(beta = 5),
    arg_name = "demo_opts"
  )

  expect_equal(merged$alpha, 1)
  expect_equal(merged$beta, 5)
  expect_error(
    ppmWass:::resolve_named_options(defaults, list(gamma = 3), arg_name = "demo_opts"),
    "demo_opts contains unknown names"
  )
  expect_error(
    ppmWass:::resolve_named_options(defaults, "bad", arg_name = "demo_opts"),
    "demo_opts must be a list"
  )
})

test_that("prepare_quant_data and spectra builders accept msdial_import", {
  tmp <- tempfile(fileext = ".txt")
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  write_msdial_fixture(tmp)

  imported <- ppmWass::read_ms_dial(tmp, show_col_types = FALSE)
  params <- ppmWass::eihrms_default_params()

  quant <- ppmWass::prepare_quant_data(
    imported,
    params,
    blank_factor = 2,
    nonzero_ratio = 0,
    write_cleaned_csv = FALSE
  )

  expect_equal(quant$sample_cols, c("S01", "S02"))
  expect_equal(quant$blank_cols, "BLK_01")
  expect_equal(quant$qc_cols, "QC_01")
  expect_equal(quant$std_cols, "STD_01")
  expect_equal(quant$signal_cols, c("S01", "S02"))
  expect_false(any(c("BLK_01", "STD_01", "Average", "Stdev") %in% names(quant$final)))

  df_spec <- ppmWass:::build_df_spec(imported, quant$final)
  expect_true(all(c("Compound_A", "Unknown_2") %in% df_spec$id))

  spectra <- ppmWass::build_spectra(imported, quant$final, params, progress = FALSE)
  expect_equal(sort(names(spectra$frag_list)), sort(df_spec$id))
})

test_that("prepare_quant_data keeps rows when all blank values are NA", {
  df <- tibble::tibble(
    `Alignment ID` = 1,
    `Average Rt(min)` = 5,
    `Average RI` = 600,
    `Metabolite name` = "Compound_A",
    `Annotation tag (VS1.0)` = "1",
    `Total score` = 95,
    `RT similarity` = 99,
    INCHIKEY = "KEY1",
    `EI spectrum` = "50:100 73:200",
    BLK_01 = NA_real_,
    S01 = 10,
    S02 = 12
  )

  params <- ppmWass::eihrms_default_params()
  quant <- ppmWass::prepare_quant_data(
    df,
    params,
    blank_factor = 2,
    nonzero_ratio = 0,
    write_cleaned_csv = FALSE
  )

  expect_equal(nrow(quant$final), 1)
  expect_equal(quant$cleaned$blank_mean, 0)
  expect_equal(quant$cleaned$S01, 10)
})

test_that("plain tibble fallback still works without metadata object", {
  tmp <- tempfile(fileext = ".txt")
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  write_msdial_fixture(tmp)

  plain <- readr::read_tsv(tmp, skip = 4, show_col_types = FALSE, name_repair = "unique_quiet")
  params <- ppmWass::eihrms_default_params()

  expect_equal(ppmWass:::get_sample_cols(plain), c("S01", "S02"))
  expect_equal(ppmWass:::get_blank_cols(plain), "BLK_01")
  expect_equal(ppmWass:::get_qc_cols(plain), "QC_01")
  expect_equal(ppmWass:::get_std_cols(plain), "STD_01")
  expect_equal(ppmWass:::get_signal_cols(plain), c("S01", "S02"))

  quant <- ppmWass::prepare_quant_data(
    plain,
    params,
    blank_factor = 2,
    nonzero_ratio = 0,
    write_cleaned_csv = FALSE
  )

  expect_true(nrow(quant$final) >= 1)
})

test_that("run_eihrms_similarity accepts msdial_import when optional deps are present", {
  skip_if_not_installed("dbscan")
  skip_if_not_installed("uwot")

  tmp <- tempfile(fileext = ".txt")
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)
  write_msdial_fixture(tmp)

  imported <- ppmWass::read_ms_dial(tmp, show_col_types = FALSE)
  params <- ppmWass::eihrms_default_params()
  params$distance_method <- "hellinger"

  res <- ppmWass::run_eihrms_similarity(
    df = imported,
    params = params,
    write_outputs = FALSE,
    progress = FALSE,
    heatmap_opts = list(threshold = 0.4, width = 9),
    umap_opts = list(seed = 42, point_size = 3),
    cluster_report_opts = list(top_pairs = 2, plot_width = 8),
    write_cluster_annotations = FALSE,
    write_cluster_reports = FALSE
  )

  expect_true(is.list(res))
  expect_true(is.list(res$quant))
  expect_true(is.list(res$spectra))
})
