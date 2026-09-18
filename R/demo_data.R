#' Generate Synthetic MS-DIAL-like Demo Data
#'
#' @param n_known Number of known compounds.
#' @param n_unknown Number of unknown compounds.
#' @param n_samples Number of sample columns.
#' @param n_blanks Number of blank columns.
#' @param n_standards Number of standard columns.
#' @param seed Random seed for reproducibility.
#' @return A tibble resembling an MS-DIAL export table.
#' @export
generate_demo_msdial <- function(
  n_known = 5,
  n_unknown = 5,
  n_samples = 6,
  n_blanks = 1,
  n_standards = 1,
  seed = 1
) {
  if (n_known < 1 || n_unknown < 0) stop("n_known must be >= 1 and n_unknown >= 0")
  if (n_samples < 1) stop("n_samples must be >= 1")

  with_seed(seed, {
    n <- n_known + n_unknown
    alignment_id <- seq_len(n)

  meta <- tibble::tibble(
    `Alignment ID` = alignment_id,
    `Average Rt(min)` = round(stats::runif(n, 2, 20), 3),
    `Average RI` = round(stats::runif(n, 400, 1200), 2),
    `Quant mass` = round(stats::runif(n, 40, 300), 5),
    `Metabolite name` = c(paste0("Compound_", seq_len(n_known)), rep("Unknown", n_unknown)),
    `Fill %` = round(stats::runif(n, 0.7, 1.0), 3),
    `Reference RT` = NA,
    `Reference RI` = NA,
    `Formula` = NA,
    `Ontology` = NA,
    `INCHIKEY` = c(paste0("KEY", seq_len(n_known)), rep(NA_character_, n_unknown)),
    `SMILES` = NA,
    `Annotation tag (VS1.0)` = c(rep("1", n_known), rep("4", n_unknown)),
    `RT/RI matched` = rep("False", n),
    `EI-MS matched` = rep("False", n),
    `Comment` = NA,
    `Manually modified for quantification` = rep("FALSE", n),
    `Manually modified for annotation` = rep("FALSE", n),
    `Total score` = NA,
    `RT similarity` = NA,
    `RI similarity` = NA,
    `Total spectrum similarity` = NA,
    `Dot product` = NA,
    `Reverse dot product` = NA,
    `Fragment presence %` = NA,
    `S/N average` = round(stats::runif(n, 50, 500), 2),
    `Spectrum reference file name` = rep("demo", n),
    `EI spectrum` = NA_character_
  )

  generate_spectrum <- function() {
    n_peaks <- sample(15:30, 1)
    mz <- sort(stats::runif(n_peaks, 35, 300))
    ints <- round(stats::runif(n_peaks, 100, 50000))
    paste(sprintf("%.5f:%d", mz, ints), collapse = " ")
  }

  meta$`EI spectrum` <- vapply(seq_len(n), function(i) generate_spectrum(), character(1))

  std_cols <- if (n_standards > 0) paste0("STD_", sprintf("%02d", seq_len(n_standards))) else character(0)
  blank_cols <- if (n_blanks > 0) paste0("BLK_", sprintf("%02d", seq_len(n_blanks))) else character(0)
  sample_cols <- paste0("S", sprintf("%02d", seq_len(n_samples)))

  base_abundance <- stats::runif(n, 0.2, 1.2)

  make_matrix <- function(cols, factor_range) {
    if (length(cols) == 0) return(NULL)
    mat <- vapply(cols, function(.x) base_abundance * stats::runif(n, factor_range[1], factor_range[2]), numeric(n))
    mat <- as.data.frame(mat)
    names(mat) <- cols
    mat
  }

  std_mat <- make_matrix(std_cols, c(0.8, 1.2))
  blk_mat <- make_matrix(blank_cols, c(0.0, 0.15))
  samp_mat <- make_matrix(sample_cols, c(0.6, 1.6))

  out <- dplyr::bind_cols(meta, std_mat, blk_mat, samp_mat)

  out$Average <- round(rowMeans(dplyr::select(out, dplyr::all_of(sample_cols))), 4)
  out$Stdev <- round(apply(dplyr::select(out, dplyr::all_of(sample_cols)), 1, stats::sd), 4)

    out
  })
}
