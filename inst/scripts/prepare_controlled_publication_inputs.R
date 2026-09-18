#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1]])
}

repo_dir <- normalizePath(get_arg("repo", getwd()), mustWork = TRUE)
recetox_msp <- normalizePath(get_arg("recetox-msp"), mustWork = TRUE)
nominal_msp <- normalizePath(get_arg("nominal-msp"), mustWork = TRUE)
mapping_csv <- normalizePath(get_arg("mapping-csv"), mustWork = TRUE)
output_dir <- get_arg("output-dir", file.path(dirname(repo_dir), "controlled_prepared"))
max_per_prefix <- as.integer(get_arg("max-per-prefix", "10"))
seed <- as.integer(get_arg("seed", "42"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
pkgload::load_all(repo_dir, quiet = TRUE)

make_params <- function(tol) {
  p <- eihrms_default_params()
  p$min_mz <- 35
  p$max_mz <- 650
  p$noise_thr <- 0.01
  p$topK <- 200L
  p$class_detection_ppm <- tol
  p$tol_ppm <- tol
  p$ot_method <- "exact"
  p$sinkhorn_epsilon <- 0.05
  p$sinkhorn_niter <- 100L
  p$wasserstein_transition_mult <- 3
  p$use_parallel <- TRUE
  p$n_cores <- 9L
  p
}

subset_spectra <- function(x, idx) {
  out <- list(
    df_spec = x$df_spec[idx, , drop = FALSE],
    frag_list = x$frag_list[idx],
    loss_list = x$loss_list[idx],
    derived_list = x$loss_list[idx],
    ri = x$ri[idx]
  )
  for (nm in c(
    "loss_typ_list", "mref_confidence", "loss_anchor_list", "loss_pair_list",
    "derived_anchor_list", "derived_pair_list", "loss_anchor_typ_list",
    "loss_pair_typ_list"
  )) {
    if (!is.null(x[[nm]])) out[[nm]] <- x[[nm]][idx]
  }
  out
}

filter_valid <- function(x) {
  valid <- vapply(x$frag_list, function(f) {
    is.matrix(f) && nrow(f) >= 2L && ncol(f) >= 2L &&
      all(is.finite(f)) && any(f[, 2] > 0)
  }, logical(1))
  subset_spectra(x, which(valid))
}

attach_source_metadata <- function(subset, full) {
  subset$raw_msp_metadata <- full$msp_metadata
  subset$msp_parser_diagnostics <- attr(
    full$msp_metadata, "msp_parser_diagnostics"
  )
  names(subset$frag_list) <- subset$df_spec$id
  names(subset$loss_list) <- subset$df_spec$id
  names(subset$derived_list) <- subset$df_spec$id
  subset
}

message("Preparing controlled RECETOX subset")
rec_full <- build_spectra_from_msp(
  recetox_msp, make_params(15), require_ri = FALSE, progress = TRUE
)
mapping <- utils::read.csv(mapping_csv, stringsAsFactors = FALSE)
mapping_valid <- mapping[!is.na(mapping$match_ik14) & mapping$match_ik14 != "", ]
mapped_idx <- match(mapping_valid$id, rec_full$df_spec$id)
if (anyNA(mapped_idx)) {
  warning("Mapped RECETOX IDs absent after parser correction: ",
          paste(mapping_valid$id[is.na(mapped_idx)], collapse = ", "))
}
mapped_idx <- mapped_idx[!is.na(mapped_idx)]

# The malformed-record repair recovers records that did not exist when the old
# index mapping was created. Include a recovered record if its own connectivity
# prefix is already in the controlled target set.
new_idx <- setdiff(seq_len(nrow(rec_full$df_spec)), match(mapping$id, rec_full$df_spec$id))
new_idx <- new_idx[!is.na(new_idx)]
new_prefix <- substr(rec_full$df_spec$inchikey[new_idx], 1, 14)
recovered_idx <- new_idx[
  !is.na(new_prefix) & new_prefix %in% unique(mapping_valid$match_ik14)
]
rec_idx <- sort(unique(c(mapped_idx, recovered_idx)))
rec_subset <- filter_valid(subset_spectra(rec_full, rec_idx))
rec_subset <- attach_source_metadata(rec_subset, rec_full)
rec_subset$controlled_selection <- list(
  source_mapping = normalizePath(mapping_csv),
  mapped_ids = rec_full$df_spec$id[mapped_idx],
  parser_recovered_ids = rec_full$df_spec$id[recovered_idx]
)

mapping_corrected <- mapping_valid[!is.na(match(mapping_valid$id, rec_full$df_spec$id)), ]
if (length(recovered_idx)) {
  recovered_rows <- data.frame(
    idx = recovered_idx,
    id = rec_full$df_spec$id[recovered_idx],
    orig_inchikey = rec_full$df_spec$inchikey[recovered_idx],
    match_ik14 = substr(rec_full$df_spec$inchikey[recovered_idx], 1, 14),
    match_source = "parser_recovered_same_prefix",
    stringsAsFactors = FALSE
  )
  mapping_corrected <- rbind(mapping_corrected, recovered_rows)
}
utils::write.csv(
  mapping_corrected,
  file.path(output_dir, "recetox_mapping_v3_parser_corrected.csv"),
  row.names = FALSE
)
saveRDS(rec_subset, file.path(output_dir, "controlled_recetox_spectra.rds"))

message("Preparing controlled NIST 2020 + Wiley-derived nominal subset")
nominal_full <- build_spectra_from_msp(
  nominal_msp, make_params(7000), require_ri = FALSE, progress = TRUE
)
df <- nominal_full$df_spec
valid_idx <- which(
  !is.na(df$inchikey) & df$inchikey != "" &
    nchar(df$inchikey) >= 14L & nchar(df$inchikey) <= 27L
)
nominal_subset <- filter_valid(subset_spectra(nominal_full, valid_idx))
prefix <- substr(nominal_subset$df_spec$inchikey, 1, 14)
set.seed(seed)
keep <- integer(0)
for (key in unique(prefix)) {
  idx <- which(prefix == key)
  if (length(idx) > max_per_prefix) idx <- sample(idx, max_per_prefix)
  keep <- c(keep, idx)
}
keep <- sort(keep)
nominal_subset <- subset_spectra(nominal_subset, keep)
nominal_subset <- attach_source_metadata(nominal_subset, nominal_full)
nominal_subset$controlled_selection <- list(
  source_label = "NIST 2020 + Wiley-derived nominal comparison set",
  max_per_connectivity_prefix = max_per_prefix,
  subsampling_seed = seed
)
saveRDS(nominal_subset, file.path(output_dir, "controlled_nominal_spectra.rds"))

qc <- rbind(
  data.frame(
    dataset = "CONTROLLED_RECETOX_HRMS",
    raw_records = nrow(rec_full$msp_metadata),
    converted_unique_names = nrow(rec_full$df_spec),
    final_spectra = nrow(rec_subset$df_spec),
    unique_full_inchikey = length(unique(rec_subset$df_spec$inchikey)),
    unique_connectivity_prefix = length(unique(substr(rec_subset$df_spec$inchikey, 1, 14))),
    recovered_parser_records = length(recovered_idx),
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = "NIST2020_WILEY_NOMINAL",
    raw_records = nrow(nominal_full$msp_metadata),
    converted_unique_names = nrow(nominal_full$df_spec),
    final_spectra = nrow(nominal_subset$df_spec),
    unique_full_inchikey = length(unique(nominal_subset$df_spec$inchikey)),
    unique_connectivity_prefix = length(unique(substr(nominal_subset$df_spec$inchikey, 1, 14))),
    recovered_parser_records = 0L,
    stringsAsFactors = FALSE
  )
)
utils::write.csv(qc, file.path(output_dir, "controlled_preparation_qc.csv"), row.names = FALSE)
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
    parameter = c(
      "timestamp_utc", "package_commit", "package_tree_dirty", "command",
      "recetox_msp", "nominal_msp", "mapping_csv",
      "controlled_hrms_hard_match_and_base_cost_ppm",
      "nominal_hard_match_and_base_cost_ppm", "transition_multiplier",
      "ot_method", "ot_estimand", "sinkhorn_epsilon", "sinkhorn_iterations",
      "nominal_max_per_prefix",
      "nominal_subsampling_seed"
    ),
    value = c(
      format(Sys.time(), tz = "UTC"), git_commit, length(git_status) > 0L,
      paste(commandArgs(), collapse = " "), recetox_msp, nominal_msp,
      mapping_csv, 15, 7000, 3, "exact",
      "unregularized_exact_transport_cost", NA, NA, max_per_prefix, seed
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "controlled_preparation_parameters.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    file = c(recetox_msp, nominal_msp, mapping_csv),
    md5 = unname(tools::md5sum(c(recetox_msp, nominal_msp, mapping_csv))),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "controlled_preparation_input_checksums.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "controlled_preparation_sessionInfo.txt"))
print(qc)
message("Controlled inputs prepared: ", output_dir)
