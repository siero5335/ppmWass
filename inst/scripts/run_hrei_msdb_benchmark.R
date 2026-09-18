#!/usr/bin/env Rscript

# Compatibility stub retained so historical command paths fail explicitly.
# The complete pre-audit implementation remains available in Git history and
# in the publication rerun bundle's immutable archive.

stop(
  paste(
    "PPMWASS_LEGACY_DRIVER_DISABLED: run_hrei_msdb_benchmark.R is a",
    "superseded pre-audit driver and is intentionally disabled because it",
    "formerly patched compute_distance_matrix() to replace nonfinite cells",
    "by the maximum finite value. Use run_main_retrieval_publication.R with",
    "--datasets=HREI-MSDB."
  ),
  call. = FALSE
)
