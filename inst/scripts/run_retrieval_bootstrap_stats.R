#!/usr/bin/env Rscript

# Compatibility stub retained so historical command paths fail explicitly.

stop(
  paste(
    "PPMWASS_LEGACY_DRIVER_DISABLED: run_retrieval_bootstrap_stats.R applies",
    "superseded optimistic-tie and query-level bootstrap rules to legacy",
    "outputs. Use run_main_retrieval_publication.R, which reports exact",
    "fractional ties and identity-cluster bootstrap inference."
  ),
  call. = FALSE
)
