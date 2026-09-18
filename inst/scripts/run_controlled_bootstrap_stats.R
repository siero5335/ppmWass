#!/usr/bin/env Rscript

# Compatibility stub retained so historical command paths fail explicitly.

stop(
  paste(
    "PPMWASS_LEGACY_DRIVER_DISABLED: run_controlled_bootstrap_stats.R",
    "reconstructs superseded controlled inputs and applies pre-audit",
    "optimistic-tie/query-bootstrap rules. Use",
    "prepare_controlled_publication_inputs.R followed by",
    "run_main_retrieval_publication.R."
  ),
  call. = FALSE
)
