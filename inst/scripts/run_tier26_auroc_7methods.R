#!/usr/bin/env Rscript

# Compatibility stub retained so historical command paths fail explicitly.
# Ordered query-to-library AUROC is now produced by the main and supplementary
# publication drivers with finite-count and orientation checks.

stop(
  paste(
    "PPMWASS_LEGACY_DRIVER_DISABLED: run_tier26_auroc_7methods.R depends on",
    "an external pre-audit helper and computes a different, legacy AUROC",
    "estimand. Use run_main_retrieval_publication.R for ordered-pair AUROC or",
    "run_perturbation_publication.R for perturbation endpoints. The legacy",
    "top-hit-correctness distance AUROC has no authoritative one-to-one",
    "successor."
  ),
  call. = FALSE
)
