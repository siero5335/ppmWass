#!/usr/bin/env Rscript

# Compatibility stub retained so historical command paths fail explicitly.
# The complete pre-audit implementation remains available in Git history and
# in the publication rerun bundle's immutable archive.

stop(
  paste(
    "PPMWASS_LEGACY_DRIVER_DISABLED: run_additional_experiments.R is a",
    "superseded pre-audit driver and is intentionally disabled because it",
    "formerly replaced nonfinite matrix cells by the maximum finite value.",
    "Use the task-specific publication drivers:",
    "run_main_retrieval_publication.R, run_secondary_ablation_publication.R,",
    "run_mixture_publication.R, run_end_to_end_reranking_publication.R,",
    "run_perturbation_publication.R, run_nonfinite_publication_audit.R, or",
    "run_supplementary_stress_publication.R as appropriate."
  ),
  call. = FALSE
)
