# ppmWass 0.0.3.9001 (experimental)

- Merge the opt-in backend into main while retaining `exact` as the default.
- Add serial runtime and numerical-agreement figures with source CSV tables,
  plus a reproducible scaling driver (serial by default; parallel is explicit).
  Deferred/interrupted scaling measurements are excluded from the figures.
- Add opt-in `ot_method = "exact_sparse"` for ppm-Wasserstein. It computes the
  same unregularized objective through sparse connected components; the
  existing `"exact"` backend remains the default and its solver is unchanged.
- Sparse execution supports alignment, fragment/loss combinations and matrix
  routes. Diagnostics report component counts and validate general-component
  transport plans. There is no automatic switch to a different backend.
- Add `run_exact_backend_comparison.R` to compare both solvers on the same
  search matrix, including timings and tolerance-aware ranking agreement.
- The experimental package now needs Rcpp and a C++ compiler for source builds.

# ppmWass 0.0.3.9000

- Fixed library-search metadata after MSP filtering and deduplication, retaining
  the selected record's metadata and hit order.
- Fixed searches against a one-entry library.
- Excluded nonfinite and RI-rejected candidates from library-search evaluation;
  queries with no surviving candidates now contribute zero rather than hits.
- Added the nominal-setting analysis drivers under `tools/` with scoped lock
  and log cleanup and regression tests. The drivers retain the
  frozen publication commit and scientific settings.

# ppmWass 0.0.3

- Added a publication-diagnostic Sinkhorn path that executes exactly the
  requested epsilon and iteration count once, with no retry or exact fallback,
  and records plan residuals, raw failures, orientation, and solver provenance.
- Kept constant-ground-cost distances mathematically exact while preserving
  any raw approximate-solver failure as a separate diagnostic observation.
- Made supplementary aggregation retain intentional inactive-setting `NA`
  groups and require complete replicate coverage for Sinkhorn retrieval
  summaries.
- Added full-library, tie-aware first-pass Top-1 and MRR baselines to the
  end-to-end reranking benchmark and its serial/parallel equivalence audit.

# ppmWass 0.0.2

- Changed the default ppm-Wasserstein backend to exact unregularized optimal
  transport through `transport::transport()`.
- Made exact-backend requests fail closed instead of silently substituting a
  finite-iteration approximate solver or Hellinger distance.
- Made matrix symmetry handling solver-aware: exact ppm-Wasserstein is safely
  mirrored, whereas explicitly requested Sinkhorn and Greenkhorn calculations
  are evaluated in both query/library orientations.
- Added publication stress tests for exact symmetry, square-versus-rectangular
  agreement, and finite-iteration Sinkhorn direction dependence.
