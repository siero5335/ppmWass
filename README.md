# ppmWass

`ppmWass` is an R package for comparing EI-HRMS fragmentation spectra from GC-HRMS experiments. It combines fragment-ion distances with a pooled derived mass-difference representation, optional typical-difference projections, derivatization-aware heuristics, clustering support, and benchmark utilities in a single workflow. The derived representation contains both differences from a heuristically estimated high-mass reference and pairwise fragment differences; it is not a pure molecular-mass-based neutral-loss spectrum.

## Features

- Multiple distance backends: Hellinger, Wasserstein, ppm-aware Wasserstein, cosine, entropy, weighted cosine, and composite distance
- Optional `approxOT` acceleration for approximate `ppm_wasserstein` solvers (`sinkhorn`, `greenkhorn`), while exact OT remains on `transport`
- Pooled anchored and pairwise mass-difference comparison, with separable channels and optional typical-difference projection
- Mref confidence weighting for the reference-anchored channel
- Automatic derivatization detection and down-weighting for TMS/TBDMS signals
- MS-DIAL and MSP import/export helpers
- Benchmark helpers for retrieval, clustering, FDR, noise robustness, and homolog-style analyses

## Installation

```r
install.packages(c("remotes", "dplyr", "tibble", "magrittr", "stringr", "readr", "tidyr"))

# Optional but recommended for additional distance methods
# - transport: exact Wasserstein / exact ppm_wasserstein backend
# - approxOT: accelerated approximate ppm_wasserstein backend
# - msentropy: entropy similarity backend
install.packages(c("transport", "approxOT", "msentropy"))

remotes::install_local(".")
```

## Minimal Workflow

The package includes a synthetic MS-DIAL-like demo table, so you can run the core pipeline without external files.

```r
library(ppmWass)

params <- eihrms_default_params()

demo_df <- generate_demo_msdial(
  n_known = 6,
  n_unknown = 4,
  n_samples = 6,
  seed = 1
)

quant <- prepare_quant_data(
  demo_df,
  params,
  write_cleaned_csv = FALSE
)

spectra <- build_spectra(
  demo_df,
  quant$final,
  params,
  progress = FALSE
)

dist_mat <- compute_distance_matrix(
  spectra$frag_list,
  spectra$loss_list,
  params,
  progress = FALSE
)

sim_res <- compute_similarity_matrices(
  spectra$frag_list,
  spectra$loss_list,
  spectra$ri,
  params
)
```

## MS-DIAL Input Workflow

```r
library(ppmWass)

params <- eihrms_default_params()

df <- read_ms_dial("path/to/msdial_export.txt")
quant <- prepare_quant_data(df, params)
spectra <- build_spectra(df, quant$final, params)

benchmark_res <- run_benchmark(
  spectra,
  methods = c("cosine", "entropy", "wasserstein", "weighted_cosine", "composite"),
  params = params
)

benchmark_res$summary
```

For an end-to-end analysis run that writes plots and tables, use:

```r
res <- run_eihrms_similarity(
  df = df,
  params = params,
  write_outputs = FALSE,
  progress = FALSE
)
```

## Library Search

```r
hits <- search_library(
  query_spectra = spectra,
  library_spectra = spectra,
  params = params,
  top_n = 5
)

head(hits)
```

## PPM-Wasserstein Backend Choice

`ppm_wasserstein` supports both exact and approximate OT backends:

- `ot_method = "exact"` uses `transport::transport()` and is the default,
  publication-grade unregularized OT backend.
- `ot_method = "sinkhorn"` uses finite-iteration `approxOT` when explicitly
  requested for sensitivity or performance work.
- `ot_method = "greenkhorn"` uses `approxOT` and can be useful as an alternative approximate solver.

Example:

```r
params <- eihrms_default_params()
params$distance_method <- "ppm_wasserstein"
params$ot_method <- "exact"
params$sinkhorn_epsilon <- 0.05
params$sinkhorn_niter <- 100L
```

The epsilon and iteration settings apply only when an approximate backend is
selected. Exact OT is intentionally kept on `transport`: it directly evaluates
the unregularized transport objective, is numerically robust for arbitrary cost
matrices, and is invariant to exchanging the two spectra. Fixed-iteration
Sinkhorn and Greenkhorn are treated as directional numerical approximations and
are therefore never mirrored from one matrix triangle.

## Parameter Sweep Example

Use `sweep_distance_params()` for manuscript-facing sensitivity checks. The
example below evaluates two ppm scales and two fragment/derived weight balances,
with repeated 80% subsampling to gauge stability.

```r
sweep_res <- sweep_distance_params(
  spectra_result = spectra,
  grid = list(
    tol_ppm = c(10, 20),
    w_frag = c(0.7, 0.5),
    w_loss = c(0.3, 0.5)
  ),
  methods = c("cosine", "entropy", "ppm_wasserstein"),
  n_replicates = 5,
  sample_frac = 0.8,
  progress = FALSE
)

head(sweep_res$summary)
```

If you leave `sample_frac = 1`, the sweep runs once per grid row on the full
dataset. That is useful for deterministic parameter scans, but not for
estimating uncertainty.

## Reproducible Benchmark Script

The repository now includes a reproducible benchmark driver at [`inst/scripts/run_benchmark_reproducible.R`](inst/scripts/run_benchmark_reproducible.R).

Run it on synthetic demo data:

```bash
Rscript inst/scripts/run_benchmark_reproducible.R --output-dir=benchmark_output
```

Run it on an MS-DIAL export:

```bash
Rscript inst/scripts/run_benchmark_reproducible.R \
  --input=path/to/msdial_export.txt \
  --output-dir=benchmark_output \
  --methods=cosine,entropy,wasserstein,weighted_cosine,composite
```

The script writes:

- `benchmark_summary.csv`
- `distance_param_sweep.csv`
- `distance_runtime_results.csv`
- `distance_runtime_summary.csv`
- `sessionInfo.txt`
- optional benchmark plots when `ggplot2` is installed

Representative synthetic outputs generated from this workflow are bundled under
[`inst/extdata/benchmark_demo`](inst/extdata/benchmark_demo)
so the table structure can be inspected without rerunning the full pipeline.

For head-to-head comparison against a local MSP library, use [`inst/scripts/run_library_head_to_head.R`](inst/scripts/run_library_head_to_head.R):

```bash
Rscript inst/scripts/run_library_head_to_head.R \
  --query=query.msp \
  --library=library.msp \
  --output-dir=library_benchmark_output \
  --methods=cosine,entropy,wasserstein,weighted_cosine,composite,ppm_wasserstein
```

This writes a per-method retrieval summary plus per-query hit tables, so it can be used directly with local NIST/MoNA-style MSP exports.

A small packaged MSP example and its benchmark summary are also included under
[`inst/extdata/library_demo`](inst/extdata/library_demo).

See the vignette [`vignettes/reproducible-benchmarks.Rmd`](vignettes/reproducible-benchmarks.Rmd)
for a walkthrough that reads these bundled artifacts and shows how to regenerate them.

## Benchmark Interpretation Notes

- `run_benchmark()` can compare retrieval and clustering behavior across multiple similarity methods with a shared parameter set.
- `tol_ppm` has method-specific semantics retained for backward compatibility. For hard-alignment comparators it is a peak-matching tolerance. For `ppm_wasserstein` it is the base ground-cost scale: displacement has nonzero cost below this value, and the cost saturates at `tol_ppm * wasserstein_transition_mult`.
- `w_loss` and other `loss` names are legacy API terminology for the pooled derived mass-difference representation.
- `sweep_distance_params()` can be used to assess sensitivity to `tol_ppm`, `w_frag`/`w_loss`, `loss_raw_beta`, and similar parameters over a grid; report the effective method-specific interpretation.
- `benchmark_distance_runtime()` can be used to generate timing tables across increasing library sizes.
- `benchmark_library_search()` benchmarks true query-vs-library retrieval when query and library spectra both carry InChIKey metadata.
- `evaluate_homolog_detection()` is most persuasive when you provide external class labels or curated homolog-series labels.
- If `compound_class` was inferred internally by this package, homolog detection should be interpreted as an internal consistency check rather than external ground truth.
- `create_homolog_ground_truth()` groups by the first 14 characters of InChIKey, which correspond to the connectivity layer. This is useful for structural grouping, but it is not a substitute for curated homolog-series annotations.
- For manuscript reporting, treat full-dataset benchmark summaries as point estimates unless you also run repeated subsampling or another resampling scheme.
- Publication reruns use the scripts whose names end in `_publication.R` and
  their pinned bundle orchestrator. Historical drivers that applied pre-audit
  parsing, tie, bootstrap, AUROC, output-layout, or nonfinite rules remain only
  as fail-closed compatibility stubs with an explicit replacement-script
  message; they must not be used to generate current manuscript results.

## Experimental LC-HRMS/MS Outer Workflow

The package includes an experimental fragment-only LC workflow that keeps
precursor mass, ion mode, adduct, retention time, and collision energy as
auditable candidate gates before Top-K spectral reranking. Start with
`read_lc_msp()`, `as_lc_spectra()`, `lc_search_plan()`, and
`search_lc_library()`.

This is an implementation scaffold, not a claim of validated LC performance.
The LC metadata contract and validation boundaries are documented in the
function references and enforced by the accompanying tests.

The publication-grade explicit-TIC boundary analysis is implemented by
`inst/scripts/run_mixture_publication.R`. It independently TIC-normalizes the
signal and interferent components, records assigned pre-merge mixing weights
and post-merge source attribution, and evaluates fragment-only and combined
fragment/derived-representation modes. The superseded
`run_additional_experiments.R` path is retained only as a fail-closed
compatibility stub because its pre-audit nonfinite handling was unsafe.

## Status

The package already supports the main algorithmic workflow for EI-HRMS similarity analysis. Current development is focused on strengthening validation, expanding direct unit tests for distance methods, and improving reproducible benchmark documentation.
## Optional sparse exact OT backend

The package provides an explicit alternative solver for the **same**
ppm-Wasserstein distance. Existing calls keep `ot_method = "exact"` and the
original dense network-simplex path. To opt in:

```r
params <- eihrms_default_params()
params$distance_method <- "ppm_wasserstein"
params$ot_method <- "exact_sparse"
D_sparse <- compute_distance_matrix(frag_list, loss_list, params)

params$ot_method <- "exact"
D_exact <- compute_distance_matrix(frag_list, loss_list, params)
stopifnot(max(abs(D_sparse - D_exact)) < 1e-8)
```

The sparse solver exploits saturation of the ground cost and solves independent
components of near-mass peak pairs. It is unregularized, supports the existing
alignment and fragment/loss options, and validates numerical transport plans
for nonanalytic components. It never silently switches to an approximate
method or the original whole-matrix solver. Dense or very small spectra can
be slower. Tiny floating-point differences can reorder exact numerical ties;
compare ties with a tolerance. Source installation now requires Rcpp and a C++
compiler. Install this package version on socket workers as well as the parent
process when using Windows parallel execution.

Compare complete query-by-library distance calculations on the same data:

```sh
Rscript inst/scripts/run_exact_backend_comparison.R --input=library.msp --output-dir=backend_comparison
```

The script records selected records, settings, input checksum, both distance
matrices, repeated timings, Top-1 tie-set agreement and rank inversions outside
a 1e-8 tolerance. Existing `run_runtime_microbenchmark.R` also accepts
`--ot-method=exact_sparse`; use separate output directories for the two runs.
Treat `exact` and `exact_sparse` as solver variants, not different scientific
distance methods. Historical publication outputs are not regenerated.

## Validation records

Completed solver comparisons, timing tables, figures, and research-specific
analysis workflows are maintained in the separate
[ppmWass-validation repository](https://github.com/siero5335/ppmWass-validation).
That repository also records exploratory analyses that were not used in the
primary publication analysis. Input spectral libraries are not bundled there.

The package retains its implementation, regression tests, synthetic examples,
and reusable benchmark drivers. `exact` remains the default; `exact_sparse`
is an opt-in solver whose runtime benefit depends on the workload.
