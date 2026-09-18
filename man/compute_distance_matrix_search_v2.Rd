\name{compute_distance_matrix_search_v2}
\alias{compute_distance_matrix_search_v2}
\title{Compute Dense Sparse-Bin Distance Matrix}
\usage{
compute_distance_matrix_search_v2(Q_set, L_set, params, use_loss = NULL,
  max_dense_cells = 1e7, progress = TRUE)
}
\arguments{
\item{Q_set, L_set}{\code{SpectraSet} objects sharing the same bin grid.}
\item{params}{Parameter list.}
\item{use_loss}{Whether to include loss spectra. \code{NULL} uses loss when both sets contain loss spectra.}
\item{max_dense_cells}{Maximum number of dense matrix cells allowed.}
\item{progress}{If \code{TRUE}, emit progress messages.}
}
\value{
Dense numeric distance matrix.
}
\description{
Experimental dense matrix API for small/medium GC benchmarks. For LC-scale
retrieval, use \code{\link{compute_distance_topk_v2}}.
}
