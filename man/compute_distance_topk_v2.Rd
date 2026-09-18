\name{compute_distance_topk_v2}
\alias{compute_distance_topk_v2}
\title{Compute Top-K Sparse-Bin Library Candidates}
\usage{
compute_distance_topk_v2(Q_set, L_set, params, top_k = 200L,
  expensive_method = NULL, expensive_top_k = NULL, block_size = 1024L,
  use_loss = NULL, return_prefilter = FALSE, progress = TRUE)
}
\arguments{
\item{Q_set, L_set}{\code{SpectraSet} objects sharing the same bin grid.}
\item{params}{Parameter list. \code{params$distance_method} is the cheap stage-1 method.}
\item{top_k}{Number of stage-1 candidates to keep per query.}
\item{expensive_method}{Optional method for rescoring stage-1 candidates.}
\item{expensive_top_k}{Number of candidates to keep after rescoring.}
\item{block_size}{Query block size.}
\item{use_loss}{Whether to include loss spectra. \code{NULL} uses loss when both sets contain loss spectra.}
\item{return_prefilter}{If \code{TRUE}, include stage-1 distances aligned to returned candidates.}
\item{progress}{If \code{TRUE}, emit block progress messages.}
}
\value{
A \code{TopKResult} object.
}
\description{
LC-scale retrieval API. Queries are processed in blocks and only the top-k
candidates per query are retained.
}
