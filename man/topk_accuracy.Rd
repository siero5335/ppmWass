\name{topk_accuracy}
\alias{topk_accuracy}
\title{Top-K Retrieval Accuracy}
\usage{
topk_accuracy(result, top_k = c(1, 3, 5, 10), ground_truth = NULL,
  exclude_self = FALSE)
}
\arguments{
\item{result}{A \code{TopKResult}.}
\item{top_k}{Integer vector of cutoffs.}
\item{ground_truth}{Optional named list or character vector of acceptable library IDs for each query.}
\item{exclude_self}{If \code{TRUE}, candidates with the same ID as the query are excluded before scoring.}
}
\value{
Named numeric vector of hit rates.
}
\description{
Computes simple ID-based hit rates from a \code{TopKResult}.
}
