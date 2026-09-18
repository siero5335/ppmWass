\name{as_spectra_set}
\alias{as_spectra_set}
\title{Build a SpectraSet}
\usage{
as_spectra_set(frag_list, loss_list = NULL, ids = NULL,
  mz_min = 35, mz_max = 650, bin_ppm = 2, smear_ppm = 15,
  smear_kind = c("tri", "rect"))
}
\arguments{
\item{frag_list}{List of two-column fragment spectra matrices.}
\item{loss_list}{Optional list of two-column loss spectra matrices.}
\item{ids}{Optional spectrum IDs. Defaults to \code{names(frag_list)}.}
\item{mz_min, mz_max}{Mass range for the ppm-spaced bin grid.}
\item{bin_ppm}{Bin spacing in ppm.}
\item{smear_ppm}{Half-width of the peak smearing kernel in ppm.}
\item{smear_kind}{Smearing kernel, either \code{"tri"} or \code{"rect"}.}
}
\value{
A \code{SpectraSet} S3 object.
}
\description{
Build a lightweight container for sparse-bin distance backends. Sparse
matrices are built lazily on first use and cached per side/weighting scheme.
}
