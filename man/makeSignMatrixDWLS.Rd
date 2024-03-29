% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_enrichment.R
\name{makeSignMatrixDWLS}
\alias{makeSignMatrixDWLS}
\title{makeSignMatrixDWLS}
\usage{
makeSignMatrixDWLS(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  reverse_log = TRUE,
  log_base = 2,
  sign_gene,
  cell_type_vector,
  cell_type = NULL
)
}
\arguments{
\item{gobject}{Giotto object of single cell}

\item{spat_unit}{spatial unit}

\item{feat_type}{feature type to use}

\item{expression_values}{expression values to use}

\item{reverse_log}{reverse a log-normalized expression matrix}

\item{log_base}{the logarithm base (default = 2)}

\item{sign_gene}{all of DE genes (signature)}

\item{cell_type_vector}{vector with cell types (length = ncol(matrix))}

\item{cell_type}{deprecated, use \code{cell_type_vector}}
}
\value{
matrix
}
\description{
Function to convert a matrix within a Giotto object into a format
 that can be used with \code{\link{runDWLSDeconv}} for deconvolution. A vector of cell types
 for parameter \code{cell_type_vector} can be created from the cell metadata (\code{\link{pDataDT}}).
}
\seealso{
\code{\link{runDWLSDeconv}}
}
