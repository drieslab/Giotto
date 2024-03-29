% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_genes.R
\name{calculate_spatial_enrichment}
\alias{calculate_spatial_enrichment}
\alias{.calc_spatial_enrichment_minimum}
\alias{.calc_spatial_enrichment_matrix}
\alias{.calc_spatial_enrichment_dt}
\title{Calculate spatial enrichment}
\usage{
.calc_spatial_enrichment_minimum(
  spatial_network,
  bin_matrix,
  adjust_method = "fdr",
  do_fisher_test = TRUE
)

.calc_spatial_enrichment_matrix(
  spatial_network,
  bin_matrix,
  adjust_method = "fdr",
  do_fisher_test = TRUE,
  do_parallel = TRUE,
  cores = NA,
  calc_hub = FALSE,
  hub_min_int = 3,
  verbose = TRUE
)

.calc_spatial_enrichment_dt(
  bin_matrix,
  spatial_network,
  calc_hub = FALSE,
  hub_min_int = 3,
  group_size = "automatic",
  do_fisher_test = TRUE,
  adjust_method = "fdr",
  cores = NA
)
}
\description{
Calculate spatial enrichment. Multiple methods are provided.
}
\section{Functions}{
\itemize{
\item \code{.calc_spatial_enrichment_minimum()}: calculate using a 'simple' and efficient for loop

\item \code{.calc_spatial_enrichment_matrix()}: calculate using 'matrix' implementation

\item \code{.calc_spatial_enrichment_dt()}: calculate using 'data.table' implementation

}}
\keyword{internal}
