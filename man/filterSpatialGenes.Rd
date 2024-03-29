% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/python_hmrf.R
\name{filterSpatialGenes}
\alias{filterSpatialGenes}
\title{filterSpatialGenes}
\usage{
filterSpatialGenes(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  spatial_genes,
  max = 2500,
  name = c("binSpect", "silhouetteRank", "silhouetteRankTest"),
  method = c("none", "elbow")
)
}
\arguments{
\item{gobject}{Giotto object}

\item{spat_unit}{spatial unit}

\item{feat_type}{feature type}

\item{spatial_genes}{input gene list}

\item{max}{max number of genes selected from spatial test}

\item{name}{name of spatial gene test for the filtering}

\item{method}{method of spatial gene selection}
}
\description{
function to filter gene list with existing spatial gene sets
}
\details{
This function filters given gene list with the gene sets of selected spatial gene test in Giotto,
also controls the total size of the gene set with given max number.
}
\keyword{external}
