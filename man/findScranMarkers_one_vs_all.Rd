% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_expression.R
\name{findScranMarkers_one_vs_all}
\alias{findScranMarkers_one_vs_all}
\title{findScranMarkers_one_vs_all}
\usage{
findScranMarkers_one_vs_all(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column,
  subset_clusters = NULL,
  pval = 0.01,
  logFC = 0.5,
  min_feats = 10,
  min_genes = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{gobject}{giotto object}

\item{spat_unit}{spatial unit}

\item{feat_type}{feature type}

\item{expression_values}{feat expression values to use}

\item{cluster_column}{clusters to use}

\item{subset_clusters}{subset of clusters to use}

\item{pval}{filter on minimal p-value}

\item{logFC}{filter on logFC}

\item{min_feats}{minimum feats to keep per cluster, overrides pval and logFC}

\item{min_genes}{deprecated, use min_feats}

\item{verbose}{be verbose}

\item{...}{additional parameters for the findMarkers function in scran}
}
\value{
data.table with marker feats
}
\description{
Identify marker feats for all clusters in a one vs all manner based on scran's implementation of findMarkers.
}
\seealso{
\code{\link{findScranMarkers}}
}
