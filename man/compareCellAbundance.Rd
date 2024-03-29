% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactivity.R
\name{compareCellAbundance}
\alias{compareCellAbundance}
\title{Compare cell types percent per polygon}
\usage{
compareCellAbundance(
  gobject,
  polygon_name = "selections",
  spat_unit = "cell",
  feat_type = "rna",
  cell_type_column = "leiden_clus",
  ...
)
}
\arguments{
\item{gobject}{A Giotto object}

\item{polygon_name}{name of polygon selections}

\item{spat_unit}{spatial unit. Default = "cell"}

\item{feat_type}{feature type. Default =  "rna"}

\item{cell_type_column}{column name within the cell metadata table to use}

\item{...}{Additional parameters passed to ComplexHeatmap::Heatmap}
}
\value{
A ComplexHeatmap::Heatmap
}
\description{
Compare cell types percent per polygon
}
