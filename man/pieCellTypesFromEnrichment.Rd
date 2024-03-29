% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_enrichment_visuals.R
\name{pieCellTypesFromEnrichment}
\alias{pieCellTypesFromEnrichment}
\title{pieCellTypesFromEnrichment}
\usage{
pieCellTypesFromEnrichment(
  gobject = NULL,
  spat_unit = NULL,
  feat_type = NULL,
  enrichment_name = "PAGE_z_score",
  title = NULL,
  save_param = list(),
  default_save_name = "cell_types_from_enrichment_pie",
  save_plot = NA,
  show_plot = NA,
  return_plot = NA
)
}
\arguments{
\item{gobject}{Giotto Object}

\item{spat_unit}{spatial unit in which the enrichment information is stored}

\item{feat_type}{feature type for which the enrichment information was calculated}

\item{enrichment_name}{name of the spatial enrichment
i.e. output from GiottoClass::list_spatial_enrichment_names()
Default value is "PAGE_Z_score"}

\item{title}{Title of the generated plot.
Default `paste0(spat_unit,"cell types (maximum", enrichment_name, ")")`}

\item{save_param}{list of saving parameters, see \code{\link{showSaveParameters}}}

\item{default_save_name}{default save name for saving, don't change, change save_name in save_param}

\item{save_plot}{logical. save the plot}

\item{show_plot}{logical. show plot}

\item{return_plot}{logical. return ggplot object}
}
\description{
pieCellTypesFromEnrichment
}
\details{
This function generates a pie chart of cell types by frequency.
These cell type resutls are based on the provided `enrichment_name`,
and will be determined by the maximum value of the z-score
or p-value for a given cell or annotation.
}
