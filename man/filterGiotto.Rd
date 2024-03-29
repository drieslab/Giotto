% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliary_giotto.R
\name{filterGiotto}
\alias{filterGiotto}
\title{filterGiotto}
\usage{
filterGiotto(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("raw", "normalized", "scaled", "custom"),
  expression_threshold = 1,
  feat_det_in_min_cells = 100,
  min_det_feats_per_cell = 100,
  spat_unit_fsub = ":all:",
  feat_type_ssub = ":all:",
  all_spat_units = NULL,
  all_feat_types = NULL,
  poly_info = NULL,
  tag_cells = FALSE,
  tag_cell_name = "tag",
  tag_feats = FALSE,
  tag_feats_name = "tag",
  verbose = TRUE
)
}
\arguments{
\item{gobject}{giotto object}

\item{spat_unit}{character. spatial unit. If more than one is provided then
the first will be filtered, the filtering results will be applied across the
other spat_units provided}

\item{feat_type}{character. feature type. If more than one is provided then
the first will be filtered, the filtering results will be applied across the
other feat_types provided.}

\item{expression_values}{expression values to use}

\item{expression_threshold}{threshold to consider a gene expressed}

\item{feat_det_in_min_cells}{minimum # of cells that need to express a feature}

\item{min_det_feats_per_cell}{minimum # of features that need to be detected in a cell}

\item{spat_unit_fsub}{character vector. (default = ':all:') limit features
to remove results to selected spat_units}

\item{feat_type_ssub}{character vector. (default = ':all:') limit cells to
remove results to selected feat_types}

\item{all_spat_units}{deprecated. Use spat_unit_fsub = ":all:"}

\item{all_feat_types}{deprecated. Use feat_type_ssub = ":all:"}

\item{poly_info}{polygon information to use}

\item{tag_cells}{tag filtered cells in metadata vs. remove cells}

\item{tag_cell_name}{column name for tagged cells in metadata}

\item{tag_feats}{tag features in metadata vs. remove features}

\item{tag_feats_name}{column name for tagged features in metadata}

\item{verbose}{verbose}
}
\value{
giotto object
}
\description{
filter Giotto object based on expression threshold
}
\details{
The function \code{\link{filterCombinations}} can be used to explore the effect of different parameter values.
Please note that this function filters data in a predefined order, features, then cells.
After filtering in this order, certain features may be left over in the metadata with a
corresponding number of cells which is less than that of the threshold value of cells,
feat_det_in_min_cells. This behavior is explained in detail here:
\url{https://github.com/drieslab/Giotto/issues/500#issuecomment-1396083446}
}
