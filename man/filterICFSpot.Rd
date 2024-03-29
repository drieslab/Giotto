% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_interaction_spot.R
\name{filterICFSpot}
\alias{filterICFSpot}
\title{filterICFSpot}
\usage{
filterICFSpot(
  icfObject,
  min_cells = 4,
  min_cells_expr_resi = 0.05,
  min_int_cells = 4,
  min_int_cells_expr_resi = 0.05,
  min_fdr = 0.5,
  min_pcc_diff = 0.05,
  min_zscore = 0.05,
  zscores_column = c("cell_type", "features"),
  direction = c("both", "up", "down")
)
}
\arguments{
\item{icfObject}{ICF (interaction changed feature) score object}

\item{min_cells}{minimum number of source cell type}

\item{min_cells_expr_resi}{minimum expression residual level for source cell type}

\item{min_int_cells}{minimum number of interacting neighbor cell type}

\item{min_int_cells_expr_resi}{minimum expression residual level for interacting neighbor cell type}

\item{min_fdr}{minimum adjusted p-value}

\item{min_pcc_diff}{minimum absolute pcc difference difference}

\item{min_zscore}{minimum z-score change}

\item{zscores_column}{calculate z-scores over cell types or features}

\item{direction}{differential expression directions to keep}
}
\value{
icfObject that contains the filtered differential feature scores
}
\description{
Filter Interaction Changed Feature scores for spots.
}
