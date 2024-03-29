% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_interaction.R
\name{findCellProximityGenes}
\alias{findCellProximityGenes}
\title{findCellProximityGenes}
\usage{
findCellProximityGenes(...)
}
\arguments{
\item{...}{
  Arguments passed on to \code{\link[=findInteractionChangedFeats]{findInteractionChangedFeats}}
  \describe{
    \item{\code{gobject}}{giotto object}
    \item{\code{feat_type}}{feature type}
    \item{\code{spat_unit}}{spatial unit}
    \item{\code{expression_values}}{expression values to use}
    \item{\code{selected_feats}}{subset of selected features (optional)}
    \item{\code{cluster_column}}{name of column to use for cell types}
    \item{\code{spatial_network_name}}{name of spatial network to use}
    \item{\code{minimum_unique_cells}}{minimum number of target cells required}
    \item{\code{minimum_unique_int_cells}}{minimum number of interacting cells required}
    \item{\code{diff_test}}{which differential expression test}
    \item{\code{mean_method}}{method to use to calculate the mean}
    \item{\code{offset}}{offset value to use when calculating log2 ratio}
    \item{\code{adjust_method}}{which method to adjust p-values}
    \item{\code{nr_permutations}}{number of permutations if diff_test = permutation}
    \item{\code{exclude_selected_cells_from_test}}{exclude interacting cells other cells}
    \item{\code{do_parallel}}{run calculations in parallel with mclapply}
    \item{\code{set_seed}}{set a seed for reproducibility}
    \item{\code{seed_number}}{seed number}
  }}
}
\description{
Identifies cell-to-cell Interaction Changed Features (ICF),
i.e. genes that are differentially expressed due to proximity to other cell types.
}
\seealso{
\code{\link{findInteractionChangedFeats}}
}
