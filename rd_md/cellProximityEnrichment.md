# `cellProximityEnrichment`

cellProximityEnrichment


## Description

Compute cell-cell interaction enrichment (observed vs expected)


## Usage

```r
cellProximityEnrichment(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  spatial_network_name = "Delaunay_network",
  cluster_column,
  number_of_simulations = 1000,
  adjust_method = c("none", "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel",
    "BY"),
  set_seed = TRUE,
  seed_number = 1234
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`spatial_network_name`     |     name of spatial network to use
`cluster_column`     |     name of column to use for clusters
`number_of_simulations`     |     number of simulations to create expected observations
`adjust_method`     |     method to adjust p.values
`set_seed`     |     use of seed
`seed_number`     |     seed number to use


## Details

Spatial proximity enrichment or depletion between pairs of cell types
 is calculated by calculating the observed over the expected frequency
 of cell-cell proximity interactions. The expected frequency is the average frequency
 calculated from a number of spatial network simulations. Each individual simulation is
 obtained by reshuffling the cell type labels of each node (cell)
 in the spatial network.


## Value

List of cell Proximity scores (CPscores) in data.table format. The first
 data.table (raw_sim_table) shows the raw observations of both the original and
 simulated networks. The second data.table (enrichm_res) shows the enrichment results.


