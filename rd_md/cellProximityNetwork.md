# `cellProximityNetwork`

cellProximityNetwork


## Description

Create network from cell-cell proximity scores


## Usage

```r
cellProximityNetwork(
  gobject,
  CPscore,
  remove_self_edges = FALSE,
  self_loop_strength = 0.1,
  color_depletion = "lightgreen",
  color_enrichment = "red",
  rescale_edge_weights = TRUE,
  edge_weight_range_depletion = c(0.1, 1),
  edge_weight_range_enrichment = c(1, 5),
  layout = c("Fruchterman", "DrL", "Kamada-Kawai"),
  only_show_enrichment_edges = F,
  edge_width_range = c(0.1, 2),
  node_size = 4,
  node_color_code = NULL,
  node_text_size = 6,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "cellProximityNetwork"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`CPscore`     |     CPscore, output from cellProximityEnrichment()
`remove_self_edges`     |     remove enrichment/depletion edges with itself
`self_loop_strength`     |     size of self-loops
`color_depletion`     |     color for depleted cell-cell interactions
`color_enrichment`     |     color for enriched cell-cell interactions
`rescale_edge_weights`     |     rescale edge weights (boolean)
`edge_weight_range_depletion`     |     numerical vector of length 2 to rescale depleted edge weights
`edge_weight_range_enrichment`     |     numerical vector of length 2 to rescale enriched edge weights
`layout`     |     layout algorithm to use to draw nodes and edges
`only_show_enrichment_edges`     |     show only the enriched pairwise scores
`edge_width_range`     |     range of edge width
`node_size`     |     size of nodes
`node_color_code`     |     color code for the nodes (e.g. cell labels)
`node_text_size`     |     size of node labels
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

This function creates a network that shows the  spatial proximity
 enrichment or depletion of cell type pairs.


## Value

igraph plot


