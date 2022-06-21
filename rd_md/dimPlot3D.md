# `dimPlot3D`

dimPlot3D


## Description

Visualize cells according to dimension reduction coordinates


## Usage

```r
dimPlot3D(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  dim_reduction_to_use = "umap",
  dim_reduction_name = "umap",
  dim1_to_use = 1,
  dim2_to_use = 2,
  dim3_to_use = 3,
  spat_enr_names = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 2,
  show_NN_network = F,
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  color_as_factor = T,
  cell_color = NULL,
  cell_color_code = NULL,
  show_cluster_center = F,
  show_center_label = T,
  center_point_size = 4,
  label_size = 4,
  edge_alpha = NULL,
  point_size = 3,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "dim3D"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`dim_reduction_to_use`     |     dimension reduction to use
`dim_reduction_name`     |     dimension reduction name
`dim1_to_use`     |     dimension to use on x-axis
`dim2_to_use`     |     dimension to use on y-axis
`dim3_to_use`     |     dimension to use on z-axis
`spat_enr_names`     |     names of spatial enrichment results to include
`select_cell_groups`     |     select subset of cells/clusters based on cell_color parameter
`select_cells`     |     select subset of cells based on cell IDs
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color of not selected cells
`other_point_size`     |     size of not selected cells
`show_NN_network`     |     show underlying NN network
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use, if show_NN_network = TRUE
`color_as_factor`     |     convert color column to factor
`cell_color`     |     color for cells (see details)
`cell_color_code`     |     named vector with colors
`show_cluster_center`     |     plot center of selected clusters
`show_center_label`     |     plot label of selected clusters
`center_point_size`     |     size of center points
`label_size`     |     size of labels
`edge_alpha`     |     column to use for alpha of the edges
`point_size`     |     size of point (cell)
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Description of parameters.


## Value

plotly


## Seealso

Other reduced dimension visualizations:
 [`dimPlot2D`](#dimplot2d) ,
 [`dimPlot`](#dimplot) ,
 [`plotPCA_2D`](#plotpca2d) ,
 [`plotPCA_3D`](#plotpca3d) ,
 [`plotPCA`](#plotpca) ,
 [`plotTSNE_2D`](#plottsne2d) ,
 [`plotTSNE_3D`](#plottsne3d) ,
 [`plotTSNE`](#plottsne) ,
 [`plotUMAP_2D`](#plotumap2d) ,
 [`plotUMAP_3D`](#plotumap3d) ,
 [`plotUMAP`](#plotumap)


