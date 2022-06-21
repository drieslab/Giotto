# `spatDimCellPlot2D`

spatDimCellPlot2D


## Description

Visualize numerical features of cells according to spatial AND dimension reduction coordinates in 2D


## Usage

```r
spatDimCellPlot2D(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  show_image = F,
  gimage = NULL,
  image_name = NULL,
  largeImage_name = NULL,
  plot_alignment = c("vertical", "horizontal"),
  spat_enr_names = NULL,
  cell_annotation_values = NULL,
  dim_reduction_to_use = "umap",
  dim_reduction_name = "umap",
  dim1_to_use = 1,
  dim2_to_use = 2,
  sdimx = "sdimx",
  sdimy = "sdimy",
  cell_color_gradient = c("blue", "white", "red"),
  gradient_midpoint = NULL,
  gradient_limits = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  dim_point_shape = c("border", "no_border"),
  dim_point_size = 1,
  dim_point_alpha = 1,
  dim_point_border_col = "black",
  dim_point_border_stroke = 0.1,
  spat_point_shape = c("border", "no_border", "voronoi"),
  spat_point_size = 1,
  spat_point_alpha = 1,
  spat_point_border_col = "black",
  spat_point_border_stroke = 0.1,
  dim_show_cluster_center = F,
  dim_show_center_label = T,
  dim_center_point_size = 4,
  dim_center_point_border_col = "black",
  dim_center_point_border_stroke = 0.1,
  dim_label_size = 4,
  dim_label_fontface = "bold",
  spat_show_cluster_center = F,
  spat_show_center_label = F,
  spat_center_point_size = 4,
  spat_center_point_border_col = "black",
  spat_center_point_border_stroke = 0.1,
  spat_label_size = 4,
  spat_label_fontface = "bold",
  show_NN_network = F,
  nn_network_to_use = "sNN",
  nn_network_name = "sNN.pca",
  dim_edge_alpha = 0.5,
  spat_show_network = F,
  spatial_network_name = "Delaunay_network",
  spat_network_color = "red",
  spat_network_alpha = 0.5,
  spat_show_grid = F,
  spatial_grid_name = "spatial_grid",
  spat_grid_color = "green",
  show_other_cells = TRUE,
  other_cell_color = "grey",
  dim_other_point_size = 0.5,
  spat_other_point_size = 0.5,
  spat_other_cells_alpha = 0.5,
  show_legend = T,
  legend_text = 8,
  legend_symbol_size = 1,
  dim_background_color = "white",
  spat_background_color = "white",
  vor_border_color = "white",
  vor_max_radius = 200,
  vor_alpha = 1,
  axis_text = 8,
  axis_title = 8,
  coord_fix_ratio = NULL,
  cow_n_col = 2,
  cow_rel_h = 1,
  cow_rel_w = 1,
  cow_align = "h",
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spatDimCellPlot2D"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`show_image`     |     show a tissue background image
`gimage`     |     a giotto image
`image_name`     |     name of a giotto image
`largeImage_name`     |     name of a giottoLargeImage
`plot_alignment`     |     direction to align plot
`spat_enr_names`     |     names of spatial enrichment results to include
`cell_annotation_values`     |     numeric cell annotation columns
`dim_reduction_to_use`     |     dimension reduction to use
`dim_reduction_name`     |     dimension reduction name
`dim1_to_use`     |     dimension to use on x-axis
`dim2_to_use`     |     dimension to use on y-axis
`sdimx`     |     = spatial dimension to use on x-axis
`sdimy`     |     = spatial dimension to use on y-axis
`cell_color_gradient`     |     vector with 3 colors for numeric data
`gradient_midpoint`     |     midpoint for color gradient
`gradient_limits`     |     vector with lower and upper limits
`select_cell_groups`     |     select subset of cells/clusters based on cell_color parameter
`select_cells`     |     select subset of cells based on cell IDs
`dim_point_shape`     |     dim reduction points with border or not (border or no_border)
`dim_point_size`     |     size of points in dim. reduction space
`dim_point_alpha`     |     transparancy of dim. reduction points
`dim_point_border_col`     |     border color of points in dim. reduction space
`dim_point_border_stroke`     |     border stroke of points in dim. reduction space
`spat_point_shape`     |     shape of points (border, no_border or voronoi)
`spat_point_size`     |     size of spatial points
`spat_point_alpha`     |     transparancy of spatial points
`spat_point_border_col`     |     border color of spatial points
`spat_point_border_stroke`     |     border stroke of spatial points
`dim_show_cluster_center`     |     show the center of each cluster
`dim_show_center_label`     |     provide a label for each cluster
`dim_center_point_size`     |     size of the center point
`dim_center_point_border_col`     |     border color of center point
`dim_center_point_border_stroke`     |     stroke size of center point
`dim_label_size`     |     size of the center label
`dim_label_fontface`     |     font of the center label
`spat_show_cluster_center`     |     show the center of each cluster
`spat_show_center_label`     |     provide a label for each cluster
`spat_center_point_size`     |     size of the spatial center points
`spat_center_point_border_col`     |     border color of the spatial center points
`spat_center_point_border_stroke`     |     stroke size of the spatial center points
`spat_label_size`     |     size of the center label
`spat_label_fontface`     |     font of the center label
`show_NN_network`     |     show underlying NN network
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`nn_network_name`     |     name of NN network to use, if show_NN_network = TRUE
`dim_edge_alpha`     |     column to use for alpha of the edges
`spat_show_network`     |     show spatial network
`spatial_network_name`     |     name of spatial network to use
`spat_network_color`     |     color of spatial network
`spat_network_alpha`     |     alpha of spatial network
`spat_show_grid`     |     show spatial grid
`spatial_grid_name`     |     name of spatial grid to use
`spat_grid_color`     |     color of spatial grid
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color of not selected cells
`dim_other_point_size`     |     size of not selected dim cells
`spat_other_point_size`     |     size of not selected spat cells
`spat_other_cells_alpha`     |     alpha of not selected spat cells
`show_legend`     |     show legend
`legend_text`     |     size of legend text
`legend_symbol_size`     |     size of legend symbols
`dim_background_color`     |     background color of points in dim. reduction space
`spat_background_color`     |     background color of spatial points
`vor_border_color`     |     border colorr for voronoi plot
`vor_max_radius`     |     maximum radius for voronoi 'cells'
`vor_alpha`     |     transparancy of voronoi 'cells'
`axis_text`     |     size of axis text
`axis_title`     |     size of axis title
`coord_fix_ratio`     |     ratio for coordinates
`cow_n_col`     |     cowplot param: how many columns
`cow_rel_h`     |     cowplot param: relative height
`cow_rel_w`     |     cowplot param: relative width
`cow_align`     |     cowplot param: how to align
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Description of parameters.


## Value

ggplot


## Seealso

Other spatial and dimension reduction cell annotation visualizations:
 [`spatDimCellPlot`](#spatdimcellplot)


