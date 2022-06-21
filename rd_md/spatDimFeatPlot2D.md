# `spatDimFeatPlot2D`

spatDimFeatPlot2D


## Description

Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode


## Usage

```r
spatDimFeatPlot2D(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  show_image = F,
  gimage = NULL,
  image_name = NULL,
  largeImage_name = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  plot_alignment = c("vertical", "horizontal"),
  feats,
  dim_reduction_to_use = "umap",
  dim_reduction_name = "umap",
  dim1_to_use = 1,
  dim2_to_use = 2,
  dim_point_shape = c("border", "no_border"),
  dim_point_size = 1,
  dim_point_alpha = 1,
  dim_point_border_col = "black",
  dim_point_border_stroke = 0.1,
  show_NN_network = F,
  show_spatial_network = F,
  dim_network_color = "gray",
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  dim_edge_alpha = NULL,
  scale_alpha_with_expression = FALSE,
  sdimx = "sdimx",
  sdimy = "sdimy",
  spatial_network_name = "Delaunay_network",
  spatial_network_color = NULL,
  show_spatial_grid = F,
  grid_color = NULL,
  spatial_grid_name = "spatial_grid",
  spat_point_shape = c("border", "no_border", "voronoi"),
  spat_point_size = 1,
  spat_point_alpha = 1,
  spat_point_border_col = "black",
  spat_point_border_stroke = 0.1,
  spat_edge_alpha = NULL,
  cell_color_gradient = c("blue", "white", "red"),
  gradient_midpoint = NULL,
  gradient_limits = NULL,
  cow_n_col = 2,
  cow_rel_h = 1,
  cow_rel_w = 1,
  cow_align = "h",
  show_legend = TRUE,
  legend_text = 10,
  dim_background_color = "white",
  spat_background_color = "white",
  vor_border_color = "white",
  vor_max_radius = 200,
  vor_alpha = 1,
  axis_text = 8,
  axis_title = 8,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spatDimFeatPlot2D"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`show_image`     |     show a tissue background image
`gimage`     |     a giotto image
`image_name`     |     name of a giotto image
`largeImage_name`     |     name of a giottoLargeImage
`expression_values`     |     feat expression values to use
`plot_alignment`     |     direction to align plot
`feats`     |     features to show
`dim_reduction_to_use`     |     dimension reduction to use
`dim_reduction_name`     |     dimension reduction name
`dim1_to_use`     |     dimension to use on x-axis
`dim2_to_use`     |     dimension to use on y-axis
`dim_point_shape`     |     dim reduction points with border or not (border or no_border)
`dim_point_size`     |     dim reduction plot: point size
`dim_point_alpha`     |     transparancy of dim. reduction points
`dim_point_border_col`     |     color of border around points
`dim_point_border_stroke`     |     stroke size of border around points
`show_NN_network`     |     show underlying NN network
`show_spatial_network`     |     show underlying spatial netwok
`dim_network_color`     |     color of NN network
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use, if show_NN_network = TRUE
`dim_edge_alpha`     |     dim reduction plot: column to use for alpha of the edges
`scale_alpha_with_expression`     |     scale expression with ggplot alpha parameter
`sdimx`     |     spatial x-axis dimension name (default = 'sdimx')
`sdimy`     |     spatial y-axis dimension name (default = 'sdimy')
`spatial_network_name`     |     name of spatial network to use
`spatial_network_color`     |     color of spatial network
`show_spatial_grid`     |     show spatial grid
`grid_color`     |     color of spatial grid
`spatial_grid_name`     |     name of spatial grid to use
`spat_point_shape`     |     spatial points with border or not (border or no_border)
`spat_point_size`     |     spatial plot: point size
`spat_point_alpha`     |     transparancy of spatial points
`spat_point_border_col`     |     color of border around points
`spat_point_border_stroke`     |     stroke size of border around points
`spat_edge_alpha`     |     edge alpha
`cell_color_gradient`     |     vector with 3 colors for numeric data
`gradient_midpoint`     |     midpoint for color gradient
`gradient_limits`     |     vector with lower and upper limits
`cow_n_col`     |     cowplot param: how many columns
`cow_rel_h`     |     cowplot param: relative height
`cow_rel_w`     |     cowplot param: relative width
`cow_align`     |     cowplot param: how to align
`show_legend`     |     show legend
`legend_text`     |     size of legend text
`dim_background_color`     |     color of plot background for dimension plot
`spat_background_color`     |     color of plot background for spatial plot
`vor_border_color`     |     border colorr for voronoi plot
`vor_max_radius`     |     maximum radius for voronoi 'cells'
`vor_alpha`     |     transparancy of voronoi 'cells'
`axis_text`     |     size of axis text
`axis_title`     |     size of axis title
`show_plot`     |     show plots
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Description of parameters.


## Value

ggplot


