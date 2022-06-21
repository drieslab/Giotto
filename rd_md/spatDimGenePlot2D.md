# `spatDimGenePlot2D`

spatDimGenePlot2D


## Description

Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode


## Usage

```r
spatDimGenePlot2D(gobject, genes, default_save_name = "spatDimGenePlot2D", ...)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`genes`     |     genes to show
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param
`...`     |      Arguments passed on to [`spatDimFeatPlot2D`](#spatdimfeatplot2d)   list("\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("show_image")), list("show a tissue background image")), "\n", "    ", list(list(list("gimage")), list("a giotto image")), "\n", "    ", list(list(list("image_name")), list("name of a giotto image")), "\n", "    ", list(list(list("largeImage_name")), list("name of a giottoLargeImage")), "\n", "    ", list(list(list("expression_values")), 
    list("feat expression values to use")), "\n", "    ", list(list(list("plot_alignment")), list("direction to align plot")), "\n", "    ", list(list(list("dim_reduction_to_use")), list("dimension reduction to use")), "\n", "    ", list(list(list("dim_reduction_name")), list("dimension reduction name")), "\n", "    ", list(list(list("dim1_to_use")), list("dimension to use on x-axis")), "\n", "    ", list(list(list("dim2_to_use")), list("dimension to use on y-axis")), "\n", "    ", list(list(list(
    "dim_point_shape")), list("dim reduction points with border or not (border or no_border)")), "\n", "    ", list(list(list("dim_point_size")), list("dim reduction plot: point size")), "\n", "    ", list(list(list("dim_point_alpha")), list("transparancy of dim. reduction points")), "\n", "    ", list(list(list("dim_point_border_col")), list("color of border around points")), "\n", "    ", list(list(list("dim_point_border_stroke")), list("stroke size of border around points")), "\n", "    ", list(
    list(list("show_NN_network")), list("show underlying NN network")), "\n", "    ", list(list(list("show_spatial_network")), list("show underlying spatial netwok")), "\n", "    ", list(list(list("nn_network_to_use")), list("type of NN network to use (kNN vs sNN)")), "\n", "    ", list(list(list("network_name")), list("name of NN network to use, if show_NN_network = TRUE")), "\n", "    ", list(list(list("dim_network_color")), list("color of NN network")), "\n", "    ", list(list(list("dim_edge_alpha")), 
    list("dim reduction plot: column to use for alpha of the edges")), "\n", "    ", list(list(list("scale_alpha_with_expression")), list("scale expression with ggplot alpha parameter")), "\n", "    ", list(list(list("sdimx")), list("spatial x-axis dimension name (default = 'sdimx')")), "\n", "    ", list(list(list("sdimy")), list("spatial y-axis dimension name (default = 'sdimy')")), "\n", "    ", list(list(list("spatial_network_name")), list("name of spatial network to use")), "\n", "    ", list(
    list(list("spatial_network_color")), list("color of spatial network")), "\n", "    ", list(list(list("show_spatial_grid")), list("show spatial grid")), "\n", "    ", list(list(list("grid_color")), list("color of spatial grid")), "\n", "    ", list(list(list("spatial_grid_name")), list("name of spatial grid to use")), "\n", "    ", list(list(list("spat_point_shape")), list("spatial points with border or not (border or no_border)")), "\n", "    ", list(list(list("spat_point_size")), list("spatial plot: point size")), 
    "\n", "    ", list(list(list("spat_point_alpha")), list("transparancy of spatial points")), "\n", "    ", list(list(list("spat_point_border_col")), list("color of border around points")), "\n", "    ", list(list(list("spat_point_border_stroke")), list("stroke size of border around points")), "\n", "    ", list(list(list("spat_edge_alpha")), list("edge alpha")), "\n", "    ", list(list(list("cell_color_gradient")), list("vector with 3 colors for numeric data")), "\n", "    ", list(list(list("gradient_midpoint")), 
        list("midpoint for color gradient")), "\n", "    ", list(list(list("gradient_limits")), list("vector with lower and upper limits")), "\n", "    ", list(list(list("show_legend")), list("show legend")), "\n", "    ", list(list(list("legend_text")), list("size of legend text")), "\n", "    ", list(list(list("dim_background_color")), list("color of plot background for dimension plot")), "\n", "    ", list(list(list("spat_background_color")), list("color of plot background for spatial plot")), 
    "\n", "    ", list(list(list("vor_border_color")), list("border colorr for voronoi plot")), "\n", "    ", list(list(list("vor_max_radius")), list("maximum radius for voronoi 'cells'")), "\n", "    ", list(list(list("vor_alpha")), list("transparancy of voronoi 'cells'")), "\n", "    ", list(list(list("axis_text")), list("size of axis text")), "\n", "    ", list(list(list("axis_title")), list("size of axis title")), "\n", "    ", list(list(list("cow_n_col")), list("cowplot param: how many columns")), 
    "\n", "    ", list(list(list("cow_rel_h")), list("cowplot param: relative height")), "\n", "    ", list(list(list("cow_rel_w")), list("cowplot param: relative width")), "\n", "    ", list(list(list("cow_align")), list("cowplot param: how to align")), "\n", "    ", list(list(list("show_plot")), list("show plots")), "\n", "    ", list(list(list("return_plot")), list("return ggplot object")), "\n", "    ", list(list(list("save_plot")), list("directly save the plot [boolean]")), "\n", "    ", list(
        list(list("save_param")), list("list of saving parameters, see ", list(list("showSaveParameters")))), "\n", "  ")


## Details

Description of parameters.


## Value

ggplot


## Seealso

[`spatDimGenePlot3D`](#spatdimgeneplot3d) 
 
 Other spatial and dimension reduction gene expression visualizations:
 [`spatDimGenePlot3D`](#spatdimgeneplot3d) ,
 [`spatDimGenePlot`](#spatdimgeneplot)


