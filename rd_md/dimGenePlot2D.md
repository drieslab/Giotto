# `dimGenePlot2D`

dimGenePlot2D


## Description

Visualize gene expression according to dimension reduction coordinates


## Usage

```r
dimGenePlot2D(gobject, genes = NULL, default_save_name = "dimGenePlot2D", ...)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`genes`     |     genes to show
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param
`...`     |      Arguments passed on to [`dimFeatPlot2D`](#dimfeatplot2d)   list("\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("expression_values")), list("gene expression values to use")), "\n", "    ", list(list(list("dim_reduction_to_use")), list("dimension reduction to use")), "\n", "    ", list(list(list("dim_reduction_name")), list("dimension reduction name")), "\n", "    ", list(list(list("dim1_to_use")), list("dimension to use on x-axis")), "\n", "    ", 
    list(list(list("dim2_to_use")), list("dimension to use on y-axis")), "\n", "    ", list(list(list("show_NN_network")), list("show underlying NN network")), "\n", "    ", list(list(list("nn_network_to_use")), list("type of NN network to use (kNN vs sNN)")), "\n", "    ", list(list(list("network_name")), list("name of NN network to use, if show_NN_network = TRUE")), "\n", "    ", list(list(list("network_color")), list("color of NN network")), "\n", "    ", list(list(list("edge_alpha")), list("column to use for alpha of the edges")), 
    "\n", "    ", list(list(list("scale_alpha_with_expression")), list("scale expression with ggplot alpha parameter")), "\n", "    ", list(list(list("point_shape")), list("point with border or not (border or no_border)")), "\n", "    ", list(list(list("point_size")), list("size of point (cell)")), "\n", "    ", list(list(list("point_alpha")), list("transparancy of points")), "\n", "    ", list(list(list("cell_color_gradient")), list("vector with 3 colors for numeric data")), "\n", "    ", list(list(
        list("gradient_midpoint")), list("midpoint for color gradient")), "\n", "    ", list(list(list("gradient_limits")), list("vector with lower and upper limits")), "\n", "    ", list(list(list("point_border_col")), list("color of border around points")), "\n", "    ", list(list(list("point_border_stroke")), list("stroke size of border around points")), "\n", "    ", list(list(list("show_legend")), list("show legend")), "\n", "    ", list(list(list("legend_text")), list("size of legend text")), 
    "\n", "    ", list(list(list("background_color")), list("color of plot background")), "\n", "    ", list(list(list("axis_text")), list("size of axis text")), "\n", "    ", list(list(list("axis_title")), list("size of axis title")), "\n", "    ", list(list(list("cow_n_col")), list("cowplot param: how many columns")), "\n", "    ", list(list(list("cow_rel_h")), list("cowplot param: relative height")), "\n", "    ", list(list(list("cow_rel_w")), list("cowplot param: relative width")), "\n", "    ", 
    list(list(list("cow_align")), list("cowplot param: how to align")), "\n", "    ", list(list(list("show_plot")), list("show plots")), "\n", "    ", list(list(list("return_plot")), list("return ggplot object")), "\n", "    ", list(list(list("save_plot")), list("directly save the plot [boolean]")), "\n", "    ", list(list(list("save_param")), list("list of saving parameters, see ", list(list("showSaveParameters")))), "\n", "  ")


## Details

Description of parameters.


## Value

ggplot


## Seealso

Other dimension reduction gene expression visualizations:
 [`dimGenePlot3D`](#dimgeneplot3d) ,
 [`dimGenePlot`](#dimgeneplot)


