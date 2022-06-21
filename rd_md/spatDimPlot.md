# `spatDimPlot`

spatDimPlot


## Description

Visualize cells according to spatial AND dimension reduction coordinates 2D


## Usage

```r
spatDimPlot(...)
```


## Arguments

Argument      |Description
------------- |----------------
`...`     |      Arguments passed on to [`spatDimPlot2D`](#spatdimplot2d)   list("\n", "    ", list(list(list("gobject")), list("giotto object")), "\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("show_image")), list("show a tissue background image")), "\n", "    ", list(list(list("gimage")), list("a giotto image")), "\n", "    ", list(list(list("image_name")), list("name of a giotto image")), "\n", "    ", list(list(list("largeImage_name")), list("name of a giottoLargeImage")), 
    "\n", "    ", list(list(list("spat_loc_name")), list("name of spatial locations")), "\n", "    ", list(list(list("plot_alignment")), list("direction to align plot")), "\n", "    ", list(list(list("dim_reduction_to_use")), list("dimension reduction to use")), "\n", "    ", list(list(list("dim_reduction_name")), list("dimension reduction name")), "\n", "    ", list(list(list("dim1_to_use")), list("dimension to use on x-axis")), "\n", "    ", list(list(list("dim2_to_use")), list("dimension to use on y-axis")), 
    "\n", "    ", list(list(list("sdimx")), list("= spatial dimension to use on x-axis")), "\n", "    ", list(list(list("sdimy")), list("= spatial dimension to use on y-axis")), "\n", "    ", list(list(list("spat_enr_names")), list("names of spatial enrichment results to include")), "\n", "    ", list(list(list("cell_color")), list("color for cells (see details)")), "\n", "    ", list(list(list("color_as_factor")), list("convert color column to factor")), "\n", "    ", list(list(list("cell_color_code")), 
        list("named vector with colors")), "\n", "    ", list(list(list("cell_color_gradient")), list("vector with 3 colors for numeric data")), "\n", "    ", list(list(list("gradient_midpoint")), list("midpoint for color gradient")), "\n", "    ", list(list(list("gradient_limits")), list("vector with lower and upper limits")), "\n", "    ", list(list(list("select_cell_groups")), list("select subset of cells/clusters based on cell_color parameter")), "\n", "    ", list(list(list("select_cells")), 
        list("select subset of cells based on cell IDs")), "\n", "    ", list(list(list("dim_point_shape")), list("point with border or not (border or no_border)")), "\n", "    ", list(list(list("dim_point_size")), list("size of points in dim. reduction space")), "\n", "    ", list(list(list("dim_point_alpha")), list("transparancy of point in dim. reduction space")), "\n", "    ", list(list(list("dim_point_border_col")), list("border color of points in dim. reduction space")), "\n", "    ", list(
        list(list("dim_point_border_stroke")), list("border stroke of points in dim. reduction space")), "\n", "    ", list(list(list("spat_point_shape")), list("shape of points (border, no_border or voronoi)")), "\n", "    ", list(list(list("spat_point_size")), list("size of spatial points")), "\n", "    ", list(list(list("spat_point_alpha")), list("transparancy of spatial points")), "\n", "    ", list(list(list("spat_point_border_col")), list("border color of spatial points")), "\n", "    ", list(
        list(list("spat_point_border_stroke")), list("border stroke of spatial points")), "\n", "    ", list(list(list("dim_show_cluster_center")), list("show the center of each cluster")), "\n", "    ", list(list(list("dim_show_center_label")), list("provide a label for each cluster")), "\n", "    ", list(list(list("dim_center_point_size")), list("size of the center point")), "\n", "    ", list(list(list("dim_center_point_border_col")), list("border color of center point")), "\n", "    ", list(list(
        list("dim_center_point_border_stroke")), list("stroke size of center point")), "\n", "    ", list(list(list("dim_label_size")), list("size of the center label")), "\n", "    ", list(list(list("dim_label_fontface")), list("font of the center label")), "\n", "    ", list(list(list("spat_show_cluster_center")), list("show the center of each cluster")), "\n", "    ", list(list(list("spat_show_center_label")), list("provide a label for each cluster")), "\n", "    ", list(list(list("spat_center_point_size")), 
        list("size of the center point")), "\n", "    ", list(list(list("spat_center_point_border_col")), list("border color of spatial center points")), "\n", "    ", list(list(list("spat_center_point_border_stroke")), list("border strike size of spatial center points")), "\n", "    ", list(list(list("spat_label_size")), list("size of the center label")), "\n", "    ", list(list(list("spat_label_fontface")), list("font of the center label")), "\n", "    ", list(list(list("show_NN_network")), list(
        "show underlying NN network")), "\n", "    ", list(list(list("nn_network_to_use")), list("type of NN network to use (kNN vs sNN)")), "\n", "    ", list(list(list("network_name")), list("name of NN network to use, if show_NN_network = TRUE")), "\n", "    ", list(list(list("nn_network_alpha")), list("column to use for alpha of the edges")), "\n", "    ", list(list(list("show_spatial_network")), list("show spatial network")), "\n", "    ", list(list(list("spat_network_name")), list("name of spatial network to use")), 
    "\n", "    ", list(list(list("spat_network_color")), list("color of spatial network")), "\n", "    ", list(list(list("spat_network_alpha")), list("alpha of spatial network")), "\n", "    ", list(list(list("show_spatial_grid")), list("show spatial grid")), "\n", "    ", list(list(list("spat_grid_name")), list("name of spatial grid to use")), "\n", "    ", list(list(list("spat_grid_color")), list("color of spatial grid")), "\n", "    ", list(list(list("show_other_cells")), list("display not selected cells")), 
    "\n", "    ", list(list(list("other_cell_color")), list("color of not selected cells")), "\n", "    ", list(list(list("dim_other_point_size")), list("size of not selected dim cells")), "\n", "    ", list(list(list("spat_other_point_size")), list("size of not selected spat cells")), "\n", "    ", list(list(list("spat_other_cells_alpha")), list("alpha of not selected spat cells")), "\n", "    ", list(list(list("dim_show_legend")), list("show legend of dimension reduction plot")), "\n", "    ", 
    list(list(list("spat_show_legend")), list("show legend of spatial plot")), "\n", "    ", list(list(list("legend_text")), list("size of legend text")), "\n", "    ", list(list(list("legend_symbol_size")), list("size of legend symbols")), "\n", "    ", list(list(list("dim_background_color")), list("background color of points in dim. reduction space")), "\n", "    ", list(list(list("spat_background_color")), list("background color of spatial points")), "\n", "    ", list(list(list("vor_border_color")), 
        list("border colorr for voronoi plot")), "\n", "    ", list(list(list("vor_max_radius")), list("maximum radius for voronoi 'cells'")), "\n", "    ", list(list(list("vor_alpha")), list("transparancy of voronoi 'cells'")), "\n", "    ", list(list(list("axis_text")), list("size of axis text")), "\n", "    ", list(list(list("axis_title")), list("size of axis title")), "\n", "    ", list(list(list("show_plot")), list("show plot")), "\n", "    ", list(list(list("return_plot")), list("return ggplot object")), 
    "\n", "    ", list(list(list("save_plot")), list("directly save the plot [boolean]")), "\n", "    ", list(list(list("save_param")), list("list of saving parameters, see ", list(list("showSaveParameters")))), "\n", "    ", list(list(list("default_save_name")), list("default save name for saving, don't change, change save_name in save_param")), "\n", "  ")


## Details

Description of parameters.


## Value

ggplot


## Seealso

[`spatDimPlot2D`](#spatdimplot2d) and [`spatDimPlot3D`](#spatdimplot3d) for 3D visualization.
 
 Other spatial and dimension reduction visualizations:
 [`spatDimPlot2D`](#spatdimplot2d) ,
 [`spatDimPlot3D`](#spatdimplot3d)


