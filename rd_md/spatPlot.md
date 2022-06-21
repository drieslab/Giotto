# `spatPlot`

spatPlot


## Description

Visualize cells according to spatial coordinates


## Usage

```r
spatPlot(...)
```


## Arguments

Argument      |Description
------------- |----------------
`...`     |      Arguments passed on to [`spatPlot2D`](#spatplot2d)   list("\n", "    ", list(list(list("gobject")), list("giotto object")), "\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("show_image")), list("show a tissue background image")), "\n", "    ", list(list(list("gimage")), list("a giotto image")), "\n", "    ", list(list(list("image_name")), list("name of a giotto image or multiple images with group_by")), "\n", "    ", list(list(list("largeImage_name")), 
    list("name of a giottoLargeImage or multiple images with group_by")), "\n", "    ", list(list(list("group_by")), list("create multiple plots based on cell annotation column")), "\n", "    ", list(list(list("group_by_subset")), list("subset the group_by factor column")), "\n", "    ", list(list(list("spat_loc_name")), list("name of spatial locations")), "\n", "    ", list(list(list("sdimx")), list("x-axis dimension name (default = 'sdimx')")), "\n", "    ", list(list(list("sdimy")), list("y-axis dimension name (default = 'sdimy')")), 
    "\n", "    ", list(list(list("spat_enr_names")), list("names of spatial enrichment results to include")), "\n", "    ", list(list(list("cell_color")), list("color for cells (see details)")), "\n", "    ", list(list(list("color_as_factor")), list("convert color column to factor")), "\n", "    ", list(list(list("cell_color_code")), list("named vector with colors")), "\n", "    ", list(list(list("cell_color_gradient")), list("vector with 3 colors for numeric data")), "\n", "    ", list(list(list(
        "gradient_midpoint")), list("midpoint for color gradient")), "\n", "    ", list(list(list("gradient_limits")), list("vector with lower and upper limits")), "\n", "    ", list(list(list("select_cell_groups")), list("select subset of cells/clusters based on cell_color parameter")), "\n", "    ", list(list(list("select_cells")), list("select subset of cells based on cell IDs")), "\n", "    ", list(list(list("point_shape")), list("shape of points (border, no_border or voronoi)")), "\n", "    ", 
    list(list(list("point_size")), list("size of point (cell)")), "\n", "    ", list(list(list("point_alpha")), list("transparancy of point")), "\n", "    ", list(list(list("point_border_col")), list("color of border around points")), "\n", "    ", list(list(list("point_border_stroke")), list("stroke size of border around points")), "\n", "    ", list(list(list("show_cluster_center")), list("plot center of selected clusters")), "\n", "    ", list(list(list("show_center_label")), list("plot label of selected clusters")), 
    "\n", "    ", list(list(list("center_point_size")), list("size of center points")), "\n", "    ", list(list(list("center_point_border_col")), list("border color of center points")), "\n", "    ", list(list(list("center_point_border_stroke")), list("border stroke size of center points")), "\n", "    ", list(list(list("label_size")), list("size of labels")), "\n", "    ", list(list(list("label_fontface")), list("font of labels")), "\n", "    ", list(list(list("show_network")), list("show underlying spatial network")), 
    "\n", "    ", list(list(list("spatial_network_name")), list("name of spatial network to use")), "\n", "    ", list(list(list("network_color")), list("color of spatial network")), "\n", "    ", list(list(list("network_alpha")), list("alpha of spatial network")), "\n", "    ", list(list(list("show_grid")), list("show spatial grid")), "\n", "    ", list(list(list("spatial_grid_name")), list("name of spatial grid to use")), "\n", "    ", list(list(list("grid_color")), list("color of spatial grid")), 
    "\n", "    ", list(list(list("show_other_cells")), list("display not selected cells")), "\n", "    ", list(list(list("other_cell_color")), list("color of not selected cells")), "\n", "    ", list(list(list("other_point_size")), list("point size of not selected cells")), "\n", "    ", list(list(list("other_cells_alpha")), list("alpha of not selected cells")), "\n", "    ", list(list(list("coord_fix_ratio")), list("fix ratio between x and y-axis (default = 1)")), "\n", "    ", list(list(list("title")), 
        list("title of plot")), "\n", "    ", list(list(list("show_legend")), list("show legend")), "\n", "    ", list(list(list("legend_text")), list("size of legend text")), "\n", "    ", list(list(list("legend_symbol_size")), list("size of legend symbols")), "\n", "    ", list(list(list("background_color")), list("color of plot background")), "\n", "    ", list(list(list("vor_border_color")), list("border colorr for voronoi plot")), "\n", "    ", list(list(list("vor_max_radius")), list("maximum radius for voronoi 'cells'")), 
    "\n", "    ", list(list(list("vor_alpha")), list("transparancy of voronoi 'cells'")), "\n", "    ", list(list(list("axis_text")), list("size of axis text")), "\n", "    ", list(list(list("axis_title")), list("size of axis title")), "\n", "    ", list(list(list("cow_n_col")), list("cowplot param: how many columns")), "\n", "    ", list(list(list("cow_rel_h")), list("cowplot param: relative height")), "\n", "    ", list(list(list("cow_rel_w")), list("cowplot param: relative width")), "\n", "    ", 
    list(list(list("cow_align")), list("cowplot param: how to align")), "\n", "    ", list(list(list("show_plot")), list("show plot")), "\n", "    ", list(list(list("return_plot")), list("return ggplot object")), "\n", "    ", list(list(list("save_plot")), list("directly save the plot [boolean]")), "\n", "    ", list(list(list("save_param")), list("list of saving parameters, see ", list(list("showSaveParameters")))), "\n", "    ", list(list(list("default_save_name")), list("default save name for saving, don't change, change save_name in save_param")), 
    "\n", "  ")


## Details

Description of parameters.


## Value

ggplot


## Seealso

[`spatPlot3D`](#spatplot3d) 
 
 Other spatial visualizations:
 [`spatPlot2D`](#spatplot2d) ,
 [`spatPlot3D`](#spatplot3d)


