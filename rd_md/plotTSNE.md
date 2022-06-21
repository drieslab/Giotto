# `plotTSNE`

plotTSNE


## Description

Short wrapper for tSNE visualization


## Usage

```r
plotTSNE(gobject, dim_reduction_name = NULL, default_save_name = "tSNE", ...)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`dim_reduction_name`     |     name of TSNE
`default_save_name`     |     default save name of TSNE plot
`...`     |      Arguments passed on to [`dimPlot2D`](#dimplot2d)   list("\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("group_by")), list("create multiple plots based on cell annotation column")), "\n", "    ", list(list(list("group_by_subset")), list("subset the group_by factor column")), "\n", "    ", list(list(list("dim1_to_use")), list("dimension to use on x-axis")), "\n", "    ", list(list(list("dim2_to_use")), list("dimension to use on y-axis")), 
    "\n", "    ", list(list(list("spat_enr_names")), list("names of spatial enrichment results to include")), "\n", "    ", list(list(list("show_NN_network")), list("show underlying NN network")), "\n", "    ", list(list(list("nn_network_to_use")), list("type of NN network to use (kNN vs sNN)")), "\n", "    ", list(list(list("network_name")), list("name of NN network to use, if show_NN_network = TRUE")), "\n", "    ", list(list(list("cell_color")), list("color for cells (see details)")), "\n", "    ", 
    list(list(list("color_as_factor")), list("convert color column to factor")), "\n", "    ", list(list(list("cell_color_code")), list("named vector with colors")), "\n", "    ", list(list(list("cell_color_gradient")), list("vector with 3 colors for numeric data")), "\n", "    ", list(list(list("gradient_midpoint")), list("midpoint for color gradient")), "\n", "    ", list(list(list("gradient_limits")), list("vector with lower and upper limits")), "\n", "    ", list(list(list("select_cell_groups")), 
        list("select subset of cells/clusters based on cell_color parameter")), "\n", "    ", list(list(list("select_cells")), list("select subset of cells based on cell IDs")), "\n", "    ", list(list(list("show_other_cells")), list("display not selected cells")), "\n", "    ", list(list(list("other_cell_color")), list("color of not selected cells")), "\n", "    ", list(list(list("other_point_size")), list("size of not selected cells")), "\n", "    ", list(list(list("show_cluster_center")), list(
        "plot center of selected clusters")), "\n", "    ", list(list(list("show_center_label")), list("plot label of selected clusters")), "\n", "    ", list(list(list("center_point_size")), list("size of center points")), "\n", "    ", list(list(list("center_point_border_col")), list("border color of center points")), "\n", "    ", list(list(list("center_point_border_stroke")), list("border stroke size of center points")), "\n", "    ", list(list(list("label_size")), list("size of labels")), "\n", 
    "    ", list(list(list("label_fontface")), list("font of labels")), "\n", "    ", list(list(list("edge_alpha")), list("column to use for alpha of the edges")), "\n", "    ", list(list(list("point_shape")), list("point with border or not (border or no_border)")), "\n", "    ", list(list(list("point_size")), list("size of point (cell)")), "\n", "    ", list(list(list("point_alpha")), list("transparancy of point")), "\n", "    ", list(list(list("point_border_col")), list("color of border around points")), 
    "\n", "    ", list(list(list("point_border_stroke")), list("stroke size of border around points")), "\n", "    ", list(list(list("title")), list("title for plot, defaults to cell_color parameter")), "\n", "    ", list(list(list("show_legend")), list("show legend")), "\n", "    ", list(list(list("legend_text")), list("size of legend text")), "\n", "    ", list(list(list("legend_symbol_size")), list("size of legend symbols")), "\n", "    ", list(list(list("background_color")), list("color of plot background")), 
    "\n", "    ", list(list(list("axis_text")), list("size of axis text")), "\n", "    ", list(list(list("axis_title")), list("size of axis title")), "\n", "    ", list(list(list("cow_n_col")), list("cowplot param: how many columns")), "\n", "    ", list(list(list("cow_rel_h")), list("cowplot param: relative height")), "\n", "    ", list(list(list("cow_rel_w")), list("cowplot param: relative width")), "\n", "    ", list(list(list("cow_align")), list("cowplot param: how to align")), "\n", "    ", 
    list(list(list("show_plot")), list("show plot")), "\n", "    ", list(list(list("return_plot")), list("return ggplot object")), "\n", "    ", list(list(list("save_plot")), list("directly save the plot [boolean]")), "\n", "    ", list(list(list("save_param")), list("list of saving parameters, see ", list(list("showSaveParameters")))), "\n", "  ")


## Details

Description of parameters, see [`dimPlot2D`](#dimplot2d) . For 3D plots see [`plotTSNE_3D`](#plottsne3d)


## Value

ggplot


## Seealso

Other reduced dimension visualizations:
 [`dimPlot2D`](#dimplot2d) ,
 [`dimPlot3D`](#dimplot3d) ,
 [`dimPlot`](#dimplot) ,
 [`plotPCA_2D`](#plotpca2d) ,
 [`plotPCA_3D`](#plotpca3d) ,
 [`plotPCA`](#plotpca) ,
 [`plotTSNE_2D`](#plottsne2d) ,
 [`plotTSNE_3D`](#plottsne3d) ,
 [`plotUMAP_2D`](#plotumap2d) ,
 [`plotUMAP_3D`](#plotumap3d) ,
 [`plotUMAP`](#plotumap)


