# `plotTSNE_3D`

plotTSNE_3D


## Description

Visualize cells according to dimension reduction coordinates


## Usage

```r
plotTSNE_3D(
  gobject,
  dim_reduction_name = "tsne",
  default_save_name = "TSNE_3D",
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`dim_reduction_name`     |     name of TSNE
`default_save_name`     |     default save name of TSNE plot
`...`     |      Arguments passed on to [`dimPlot3D`](#dimplot3d)   list("\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("dim1_to_use")), list("dimension to use on x-axis")), "\n", "    ", list(list(list("dim2_to_use")), list("dimension to use on y-axis")), "\n", "    ", list(list(list("dim3_to_use")), list("dimension to use on z-axis")), "\n", "    ", list(list(list("spat_enr_names")), list("names of spatial enrichment results to include")), "\n", "    ", 
    list(list(list("show_NN_network")), list("show underlying NN network")), "\n", "    ", list(list(list("nn_network_to_use")), list("type of NN network to use (kNN vs sNN)")), "\n", "    ", list(list(list("network_name")), list("name of NN network to use, if show_NN_network = TRUE")), "\n", "    ", list(list(list("cell_color")), list("color for cells (see details)")), "\n", "    ", list(list(list("color_as_factor")), list("convert color column to factor")), "\n", "    ", list(list(list("cell_color_code")), 
        list("named vector with colors")), "\n", "    ", list(list(list("select_cell_groups")), list("select subset of cells/clusters based on cell_color parameter")), "\n", "    ", list(list(list("select_cells")), list("select subset of cells based on cell IDs")), "\n", "    ", list(list(list("show_other_cells")), list("display not selected cells")), "\n", "    ", list(list(list("other_cell_color")), list("color of not selected cells")), "\n", "    ", list(list(list("other_point_size")), list("size of not selected cells")), 
    "\n", "    ", list(list(list("show_cluster_center")), list("plot center of selected clusters")), "\n", "    ", list(list(list("show_center_label")), list("plot label of selected clusters")), "\n", "    ", list(list(list("center_point_size")), list("size of center points")), "\n", "    ", list(list(list("label_size")), list("size of labels")), "\n", "    ", list(list(list("edge_alpha")), list("column to use for alpha of the edges")), "\n", "    ", list(list(list("point_size")), list("size of point (cell)")), 
    "\n", "    ", list(list(list("show_plot")), list("show plot")), "\n", "    ", list(list(list("return_plot")), list("return ggplot object")), "\n", "    ", list(list(list("save_plot")), list("directly save the plot [boolean]")), "\n", "    ", list(list(list("save_param")), list("list of saving parameters, see ", list(list("showSaveParameters")))), "\n", "  ")


## Details

Description of parameters.


## Value

plotly


## Seealso

Other reduced dimension visualizations:
 [`dimPlot2D`](#dimplot2d) ,
 [`dimPlot3D`](#dimplot3d) ,
 [`dimPlot`](#dimplot) ,
 [`plotPCA_2D`](#plotpca2d) ,
 [`plotPCA_3D`](#plotpca3d) ,
 [`plotPCA`](#plotpca) ,
 [`plotTSNE_2D`](#plottsne2d) ,
 [`plotTSNE`](#plottsne) ,
 [`plotUMAP_2D`](#plotumap2d) ,
 [`plotUMAP_3D`](#plotumap3d) ,
 [`plotUMAP`](#plotumap)


