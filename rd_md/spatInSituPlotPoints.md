# `spatInSituPlotPoints`

spatInSituPlotPoints


## Description

Function to plot multiple features for multiple modalities at the spatial in situ level


## Usage

```r
spatInSituPlotPoints(
  gobject,
  show_image = F,
  gimage = NULL,
  image_name = NULL,
  largeImage_name = NULL,
  spat_unit = NULL,
  spat_loc_name = NULL,
  feats = NULL,
  feat_type = "rna",
  feats_color_code = NULL,
  feat_shape_code = NULL,
  sdimx = "x",
  sdimy = "y",
  point_size = 1.5,
  show_polygon = TRUE,
  use_overlap = TRUE,
  polygon_feat_type = "cell",
  polygon_color = "black",
  polygon_bg_color = "black",
  polygon_fill = NULL,
  polygon_fill_gradient = c("blue", "white", "red"),
  polygon_fill_gradient_midpoint = NULL,
  polygon_fill_as_factor = NULL,
  polygon_fill_code = NULL,
  polygon_alpha = 0.5,
  polygon_line_size = 2,
  axis_text = 8,
  axis_title = 8,
  legend_text = 6,
  coord_fix_ratio = NULL,
  background_color = "black",
  show_legend = TRUE,
  plot_method = c("ggplot", "scattermore", "scattermost"),
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spatInSituPlotPoints"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`show_image`     |     show a tissue background image
`gimage`     |     a giotto image
`image_name`     |     name of a giotto image
`largeImage_name`     |     name of a giottoLargeImage
`spat_unit`     |     spatial unit
`spat_loc_name`     |     name of spatial locations
`feats`     |     features to plot
`feat_type`     |     feature types of the feats
`feats_color_code`     |     code to color the provided features
`feat_shape_code`     |     code to shape the provided feature types
`sdimx`     |     spatial dimension x
`sdimy`     |     spatial dimension y
`point_size`     |     size of the points
`show_polygon`     |     overlay polygon information (e.g. cell shape)
`use_overlap`     |     use polygon and feature coordinates overlap results
`polygon_feat_type`     |     feature type associated with polygon information
`polygon_color`     |     color for polygon border
`polygon_bg_color`     |     color for polygon background (overruled by polygon_fill)
`polygon_fill`     |     fill color or column for polygon
`polygon_fill_gradient`     |     polygon fill gradient colors given in order from low to high
`polygon_fill_gradient_midpoint`     |     value to set as gradient midpoint (optional). If left as `NULL` , the median value detected will be chosen
`polygon_fill_as_factor`     |     is fill color a factor
`polygon_fill_code`     |     code to color the fill column
`polygon_alpha`     |     alpha of polygon
`polygon_line_size`     |     line width of the polygon's outline
`axis_text`     |     axis text size
`axis_title`     |     title text size
`legend_text`     |     legend text size
`coord_fix_ratio`     |     fix ratio of coordinates
`background_color`     |     background color
`show_legend`     |     show legend
`plot_method`     |     method to plot points
`show_plot`     |     show plots
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

TODO


## Value

ggplot


## Seealso

Other In Situ visualizations:
 [`spatInSituPlotDensity`](#spatinsituplotdensity) ,
 [`spatInSituPlotHex`](#spatinsituplothex)


