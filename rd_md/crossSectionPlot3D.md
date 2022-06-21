# `crossSectionPlot3D`

crossSectionPlot3D


## Description

Visualize cells in a virtual cross section according to spatial coordinates


## Usage

```r
crossSectionPlot3D(
  gobject,
  crossSection_obj = NULL,
  name = NULL,
  spatial_network_name = "Delaunay_network",
  show_other_cells = T,
  other_cell_color = alpha("lightgrey", 0),
  default_save_name = "crossSection3D",
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`crossSection_obj`     |     cross section object as alternative input. default = NULL.
`name`     |     name of virtual cross section to use
`spatial_network_name`     |     name of spatial network to use
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color of cells outside the cross section. default = transparent.
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param
`...`     |     parameters for spatPlot3D


## Details

Description of parameters.


## Value

ggplot


