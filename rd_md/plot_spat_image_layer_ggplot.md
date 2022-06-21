# `plot_spat_image_layer_ggplot`

plot_spat_image_layer_ggplot


## Description

create background image in ggplot


## Usage

```r
plot_spat_image_layer_ggplot(
  ggplot,
  gobject,
  gimage,
  feat_type = NULL,
  spat_unit = NULL,
  spat_loc_name = NULL,
  sdimx = NULL,
  sdimy = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`gimage`     |     a giotto image or a list/vector of giotto images
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`spat_loc_name`     |     name for spatial locations
`sdimx`     |     x-axis dimension name (default = 'sdimx')
`sdimy`     |     y-axis dimension name (default = 'sdimy')


## Value

ggplot


