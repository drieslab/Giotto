# `plot_spat_scatterpie_layer_ggplot`

plot_spat_scatterpie_layer_ggplot


## Description

create scatterpie in ggplot


## Usage

```r
plot_spat_scatterpie_layer_ggplot(
  ggobject,
  sdimx = "sdimx",
  sdimy = "sdimy",
  spatial_locations = NULL,
  spatial_enrichment = NULL,
  radius = 10,
  color = NA,
  alpha = 1,
  cell_color_code = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`sdimx`     |     x-axis dimension name (default = 'sdimx')
`sdimy`     |     y-axis dimension name (default = 'sdimy')
`spatial_locations`     |     spatial locations
`spatial_enrichment`     |     spatial enrichment results
`radius`     |     radius of scatterpie
`color`     |     color of lines within pie chart
`cell_color_code`     |     color code for the cell types
`gobject`     |     giotto object


## Value

ggplot


