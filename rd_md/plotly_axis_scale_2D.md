# `plotly_axis_scale_2D`

plotly_axis_scale_2D


## Description

adjust the axis scale in 2D plotly plot


## Usage

```r
plotly_axis_scale_2D(
  cell_locations,
  sdimx = NULL,
  sdimy = NULL,
  mode = c("cube", "real", "custom"),
  custom_ratio = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`cell_locations`     |     spatial_loc in giotto object
`sdimx`     |     x axis of cell spatial location
`sdimy`     |     y axis of cell spatial location
`mode`     |     axis adjustment mode
`custom_ratio`     |     set the ratio artificially


## Value

edges in spatial grid as data.table()


