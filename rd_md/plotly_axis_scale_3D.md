# `plotly_axis_scale_3D`

plotly_axis_scale_3D


## Description

adjust the axis scale in 3D plotly plot


## Usage

```r
plotly_axis_scale_3D(
  cell_locations,
  sdimx = NULL,
  sdimy = NULL,
  sdimz = NULL,
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
`sdimz`     |     z axis of cell spatial location
`mode`     |     axis adjustment mode
`custom_ratio`     |     set the ratio artificially


## Value

edges in spatial grid as data.table()


