# `createSpatialGrid`

createSpatialGrid


## Description

Create a spatial grid using the default method


## Usage

```r
createSpatialGrid(
  gobject,
  spat_unit = NULL,
  spat_loc_name = "raw",
  name = NULL,
  method = c("default"),
  sdimx_stepsize = NULL,
  sdimy_stepsize = NULL,
  sdimz_stepsize = NULL,
  minimum_padding = 1,
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`spat_loc_name`     |     spatial location name
`name`     |     name for spatial grid
`method`     |     method to create a spatial grid
`sdimx_stepsize`     |     stepsize along the x-axis
`sdimy_stepsize`     |     stepsize along the y-axis
`sdimz_stepsize`     |     stepsize along the z-axis
`minimum_padding`     |     minimum padding on the edges
`return_gobject`     |     boolean: return giotto object (default = TRUE)


## Details

Creates a spatial grid with defined x, y (and z) dimensions.
 The dimension units are based on the provided spatial location units.
   

*  default method:  list(list("createSpatialDefaultGrid"))


## Value

giotto object with updated spatial grid slot


