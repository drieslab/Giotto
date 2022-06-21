# `joinGiottoObjects`

joinGiottoObjects


## Description

Function to join multiple giotto objects together


## Usage

```r
joinGiottoObjects(
  gobject_list,
  gobject_names = NULL,
  join_method = c("shift", "z_stack"),
  z_vals = 1000,
  x_shift = NULL,
  y_shift = NULL,
  x_padding = 0,
  y_padding = 0,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject_list`     |     list of giotto objects
`gobject_names`     |     unique giotto names for each giotto object
`join_method`     |     method to join giotto objects
`z_vals`     |     distance(s) along z-axis if method is z-stack
`x_shift`     |     shift along x-axis if method is shift
`y_shift`     |     shift along y-axis if method is shift
`x_padding`     |     padding between datasets/images if method is shift
`y_padding`     |     padding between datasets/images if method is shift
`verbose`     |     be verbose


## Value

giotto object


