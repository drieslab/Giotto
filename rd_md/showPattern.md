# `showPattern`

showPattern


## Description

show patterns for 2D spatial data


## Usage

```r
showPattern(gobject, spatPatObj, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spatPatObj`     |     Output from detectSpatialPatterns
`...`     |      Arguments passed on to [`showPattern2D`](#showpattern2d)   list("\n", "    ", list(list(list("dimension")), list("dimension to plot")), "\n", "    ", list(list(list("trim")), list("Trim ends of the PC values.")), "\n", "    ", list(list(list("background_color")), list("background color for plot")), "\n", "    ", list(list(list("grid_border_color")), list("color for grid")), "\n", "    ", list(list(list("show_legend")), list("show legend of ggplot")), "\n", "    ", list(list(list("point_size")), list("size of points")), "\n", "    ", list(list(list("show_plot")), 
    list("show plot")), "\n", "    ", list(list(list("return_plot")), list("return ggplot object")), "\n", "    ", list(list(list("save_plot")), list("directly save the plot [boolean]")), "\n", "    ", list(list(list("save_param")), list("list of saving parameters, see ", list(list("showSaveParameters")))), "\n", "    ", list(list(list("default_save_name")), list("default save name for saving, don't change, change save_name in save_param")), "\n", "  ")


## Value

ggplot


## Seealso

[`showPattern2D`](#showpattern2d)


