# `plotInteractivePolygons`

Select image regions by plotting interactive polygons


## Description

Plot interactive polygons on an image and retrieve the polygons coordinates.


## Usage

```r
plotInteractivePolygons(x, width = "auto", height = "auto", ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     A `ggplot` or `rast` plot object to draw polygons on
`width, height`     |     An integer, defining the width/height in pixels.
`...`     |     Graphical parameters passed on to `polygon` or `geom_point`.


## Value

A `data.table` containing x,y coordinates from the plotted polygons.


## Examples

```r
# Using a ggplot2 plot
library(ggplot2)
df <- data.frame(x = 1:5, y = 1:5)
my_plot <- ggplot(df, aes(x,y)) + geom_point()
plotInteractivePolygons(my_plot)

# Using a terra rast image
library(terra)
r = rast(system.file("ex/elev.tif", package="terra"))
plotInteractivePolygons(r)
plotInteractivePolygons(r, border = "red", lwd = 2)

# Using an image contained in Giotto object
library(Giotto)
my_spatPlot <- spatPlot2D(gobject = my_giotto_object,
show_image = TRUE,
point_alpha = 0.75,
save_plot = FALSE)
plotInteractivePolygons(my_spatPlot, height = 500)
my_polygon_coordinates <- plotInteractivePolygons(my_spatPlot, height = 500)
```


