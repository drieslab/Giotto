# `plotly_network`

plotly_network


## Description

provide network segment to draw in 3D plot_ly()


## Usage

```r
plotly_network(
  network,
  x = "sdimx_begin",
  y = "sdimy_begin",
  z = "sdimz_begin",
  x_end = "sdimx_end",
  y_end = "sdimy_end",
  z_end = "sdimz_end"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     network in giotto object


## Value

edges in network as data.table()


