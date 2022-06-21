# `createGiottoInstructions`

Create instructions for giotto functions


## Description

Function to set global instructions for giotto functions


## Usage

```r
createGiottoInstructions(
  python_path = NULL,
  show_plot = NULL,
  return_plot = NULL,
  save_plot = NULL,
  save_dir = NULL,
  plot_format = NULL,
  dpi = NULL,
  units = NULL,
  height = NULL,
  width = NULL,
  is_docker = FALSE,
  plot_count = 0,
  fiji_path = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`python_path`     |     path to python binary to use
`show_plot`     |     print plot to console, default = TRUE
`return_plot`     |     return plot as object, default = TRUE
`save_plot`     |     automatically save plot, dafault = FALSE
`save_dir`     |     path to directory where to save plots
`plot_format`     |     format of plots (defaults to png)
`dpi`     |     resolution for raster images
`units`     |     units of format (defaults to in)
`height`     |     height of plots
`width`     |     width of  plots
`is_docker`     |     using docker implementation of Giotto (defaults to FALSE)
`plot_count`     |     [global option] start count for creating automatic unique plots
`fiji_path`     |     path to fiji executable


## Value

named vector with giotto instructions


## Seealso

More online information can be found here [https://rubd.github.io/Giotto_site/articles/instructions_and_plotting.html](https://rubd.github.io/Giotto_site/articles/instructions_and_plotting.html)


