# `all_plots_save_function`

all_plots_save_function


## Description

Function to automatically save plots to directory of interest


## Usage

```r
all_plots_save_function(
  gobject,
  plot_object,
  save_dir = NULL,
  save_folder = NULL,
  save_name = NULL,
  default_save_name = "giotto_plot",
  save_format = NULL,
  show_saved_plot = F,
  ncol = 1,
  nrow = 1,
  scale = 1,
  base_width = NULL,
  base_height = NULL,
  base_aspect_ratio = NULL,
  units = NULL,
  dpi = NULL,
  limitsize = TRUE,
  plot_count = NULL,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`plot_object`     |     object to plot
`save_dir`     |     directory to save to
`save_folder`     |     folder in save_dir to save to
`save_name`     |     name of plot
`default_save_name`     |     default name to save a plot
`save_format`     |     format (e.g. png, tiff, pdf, ...)
`show_saved_plot`     |     load & display the saved plot
`ncol`     |     number of columns
`nrow`     |     number of rows
`scale`     |     scale
`base_width`     |     width
`base_height`     |     height
`base_aspect_ratio`     |     aspect ratio
`units`     |     units
`dpi`     |     Plot resolution
`limitsize`     |     When TRUE (the default), ggsave will not save images larger than 50x50 inches, to prevent the common error of specifying dimensions in pixels.
`plot_count`     |     count number for plot
`list()`     |     additional parameters to ggplot_save_function or general_save_function


## Seealso

[`general_save_function`](#generalsavefunction)


