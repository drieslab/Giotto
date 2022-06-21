# `spatPlot_2D_plotly`

spatPlot_2D_plotly


## Description

Visualize cells at their 2D spatial locations with plotly


## Usage

```r
spatPlot_2D_plotly(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = NULL,
  sdimx = NULL,
  sdimy = NULL,
  spat_enr_names = NULL,
  point_size = 3,
  cell_color = NULL,
  cell_color_code = NULL,
  color_as_factor = T,
  select_cell_groups = NULL,
  select_cells = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 0.5,
  show_network = FALSE,
  spatial_network_name = "spatial_network",
  network_color = "lightgray",
  network_alpha = 1,
  other_cell_alpha = 0.5,
  show_grid = FALSE,
  spatial_grid_name = "spatial_grid",
  grid_color = NULL,
  grid_alpha = 1,
  show_legend = T,
  axis_scale = c("cube", "real", "custom"),
  custom_ratio = NULL,
  x_ticks = NULL,
  y_ticks = NULL,
  show_plot = F
)
```


## Value

plotly object


