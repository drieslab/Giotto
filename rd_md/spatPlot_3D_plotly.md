# `spatPlot_3D_plotly`

spatPlot_3D_plotly


## Description

Visualize cells at their 3D spatial locations with plotly


## Usage

```r
spatPlot_3D_plotly(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = NULL,
  sdimx = NULL,
  sdimy = NULL,
  sdimz = NULL,
  spat_enr_names = NULL,
  point_size = 3,
  cell_color = NULL,
  cell_color_code = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 0.5,
  show_network = FALSE,
  spatial_network_name = "spatial_network",
  network_color = NULL,
  network_alpha = 1,
  other_cell_alpha = 0.5,
  show_grid = FALSE,
  spatial_grid_name = "spatial_grid",
  title = "",
  show_legend = TRUE,
  axis_scale = c("cube", "real", "custom"),
  custom_ratio = NULL,
  x_ticks = NULL,
  y_ticks = NULL,
  z_ticks = NULL,
  show_plot = FALSE
)
```


## Value

plotly object


