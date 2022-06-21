# `plot_feature_points_layer`

plot_feature_points_layer


## Description

low level function to plot a points at the spatial in situ level


## Usage

```r
plot_feature_points_layer(
  ggobject,
  spatial_feat_info,
  feats,
  feats_color_code = NULL,
  feat_shape_code = NULL,
  sdimx = "x",
  sdimy = "y",
  color = "feat_ID",
  shape = "feat",
  point_size = 1.5,
  show_legend = TRUE,
  plot_method = c("ggplot", "scattermore", "scattermost")
)
```


## Details

This function can plot multiple features over multiple modalities. These plots can get very big very fast.


## Value

ggplot


