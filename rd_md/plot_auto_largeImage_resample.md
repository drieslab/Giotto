# `plot_auto_largeImage_resample`

plot_auto_largeImage_resample


## Description

resamples largeImage for plotting


## Usage

```r
plot_auto_largeImage_resample(
  gobject,
  giottoLargeImage = NULL,
  largeImage_name = NULL,
  spat_unit = NULL,
  spat_loc_name = NULL,
  include_image_in_border = TRUE,
  flex_resample = TRUE,
  max_crop = 1e+08,
  max_resample_scale = 100
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     gobject containing largeImage
`giottoLargeImage`     |     giottoLargeImage
`largeImage_name`     |     name of largeImage
`spat_unit`     |     spatial unit
`spat_loc_name`     |     name of spatial locations to plot
`include_image_in_border`     |     expand the extent sampled to also show image in border regions not included in spatlocs
`flex_resample`     |     decide to crop then resample or resample then crop depending on plotting ROI dimensions
`max_crop`     |     maximum crop size allowed for crop THEN resample
`max_resampl_scale`     |     maximum cells allowed to resample to compensate for decreased resolution after subsequent cropping


## Value

a giottoLargeImage cropped and resampled properly for plotting


