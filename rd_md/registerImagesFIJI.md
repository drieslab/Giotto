# `registerImagesFIJI`

registerImagesFIJI


## Description

Wrapper function for Register Virtual Stack Slices plugin in FIJI


## Usage

```r
registerImagesFIJI(
  source_img_dir,
  output_img_dir,
  transforms_save_dir,
  ref_img_name,
  init_gauss_blur = 1.6,
  steps_per_scale_octave = 3,
  min_img_size = 64,
  max_img_size = 1024,
  feat_desc_size = 8,
  feat_desc_orient_bins = 8,
  closest_next_closest_ratio = 0.92,
  max_align_err = 25,
  inlier_ratio = 0.05,
  headless = FALSE,
  batch = TRUE,
  MinMem = MaxMem,
  MaxMem = 2500,
  IncrementalGC = TRUE,
  Threads = NULL,
  fijiArgs = NULL,
  javaArgs = NULL,
  ijArgs = NULL,
  jython = FALSE,
  fijiPath = fiji(),
  DryRun = FALSE
)
```


## Arguments

Argument      |Description
------------- |----------------
`source_img_dir`     |     Folder containing images to be registered
`output_img_dir`     |     Folder to save registered images to
`transforms_save_dir`     |     (jython implementation only) Folder to save transforms to
`ref_img_name`     |     (jython implementation only) File name of reference image for the registration
`init_gauss_blur`     |     Point detector option: initial image blurring
`steps_per_scale_octave`     |     Point detector option
`min_img_size`     |     Point detector option
`max_img_size`     |     Point detector option
`feat_desc_size`     |     Feature descriptor option
`feat_desc_orient_bins`     |     Feature descriptor option
`closest_next_closest_ratio`     |     Feature descriptor option
`max_align_err`     |     Geometric consensus filter option
`inlier_ratio`     |     Geometric consensus filter option
`headless`     |     Whether to have ImageJ/Fiji running headless #TODO
`batch`     |     Use batch mode #TODO
`MinMem, MaxMem`     |     Memory limits
`IncrementalGC`     |     Whether to use incremental garbage collection
`Threads`     |     Number of threads
`fijiArgs`     |     Arguments for ImageJ/FIJI
`javaArgs`     |     Arguments for Java
`ijArgs`     |     Arguments for ImageJ
`jython`     |     Use jython wrapper script
`fijiPath`     |     Path to fiji executable (can be set by `options(giotto.fiji="/some/path")` )
`DryRun`     |     Whether to return the command to be run rather than actually executing it.


## Details

This function was adapted from runFijiMacro function in jimpipeline by jefferislab


## Value

list of registered giotto objects where the registered images and spatial locations


