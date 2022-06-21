# `registerGiottoObjectListRvision`

registerGiottoObjectListRvision


## Description

Function to spatially align gobject data based on Rvision image registration.


## Usage

```r
registerGiottoObjectListRvision(
  gobject_list = gobject_list,
  image_list = NULL,
  save_dir = NULL,
  spatloc_unreg = NULL,
  spatloc_reg_name = "raw",
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject_list`     |     list of gobjects to register
`image_list`     |     Filepaths to unregistered images
`save_dir`     |     (Optional) If given, save registered images to this directory
`spatloc_unreg`     |     spatial locations to use
`spatloc_reg_name`     |     name for registered spatial locations to. Defaults to replacement of spat_unreg (optional)
`verbose`     |     be verbose


## Value

list of registered giotto objects where the registered images and spatial locations


