# `registerGiottoObjectListFiji`

registerGiottoObjectListFiji


## Description

Function to spatially align gobject data based on FIJI image registration.


## Usage

```r
registerGiottoObjectListFiji(
  gobject_list,
  spat_unit = NULL,
  image_unreg = "image",
  image_reg_name = "image",
  image_replace_name = "unregistered",
  registered_images = NULL,
  spatloc_unreg = "raw",
  spatloc_reg_name = "raw",
  spatloc_replace_name = "unregistered",
  xml_files,
  scale_factor = NULL,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject_list`     |     list of gobjects to register
`spat_unit`     |     spatial unit
`image_unreg`     |     name of original unregistered images. Defaults to 'image' (optional)
`image_reg_name`     |     arbitrary name for registered images to occupy. Defaults to replacement of 'image' (optional)
`image_replace_name`     |     arbitrary name for any images replaced due to image_reg_name argument (optional)
`registered_images`     |     registered images output by FIJI register_virtual_stack_slices
`spatloc_unreg`     |     spatial locations to use. Defaults to 'raw' (optional)
`spatloc_reg_name`     |     name for registered spatial locations. Defaults to replacement of 'raw' (optional)
`spatloc_replace_name`     |     arbitrary name for any spatial locations replaced due to spatloc_reg_name argument (optional)
`xml_files`     |     atomic vector of filepaths to xml outputs from FIJI register_virtual_stack_slices
`scale_factor`     |     vector of scaling factors of images used in registration vs spatlocs
`verbose`     |     be verbose


## Value

list of registered giotto objects where the registered images and spatial locations


