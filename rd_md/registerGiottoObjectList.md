# `registerGiottoObjectList`

registerGiottoObjectList


## Description

Wrapper function for registerGiottoObjectListFiji and registerGiottoObjectListRvision


## Usage

```r
registerGiottoObjectList(
  gobject_list,
  spat_unit = NULL,
  method = c("fiji", "rvision"),
  image_unreg = "image",
  image_reg_name = "image",
  image_list = NULL,
  save_dir = NULL,
  spatloc_unreg = "raw",
  spatloc_reg_name = "raw",
  fiji_xml_files,
  fiji_registered_images,
  scale_factor = NULL,
  allow_rvision_autoscale = TRUE,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject_list`     |     List of gobjects to register
`spat_unit`     |     spatial unit
`method`     |     Method used to align gobjects. Current options are either using FIJI register_virtual_stack_slices output or rvision
`image_unreg`     |     Gobject image slot to use. Defaults to 'image' (optional)
`image_reg_name`     |     Arbitrary image slot name for registered images to occupy. Defaults to replacement of 'image' slot (optional)
`image_list`     |     RVISION - under construction
`save_dir`     |     RVISION - under construction
`spatloc_unreg`     |     Unregistered spatial locations to align. Defaults to 'raw' slot (optional)
`spatloc_reg_name`     |     Arbitrary name for registered spatial locations. Defaults to replacement of 'raw' slot (optional)
`fiji_xml_files`     |     Filepaths to FIJI registration XML outputs
`fiji_registered_images`     |     Registered images output by FIJI register_virtual_stack_slices
`scale_factor`     |     Scaling to be applied to spatial coordinates
`allow_rvision_autoscale`     |     Whether or not to allow rvision to automatically scale the images when performing image registration
`verbose`     |     Be verbose


## Value

List of registered giotto objects where the registered images and spatial locations


