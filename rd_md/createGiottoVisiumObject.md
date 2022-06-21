# `createGiottoVisiumObject`

Create a giotto object from 10x visium data


## Description

creates Giotto object directly from a 10X visium folder


## Usage

```r
createGiottoVisiumObject(
  visium_dir = NULL,
  expr_data = c("raw", "filter"),
  gene_column_index = 1,
  h5_visium_path = NULL,
  h5_gene_ids = c("symbols", "ensembl"),
  h5_tissue_positions_path = NULL,
  h5_image_png_path = NULL,
  h5_json_scalefactors_path = NULL,
  png_name = NULL,
  do_manual_adj = FALSE,
  xmax_adj = 0,
  xmin_adj = 0,
  ymax_adj = 0,
  ymin_adj = 0,
  instructions = NULL,
  cores = NA,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`visium_dir`     |     path to the 10X visium directory [required]
`expr_data`     |     raw or filtered data (see details)
`gene_column_index`     |     which column index to select (see details)
`h5_visium_path`     |     path to visium 10X .h5 file
`h5_gene_ids`     |     gene names as symbols (default) or ensemble gene ids
`h5_tissue_positions_path`     |     path to tissue locations (.csv file)
`h5_image_png_path`     |     path to tissue .png file (optional)
`h5_json_scalefactors_path`     |     path to .json scalefactors (optional)
`png_name`     |     select name of png to use (see details)
`do_manual_adj`     |     flag to use manual adj values instead of automatic image alignment
`xmax_adj`     |     adjustment of the maximum x-value to align the image
`xmin_adj`     |     adjustment of the minimum x-value to align the image
`ymax_adj`     |     adjustment of the maximum y-value to align the image
`ymin_adj`     |     adjustment of the minimum y-value to align the image
`instructions`     |     list of instructions or output result from [`createGiottoInstructions`](#creategiottoinstructions)
`cores`     |     how many cores or threads to use to read data if paths are provided
`verbose`     |     be verbose


## Details

If starting from a Visium 10X directory:
   

*  expr_data: raw will take expression data from raw_feature_bc_matrix and filter from filtered_feature_bc_matrix   

*  gene_column_index: which gene identifiers (names) to use if there are multiple columns (e.g. ensemble and gene symbol)   

*  png_name: by default the first png will be selected, provide the png name to override this (e.g. myimage.png)   

*  the file scalefactors_json.json will be detected automaticated and used to attempt to align the data  
 
 If starting from a Visium 10X .h5 file
   

*  h5_visium_path: full path to .h5 file: /your/path/to/visium_file.h5   

*  h5_tissue_positions_path: full path to spatial locations file: /you/path/to/tissue_positions_list.csv   

*  h5_image_png_path: full path to png: /your/path/to/images/tissue_lowres_image.png   

*  h5_json_scalefactors_path: full path to .json file: /your/path/to/scalefactors_json.json


## Value

giotto object


