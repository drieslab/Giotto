# `write_giotto_viewer_dim_reduction`

write_giotto_viewer_dim_reduction


## Description

write out dimensional reduction data from a giotto object for the Viewer


## Usage

```r
write_giotto_viewer_dim_reduction(
  dim_reduction_cell,
  dim_red = NULL,
  dim_red_name = NULL,
  dim_red_rounding = NULL,
  dim_red_rescale = c(-20, 20),
  output_directory = getwd()
)
```


## Arguments

Argument      |Description
------------- |----------------
`dim_reduction_cell`     |     dimension reduction slot from giotto object
`dim_red`     |     high level name of dimension reduction
`dim_red_name`     |     specific name of dimension reduction to use
`dim_red_rounding`     |     numerical indicating how to round the coordinates
`dim_red_rescale`     |     numericals to rescale the coordinates
`output_directory`     |     directory where to save the files


## Value

write a .txt and .annot file for the selection annotation


