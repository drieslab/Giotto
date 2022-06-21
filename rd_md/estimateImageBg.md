# `estimateImageBg`

estimateImageBg


## Description

helps to estimate which color is the background color of your plot


## Usage

```r
estimateImageBg(mg_object, top_color_range = 1:50)
```


## Arguments

Argument      |Description
------------- |----------------
`mg_object`     |     magick image or Giotto image object
`top_color_range`     |     top possible background colors to return


## Value

vector of pixel color frequencies and an associated barplot


