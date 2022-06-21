# `changeImageBg`

changeImageBg


## Description

Function to change the background color of a magick image plot to another color


## Usage

```r
changeImageBg(
  mg_object,
  bg_color,
  perc_range = 10,
  new_color = "#FFFFFF",
  new_name = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`mg_object`     |     magick image or giotto image object
`bg_color`     |     estimated current background color
`perc_range`     |     range around estimated background color to include (percentage)
`new_color`     |     new background color
`new_name`     |     change name of Giotto image


## Value

magick image or giotto image object with updated background color


