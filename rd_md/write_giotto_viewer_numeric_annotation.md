# `write_giotto_viewer_numeric_annotation`

write_giotto_viewer_numeric_annotation


## Description

write out numeric annotation data from a giotto object for the Viewer


## Usage

```r
write_giotto_viewer_numeric_annotation(
  annotation,
  annot_name = "test",
  output_directory = getwd()
)
```


## Arguments

Argument      |Description
------------- |----------------
`annotation`     |     annotation from the data.table from giotto object
`annot_name`     |     name of the annotation
`output_directory`     |     directory where to save the files


## Value

write a .txt and .annot file for the selection annotation


