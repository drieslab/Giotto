# `annotateGiotto`

Annotate giotto clustering


## Description

Converts cluster results into a user provided annotation.


## Usage

```r
annotateGiotto(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  annotation_vector = NULL,
  cluster_column = NULL,
  name = "cell_types"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`annotation_vector`     |     named annotation vector (names = cluster ids)
`cluster_column`     |     cluster column to convert to annotation names
`name`     |     new name for annotation column


## Details

You need to specifify which (cluster) column you want to annotate
 and you need to provide an annotation vector like this:
   

*  1. identify the cell type of each cluster   

*  2. create a vector of these cell types, e.g. cell_types =  c('T-cell', 'B-cell', 'Stromal')   

*  3. provide original cluster names to previous vector, e.g. names(cell_types) = c(2, 1, 3)


## Value

giotto object


