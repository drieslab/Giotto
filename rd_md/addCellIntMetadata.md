# `addCellIntMetadata`

addCellIntMetadata


## Description

Creates an additional metadata column with information about interacting and non-interacting cell types of the
 selected cell-cell interaction.


## Usage

```r
addCellIntMetadata(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  spatial_network = "spatial_network",
  cluster_column,
  cell_interaction,
  name = "select_int",
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`spatial_network`     |     name of spatial network to use
`cluster_column`     |     column of cell types
`cell_interaction`     |     cell-cell interaction to use
`name`     |     name for the new metadata column
`return_gobject`     |     return an updated giotto object


## Details

This function will create an additional metadata column which selects interacting cell types for a specific cell-cell
 interaction. For example, if you want to color interacting astrocytes and oligodendrocytes it will create a new metadata column with
 the values "select_astrocytes", "select_oligodendrocytes", "other_astrocytes", "other_oligodendroyctes" and "other". Where "other" is all
 other cell types found within the selected cell type column.


## Value

Giotto object


