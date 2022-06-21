# `addCellMetadata`

Add cell metadata


## Description

Adds cell metadata to the giotto object


## Usage

```r
addCellMetadata(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  new_metadata,
  vector_name = NULL,
  by_column = FALSE,
  column_cell_ID = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`new_metadata`     |     new cell metadata to use (data.table, data.frame, ...)
`vector_name`     |     (optional) custom name if you provide a single vector
`by_column`     |     merge metadata based on cell_ID column in [`pDataDT`](#pdatadt) (default = FALSE)
`column_cell_ID`     |     column name of new metadata to use if by_column = TRUE


## Details

You can add additional cell metadata in two manners:
   

*  list("1. Provide a data.table or data.frame with cell annotations in the same order as the ", list("cell_ID"), " column in pDataDT(gobject) ")   

*  list("2. Provide a data.table or data.frame with cell annotations and specify which column contains the cell IDs, these cell IDs need to match with the ", list("cell_ID"), " column in pDataDT(gobject)")


## Value

giotto object


