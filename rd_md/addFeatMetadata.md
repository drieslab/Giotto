# `addFeatMetadata`

Add feature metadata


## Description

Adds feature metadata to the giotto object


## Usage

```r
addFeatMetadata(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  new_metadata,
  by_column = F,
  column_feat_ID = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`new_metadata`     |     new metadata to use
`by_column`     |     merge metadata based on feat_ID column in [`fDataDT`](#fdatadt)
`column_feat_ID`     |     column name of new metadata to use if by_column = TRUE


## Details

You can add additional feature metadata in two manners: list() 
 1. Provide a data.table or data.frame with feature annotations in the same order as the list("feat_ID") column in fDataDT(gobject) list() 
 2. Provide a data.table or data.frame with feature annotations and specify which column contains the feature IDs,
 these feature IDs need to match with the list("feat_ID") column in fDataDT(gobject)


## Value

giotto object


