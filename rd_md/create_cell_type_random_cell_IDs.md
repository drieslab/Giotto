# `create_cell_type_random_cell_IDs`

create_cell_type_random_cell_IDs


## Description

creates randomized cell ids within a selection of cell types


## Usage

```r
create_cell_type_random_cell_IDs(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  cluster_column = "cell_types",
  needed_cell_types,
  set_seed = FALSE,
  seed_number = 1234
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object to use
`feat_type`     |     feature type
`cluster_column`     |     cluster column with cell type information
`needed_cell_types`     |     vector of cell type names for which a random id will be found
`set_seed`     |     set a seed for reproducibility
`seed_number`     |     seed number


## Value

list of randomly sampled cell ids with same cell type composition


