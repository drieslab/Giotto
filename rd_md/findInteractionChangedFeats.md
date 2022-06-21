# `findInteractionChangedFeats`

findInteractionChangedFeats


## Description

Identifies cell-to-cell Interaction Changed Features (ICF),
 i.e. features that are differentially expressed due to proximity to other cell types.#'


## Usage

```r
findInteractionChangedFeats(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  expression_values = "normalized",
  selected_feats = NULL,
  cluster_column,
  spatial_network_name = "Delaunay_network",
  minimum_unique_cells = 1,
  minimum_unique_int_cells = 1,
  diff_test = c("permutation", "limma", "t.test", "wilcox"),
  mean_method = c("arithmic", "geometric"),
  offset = 0.1,
  adjust_method = c("bonferroni", "BH", "holm", "hochberg", "hommel", "BY", "fdr",
    "none"),
  nr_permutations = 1000,
  exclude_selected_cells_from_test = T,
  do_parallel = TRUE,
  set_seed = TRUE,
  seed_number = 1234
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`expression_values`     |     expression values to use
`selected_feats`     |     subset of selected features (optional)
`cluster_column`     |     name of column to use for cell types
`spatial_network_name`     |     name of spatial network to use
`minimum_unique_cells`     |     minimum number of target cells required
`minimum_unique_int_cells`     |     minimum number of interacting cells required
`diff_test`     |     which differential expression test
`mean_method`     |     method to use to calculate the mean
`offset`     |     offset value to use when calculating log2 ratio
`adjust_method`     |     which method to adjust p-values
`nr_permutations`     |     number of permutations if diff_test = permutation
`exclude_selected_cells_from_test`     |     exclude interacting cells other cells
`do_parallel`     |     run calculations in parallel with mclapply
`set_seed`     |     set a seed for reproducibility
`seed_number`     |     seed number


## Details

Function to calculate if genes are differentially expressed in cell types
 when they interact (approximated by physical proximity) with other cell types.
 The results data.table in the cpgObject contains - at least - the following columns:
   

*  features:  All or selected list of tested features   

*  sel:  average feature expression in the interacting cells from the target cell type    

*  other:  average feature expression in the NOT-interacting cells from the target cell type    

*  log2fc:  log2 fold-change between sel and other   

*  diff:  spatial expression difference between sel and other   

*  p.value:  associated p-value   

*  p.adj:  adjusted p-value   

*  cell_type:  target cell type   

*  int_cell_type:  interacting cell type   

*  nr_select:  number of cells for selected target cell type   

*  int_nr_select:  number of cells for interacting cell type   

*  nr_other:  number of other cells of selected target cell type   

*  int_nr_other:  number of other cells for interacting cell type   

*  unif_int:  cell-cell interaction


## Value

cpgObject that contains the Interaction Changed differential gene scores


