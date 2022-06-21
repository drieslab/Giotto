# `findCPG`

findCPG


## Description

Identifies cell-to-cell Interaction Changed Genes (ICG),
 i.e. genes that are differentially expressed due to proximity to other cell types.


## Usage

```r
findCPG(...)
```


## Arguments

Argument      |Description
------------- |----------------
`...`     |      Arguments passed on to [`findICG`](#findicg)   list("\n", "    ", list(list(list("gobject")), list("giotto object")), "\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("expression_values")), list("expression values to use")), "\n", "    ", list(list(list("selected_genes")), list("subset of selected genes (optional)")), "\n", "    ", list(list(list("cluster_column")), list("name of column to use for cell types")), "\n", "    ", list(
    list(list("spatial_network_name")), list("name of spatial network to use")), "\n", "    ", list(list(list("minimum_unique_cells")), list("minimum number of target cells required")), "\n", "    ", list(list(list("minimum_unique_int_cells")), list("minimum number of interacting cells required")), "\n", "    ", list(list(list("diff_test")), list("which differential expression test")), "\n", "    ", list(list(list("mean_method")), list("method to use to calculate the mean")), "\n", "    ", list(list(
    list("offset")), list("offset value to use when calculating log2 ratio")), "\n", "    ", list(list(list("adjust_method")), list("which method to adjust p-values")), "\n", "    ", list(list(list("nr_permutations")), list("number of permutations if diff_test = permutation")), "\n", "    ", list(list(list("exclude_selected_cells_from_test")), list("exclude interacting cells other cells")), "\n", "    ", list(list(list("do_parallel")), list("run calculations in parallel with mclapply")), "\n", "    ", 
    list(list(list("set_seed")), list("set a seed for reproducibility")), "\n", "    ", list(list(list("seed_number")), list("seed number")), "\n", "  ")


## Seealso

[`findICG`](#findicg)


