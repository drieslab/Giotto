# `createSpatialEnrich`

createSpatialEnrich


## Description

Function to calculate gene signature enrichment scores per spatial position using an enrichment test.


## Usage

```r
createSpatialEnrich(...)
```


## Arguments

Argument      |Description
------------- |----------------
`...`     |      Arguments passed on to [`runSpatialEnrich`](#runspatialenrich)   list("\n", "    ", list(list(list("gobject")), list("Giotto object")), "\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("enrich_method")), list("method for gene signature enrichment calculation")), "\n", "    ", list(list(list("sign_matrix")), list("Matrix of signature genes for each cell type / process")), "\n", "    ", list(list(list("expression_values")), list("expression values to use")), 
    "\n", "    ", list(list(list("reverse_log_scale")), list("reverse expression values from log scale")), "\n", "    ", list(list(list("min_overlap_genes")), list("minimum number of overlapping genes in sign_matrix required to calculate enrichment (PAGE)")), "\n", "    ", list(list(list("logbase")), list("log base to use if reverse_log_scale = TRUE")), "\n", "    ", list(list(list("p_value")), list("calculate p-value (default = FALSE)")), "\n", "    ", list(list(list("n_times")), list("(page/rank) number of permutation iterations to calculate p-value")), 
    "\n", "    ", list(list(list("rbp_p")), list("(rank) fractional binarization threshold (default = 0.99)")), "\n", "    ", list(list(list("num_agg")), list("(rank) number of top genes to aggregate (default = 100)")), "\n", "    ", list(list(list("max_block")), list("number of lines to process together (default = 20e6)")), "\n", "    ", list(list(list("top_percentage")), list("(hyper) percentage of cells that will be considered to have gene expression with matrix binarization")), "\n", "    ", 
    list(list(list("output_enrichment")), list("how to return enrichment output")), "\n", "    ", list(list(list("name")), list("to give to spatial enrichment results, default = PAGE")), "\n", "    ", list(list(list("verbose")), list("be verbose")), "\n", "    ", list(list(list("return_gobject")), list("return giotto object")), "\n", "  ")


## Seealso

[`runSpatialEnrich`](#runspatialenrich)


