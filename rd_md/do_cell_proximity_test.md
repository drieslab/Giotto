# `do_cell_proximity_test`

Do cell proximity test


## Description

Performs a selected differential test on subsets of a matrix


## Usage

```r
do_cell_proximity_test(
  expr_values,
  select_ind,
  other_ind,
  diff_test = c("permutation", "limma", "t.test", "wilcox"),
  mean_method = c("arithmic", "geometric"),
  offset = 0.1,
  n_perm = 100,
  adjust_method = c("bonferroni", "BH", "holm", "hochberg", "hommel", "BY", "fdr",
    "none"),
  set_seed = TRUE,
  seed_number = 1234
)
```


