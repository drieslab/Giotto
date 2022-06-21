# `do_permuttest_random`

do_permuttest_random


## Description

calculate random values
 
 Performs permutation test on subsets of a matrix


## Usage

```r
do_permuttest_random(
  expr_values,
  select_ind,
  other_ind,
  name = "perm_1",
  mean_method,
  offset = 0.1,
  set_seed = TRUE,
  seed_number = 1234
)
do_permuttest(
  expr_values,
  select_ind,
  other_ind,
  n_perm = 1000,
  adjust_method = "fdr",
  mean_method,
  offset = 0.1,
  set_seed = TRUE,
  seed_number = 1234
)
```


