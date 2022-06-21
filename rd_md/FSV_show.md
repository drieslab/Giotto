# `FSV_show`

FSV_show


## Description

Visualize spatial varible genes caculated by spatial_DE


## Usage

```r
FSV_show(
  results,
  ms_results = NULL,
  size = c(4, 2, 1),
  color = c("blue", "green", "red"),
  sig_alpha = 0.5,
  unsig_alpha = 0.5
)
```


## Arguments

Argument      |Description
------------- |----------------
`results`     |     results caculated by spatial_DE
`ms_results`     |     ms_results caculated by spatial_DE
`size`     |     indicate different levels of qval
`color`     |     indicate different SV features
`sig_alpha`     |     transparency of significant genes
`unsig_alpha`     |     transparency of unsignificant genes


## Details

Description of parameters.


## Value

ggplot object


