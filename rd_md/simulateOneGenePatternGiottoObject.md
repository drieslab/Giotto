# `simulateOneGenePatternGiottoObject`

simulateOneGenePatternGiottoObject


## Description

Create a simulated spatial pattern for one selected gnee


## Usage

```r
simulateOneGenePatternGiottoObject(
  gobject,
  pattern_name = "pattern",
  pattern_cell_ids = NULL,
  gene_name = NULL,
  spatial_prob = 0.95,
  gradient_direction = NULL,
  show_pattern = TRUE,
  pattern_colors = c(`in` = "green", out = "red"),
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`pattern_name`     |     name of spatial pattern
`pattern_cell_ids`     |     cell ids that make up the spatial pattern
`gene_name`     |     selected gene
`spatial_prob`     |     probability for a high expressing gene value to be part of the spatial pattern
`gradient_direction`     |     direction of gradient
`show_pattern`     |     show the discrete spatial pattern
`pattern_colors`     |     2 color vector for the spatial pattern
`list()`     |     additional parameters for (re-)normalizing


## Value

Reprocessed Giotto object for which one gene has a forced spatial pattern


