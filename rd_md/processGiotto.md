# `processGiotto`

processGiotto


## Description

Wrapper for the different Giotto object processing functions


## Usage

```r
processGiotto(
  gobject,
  filter_params = list(),
  norm_params = list(),
  stat_params = list(),
  adjust_params = list(),
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`filter_params`     |     additional parameters to filterGiotto
`norm_params`     |     additional parameters to normalizeGiotto
`stat_params`     |     additional parameters to addStatistics
`adjust_params`     |     additional parameters to adjustGiottoMatrix
`verbose`     |     be verbose (default is TRUE)


## Details

See [`filterGiotto`](#filtergiotto) , [`normalizeGiotto`](#normalizegiotto) ,
 [`addStatistics`](#addstatistics) , and [`adjustGiottoMatrix`](#adjustgiottomatrix) . For more
 information about the different parameters in each step. If you do not provide
 them it will use the default values.


## Value

giotto object


