# `addHMRF`

addHMRF


## Description

Add selected results from doHMRF to the giotto object


## Usage

```r
addHMRF(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  HMRFoutput,
  k = NULL,
  betas_to_add = NULL,
  hmrf_name = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`HMRFoutput`     |     HMRF output from doHMRF()
`k`     |     number of domains
`betas_to_add`     |     results from different betas that you want to add
`hmrf_name`     |     specify a custom name


## Value

giotto object


