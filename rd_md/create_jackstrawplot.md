# `create_jackstrawplot`

create_jackstrawplot


## Description

create jackstrawplot with ggplot


## Usage

```r
create_jackstrawplot(
  jackstraw_data,
  ncp = 20,
  ylim = c(0, 1),
  threshold = 0.01
)
```


## Arguments

Argument      |Description
------------- |----------------
`jackstraw_data`     |     result from jackstraw function
`ncp`     |     number of principal components to calculate
`ylim`     |     y-axis limits on jackstraw plot
`threshold`     |     p.value threshold to call a PC significant


## Value

ggplot


