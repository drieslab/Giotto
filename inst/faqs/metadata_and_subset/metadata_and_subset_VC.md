
<!-- mouse_cortex_example.md is generated from mouse_cortex_example.Rmd Please edit that file -->

``` r
library(Giotto)
```

    #> Number of cells removed:  0  out of  913 
    #> Number of genes removed:  0  out of  10000

Start from visual cortext Giotto object:

``` r
visPlot(gobject = VC_test)
#> first and second dimenion need to be defined, default is first 2
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="75%" style="display: block; margin: auto;" />

### 1\. add annotation

### 2\. subset Giotto

  - subset first 500 cells

<!-- end list -->

``` r
cell_metadata = pDataDT(VC_test)
first_500_cells = cell_metadata[1:500][['cell_ID']]

VC_subset = subsetGiotto(VC_test, cell_ids = first_500_cells)

visPlot(gobject = VC_subset)
#> first and second dimenion need to be defined, default is first 2
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="75%" style="display: block; margin: auto;" />
