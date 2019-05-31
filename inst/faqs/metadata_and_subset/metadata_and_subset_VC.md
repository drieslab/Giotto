
<!-- mouse_cortex_example.md is generated from mouse_cortex_example.Rmd Please edit that file -->

Start from Giotto object from the Mouse Visual Cortex and SVZ:

``` r
visPlot(gobject = VC_test, point_size = 2)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="75%" style="display: block; margin: auto;" />

### 1\. add annotation

``` r
# cell metadata before 
print(pDataDT(VC_test))
#>       cell_ID nr_genes perc_genes total_expr
#>   1:   cell_1     5948      59.48  2367.0630
#>   2:   cell_2     2945      29.45  2118.8796
#>   3:   cell_3     3953      39.53  2243.4318
#>   4:   cell_4     1429      14.29  1721.6024
#>   5:   cell_5     3189      31.89  2148.4464
#>  ---                                        
#> 909: cell_909     1747      17.47  1820.9809
#> 910: cell_910     1292      12.92  1668.5748
#> 911: cell_911      300       3.00   861.8407
#> 912: cell_912     1267      12.67  1650.0913
#> 913: cell_913      179       1.79   632.5006

# add metadata by column name
cortex_annotation = fread(system.file("extdata", "cortex_annotation.txt", package = "Giotto"))
VC_test = addCellMetadata(gobject = VC_test, new_metadata = cortex_annotation, by_column = T, column_cell_ID = 'uniq_ID')

# cell metadata after
cell_metadata = pDataDT(VC_test)
print(cell_metadata)
#>       cell_ID nr_genes perc_genes total_expr Field_of_View
#>   1:   cell_1     5948      59.48  2367.0630             0
#>   2:   cell_2     2945      29.45  2118.8796             0
#>   3:   cell_3     3953      39.53  2243.4318             0
#>   4:   cell_4     1429      14.29  1721.6024             0
#>   5:   cell_5     3189      31.89  2148.4464             0
#>  ---                                                      
#> 909: cell_909     1747      17.47  1820.9809             6
#> 910: cell_910     1292      12.92  1668.5748             6
#> 911: cell_911      300       3.00   861.8407             6
#> 912: cell_912     1267      12.67  1650.0913             6
#> 913: cell_913      179       1.79   632.5006             6
#>             cell_types
#>   1: Excitatory neuron
#>   2: Excitatory neuron
#>   3: Excitatory neuron
#>   4:       Interneuron
#>   5: Excitatory neuron
#>  ---                  
#> 909:    Choroid Plexus
#> 910:    Choroid Plexus
#> 911:       Endothelial
#> 912:    Choroid Plexus
#> 913:         Microglia
```

### 2\. subset Giotto

  - subset first 500 cells

<!-- end list -->

``` r
first_500_cells = cell_metadata[1:500][['cell_ID']]
VC_subset = subsetGiotto(VC_test, cell_ids = first_500_cells)
visPlot(gobject = VC_subset, point_size = 2, cell_color = 'cell_types')
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="75%" style="display: block; margin: auto;" />

  - subset for SVZ cells (located in field 5)

<!-- end list -->

``` r
SVZ_cells = cell_metadata[Field_of_View == 5][['cell_ID']]
VC_SVZ = subsetGiotto(VC_test, cell_ids = SVZ_cells)
visPlot(gobject = VC_SVZ, point_size = 4, cell_color = 'cell_types')
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="75%" style="display: block; margin: auto;" />

  - subset for all interneurons

<!-- end list -->

``` r
interneuron_cells = cell_metadata[cell_types == 'Interneuron'][['cell_ID']]
VC_interneuron = subsetGiotto(VC_test, cell_ids = interneuron_cells)
visPlot(gobject = VC_interneuron, point_size = 2, cell_color = 'cell_types')
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="75%" style="display: block; margin: auto;" />
