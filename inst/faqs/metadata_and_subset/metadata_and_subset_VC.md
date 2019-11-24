
<!-- mouse_cortex_example.md is generated from mouse_cortex_example.Rmd Please edit that file -->

Start from Giotto object from the Mouse Visual Cortex and SVZ:

``` r
spatPlot2D(gobject = VC_test, point_size = 2, return_plot = F)
#> create 2D plot with ggplot
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="75%" style="display: block; margin: auto;" />

### 1\. add annotation

``` r
# cell metadata before 
print(pDataDT(VC_test))
#>       cell_ID nr_genes perc_genes total_expr
#>   1:   cell_1     5948      59.48  5580.4220
#>   2:   cell_2     2945      29.45  4512.3119
#>   3:   cell_3     3953      39.53  5005.0044
#>   4:   cell_4     1429      14.29  3252.5006
#>   5:   cell_5     3189      31.89  4632.9810
#>  ---                                        
#> 909: cell_909     1747      17.47  3563.6841
#> 910: cell_910     1292      12.92  3099.3098
#> 911: cell_911      300       3.00  1294.1751
#> 912: cell_912     1267      12.67  3054.7134
#> 913: cell_913      179       1.79   900.0755

# add metadata by column name
cortex_annotation = fread(system.file("extdata", "cortex_annotation.txt", package = "Giotto"))
VC_test = addCellMetadata(gobject = VC_test, new_metadata = cortex_annotation, by_column = T, column_cell_ID = 'uniq_ID')

# cell metadata after
cell_metadata = pDataDT(VC_test)
print(cell_metadata)
#>       cell_ID nr_genes perc_genes total_expr Field_of_View
#>   1:   cell_1     5948      59.48  5580.4220             0
#>   2:   cell_2     2945      29.45  4512.3119             0
#>   3:   cell_3     3953      39.53  5005.0044             0
#>   4:   cell_4     1429      14.29  3252.5006             0
#>   5:   cell_5     3189      31.89  4632.9810             0
#>  ---                                                      
#> 909: cell_909     1747      17.47  3563.6841             6
#> 910: cell_910     1292      12.92  3099.3098             6
#> 911: cell_911      300       3.00  1294.1751             6
#> 912: cell_912     1267      12.67  3054.7134             6
#> 913: cell_913      179       1.79   900.0755             6
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
spatPlot2D(gobject = VC_subset, point_size = 2, cell_color = 'cell_types', return_plot = F)
#> create 2D plot with ggplot
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="75%" style="display: block; margin: auto;" />

  - subset for SVZ cells (located in field 5)

<!-- end list -->

``` r
SVZ_cells = cell_metadata[Field_of_View == 5][['cell_ID']]
VC_SVZ = subsetGiotto(VC_test, cell_ids = SVZ_cells)
spatPlot2D(gobject = VC_SVZ, point_size = 4, cell_color = 'cell_types', return_plot = F)
#> create 2D plot with ggplot
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="75%" style="display: block; margin: auto;" />

  - subset for all interneurons

<!-- end list -->

``` r
interneuron_cells = cell_metadata[cell_types == 'Interneuron'][['cell_ID']]
VC_interneuron = subsetGiotto(VC_test, cell_ids = interneuron_cells)
spatPlot2D(gobject = VC_interneuron, point_size = 2, cell_color = 'cell_types', return_plot = F)
#> create 2D plot with ggplot
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="75%" style="display: block; margin: auto;" />
