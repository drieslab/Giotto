
<!-- mouse_cortex_example.md is generated from mouse_cortex_example.Rmd Please edit that file -->

Start from Giotto object from the Mouse Visual Cortex and SVZ:

``` r
spatPlot(gobject = VC_test, point_size = 2, return_plot = F)
```

<img src="../inst/faqs/metadata_and_subset/README-unnamed-chunk-5-1.png" />

### 1\. add annotation

``` r
# cell metadata before 
print(pDataDT(VC_test))

# add metadata by column name
cortex_annotation = fread(system.file("extdata", "cortex_annotation.txt", package = "Giotto"))
VC_test = addCellMetadata(gobject = VC_test, new_metadata = cortex_annotation, by_column = T, column_cell_ID = 'uniq_ID')

# cell metadata after
cell_metadata = pDataDT(VC_test)
print(cell_metadata)
```

### 2\. subset Giotto

  - subset first 500 cells

<!-- end list -->

``` r
first_500_cells = cell_metadata[1:500][['cell_ID']]
VC_subset = subsetGiotto(VC_test, cell_ids = first_500_cells)
spatPlot2D(gobject = VC_subset, point_size = 2, cell_color = 'cell_types', return_plot = F)
```

<center>

![](../inst/faqs/metadata_and_subset/README-unnamed-chunk-8-1.png)

</center>

  - subset for SVZ cells (located in field 5)

<!-- end list -->

``` r
SVZ_cells = cell_metadata[Field_of_View == 5][['cell_ID']]
VC_SVZ = subsetGiotto(VC_test, cell_ids = SVZ_cells)
spatPlot2D(gobject = VC_SVZ, point_size = 4, cell_color = 'cell_types', return_plot = F)
```

  - subset for all interneurons

<!-- end list -->

``` r
interneuron_cells = cell_metadata[cell_types == 'Interneuron'][['cell_ID']]
VC_interneuron = subsetGiotto(VC_test, cell_ids = interneuron_cells)
spatPlot2D(gobject = VC_interneuron, point_size = 2, cell_color = 'cell_types', return_plot = F)
```
