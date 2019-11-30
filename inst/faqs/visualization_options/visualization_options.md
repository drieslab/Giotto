
## Different ways to visualize your spatial data?

Giotto provides ways to visualize your spatial data:

##### 1\. Functions specific to a task or analysis:

  - **plotHeatmap**
  - **plotMetadataHeatmap**
  - **violinPlot**
  - check out the different dataset examples to see the wide range of
    visualization options

##### 2\. Functions that use the spatial coordinates of each cell or spot, they usually start with **spat**XXX

  - **spatPlot**, to plot cell annotation information (e.g. cluster or
    cell type)
  - **spatGenePlot**, to overlay (multiple) gene expression values
  - **spatCellPlot**, to overlay (multiple) numerical cell values
    (e.g. spatial enrichment)

##### 3\. Functions that use the dimension reduction coordinates of each cell or spot, they usually start with **dim**XXX

  - **dimPlot**, to plot cell annotation information (e.g. cluster or
    cell type)
  - **dimGenePlot**, to overlay (multiple) gene expression values
  - **dimCellPlot**, to overlay (multiple) numerical cell values
    (e.g. spatial enrichment)  
  - **plotPCA**, **plotUMAP** and **plotTSNE** are shortcuts for
    **dimPlot**

##### 4\. Functions for co-visualization that combine **spat** (2) and **dim** (3)

  - **spatDimPlot**, to plot cell annotation information (e.g. cluster
    or cell type)
  - **spatDimGenePlot**, to overlay (multiple) gene expression values
  - **spatDimCellPlot**, to overlay (multiple) numerical cell values
    (e.g. spatial enrichment)

##### 5\. Both in 2D and 3D

Most functions both have a 2D and 3D version, like **spatDimPlot2D** and
**spatDimPlot3D**. In those cases the **spatDimPlot2D** is the same as
**spatDimPlot**. So only in case you want to plot your spatial or
dimension data in 3D, you need to specifically say so.

##### 6\. ways to save plot

Hypothetical example showing a number of options

``` r
# 1. standard R way
pl = spatPlot(mygobject, cell_color = 'cell_types')
pdf(file = 'path/to/save/to/plot.pdf')
print(pl)
dev.off()

# 2. indicate to save plot, this will save the plot according to Giotto instructions file
# If the instruction file is not provided in the beginning, it uses the defaults (e.g. working directory)
spatPlot(mygobject, cell_color = 'cell_types', save_plot = TRUE)

# 3. indicate to save plot and specifiy specific saving parameters by providing a list to save_param
# they will overrule the giotto instructions
spatPlot(mygobject, cell_color = 'cell_types', save_plot = TRUE,
         save_param = list(save_folder = 'my_subfolder', save_name = 'my_name', save_format = 'png', units = 'in'))

# 4. don't save or return plot, but just view plot
# defaults are: save_plot = F, return_plot = T and show_plot = T
# this can be changed in the instructions file or at each specific plotting function
spatPlot(mygobject, cell_color = 'cell_types', save_plot = F, return_plot = F show_plot = T)
```
