
# How to visualize and save plots in Giotto?

Each Giotto function that creates a plot as 3 important parameters:  
\- **show\_plot**: to print the plot to the console, default is TRUE  
\- **return\_plot**: to return the plot as an object, default is TRUE  
\- **save\_plot**: to automatically save the plot, default is FALSE
(only for Unix based systems\!)

These parameters are automatically provided by **createGiottoObject**,
but can also be explicitely provided using **createGiottoInstructions**
or a named list, which can be given to **createGiottoObject**.

Besides those 3 parameters **createGiottoInstructions** also allows to
provide other general Giotto parameters, such as your python path and
other information for automatically saving a plot, like size, plotting
format, etc.

In total there are 4 functions to work with setting instructions:  
\- **createGiottoInstructions**: creates instructions that can be
provided to createGiottoObject  
\- **showGiottoInstructions**: to view the instructions of a Giotto
object  
\- **changeGiottoInstructions**: to replace 1 or more of the instruction
parameters (e.g. plotting format)  
\- **replaceGiottoInstructions**: to replace all instructions with new
instructions (e.g after subsetting)  
In addition, for each plot the parameters can be *manually overwritten*
within the plotting function itself.

Here we use the osmFISH as an example:

### 1\. create Giotto instructions and Giotto object

  - we will override automatically showing plots by setting show\_plot =
    FALSE
  - we keep return\_plot = TRUE, to store plots and modify them if
    necessary  
  - we set save\_plot = TRUE to automatically saving each plot

We also specify some other plotting parameters for the automatic saving
functionality, such as the save\_dir whose default is the current
working directory.

``` r
library(Giotto)

# create instructions for your python path, how to view your plots and
# parameters to save your plot if wanted
my_python_path = "/your/python/path/"
results_folder = '/path/to/your/results/'
instrs = createGiottoInstructions(python_path = my_python_path,
                                  show_plot = FALSE,  
                                  return_plot = TRUE,
                                  save_plot = TRUE,
                                  save_dir = results_folder,
                                  plot_format = 'png',
                                  dpi = 200,
                                  height = 9,
                                  width = 9)
```

### 2\. create Giotto object

  - provide an expression matrix  
  - provide the cell locations
  - provide the previously generated instructions file (named list)

<!-- end list -->

``` r
## PREPARE osmFISH DATA ####
data_dir = '/directory/for/osmFISH_data/'
osm_exprs = read.table(file = paste0(data_dir,'/','osmFISH_prep_expression.txt'))
osm_locs = read.table(file = paste0(data_dir,'/','osmFISH_prep_cell_coordinates.txt'))
osm_locs = osm_locs[rownames(osm_locs) %in% colnames(osm_exprs),]

## CREATE GIOTTO OBJECT with expression data, location data and instructions
osm_test <- createGiottoObject(raw_exprs = osm_exprs,
                               spatial_locs = osm_locs,
                               instructions = instrs)

## add provided annotation and filter + normalize
metadata = fread(file = paste0(data_dir,'/','osmFISH_prep_cell_metadata.txt'))
osm_test = addCellMetadata(osm_test, new_metadata = metadata,
                           by_column = T, column_cell_ID = 'CellID')
osm_test <- filterGiotto(gobject = osm_test)
osm_test <- normalizeGiotto(gobject = osm_test)
```

### 2\. work with Giotto instructions

``` r
# show the provided Giotto instructions
showGiottoInstructions(osm_test)
# change a previously set parameter, e.g. change dpi = 200 to dpi = 300
osm_test = changeGiottoInstructions(osm_test, param = 'dpi', new_value = 300)
```

### 3\. Different ways of saving a plot

Here we will show a couple of ways to save plots.

  - default way according to Giotto instructions, if save\_plot = TRUE
    (1)
  - changing default parameters by by providing a named list to
    save\_param parameter (2 & 3)
  - block automatic saving, but modify the created ggplot object and
    save that one (4)

<!-- end list -->

``` r
# 1. default instructions from Giotto object
spatPlot2D(gobject = osm_test, cell_color = 'ClusterName')

# 2. overwrite save_name instruction by providing a named list for save_param
# save_param takes all parameters of all_plots_save_function
spatPlot2D(gobject = osm_test, cell_color = 'ClusterName', save_param = list(save_name = 'myplot'))

# 3. overwrite save_name instruction and add specific subfolder
# save_folder creates a specific subfolder in the provided directory
spatPlot2D(gobject = osm_test, cell_color = 'ClusterName', point_size = 1.5,
           save_param = list(save_folder = '2_Gobject', save_name = 'original_clusters', units = 'in', base_height = 6, base_width = 6))

# 4. do not save, but return as object, modify and save
# for example to create a black background
mypl = spatPlot2D(gobject = osm_test, cell_color = 'ClusterName', save_plot = F)
mypl = mypl + theme(panel.background = element_rect(fill ='black'),
                    panel.grid = element_blank())
mypl = mypl + guides(fill = guide_legend(override.aes = list(size=5)))
mypl

gobject_folder = paste0(results_folder,'/','2_Gobject/')
if(!file.exists(gobject_folder)) dir.create(gobject_folder, recursive = T)

cowplot::save_plot(plot = mypl,
                   filename = 'original_clusters_black.png', path = gobject_folder,
                   device = png(),
                   dpi = 300, base_height = 10, base_width = 10)
```

![](./figures/original_clusters_black.png)
