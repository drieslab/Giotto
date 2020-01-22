
  - [Giotto global instructions](#giotto-global-instructions)
  - [Data input](#data-input)
  - [part 1: Create Giotto object & process
    data](#part-1-create-giotto-object-process-data)
  - [part 2: dimension reduction](#part-2-dimension-reduction)
  - [part 3: cluster](#part-3-cluster)
  - [part 4: co-visualize](#part-4-co-visualize)
  - [part 5: differential expression](#part-5-differential-expression)
  - [part 6: cell-type annotation](#part-6-cell-type-annotation)
  - [part 7: spatial grid](#part-7-spatial-grid)
  - [part 8: spatial network](#part-8-spatial-network)
  - [part 9: spatial genes](#part-9-spatial-genes)
  - [part 10: HMRF domains](#part-10-hmrf-domains)
  - [part 11: Cell-cell preferential
    proximity](#part-11-cell-cell-preferential-proximity)
  - [Export and create Giotto Viewer](#export-and-create-giotto-viewer)

<!-- mouse_cortex_1_simple.md is generated from mouse_cortex_1_simple.Rmd Please edit that file -->

### Giotto global instructions

``` r
# this example works with Giotto v.0.1.3
library(Giotto)

## create instructions
## instructions allow us to automatically save all plots into a chosen results folder
## Here we will not automatically save plots, for an example see the visium kidney dataset
## We will only set the python path which is needed for certain analyses
my_python_path = "/Users/rubendries/Bin/anaconda3/envs/py36/bin/pythonw"
results_folder = '/path/to/Visium/Brain/191204_results//'
instrs = createGiottoInstructions(python_path = my_python_path)
```

### Data input

[10X genomics](https://www.10xgenomics.com/spatial-transcriptomics/)
recently launched a new platform to obtain spatial expression data using
a Visium Spatial Gene Expression slide.

![](./visium_technology.png)

``` r
## expression and cell location
## expression data
data_path = '/path/to/Visium_data/Brain_data/raw_feature_bc_matrix//'
raw_matrix = get10Xmatrix(path_to_data = data_path)
library("biomaRt") # convert ensembl to gene names
raw_matrix = convertEnsemblToGeneSymbol(matrix = raw_matrix, species = 'mouse')

## spatial results
spatial_results = fread('/path/to/Visium_data/Brain_data/spatial/tissue_positions_list.csv')
spatial_results = spatial_results[match(colnames(raw_matrix), V1)]
colnames(spatial_results) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
```

-----

### part 1: Create Giotto object & process data

<details>

<summary>Expand</summary>  

``` r
## create
## we need to reverse the column pixel column (col_pxl) to get the same .jpg image as provided by 10X
visium_brain <- createGiottoObject(raw_exprs = raw_matrix,
                                    spatial_locs = spatial_results[,.(row_pxl,-col_pxl)],
                                    instructions = instrs,
                                    cell_metadata = spatial_results[,.(in_tissue, array_row, array_col)])

## check metadata
pDataDT(visium_brain)

## compare in tissue with provided jpg
spatPlot2D(gobject = visium_brain,  point_size = 2,
           cell_color = 'in_tissue', cell_color_code = c('0' = 'lightgrey', '1' = 'blue'))

## subset on spots that were covered by brain tissue
metadata = pDataDT(visium_brain)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_brain = subsetGiotto(visium_brain, cell_ids = in_tissue_barcodes)

## filter genes and cells
visium_brain <- filterGiotto(gobject = visium_brain,
                              expression_threshold = 1,
                              gene_det_in_min_cells = 50,
                              min_det_genes_per_cell = 1000,
                              expression_values = c('raw'),
                              verbose = T)

## normalize
visium_brain <- normalizeGiotto(gobject = visium_brain, scalefactor = 6000, verbose = T)

## add gene & cell statistics
visium_brain <- addStatistics(gobject = visium_brain)

## visualize
# location of spots
spatPlot2D(gobject = visium_brain,  point_size = 2)

# number of genes per spot
spatPlot2D(gobject = visium_brain, cell_color = 'nr_genes', color_as_factor = F,  point_size = 2)
```

High resolution png from original tissue.  
![](./mouse_brain_highres.png)

Spots labeled according to whether they were covered by tissue or not:  
![](./figures/1_in_tissue.png)

Spots after subsetting and filtering:  
![](./figures/1_spatial_locations.png)

Overlay with number of genes detected per spot:  
![](./figures/1_nr_genes.png)

</details>

### part 2: dimension reduction

<details>

<summary>Expand</summary>  

``` r
## highly variable genes (HVG)
visium_brain <- calculateHVG(gobject = visium_brain)

## select genes based on HVG and gene statistics, both found in gene metadata
gene_metadata = fDataDT(visium_brain)
featgenes = gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID

## run PCA on expression values (default)
visium_brain <- runPCA(gobject = visium_brain, genes_to_use = featgenes, scale_unit = F)
# significant PCs
signPCA(visium_brain, genes_to_use = featgenes, scale_unit = F)
# plot PCA
plotPCA(gobject = visium_brain)

## run UMAP 
visium_brain <- runUMAP(visium_brain, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_brain)

## run tSNE
visium_brain <- runtSNE(visium_brain, dimensions_to_use = 1:10)
plotTSNE(gobject = visium_brain)
```

highly variable genes:  
![](./figures/2_HVGplot.png)

screeplot to determine number of Principal Components to keep:  
![](./figures/2_screeplot.png)

PCA:  
![](./figures/2_PCA_reduction.png)

UMAP:  
![](./figures/2_UMAP_reduction.png)

tSNE:  
![](./figures/2_tSNE_reduction.png) \*\*\*

</details>

### part 3: cluster

<details>

<summary>Expand</summary>  

``` r
## sNN network (default)
visium_brain <- createNearestNetwork(gobject = visium_brain, dimensions_to_use = 1:10, k = 15)

## Leiden clustering
visium_brain <- doLeidenCluster(gobject = visium_brain, resolution = 0.4, n_iterations = 1000)

# default cluster result name from doLeidenCluster = 'leiden_clus'
plotUMAP(gobject = visium_brain, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2)
```

Leiden clustering:  
![](./figures/3_UMAP_leiden.png)

-----

</details>

### part 4: co-visualize

<details>

<summary>Expand</summary>  

``` r
# leiden clustering results
spatDimPlot(gobject = visium_brain, cell_color = 'leiden_clus',
            dim_point_size = 1.5, spat_point_size = 1.5)

# number of genes detected per spot
spatDimPlot(gobject = visium_brain, cell_color = 'nr_genes', color_as_factor = F,
            dim_point_size = 1.5, spat_point_size = 1.5)

# zoom-in on Dentate Gyrus by subsetting giotto object based on spatial coordinates/locations
DG_subset = subsetGiottoLocs(visium_brain, x_max = 6500, x_min = 3000, y_max = -2500, y_min = -5500, return_gobject = T)
spatDimPlot(gobject = DG_subset, cell_color = 'leiden_clus', point_size = 5)
```

Co-visualzation: ![](./figures/4_covis_leiden.png)

Co-visualzation overlaid with number of genes detected:  
![](./figures/4_nr_genes.png)

Zoom-in on Dentate Gyrus:  
![](./figures/4_zoom_dentate_gyrus.png)

-----

</details>

### part 5: differential expression

<details>

<summary>Expand</summary>  

``` r
## gini ##
## ---- ##
gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_brain,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_genes = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)

# violinplot
topgenes_gini = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes
violinPlot(visium_brain, genes = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right')

# cluster heatmap
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes
my_cluster_order = c(5, 13, 7, 2, 1, 10, 14, 6, 12, 9, 3, 4 , 8, 11, 15)
plotMetaDataHeatmap(visium_brain, selected_genes = topgenes_gini, custom_cluster_order = my_cluster_order,
                    metadata_cols = c('leiden_clus'), x_text_size = 10, y_text_size = 10)

# umap plots
dimGenePlot2D(visium_brain, expression_values = 'scaled',
              genes = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes,
              cow_n_col = 3, point_size = 1,
              genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0)



## scran ##
## ----- ##
scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_brain,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')

# violinplot
topgenes_scran = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes
violinPlot(visium_brain, genes = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right')

# cluster heatmap
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes
plotMetaDataHeatmap(visium_brain, selected_genes = topgenes_scran, custom_cluster_order = my_cluster_order,
                    metadata_cols = c('leiden_clus'))

# umap plots
dimGenePlot2D(visium_brain, expression_values = 'scaled',
              genes = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes,
              cow_n_col = 3, point_size = 1,
              genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0)
```

Gini: - violinplot: ![](./figures/5_violinplot_gini.png)

  - Heatmap clusters: ![](./figures/5_metaheatmap_gini.png)

  - UMAPs: ![](./figures/5_gini_umap.png)

Scran: - violinplot: ![](./figures/5_violinplot_scran.png)

  - Heatmap clusters: ![](./figures/5_metaheatmap_scran.png)

  - UMAPs: ![](./figures/5_scran_umap.png)

-----

</details>

### part 6: cell-type annotation

<details>

<summary>Expand</summary>  

Visium spatial transcriptomics does not provide single-cell resolution,
making cell type annotation a harder problem. Giotto provides 3 ways to
calculate enrichment of specific cell-type signature gene list:  
\- PAGE  
\- rank  
\- hypergeometric test

To generate the cell-type specific gene lists for the mouse brain data
we used cell-type specific gene sets as identified in [Zeisel, A. et
al. Molecular Architecture of the Mouse Nervous
System](https://www.mousebrain.org)

![](./clusters_Zeisel.png)

``` r

# known markers for different mouse brain cell types:
# Zeisel, A. et al. Molecular Architecture of the Mouse Nervous System. Cell 174, 999-1014.e22 (2018).

## cell type signatures ##
## combination of all marker genes identified in Zeisel et al
brain_sc_markers = fread('/path/to/Visium_data/Brain_data/sig_matrix.txt')
sig_matrix = as.matrix(brain_sc_markers[,-1]); rownames(sig_matrix) = brain_sc_markers$Event
  
## enrichment tests 
visium_brain = createSpatialEnrich(visium_brain, sign_matrix = sig_matrix, enrich_method = 'PAGE') #default = 'PAGE'

## heatmap of enrichment versus annotation (e.g. clustering result)
cell_types = colnames(sig_matrix)
plotMetaDataCellsHeatmap(gobject = visium_brain,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types,
                         spat_enr_names = 'PAGE',x_text_size = 8, y_text_size = 8,
                         save_param = c(save_folder = '7_annotation', save_name = 'heatmap_PAGE',
                                        base_width = 12, base_height = 8, units = 'cm'))

## spatial enrichment results for all cell type signatures ##
cell_types_subset = colnames(sig_matrix)[1:10]
spatCellPlot(gobject = visium_brain, spat_enr_names = 'PAGE',
             cell_annotation_values = cell_types_subset,
             cow_n_col = 4,coord_fix_ratio = NULL, point_size = 0.75,
             save_param = c(save_folder = '7_annotation', save_name = 'PAGE_spatplot_1_10',
                            base_width = 10, base_height = 6))

cell_types_subset = colnames(sig_matrix)[11:20]
spatCellPlot(gobject = visium_brain, spat_enr_names = 'PAGE',
             cell_annotation_values = cell_types_subset,
             cow_n_col = 4,coord_fix_ratio = NULL, point_size = 0.75,
             save_param = c(save_folder = '7_annotation', save_name = 'PAGE_spatplot_11_20',
                            base_width = 10, base_height = 6))


## spatial and dimension reduction visualization with
spatDimCellPlot(gobject = visium_brain, spat_enr_names = 'PAGE',
                cell_annotation_values = c('Cortex_hippocampus', 'Granule_neurons', 'di_mesencephalon_1', 'Oligo_dendrocyte', 'Vascular'),
                cow_n_col = 1, spat_point_size = 1, plot_alignment = 'horizontal')


## visualize individual spatial enrichments
spatDimPlot(gobject = visium_brain,
            spat_enr_names = 'PAGE',
            cell_color = 'Cortex_hippocampus', color_as_factor = F,
            spat_show_legend = T, dim_show_legend = T,
            gradient_midpoint = 0, 
            dim_point_size = 1.5, spat_point_size = 1.5)

spatDimPlot(gobject = visium_brain,
            spat_enr_names = 'PAGE',
            cell_color = 'Granule_neurons', color_as_factor = F,
            spat_show_legend = T, dim_show_legend = T,
            gradient_midpoint = 0, 
            dim_point_size = 1.5, spat_point_size = 1.5
```

Heatmap:

![](./figures/6_heatmap_PAGE.png)

Spatial enrichment plots for all cell types/clusters:

![](./figures/6_PAGE_spatplot_1_10.png)

![](./figures/6_PAGE_spatplot_11_20.png)

Co-visualization for selected subset:

![](./figures/6_PAGE_spatdimplot.png)

Cortex hippocampus specific plot:  
![](./figures/6_PAGE_spatdimplot_cortex_hc.png)

Dentate gyrus specific plot:  
![](./figures/6_PAGE_spatdimplot_gc.png)

-----

</details>

### part 7: spatial grid

<details>

<summary>Expand</summary>  

``` r
# create spatial grid
visium_brain <- createSpatialGrid(gobject = visium_brain,
                                   sdimx_stepsize = 400,
                                   sdimy_stepsize = 400,
                                   minimum_padding = 0)
spatPlot(visium_brain, cell_color = 'leiden_clus', show_grid = T,
         grid_color = 'red', spatial_grid_name = 'spatial_grid')


### spatial patterns ###
pattern_osm = detectSpatialPatterns(gobject = visium_brain, 
                                    spatial_grid_name = 'spatial_grid',
                                    min_cells_per_grid = 3, 
                                    scale_unit = T, 
                                    PC_zscore = 0.5, 
                                    show_plot = T)

# dimension 1
PC_dim = 1
showPattern2D(visium_brain, pattern_osm, dimension = PC_dim, point_size = 4)
showPatternGenes(visium_brain, pattern_osm, dimension = PC_dim)

# dimension 2
PC_dim = 2
showPattern2D(visium_brain, pattern_osm, dimension = PC_dim, point_size = 4)
showPatternGenes(visium_brain, pattern_osm, dimension = PC_dim)

# dimension 3
PC_dim = 3
showPattern2D(visium_brain, pattern_osm, dimension = PC_dim, point_size = 4)
showPatternGenes(visium_brain, pattern_osm, dimension = PC_dim)

view_pattern_genes = selectPatternGenes(pattern_osm, return_top_selection = TRUE)
```

![](./figures/7_grid.png)

Dimension 1: ![](./figures/7_pattern1_PCA.png)
![](./figures/7_pattern1_genes.png)

Dimension 2: ![](./figures/7_pattern2_PCA.png)

![](./figures/7_pattern2_genes.png)

Dimension 2: ![](./figures/7_pattern3_PCA.png)

![](./figures/7_pattern3_genes.png)

-----

</details>

### part 8: spatial network

<details>

<summary>Expand</summary>  

``` r
# create spatial network
visium_brain <- createSpatialNetwork(gobject = visium_brain, k = 5, maximum_distance = 400)
spatPlot(gobject = visium_brain, show_network = T, point_size = 1,
         network_color = 'blue', spatial_network_name = 'spatial_network'))
```

![](./figures/8_spatial_network_k5.png)

-----

</details>

### part 9: spatial genes

<details>

<summary>Expand</summary>  

``` r
## kmeans binarization
kmtest = binGetSpatialGenes(visium_brain, bin_method = 'kmeans',
                            do_fisher_test = T, community_expectation = 5,
                            spatial_network_name = 'spatial_network', verbose = T)
spatGenePlot(visium_brain, expression_values = 'scaled',
             genes = kmtest$genes[1:6], cow_n_col = 2, point_size = 1,
             genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0)

## rank binarization
ranktest = binGetSpatialGenes(visium_brain, bin_method = 'rank',
                              do_fisher_test = T, community_expectation = 5,
                              spatial_network_name = 'spatial_network', verbose = T)
spatGenePlot(visium_brain, expression_values = 'scaled',
             genes = ranktest$genes[1:6], cow_n_col = 2, point_size = 1,
             genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0)

## distance
spatial_genes = calculate_spatial_genes_python(gobject = visium_brain,
                                               expression_values = 'scaled',
                                               rbp_p=0.95, examine_top=0.3)
spatGenePlot(visium_brain, expression_values = 'scaled',
             genes = spatial_genes$genes[1:6], cow_n_col = 2, point_size = 1,
             genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0)
```

Spatial genes:  
\- kmeans ![](./figures/9_spatial_genes_km.png)

  - rank ![](./figures/9_spatial_genes_rank.png)

  - distance  
    ![](./figures/9_spatial_genes.png)

-----

</details>

### part 10: HMRF domains

<details>

<summary>Expand</summary>  

``` r
# spatial genes
my_spatial_genes <- spatial_genes[1:100]$genes

# do HMRF with different betas
hmrf_folder = paste0(results_folder,'/','11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

HMRF_spatial_genes = doHMRF(gobject = visium_brain, expression_values = 'scaled',
                            spatial_genes = my_spatial_genes,
                            k = 12,
                            betas = c(0, 0.5, 6), 
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k12_scaled'))

## view results of HMRF
for(i in seq(0, 3, by = 0.5)) {
  viewHMRFresults2D(gobject = visium_brain,
                    HMRFoutput = HMRF_spatial_genes,
                    k = 12, betas_to_view = i,
                    point_size = 2)
}


## alternative way to add all HMRF results 
#results = writeHMRFresults(gobject = ST_test,
#                           HMRFoutput = HMRF_spatial_genes,
#                           k = 12, betas_to_view = seq(0, 3, by = 0.5))
#ST_test = addCellMetadata(ST_test, new_metadata = results, by_column = T, column_cell_ID = 'cell_ID')


## add HMRF of interest to giotto object
visium_brain = addHMRF(gobject = visium_brain,
                        HMRFoutput = HMRF_spatial_genes,
                        k = 12, betas_to_add = c(0, 0.5),
                        hmrf_name = 'HMRF')

## visualize
# b = 0
spatPlot(gobject = visium_brain, cell_color = 'HMRF_k12_b.0', point_size = 2)

# b = 1
spatPlot(gobject = visium_brain, cell_color = 'HMRF_k12_b.0.5', point_size = 2)
```

HMRF:  
b = 0  
![](./figures/10_HMRF_k12_b.0.png)

b = 1 results in smoother spatial domains:  
![](./figures/10_HMRF_k12_b.0.5.png)

-----

</details>

### part 11: Cell-cell preferential proximity

<details>

<summary>Expand</summary>  

![cell-cell](./cell_cell_neighbors.png)

``` r
## calculate frequently seen proximities
cell_proximities = cellProximityEnrichment(gobject = visium_brain,
                                           cluster_column = 'leiden_clus',
                                           spatial_network_name = 'spatial_network',
                                           number_of_simulations = 1000)

## barplot
cellProximityBarplot(gobject = visium_brain, CPscore = cell_proximities, min_orig_ints = 5, min_sim_ints = 5)

## heatmap
cellProximityHeatmap(gobject = visium_brain, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'))
## network
cellProximityNetwork(gobject = visium_brain, CPscore = cell_proximities, remove_self_edges = T, only_show_enrichment_edges = F)

## visualization
spec_interaction = "5--13"
cellProximitySpatPlot2D(gobject = visium_brain,
                        interaction_name = spec_interaction, cell_color_code = c('5' = 'orange', '13' = 'blue'),
                        cluster_column = 'leiden_clus', show_network = T,
                        cell_color = 'leiden_clus', coord_fix_ratio = 0.5,
                        point_size_select = 1.5, point_size_other = 1)
```

barplot:  
![](./figures/11_barplot_cell_cell_enrichment.png)

heatmap:  
![](./figures/11_heatmap_cell_cell_enrichment.png)

network:  
![](./figures/11_network_cell_cell_enrichment.png)

selected enrichment:  
![](./figures/11_selected_enrichment.png)

-----

</details>

### Export and create Giotto Viewer

<details>

<summary>Expand</summary>  

##### 1\. Export Giotto results to a specificied directory

  - export spot/cell annotations  
  - export dimension reduction coordinates (umap, tsne, …)  
  - export expression data

This function will create a directory that, together with the 10X
provided .tiff file, can be used to create an interactive Giotto Viewer

``` r
# select annotations, reductions and expression values to view in Giotto Viewer
viewer_folder = paste0(results_folder, '/', 'mouse_visium_brain_viewer')

exportGiottoViewer(gobject = visium_brain,
                   output_directory = viewer_folder,
                   spat_enr_names = 'PAGE', 
                   factor_annotations = c('in_tissue',
                                          'leiden_clus',
                                          'HMRF_k12_b.1'),
                   numeric_annotations = c('nr_genes',
                                           'Granule_neurons'),
                   dim_reductions = c('tsne', 'umap'),
                   dim_reduction_names = c('tsne', 'umap'),
                   expression_values = 'scaled',
                   expression_rounding = 2,
                   overwrite_dir = T)
```

##### 2\. Create a Giotto Viewer within the specified directory

After installing the Giotto Viewer requirements (see README), you can
create a viewer directly from the **terminal** using the following
commands:

``` bash
#===========================================
# 10X Genomics Visium specific instructions
#===========================================

# first go the working directory where you exported the Giotto results
cd viewer_folder #should be the directory where files are exported to by giotto
cp ~/Downloads/V1_Adult_Mouse_Brain_image.tif . # copy the .tiff image provided by 10X Visium here to this directory

## STEP 1: ##
# create step1 json file
~/.local/bin/giotto_setup_image --require-stitch=n --image=y --image-multi-channel=n --segmentation=n --multi-fov=n --output-json=step1.json

# automatically fill in image dimension in the step1 json file
~/.local/bin/giotto_step1_modify_json --add-image V1_Adult_Mouse_Brain_image.tif --input step1.json --output step1.json

# do the step1 actions
~/.local/bin/smfish_step1_setup -c step1.json


## STEP 2: ##
# create step2 json file
~/.local/bin/giotto_setup_viewer --num-panel=2 --input-preprocess-json=step1.json --panel-1=PanelPhysical10X --panel-2=PanelTsne --output-json=step2.json --input-annotation-list=annotation_list.txt

# do the step2 actions
~/.local/bin/smfish_read_config -c step2.json -o test.dec6.js -p test.dec6.html -q test.dec6.css

# copy extra js and css folders
~/.local/bin/giotto_copy_js_css --output .

# LAUNCH viewer
python3 -m http.server
# the viewer is now launched and can be seen in a browser (e.g. Chrome) using the URL localhost:8000/test.dec6.html
```

##### 3\. Update the Giotto Viewer

If you want to update the Viewer with new results (e.g. new
annotations), you first need to call the **exportGiottoViewer** function
again with the updated results.

Then follow these instructions:

``` bash
#=================================
# Updating annotations
#=================================

cd viewer_folder # should be the directory where files are exported to by giotto

# no need to re-do step 1

# re-create step2 json file
~/.local/bin/giotto_setup_viewer --num-panel=2 --input-preprocess-json=step1.json --panel-1=PanelPhysical10X --panel-2=PanelTsne --output-json=step2.json --input-annotation-list=annotation_list.txt

# re-do the step2 actions
~/.local/bin/smfish_read_config -c step2.json -o test.dec6.js -p test.dec6.html -q test.dec6.css

# no need to restart http.server
# go back to localhost:8000/test.dec6.html (refresh the page) (recommend incognito mode in Chrome)
```

##### 4\. Customize the viewer

If you want to customize the viewer output results, then modify the
step2.json after the create step. Then re-do step2 actions. Things that
can be modified at this stage:  
\- map height  
\- default annotation of each panel  
\- umap or tsne in the PanelTsne

To customize number of panels, what goes in each panel (PanelTsne,
PanelPhysical, PanelPhysical10X, or PanelPhysicalSimple), and add
annotations, please re-run the step2.json creation step because these
changes are major.

-----

</details>
