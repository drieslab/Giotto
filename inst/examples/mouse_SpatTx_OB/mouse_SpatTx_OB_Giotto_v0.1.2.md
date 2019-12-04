
<!-- mouse_cortex_1_simple.md is generated from mouse_cortex_1_simple.Rmd Please edit that file -->

### Giotto global instructions

``` r
# this example works with Giotto v.0.1.2
library(Giotto)

## create instructions
my_python_path = "/Users/rubendries/Bin/anaconda3/envs/py36/bin/pythonw"
results_folder = '/Path/to/Results/SpatTx_OB_results/'
instrs = createGiottoInstructions(python_path = my_python_path,
                                  show_plot = F, return_plot = T, save_plot = T,
                                  save_dir = results_folder,
                                  plot_format = 'png',
                                  dpi = 300, height = 9, width = 9)
```

### Data input

[Stahl et al.](https://science.sciencemag.org/content/353/6294/78) used
immobilized reverse-transcription oligo(dT) primers on glass glides to
profile the spatial expression of the mouse olfactory bulb.

![](./spatial_transcriptomics_summary.png) .

``` r
## select the directory where you have saved the Spatial Transcriptomics data
data_folder = '/path/to/ST_data/'
expr = read.table(paste0(data_folder,'/','Rep11_MOB_0_expr.txt'))
cell_locs = fread(paste0(data_folder,'/','Rep11_MOB_0_location.txt'))
```

-----

### 1\. Create Giotto object & process data

<details>

<summary>Expand</summary>  

``` r
## create
ST_test <- createGiottoObject(raw_exprs = expr, spatial_locs = cell_locs[,-1], instructions = instrs)
showGiottoInstructions(ST_test)

## filter
ST_test <- filterGiotto(gobject = ST_test,
                         expression_threshold = 1,
                         gene_det_in_min_cells = 1,
                         min_det_genes_per_cell = 100,
                         expression_values = c('raw'),
                         verbose = T)

## normalize
ST_test <- normalizeGiotto(gobject = ST_test)
ST_test <- addStatistics(gobject = ST_test)

## visualize
spatPlot(gobject = ST_test,
           save_param = list(save_folder = '2_Gobject', save_name = 'spatial_locations'))
```

![](./figures/1_spatial_locations.png)

</details>

### 2\. dimension reduction

<details>

<summary>Expand</summary>  

``` r
## highly variable genes (HVG)
ST_test <- calculateHVG(gobject = ST_test)
gene_metadata = fDataDT(ST_test)
featgenes = gene_metadata[hvg == 'yes' & perc_cells > 4 & mean_expr_det > 0.25]$gene_ID

## run PCA on expression values (default)
ST_test <- runPCA(gobject = ST_test, genes_to_use = featgenes, expression_values = 'scaled', scale_unit = F)
signPCA(gobject = ST_test, expression_values = 'scaled', scale_unit = F, scree_ylim = c(0,1),
        save_param = list(save_folder = '3_DimRed', save_name = 'screeplot'))
plotPCA(ST_test,
           save_param = list(save_folder = '3_DimRed', save_name = 'PCA_reduction'))

## run UMAP and tSNE on PCA space (default)
ST_test <- runUMAP(ST_test, dimensions_to_use = 1:10, expression_values = 'scaled', n_threads = 2)
plotUMAP(gobject = ST_test,
            save_param = list(save_folder = '3_DimRed', save_name = 'UMAP_reduction'))
```

![](./figures/2_screeplot.png)

![](./figures/2_PCA_reduction.png)

![](./figures/2_UMAP_reduction.png)

-----

</details>

### 3\. cluster

<details>

<summary>Expand</summary>  

``` r

## Leiden clustering
ST_test <- createNearestNetwork(gobject = ST_test, dimensions_to_use = 1:10, k = 10)
ST_test <- doLeidenCluster(gobject = ST_test, resolution = 0.2, n_iterations = 200)
plotUMAP(gobject = ST_test, cell_color = 'leiden_clus', point_size = 3,
            save_param = list(save_folder = '4_Cluster', save_name = 'UMAP_leiden'))
```

![](./figures/3_UMAP_leiden.png)

-----

</details>

### 4\. co-visualize

<details>

<summary>Expand</summary>  

``` r
spatDimPlot(gobject = ST_test, cell_color = 'leiden_clus', 
            dim_point_size = 2, spat_point_size = 6,
            save_param = list(save_folder = '5_Covisuals', save_name = 'covisual_leiden'))
```

Co-visualzation: ![](./figures/4_covisual_leiden.png)

-----

</details>

### 5\. differential expression

<details>

<summary>Expand</summary>  

``` r
## gini ##
## very specific to a cluster, but not necessarily expressed in all cells of that cluster
gini_markers_subclusters = findMarkers_one_vs_all(gobject = ST_test,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_genes = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes

# violinplot
violinPlot(ST_test, genes = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'top',
           save_param = c(save_folder = '6_DEG', save_name = 'violinplot_gini',  base_width = 5, base_height = 10))

# cluster heatmap
plotMetaDataHeatmap(ST_test, selected_genes = topgenes_gini,metadata_cols = c('leiden_clus'),
                    save_param = c(save_folder = '6_DEG', save_name = 'metaheatmap_gini'))


# visualize genes
topgenes_gini = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes
spatDimGenePlot(ST_test, expression_values = 'scaled',
                  genes = topgenes_gini,
                  plot_alignment = 'horizontal', cow_n_col = 1, point_size = 2,
                  genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0,
                  save_param = c(save_folder = '6_DEG', save_name = 'genes_gini', base_width = 6, base_height = 14))



## scran ##
scran_markers_subclusters = findMarkers_one_vs_all(gobject = ST_test,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster_ID']$gene_ID

# violinplot
violinPlot(ST_test, genes = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'top',
           save_param = c(save_folder = '6_DEG', save_name = 'violinplot_scran',  base_width = 5, base_height = 10))

# cluster heatmap
plotMetaDataHeatmap(ST_test, selected_genes = topgenes_scran, metadata_cols = c('leiden_clus'),
                    save_param = c(save_folder = '6_DEG', save_name = 'metaheatmap_scran'))

# visualize genes
topgenes_scran = scran_markers_subclusters[, head(.SD, 1), by = 'cluster_ID']$gene_ID
spatDimGenePlot(ST_test, expression_values = 'scaled',
                  genes = topgenes_scran,
                  plot_alignment = 'horizontal', cow_n_col = 1, point_size = 2,
                  genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0,
                  save_param = c(save_folder = '6_DEG', save_name = 'genes_scran', base_width = 6, base_height = 14))
```

Gini: - violinplot: ![](./figures/5_violinplot_gini.png)

  - Heatmap clusters: ![](./figures/5_metaheatmap_gini.png)

  - Genes: ![](./figures/5_genes_gini.png)

Scran: - violinplot: ![](./figures/5_violinplot_scran.png)

  - Heatmap clusters: ![](./figures/5_metaheatmap_scran.png)

  - Genes: ![](./figures/5_genes_scran.png)

-----

</details>

### 6\. cell-type annotation

<details>

<summary>Expand</summary>  

Spatial transcriptomics does not provide single-cell resolution, making
cell type annotation a harder problem. Giotto provides 3 ways to
calculate enrichment of specific cell-type signature gene list:  
\- PAGE  
\- rank  
\- hypergeometric test

To generate the cell-type specific gene lists for the olfactory bulb
(OB) data we reanalyzed the paper from Zeisel et al, associated with the
[mouse brain atlas](http://mousebrain.org/), and identified cell-type
specific genes for cells from the OB.

![paper](./Zeisel_paper.png)

``` r

## cell type identification based on individual marker genes is hard 
# known markers for different interneuron subtypes
spatDimGenePlot(ST_test, expression_values = 'scaled',
                genes = c('Vip', 'Camk4', 'Th', 'Igfbpl1'),
                plot_alignment = 'vertical', cow_n_col = 4, point_size = 3,
                genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0,
                save_param = c(save_folder = '7_annotation', save_name = 'interneuron_genes', base_width = 13, base_height = 5))




## cell type identification based on signatures from single-cell RNAseq
## for PAGE ##
sig_matrix = fread('/path/to/sig_matrix_PAGE.txt')
sig_matrix = Giotto:::dt_to_matrix(sig_matrix)

## example to make PAGE signature matrix from list of 2 signature genesets
OBDOP1_sig = sig_matrix[,'OBDOP1']; OBDOP1_sig_genes = names(OBDOP1_sig[OBDOP1_sig == 1])
ACOB_sig = sig_matrix[,'ACOB'];ACOB_sig_genes = names(ACOB_sig[ACOB_sig == 1])
small_sign_matrix = convertSignListToMatrix(sign_names = c('OBDOP1', 'ACOB'),
                               sign_list = list(OBDOP1_sig_genes, ACOB_sig_genes))


## for rank ##
sig_rank_matrix = fread('/path/to/sig_matrix_rank.txt')
sig_rank_matrix = Giotto:::dt_to_matrix(sig_rank_matrix)

## enrichment tests 
ST_test = createSpatialEnrich(ST_test, sign_matrix = sig_matrix) #default = 'PAGE' for method and name
ST_test = createSpatialEnrich(ST_test, sign_matrix = sig_matrix, output_enrichment = 'zscore', name = 'PAGEz') 
ST_test = createSpatialEnrich(ST_test, sign_matrix = sig_rank_matrix, enrich_method = 'rank', name = 'rank')


## heatmap
value_columns = c('ACOB', 'OBDOP1', 'OBDOP2-OBINH123', 'OBINH4', 'OBINH5', 'OBNBL12', 'OBNBL3', 'OBNBL45', 'OEC')
meta_columns = c('leiden_clus')

plotMetaDataCellsHeatmap(gobject = ST_test,
                         metadata_cols = 'leiden_clus',
                         value_cols = value_columns,
                         spat_enr_names = 'PAGE',
                         save_param = c(save_folder = '7_annotation', save_name = 'heatmap_PAGE',
                                        base_width = 4, base_height = 4))

plotMetaDataCellsHeatmap(gobject = ST_test,
                         metadata_cols = 'leiden_clus',
                         value_cols = value_columns,
                         spat_enr_names = 'rank',
                         save_param = c(save_folder = '7_annotation', save_name = 'heatmap_rank',
                                        base_width = 4, base_height = 4))

## visualize individual enrichments
spatDimPlot(gobject = ST_test,
            spat_enr_names = 'PAGE',
            cell_color = 'OEC', color_as_factor = F,
            spat_show_legend = T, dim_show_legend = T,
            gradient_midpoint = 3, 
            dim_point_size = 2, spat_point_size = 4, save_plot = F)

spatDimPlot(gobject = ST_test,
            spat_enr_names = 'PAGE',
            cell_color = 'OBINH4', color_as_factor = F,
            spat_show_legend = T, dim_show_legend = T,
            gradient_midpoint = 1, 
            dim_point_size = 2, spat_point_size = 4, save_plot = F)

spatDimPlot(gobject = ST_test,
            spat_enr_names = 'PAGE',
            cell_color = 'OBNBL3', color_as_factor = F,
            spat_show_legend = T, dim_show_legend = T,
            gradient_midpoint = 1, gradient_limits = c(0,4), 
            dim_point_size = 2, spat_point_size = 4, save_plot = F)



## multiple value columns with spatPlot ##
value_columns = c('ACOB', 'OBDOP1', 'OBDOP2-OBINH123', 'OBINH4', 'OBINH5', 'OBNBL12', 'OBNBL3', 'OBNBL45', 'OEC')

spatCellPlot(gobject = ST_test, spat_enr_names = 'PAGE',
             cell_annotation_values = value_columns,
             cow_n_col = 3,coord_fix_ratio = NULL,
             save_param = c(save_folder = '7_annotation', save_name = 'PAGE_spatplot',
                            base_width = 10, base_height = 6))

spatCellPlot(gobject = ST_test, spat_enr_names = 'PAGEz',
             cell_annotation_values = value_columns,
             cow_n_col = 3, coord_fix_ratio = NULL,
             gradient_limits = c(-2,2),
             save_param = c(save_folder = '7_annotation', save_name = 'PAGEz_spatplot',
                            base_width = 10, base_height = 6))

spatCellPlot(gobject = ST_test, spat_enr_names = 'rank',
             cell_annotation_values = value_columns,
             cow_n_col = 3, coord_fix_ratio = NULL,
             save_param = c(save_folder = '7_annotation', save_name = 'rank_spatplot',
                            base_width = 10, base_height = 6))


## multiple value columns with spatDimPlot ##
spatDimCellPlot(gobject = ST_test, spat_enr_names = 'PAGE',
                cell_annotation_values = value_columns[1:4],
                cow_n_col = 1, spat_point_size = 3, plot_alignment = 'horizontal',
                save_param = c(save_folder = '7_annotation', save_name = 'PAGE_spatdimplot',
                               base_width = 6, base_height = 10))
```

Markers for interneuron genes: ![](./figures/6_interneuron_genes.png)

Heatmap:

PAGE

![](./figures/6_heatmap_PAGE.png)

rank

![](./figures/6_heatmap_rank.png)

Spatial enrichment plots for all cell types:

PAGE enrichment:

![](./figures/6_PAGE_spatplot.png)

PAGE enrichment z-scores:

![](./figures/6_PAGEz_spatplot.png)

rank enrichment:

![](./figures/6_rank_spatplot.png)

Spatial and dimension reduction PAGE enrichment plots for first 4 cell
types:

![](./figures/6_PAGE_spatdimplot.png)

-----

</details>

### 7\. spatial grid

<details>

<summary>Expand</summary>  

``` r
## create spatial grid
ST_test <- createSpatialGrid(gobject = ST_test,
                              sdimx_stepsize = 2,
                              sdimy_stepsize = 2,
                              minimum_padding = 0)
spatPlot(ST_test, cell_color = 'leiden_clus', show_grid = T,
           grid_color = 'lightblue', spatial_grid_name = 'spatial_grid', 
           save_param = c(save_folder = '8_grid', save_name = 'grid'))

## identify spatial patterns with spatial grid #
pattern_osm = detectSpatialPatterns(gobject = ST_test, 
                                    spatial_grid_name = 'spatial_grid',
                                    min_cells_per_grid = 2, 
                                    scale_unit = T, 
                                    PC_zscore = 1, 
                                    show_plot = T)

# dimension 1
showPattern2D(ST_test, pattern_osm, dimension = 1, point_size = 4,
              save_param = c(save_folder = '8_grid', save_name = 'pattern1_PCA'))
showPatternGenes(ST_test, pattern_osm, dimension = 1,
                 save_param = c(save_folder = '8_grid', save_name = 'pattern1_genes'))

# dimension 2
showPattern2D(ST_test, pattern_osm, dimension = 2, point_size = 4,
              save_param = c(save_folder = '8_grid', save_name = 'pattern2_PCA'))
showPatternGenes(ST_test, pattern_osm, dimension = 2,
                 save_param = c(save_folder = '8_grid', save_name = 'pattern2_genes'))

view_pattern_genes = selectPatternGenes(pattern_osm, return_top_selection = TRUE)
```

![](./figures/7_grid.png)

Dimension 1: ![](./figures/7_pattern1_PCA.png)
![](./figures/7_pattern1_genes.png)

Dimension 2: ![](./figures/7_pattern2_PCA.png)

![](./figures/7_pattern2_genes.png)

-----

</details>

### 8\. spatial network

<details>

<summary>Expand</summary>  

``` r
ST_test <- createSpatialNetwork(gobject = ST_test, k = 5)
spatPlot(gobject = ST_test, show_network = T,
           network_color = 'blue', spatial_network_name = 'spatial_network',
           save_param = c(save_name = 'spatial_network_k5', save_folder = '9_spatial_network'))
```

![](./figures/8_spatial_network_k5.png)

-----

</details>

### 9\. spatial genes

<details>

<summary>Expand</summary>  

``` r
## kmeans binarization
kmtest = binGetSpatialGenes(ST_test, bin_method = 'kmeans',
                            do_fisher_test = T, community_expectation = 5,
                            spatial_network_name = 'spatial_network', verbose = T)

spatGenePlot(ST_test, expression_values = 'scaled',
               genes = kmtest$genes[1:6], cow_n_col = 2, point_size = 3,
               genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0,
               save_param = c(save_name = 'spatial_genes_km', save_folder = '10_spatial_genes'))

## rank binarization
ranktest = binGetSpatialGenes(ST_test, bin_method = 'rank',
                              do_fisher_test = T, community_expectation = 5,
                              spatial_network_name = 'spatial_network', verbose = T)

spatGenePlot(ST_test, expression_values = 'scaled',
               genes = ranktest$genes[1:6], cow_n_col = 2, point_size = 3,
               genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue', midpoint = 0,
               save_param = c(save_name = 'spatial_genes_rank', save_folder = '10_spatial_genes'))
```

Spatial genes: - kmeans ![](./figures/9_spatial_genes_km.png)

  - rank ![](./figures/9_spatial_genes_rank.png)

-----

</details>

### 10\. HMRF domains

<details>

<summary>Expand</summary>  

not available for this specific dataset

-----

</details>

### 11\. Cell-cell preferential proximity

<details>

<summary>Expand</summary>  

![cell-cell](./cell_cell_neighbors.png)

``` r
## calculate frequently seen proximities
cell_proximities = cellProximityEnrichment(gobject = ST_test,
                                           cluster_column = 'leiden_clus',
                                           spatial_network_name = 'spatial_network',
                                           number_of_simulations = 1000)

## barplot
cellProximityBarplot(gobject = ST_test, CPscore = cell_proximities, min_orig_ints = 5, min_sim_ints = 5, 
                     save_param = c(save_name = 'barplot_cell_cell_enrichment', save_folder = '12_cell_proxim'))
## heatmap
cellProximityHeatmap(gobject = ST_test, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),
                     save_param = c(save_name = 'heatmap_cell_cell_enrichment', save_folder = '12_cell_proxim', unit = 'in'))
## network
cellProximityNetwork(gobject = ST_test, CPscore = cell_proximities, remove_self_edges = T, only_show_enrichment_edges = F,
                     save_param = c(save_name = 'network_cell_cell_enrichment', save_folder = '12_cell_proxim'))

## visualization
spec_interaction = "2--6"
cellProximitySpatPlot2D(gobject = ST_test,
                        interaction_name = spec_interaction,
                        cluster_column = 'leiden_clus', show_network = T,
                        cell_color = 'leiden_clus', coord_fix_ratio = 0.5,
                        point_size_select = 4, point_size_other = 2,
                        save_param = c(save_name = 'selected_enrichment', save_folder = '12_cell_proxim'))
```

barplot:  
![](./figures/10_barplot_cell_cell_enrichment.png)

heatmap:  
![](./figures/10_heatmap_cell_cell_enrichment.png)

network:  
![](./figures/10_network_cell_cell_enrichment.png)

selected enrichment:  
![](./figures/10_selected_enrichment.png)

-----

</details>
