

library(Giotto)

# createGiottoInstructions(python_path = '/your/path')

## Giotto 0.3.5 ##
## mini-test Visium Brain Giotto 0.3.5 ##

## !! change this directory path !!:
#temp_dir = '/path/to/your/temporary/directory/results'
temp_dir = getwd()

## 1. giotto object ####
expr_path = system.file("extdata", "starmap_expr.txt", package = 'Giotto')
loc_path = system.file("extdata", "starmap_cell_loc.txt", package = 'Giotto')

# default
star_small <- createGiottoObject(raw_exprs = expr_path, spatial_locs = loc_path)

## 2. processing steps ####
filterDistributions(star_small, detection = 'genes')
filterDistributions(star_small, detection = 'cells')
filterCombinations(star_small,
                   expression_thresholds = c(1),
                   gene_det_in_min_cells = c(50, 100, 200),
                   min_det_genes_per_cell = c(20, 28, 28))

star_small <- filterGiotto(gobject = star_small,
                           expression_threshold = 1,
                           gene_det_in_min_cells = 50,
                           min_det_genes_per_cell = 20,
                           expression_values = c('raw'),
                           verbose = T)
star_small <- normalizeGiotto(gobject = star_small, scalefactor = 6000, verbose = T)
star_small <- addStatistics(gobject = star_small)


## 3. dimension reduction ####
star_small <- runPCA(gobject = star_small, method = 'factominer')
screePlot(star_small, ncp = 30)
plotPCA(gobject = star_small)

# 2D umap
star_small <- runUMAP(star_small, dimensions_to_use = 1:8)
plotUMAP(gobject = star_small)

# 3D umap
star_small <- runUMAP(star_small, dimensions_to_use = 1:8, name = '3D_umap', n_components = 3)
plotUMAP_3D(gobject = star_small, dim_reduction_name = '3D_umap')

# 2D tsne
star_small <- runtSNE(star_small, dimensions_to_use = 1:8)
plotTSNE(gobject = star_small)

## 4. clustering ####
star_small <- createNearestNetwork(gobject = star_small, dimensions_to_use = 1:8, k = 25)
star_small <- doLeidenCluster(gobject = star_small, resolution = 0.5, n_iterations = 1000)

# 2D umap
plotUMAP(gobject = star_small, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)

# 3D umap
plotUMAP_3D(gobject = star_small, dim_reduction_name = '3D_umap', cell_color = 'leiden_clus')


## 5. co-visualize ####

# 2D
spatDimPlot(gobject = star_small, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

# 3D
spatDimPlot3D(gobject = star_small, cell_color = 'leiden_clus', dim_reduction_name = '3D_umap')



## 6. differential expression ####
markers = findMarkers_one_vs_all(gobject = star_small,
                                 method = 'gini',
                                 expression_values = 'normalized',
                                 cluster_column = 'leiden_clus')

# violinplot
topgenes = markers[, head(.SD, 2), by = 'cluster']$genes
violinPlot(star_small, genes = topgenes, cluster_column = 'leiden_clus')

# genes heatmap
plotHeatmap(star_small, genes = star_small@gene_ID, cluster_column = 'leiden_clus',
            legend_nrows = 2, expression_values = 'scaled',
            cluster_order = 'correlation', gene_order = 'correlation')

# cluster heatmap
plotMetaDataHeatmap(star_small, expression_values = 'scaled', metadata_cols = c('leiden_clus'))



## 7. gene expression
dimGenePlot3D(star_small,
              dim_reduction_name = '3D_umap',
              expression_values = 'scaled',
              genes = "Pcp4",
              genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue')

spatGenePlot3D(star_small,
               expression_values = 'scaled',
               genes = "Pcp4",
               show_other_cells = F,
               genes_high_color = 'red', genes_mid_color = 'white', genes_low_color = 'darkblue')




## 8. cross section
# install.packages('geometry') # necessary for 3D delaunay network
star_small <- createSpatialNetwork(gobject = star_small, delaunay_method = 'delaunayn_geometry')

star_small = createCrossSection(star_small,
                                method="equation",
                                equation=c(0,1,0,600),
                                extend_ratio = 0.6)

# show cross section
insertCrossSectionSpatPlot3D(star_small, cell_color = 'leiden_clus',
                             axis_scale = 'cube',
                             point_size = 2)

insertCrossSectionGenePlot3D(star_small, expression_values = 'scaled',
                             axis_scale = "cube",
                             genes = "Slc17a7")


# for cell annotation
crossSectionPlot(star_small,
                 point_size = 2, point_shape = "border",
                 cell_color = "leiden_clus")

crossSectionPlot3D(star_small,
                   point_size = 2, cell_color = "leiden_clus",
                   axis_scale = "cube")


# for gene expression
crossSectionGenePlot(star_small,
                     genes = "Slc17a7",
                     point_size = 2,
                     point_shape = "border",
                     cow_n_col = 1.5,
                     expression_values = 'scaled')

crossSectionGenePlot3D(star_small,
                       point_size = 2,
                       genes = c("Slc17a7"),
                       expression_values = 'scaled')






