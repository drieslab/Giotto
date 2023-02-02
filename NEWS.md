# Giotto Suite 3.2.0 (2023-02-02)

## Breaking Changes
- Removed support for deprecated nesting in `@nn_network` slot
- `createSpatialNetwork()` will now output a `spatialNetworkObj` by default when `return_gobject = FALSE`. It is possible to change this back to the data.table output by setting `output = 'data.table'`
- Set incomplete classes in classes.R as virtual to prevent their instantiation
- Removed `createGiottoCosMxObject()` `aggregate` and `all` workflows until they are updated

## Added
- New `gefToGiotto()` interoperability function to convert gef object from Stereo-seq to giotto
- New `giottoToAnnData()` interoperability function to convert giotto object to squidpy flavor anndata .h5ad file(s)
- New `giottoToSpatialExperiment()` and `spatialExperimentToGiotto()` to convert between Giotto and SpatialExperiment
- New `spatialAutoCorLocal()` and `spatialAutoCorGlobal()` functions to find spatial autocorrelations from expression and cell metadata information
- New `createSpatialWeightMatrix()` function to generate spatial weight matrix from spatial networks for autocorrelation
- Add spatial_interaction_spot.R with functions adapted from master branch for working with the Giotto suite object.
- New exported accessors for slots (experimental)
- Add `coord_fix_ratio` param to `spatFeatPlot2D()` and `spatFeatPlot2D_single()`
- Add `order` parameter to `dimFeatPlot2D` and `spatDimFeatPlot2d` to plot and order cells according to the levels of the selected feature ([#477](https://github.com/drieslab/Giotto/issues/477))
- Add `plot()` method for `spatialNetworkObj`
- New `set_row_order_dt()` internal for setting `data.table` to a specific row order by reference
- New `fread_colmatch()` internal for fread loading a subset of rows based on matches in a specified column
- Add missing `create_nn_net_obj()` internal constructor function for S4 `nnNetObj`
- Add `id_col`, `x_col`, `y_col` params to `polyStamp()` to make stamp location input more flexible
- Add `optional` and `custom_msg` params to `package_check()`
- New `wrap()` and `vect()` generics for `giotto`, `giottoPoints`, and `giottoPolygons`
- New `rotate()`, `t()`, and `spatShift` generics for giotto subobject spatial manipulation
- New `spatIDs()` and `featIDs()` generics
- New `objName()` and `objName` generics for setting the names of relevant S4 subobjects
- New `rbind()` generic to append `giottoPolygon` objects
- Add packages `exactextractr` and `sf` to "suggests" packages
- Add package `progressr` to "imports" packages

## Changes
- Move giotto object method-specific creation functions from `giotto.R` to `convenience.R`
- Update `addFeatMetadata()` to handle replacement of existing columns
- Update `show()` method for `giotto`
- Update `show()` method for `spatEnrObj`
- Deprecate older snake_case accessors
- Deprecate `polygon_feat_names` param in favor of `z_indices` in `readPolygonFilesVizgenHDF5()`
- Deprecate `xy_translate_spatial_locations()` in favor of `shift_spatial_locations()`
- Optimize `readPolygonFilesVizgen()`
- Fix bug in `replaceGiottoInstructions()` where instructions with more slots than previous are not allowed
- Fix bug in `loadGiotto()` that prevents proper parsing of filenames when spat_unit or feat_type contains '_' characters
- Fix `loadGiotto()` loss of over-allocation for data.tables-based objects after loading from disk


# Giotto Suite 3.1.0 (2202-12-01)  
  

## Added
- New `initialize()` generic that calls `setalloccol()` for data.table-based S4 subobjects to allow setting by reference
- New `spatUnit`, `spatUnit<-`, `featType`, and `featType<-` feat type generics for S4 subobjects for setting the relevant slots
- Add `hexVertices()` to polygon shape array generation functionality

## Changes  
- Update `createGiottoCosMxObject()` for 3.0 and modularization of functions. 'subcellular' workflow has been tested to work along with an updated tutorial.
- Update grid plotting behavior to set a default number columns to use based on number of elements to plot. Can be overridden by explicitly providing input to `cow_n_col` param
- Fix bug in `annotateGiotto()` after 3.0 update ([#433](https://github.com/drieslab/Giotto/issues/433#issuecomment-1324211224))
- Fix bug in `joinGiottoObjects()` metadata processing
- Update seed setting behavior in [dimension_reduction.R](https://github.com/drieslab/Giotto/blob/suite/R/dimension_reduction.R) and [clustering.R](https://github.com/drieslab/Giotto/blob/suite/R/clustering.R)





# Giotto Suite 3.0.1 (2022-11-20)

## Added
- New system color support detection (based on crayon package logic)
- Add ability to turn off colored text in `show` functions with `options("giotto.color_show" = FALSE)`

## Changes
- Fix bug in `extract_polygon_list()` ([#433](https://github.com/drieslab/Giotto/issues/433#issuecomment-1321221382))
- Update Unicode character printing with `show` functions for Latin1 systems





# Giotto Suite 3.0.0 (2022-11-18)

## Breaking Changes
- S4 subobjects framework will require giotto objects to be remade

## Added
- New `createGiottoXeniumObject()` for loading 10x Xenium data
- New S4 subobjects. Details can be found in [classes.R](https://github.com/drieslab/Giotto/blob/suite/R/classes.R)
- New basic generics for S4 subobjects. Mainly the use of `[]` and `[]<-` to get or set information into the main data slot
- New `@provenance` slot in S4 subobjects to track provenance of aggregated information (z_layers used for example)
- New `calculateOverlapPolygonImages()` for calculating overlapped intensities from image-based information (e.g. IMC, IF, MIBI, ...) and polygon data (e.g. cell)
- New `overlapImagesToMatrix()` converts intensity-polygon overlap info into an expression matrix (e.g. cell by protein)
- New `aggregateStacks()` set of functions work with multiple subcellular layers when generating aggregated expression matrices

## Changes
- Update `setter` functions to read the `@spat_unit` and `@feat_type` slots of subobjects to determine nesting
- Update of `show` functions to display color coded nesting names and tree structure





# Giotto Suite 2.1.0 (2202-11-09)

## Breaking Changes
- Update of python version to **3.10.2** [details](https://giottosuite.readthedocs.io/en/latest/additionalinformation.html#giotto-suite-2-1-0-2202-11-09)

## Added
- New `anndataToGiotto()` to convert scanpy anndata to Giotto






# Giotto Suite 2.0.0.998

## Added
- New `GiottoData` package to work with spatial datasets associated with Giotto
	- Stores the minidatasets: preprocessed giotto objects that are ready to be used in any function
	- Moved: `getSpatialDataset()` and `loadGiottoMini()` functions to this package
- New `saveGiotto()` and `loadGiotto()` for preserving memory-pointer based objects. In [general_help.R](https://github.com/drieslab/Giotto/blob/suite/R/general_help.R)
  - It saves a Giotto object into a folder using a specific structure. Essentially a wrapper around `saveRDS()` that also works with spatVector and spatRaster pointers.
- New `plotInteractivePolygon()` for plot-interactive polygonal selection of points.
- New polygon shape array creation through `polyStamp()`, `circleVertices`, `rectVertices`. In [giotto_structures.R](https://github.com/drieslab/Giotto/blob/suite/R/giotto_structures.R)
- Add accessor functions `get_CellMetadata` (alias of `pDataDT()`), `set_CellMetadata`, `get_FeatMetadata` (alias of `fDataDT()`), `set_FeatMetadata`. See [accessors.R](https://github.com/drieslab/Giotto/blob/suite/R/accessors.R)
- New `filterDistributions()` to generate histogram plots from expression statistics

## Changes
- Deprecate `plotInteractionChangedGenes()` ,`plotICG()`, `plotCPG()` in favor of `plotInteractionChangedFeatures()` and `plotICF()` and `plotCPF()`
- Deprecate `plotCellProximityGenes()`, in favor of `plotCellProximityFeatures()`
- Deprecate `plotCombineInteractionChangedGenes()`, `plotCombineICG()`, `plotCombineCPG()` in favor of `plotCombineInteractionChangedFeatures()` and `plotCombineICF()`
- Deprecate `findInteractionChangedGenes()`, `findICG()`, `findCPG()` in favor of `findInteractionChangedFeats()` and `findICF`
- Deprecate `filterInteractionChangedGenes()`, `filterICG()`, `filterCPG()` in favor of `filterInteractionChangedFeats()` and `filterICF()`
- Deprecate `combineInteractionChangedGenes()`, `combineICG()`, `combineCPG()` in favor of `combineInteractionChangedFeats()` and `combineICF()`
- Deprecate `combineCellProximityGenes_per_interaction()` in favor of `combineCellProximityFeatures_per_interaction()`

## Breaking Changes
- ICF output internal object structure names have changed to use feats instead of genes












