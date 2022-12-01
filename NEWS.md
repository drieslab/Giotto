# Giotto Suite 3.1.0 (2202-12-01)  
  

## Added
- Add `initialize()` generic that calls `setalloccol()` for data.table-based S4 subobjects to allow setting by reference
- Add `spatUnit`, `spatUnit<-`, `featType`, and `featType<-` feat type generics for S4 subobjects for setting the relevant slots

## Changes  
- Update `createGiottoCosMxObject()` for 3.0 and modularization of functions. 'subcellular' workflow has been tested to work along with an updated tutorial.
- Update grid plotting behavior to set a default number columns to use based on number of elements to plot. Can be overridden by explicitly providing input to `cow_n_col` param
- Fix bug in `annotateGiotto()` after 3.0 update ([#455](https://github.com/drieslab/Giotto/issues/433#issuecomment-1324211224))
- Update seed setting behavior


# Giotto Suite 3.0.1 (2022-11-20)

## Added
- Add system color support detection (based on crayon package logic)
- Add ability to turn off colored text in `show` functions with `options("giotto.color_show" = FALSE)`

## Changes
- Fix bug in `extract_polygon_list()` ([#455](https://github.com/drieslab/Giotto/issues/433#issuecomment-1321221382))
- Update Unicode character printing with `show` functions for Latin1 systems


# Giotto Suite 3.0.0 (2022-11-18)

## Added
- New S4 subobjects. Details can be found in [classes.R](https://github.com/drieslab/Giotto/blob/suite_dev/R/classes.R)
- New basic generics for S4 subobjects. Mainly the use of `[]` and `[]<-` to get or set information into the main data slot
- New `@provenance` slot in S4 subobjects to track provenance of aggregated information (z_layers used for example)


## Changes
- Update of python environment to **3.10**
- Update `anndataToGiotto()`
- Update `setter` functions to read the `@spat_unit` and `@feat_type` slots of subobjects to determine nesting
- Update of `show` functions to display color coded nesting names and tree structure
