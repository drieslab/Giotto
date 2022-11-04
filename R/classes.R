# MISC ####
## * Define class unions ####

#' @title NULL or char class union
#' @description class to allow either NULL or character
#' @keywords internal
setClassUnion('nullOrChar', c('NULL', 'character'))

#' @title NULL or list class union
#' @description class to allow either NULL or list
#' @keywords internal
setClassUnion('nullOrList', c('NULL', 'list'))

#' @title NULL or data.table class union
#' @description class to allow either NULL or data.table
#' @keywords internal
setClassUnion('nullOrDatatable', c('NULL', 'data.table'))

## * Define external classes ####

#' data.table S4 class for method dispatch
#' @name data.table-class
#' @aliases data.table
#' @family data.table
#'
#' @exportClass data.table
setOldClass('data.table')






# BASIC CLASSES ####

# ** exprData Class ####
#' Basic class for classes with expression information
#'
setClass('exprData',
         representation = list(exprMat = 'ANY',
                               sparse = 'logical'),
         prototype = prototype(exprMat = NULL,
                               sparse = NA))



# ** coordData Class ####
#' Basic class for classes with coordinate information
#'
#' coordDataDT is the specific flavor that deals with objects where the coordinate
#' information is stored within data.table objects and should work similarly to
#' data.table when interacting with some basic generic operators for data
#' retreival and setting.
setClass('coordDataDT',
         representation = list(coordinates = 'data.table'),
         prototype = prototype(coordinates = data.table::data.table()))

# setClass('coordDataMT',
#          representation = list(coordinates = 'matrix'),
#          prototype = prototype(coordinates = matrix()))


# ** metaData Class ####
#' Basic class for classes with metadata information
#'
#' Classes that inherit from this class will contain a metadata slot that stores
#' information in a data.table and should work similarly to data.table when interacting
#' with some basic generic operators for data retrieval and setting
setClass('metaData',
         representation = list(metaDT = 'data.table',
                               col_desc = 'character'),
         prototype = prototype(metaDT = data.table::data.table(),
                               col_desc = NA_character_))


# ** nnData ####
setClass('nnData',
         representation = list(nn_type = 'character',
                               igraph = 'ANY'),
         prototype = prototype(nn_type = NA_character_,
                               igraph = NULL))




# ** provData Class ####
#' Basic class for classes with provenance information.
#'
#' This kind of information is necesssary when generating data that is aggregated
#' from multiple original sources of raw information. This could refer to situations
#' such as when producing cellxfeature expression matrices from subcellular transcript
#' information and polygons that are provided as multiple z layers. Provenance
#' is Giotto's method of mapping this aggregated information back to the original
#' z layers that were used in its generation.
setClass('provData',
         representation = list(provenance = 'ANY'),
         prototype = prototype(provenance = NULL))



# ** spatData Class ####
#' Basic class for classes with spatial information
#'
#' Classes that inherit from this class will contain a spat_unit slot that describes
#' which spatial unit the data belongs to. This is most relevant to aggregated information.
#' Subcellular information such as poly data in \code{spatial_info} slot essentially define their
#' own spatial units. Within slots that deal with classes that contain spatData,
#' there is a nesting structure that first nests by spatial unit.
#'
setClass('spatData',
         representation = list(spat_unit = 'character'), # not allowed to be NULL
         prototype = prototype(spat_unit = NA_character_))



# ** featData Class ####
#' @title Basic class for classes with feature information
#'
#' @description
#' Features in Giotto are a blanket term for any features that are detected, covering
#' modalities such as, but not limited to rna, protein, ATAC, and even QC probes.
#' Classes that inherit from this class will contain a feat_type slot that describes
#' which feature type the data is. Within slots that deal with classes that contain
#' featData, there is a nesting structure that usually first nests by spatial unit
#' and then by feature type
#'
setClass('featData',
         representation = list(feat_type = 'character'), # not allowed to be NULL
         prototype = prototype(feat_type = NA_character_))



# ** miscData Class ####
#' @title Basic class for additional miscellaneous information
#'
#' @description
#' Classes (such as dimObj) that can hold information from multiple types of methods
#' use the misc slot to hold additional information specific to each method.
#' Information may be stored within as S3 structures.
setClass('miscData',
         representation = list(misc = 'ANY'),
         prototype = prototype(misc = NULL))








# SUPER CLASSES ####

# ** spatFeatData ####
#' Superclass for classes that contain both spatial and feature data
setClass('spatFeatData',
         contains = c('spatData', 'featData'))









# CORE ####

## Giotto class ####


## * Check ####
# Giotto class

#TODO
# @title Check Giotto Object
# @name check_giottoObj
# @description check function for S4 giotto object
# @param object giotto object to check
# @keywords internal
# checkGiottoObj = function(object) {
#
#   # !!! (validity check still under construction) !!!
#   # Add check class definition when finished
#
#   errors = character()
#   ch = box_chars()
#
#   gDataSlots = list(expression = 'expression',
#                     expression_feat = 'expression_feat',
#                     spatial_locs = 'spatial_locs',
#                     spatial_info = 'spatial_info',
#                     feat_info = 'feat_info',
#                     cell_metadata = 'cell_metadata',
#                     feat_metadata = 'feat_metadata',
#                     cell_ID = 'cell_ID',
#                     feat_ID = 'feat_ID',
#                     spatial_network = 'spatial_network',
#                     spatial_grid = 'spatial_grid',
#                     spatial_enrichment = 'spatial_enrichment',
#                     dimension_reduction = 'dimension_reduction',
#                     nn_network = 'nn_network')
#
#
#   if(!is.null(slot(object, 'expression'))) {
#     avail_expr = list_expression(object)
#     for(expr in seq(nrow(avail_expr)))
#     if(!inherits(slot(object, 'expression')[[avail_expr[[expr]]]][[avail_expr[[expr]]]][[avail_expr[[expr]]]], 'data.table')) {
#       TODO
#     }
#   }
#
#
#
#
#
#   # Check for nesting issues (due to updates) - to be deprecated
#   available_nn = list_nearest_networks(object)
#   available_grids = list_spatial_grids(object)
#
#
#
#
#   ## Match spatial units and feature types ##
#
#   # Find existing spat_units in source data
#   uniqueSrcSU = list_expression(object)$spat_unit
#   uniqueSrcSU = unique(uniqueSrcSU, list_spatial_info_names(object))
#
#   # Find existing feat_types in source data
#   uniqueSrcFT = list_expression(object)$feat_type
#   uniqueSrcFT = unique(uniqueSrcFT, list_feature_info_names(object))
#
#   # Find existing spat_units and feat_types within other data slots
#   uniqueDatSU = lapply(gDataSlots, function(x) list_giotto_data(gobject = object, slot = x)$spat_unit)
#   uniqueDatFT = lapply(gDataSlots, function(x) list_giotto_data(gobject = object, slot = x)$feat_type)
#
#   # # If any spat_units exist, ensure that associated objects have identical cell_IDs
#   # if(length(uniqueSU) > 0) {
#   #   for(SU in seq_along(uniqueSU)) {
#   #
#   #     all(lapply(gDataSlots, function(x) {
#   #       list_giotto_data(gobject = object,
#   #                        slot = x)
#   #     }))
#   #   }
#   # }
#
#
#   # cell_ID slot
#   slot_cell_ID = slot(object, 'cell_ID')
#
#   ## check nesting depth
#   if(depth(slot_cell_ID) != 1) {
#     msg = paste0('Nesting depth of cell_ID is ', depth(slot_cell_ID), '.
#         Should be 1. Ex:
#         \t.
#         \t', ch$b,'spat_unit
#         \t', ch$s, ch$b,'cell IDs\n')
#     errors = c(errors, msg)
#   }
#
#   for(spat_unit in names(slot_cell_ID)) {
#     if(!inherits(slot_cell_ID[[spat_unit]], 'character')) {
#       msg = paste0('cell_ID slot info for spat_unit "', spat_unit,'" is of class "',
#                    class(slot_cell_ID[[spat_unit]]), '".\n    Should be "character"\n')
#       errors = c(errors, msg)
#     }
#   }
#
#
#   if(length(errors) == 0) TRUE else errors
# }




##### * Definition ####
# Giotto class

#' @title S4 giotto Class
#' @description Framework of giotto object to store and work with spatial expression data
#' @concept giotto object
#' @slot expression expression information
#' @slot expression_feat available features (e.g. rna, protein, ...)
#' @slot spatial_locs spatial location coordinates for cells/spots/grids
#' @slot spatial_info information about spatial units (Giotto spatVector)
#' @slot cell_metadata metadata for cells
#' @slot feat_metadata metadata for available features
#' @slot feat_info information about features (Giotto spatVector)
#' @slot cell_ID unique cell IDs
#' @slot feat_ID unique feature IDs for all features or modalities
#' @slot spatial_network spatial network in data.table/data.frame format
#' @slot spatial_grid spatial grid in data.table/data.frame format
#' @slot spatial_enrichment slot to save spatial enrichment-like results
#' @slot dimension_reduction slot to save dimension reduction coordinates
#' @slot nn_network nearest neighbor network in igraph format
#' @slot images slot to store giotto images
#' @slot largeImages slot to store giottoLargeImage objects
#' @slot parameters slot to save parameters that have been used
#' @slot instructions slot for global function instructions
#' @slot offset_file offset file used to stitch together image fields
#' @slot OS_platform Operating System to run Giotto analysis on
#' @slot join_info information about joined Giotto objects
#' @details
#' [\strong{expression}] There are several ways to provide expression information:
#'
#' [\strong{expression_feat}] The different features or modalities such as rna, protein, metabolites, ...
#' that are provided in the expression slot.
#'
#'
#' @export
giotto <- setClass(
  "giotto",
  slots = c(
    expression = "nullOrList",
    expression_feat = "ANY",
    spatial_locs = "ANY",
    spatial_info = "ANY",
    cell_metadata = "ANY",
    feat_metadata = "ANY",
    feat_info = "ANY",
    cell_ID = "ANY",
    feat_ID = "ANY",
    spatial_network = "ANY",
    spatial_grid = "ANY",
    spatial_enrichment = "ANY",
    dimension_reduction = 'ANY',
    nn_network = "ANY",
    images = "ANY",
    largeImages = "ANY",
    parameters = "ANY",
    instructions = "ANY",
    offset_file = "ANY",
    OS_platform = "ANY",
    join_info = "ANY"

  ),

  prototype = list(
    expression = NULL,
    expression_feat = NULL,
    spatial_locs = NULL,
    spatial_info = NULL,
    cell_metadata = NULL,
    feat_metadata = NULL,
    feat_info = NULL,
    cell_ID = NULL,
    feat_ID = NULL,
    spatial_network = NULL,
    spatial_grid = NULL,
    spatial_enrichment = NULL,
    dimension_reduction = NULL,
    nn_network = NULL,
    images = NULL,
    largeImages = NULL,
    parameters = NULL,
    instructions = NULL,
    offset_file = NULL,
    OS_platform = NULL,
    join_info = NULL
  )
)



##### * Show ####
# Giotto class

#' show method for giotto class
#' @param object giotto object
#' @aliases show,giotto-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giotto",
  definition = function(object) {

    cat("An object of class",  class(object), "\n")

    for(spat_unit in names(object@expression_feat)) {
      cat('spatial units = ', spat_unit, '\n')
      for(feat_type in unique(object@expression_feat)) {
        cat("features = ", feat_type, "\n")

        cat(
          nrow(x = object@expression[[spat_unit]][[feat_type]][['raw']]),
          "features across",
          ncol(x = object@expression[[spat_unit]][[feat_type]][['raw']]),
          "samples.\n \n"
        )
      }
    }


    cat('Steps and parameters used: \n \n')
    print(object@parameters)
    invisible(x = NULL)
  }
)












# EXPRESSION ####

## exprObj Class ####

## * Check ####
# exprObj Class

#' @title Check exprObj
#' @name check_expr_obj
#' @description Check function for S4 exprObj
#' @param object S4 exprObj to check
#' @keywords internal
check_expr_obj = function(object) {
  errors = character()

  # Check for expr info
  if(is.null(slot(object, 'exprMat'))) {
    obj_info = paste0('exprObj ',
                      'spat_unit "', slot(object, 'spat_unit'), '", ',
                      'feat_type "', slot(object, 'feat_type'), '", ',
                      'name "', slot(object, 'name'), '": \n')

    msg = paste0(obj_info, 'No expression information found.\n')
    errors = c(errors, msg)
  }

  # TODO Check expr matrix class
  # if(!inherits(slot(object, 'exprMat'), c('dgeMatrix', 'dgCMatrix'))) {
  #   msg = paste0(obj_info, 'Expression matrix should be provided as either dgCmatrix (sparse)',
  #                ' or dgeMatrix (dense) from the Matrix package.\n')
  #   errors = c(errors, msg)
  # }
  #
  # if(inherits(slot(object, 'exprMat'), 'dgCMatrix')) {
  #   if(!isTRUE(slot(object, 'sparse'))) {
  #     msg = paste0(obj_info, 'Contains a dgCmatrix. Slot sparse should be TRUE')
  #     errors = c(errors, msg)
  #   }
  # }
  #
  # if(inherits(slot(object, 'exprMat'), 'dgeMatrix')) {
  #   if(isTRUE(slot(object, 'sparse'))) {
  #     msg = paste0(obj_info, 'Contains a dgematrix. Slot sparse should be FALSE')
  #     errors = c(errors, msg)
  #   }
  # }

  if(length(errors) == 0) TRUE else errors
}



## * Definition ####
# exprObj Class

#' @title S4 exprObj
#' @description Framework to store aggregated expression information
#' @slot name name of exprObj
#' @slot expression matrix of expression information
#' @slot spat_unit spatial unit of expression (e.g. 'cell')
#' @slot feat_type feature type of expression (e.g. 'rna', 'protein')
#' @slot provenance origin data of expression information (if applicable)
#' @slot misc misc
#' @export
setClass('exprObj',
         contains = c('exprData', 'spatFeatData', 'provData', 'miscData'),
         slots = c(name = 'nullOrChar'),
         prototype = list(name = NULL),
         validity = check_expr_obj)

## * Show ####
# exprObj Class

#' show method for exprObj class
#' @param object expression object
#' @aliases show,exprObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('exprObj'), function(object) {

    # Print if data is sparse
    cat("An object of class",  class(object), "\n\n")
    if(!is.null(slot(object, 'sparse'))) {
      if(isTRUE(slot(object, 'sparse'))) {
        cat('Contains sparse matrix with aggregated expression information\n')
      } else if(isFALSE(slot(object, 'sparse'))) {
        cat('Contains dense matrix with aggregated expression information\n')
      }
    }

    # print spat/feat and provenance info
    cat(paste0('for spatial unit: "', slot(object, 'spat_unit'), '" and feature type: "', slot(object, 'feat_type'),'" \n'))
    if(!is.null(slot(object, 'provenance'))) cat('  Provenance: ', unlist(slot(object, 'provenance')),'\n')

    # preview matrix
    print(slot(object, 'exprMat'))

    cat('\n')

  })





# METADATA ####


## cellMetaObj class ####

# * Check ####
#' @title Check cell metadata object
#' @name check_cell_meta_obj
#' @description Function to check S4 cellMetaObj
#' @param object S4 cellMetaObj to check
#' @keywords internal
check_cell_meta_obj = function(object) {

  errors = character()

  if(!'cell_ID' %in% colnames(object@metaDT)) {
    msg = 'No "cell_ID" column found.'
    errors = c(errors, msg)
  } else {

    if(!is.character(object@metaDT[['cell_ID']])) {
      msg = '"cell_ID" column must be of class character.'
    }

  }
  if(length(errors) == 0) TRUE else errors
}

# * Definition ####
#' @title S4 cellMetaObj
#' @description Framework to store cell metadata
#' @slot metadata metadata info
#' @slot col_desc (optional) character vector describing columns of the metadata
#' @slot spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @slot feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @slot provenance origin data of aggregated expression information (if applicable)
#' @slot misc misc
#' @export
setClass('cellMetaObj',
         contains = c('metaData', 'spatFeatData', 'provData'),
         validity = check_cell_meta_obj)


setMethod('show', signature('cellMetaObj'), function(object) {

  cat('An object of class', class(object), '\n')
  # TODO print spat/feat/prov
  if(!is.null(slot(object, 'metaDT'))) print(slot(object, 'metaDT'))

})


## featMetaObj class ####

# * Check ####
#' @title Check feature metadata object
#' @name check_feat_meta_obj
#' @description Function to check S4 featMetaObj
#' @param object S4 featMetaObj to check
#' @keywords internal
check_feat_meta_obj = function(object) {

  errors = character()

  if(!'feat_ID' %in% colnames(object@metaDT)) {
    msg = 'No "feat_ID" column found.'
    errors = c(errors, msg)
  } else {

    if(!is.character(object@metaDT[['feat_ID']])) {
      msg = '"feat_ID" column must be of class character.'
    }

  }
  if(length(errors) == 0) TRUE else errors
}

# * Definition ####
#' @title S4 featMetaObj
#' @description Framework to store feature metadata
#' @slot metadata metadata info
#' @slot col_desc (optional) character vector describing columns of the metadata
#' @slot spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @slot feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @slot provenance origin data of aggregated expression information (if applicable)
#' @slot misc misc
#' @export
setClass('featMetaObj',
         contains = c('metaData', 'spatFeatData', 'provData'),
         validity = check_feat_meta_obj)


setMethod('show', signature('featMetaObj'), function(object) {

  cat('An object of class', class(object), '\n')
  # TODO print spat/feat/prov
  if(!is.null(slot(object, 'metaDT'))) print(slot(object, 'metaDT'))

})


# DIMENSION REDUCTION ####

## dimObj Class ####



##### * Check #####
# dimObj Class

#' @title Check dimOjb
#' @name check_dim_obj
#' @description check function for S4 dimObj
#' @param object S4 dimObj to check
#' @keywords internal
check_dim_obj = function(object) {
  errors = character()
  length_reduction_method = length(object@reduction_method)
  if(length_reduction_method > 1) {
    msg = paste0('reduction_method is length ', length_reduction_method, '. Should be 1')
    errors = c(errors, msg)
  }

  if(length_reduction_method == 0) {
    msg = 'A reduction_method must be given'
    errors = c(errors, msg)
  }

  lastCols = tail(colnames(object@coordinates),2)
  col_dims = all(grepl(pattern = 'Dim.', x = lastCols))
  if(!isTRUE(col_dims)) {
    msg = 'Dim reduction coordinates should be provided with dimensions ("Dim.#") as columns and samples as rows\n'
    errors = c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}



## * Definition ####
# dimObj Class

#' @title S4 dimObj Class
#' @description Framework to store dimension reduction information
#' @slot name name of dimObject
#' @slot feat_type feature type of data
#' @slot spat_unit spatial unit of data
#' @slot provenance origin of aggregated information (if applicable)
#' @slot reduction_method method used to generate dimension reduction
#' @slot coordinates embedding coordinates
#' @slot misc method-specific additional outputs
#' @export
setClass('dimObj',
         contains = c('spatFeatData'),
         slots = c(name = 'character',
                   reduction_method = 'character',
                   coordinates = 'ANY',
                   misc = 'ANY'),
         prototype = list(name = NA_character_,
                          reduction_method = NA_character_,
                          coordinates = NULL,
                          misc = NULL),
         validity = check_dim_obj)



##### * Show ####
# dimObj Class

#' show method for dimObj class
#' @param object dimension reduction object
#' @aliases show,dimObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('dimObj'), function(object) {

    cat("An object of class",  class(object), "\n")
    if(!is.null(object@reduction_method)) cat('--| Contains dimension reduction generated with:', object@reduction_method, '\n')
    if(!is.null(object@feat_type) & !is.null(object@spat_unit)) {
      cat('----| for feat_type:', object@feat_type, '\n')
      cat('----|     spat_unit:', object@spat_unit, '\n\n')
    }

    if(!is.null(object@coordinates)) cat('  ', ncol(object@coordinates), 'dimensions for', nrow(object@coordinates),'data points\n\n')

    if(!is.null(object@misc)) {
      cat('Additional included info:\n')
      print(names(object@misc))
      cat('\n')
    }

  })


## * Additional functions ####
# dimObj Class

#' @title Dimension reductions
#' @name S3toS4dimObj
#' @description Convert S3 dimObj to S4
#' @param object S3 dimObj
#' @keywords internal
S3toS4dimObj = function(object) {
  if(!isS4(object)) {
    object = new('dimObj',
                 name = object$name,
                 feat_type = object$feat_type,
                 spat_unit = object$spat_unit,
                 reduction_method = object$reduction_method,
                 coordinates = object$coordinates,
                 misc = object$misc)
  }
  object
}





## nnNetObj ####

setClass('nnNetObj',
         contains = c('nnData', 'spatFeatData', 'provData', 'miscData'),
         representation = list(name = 'character'),
         prototype = prototype(name = NA_character_))








# SPATIAL ####



## spatLocsObj Class ####

## * check ####
# spatLocsObj Class

#' @title Check spatLocsObj
#' @name check_spat_locs_obj
#' @description Check function for S4 spatLocsObj
#' @param object S4 spatLocsObj to check
#' @keywords internal
check_spat_locs_obj = function(object) {
  errors = character()

  if(!'sdimx' %in% colnames(slot(object, 'coordinates'))) {
    msg = 'Column "sdimx" for x spatial location was not found'
    errors = c(errors, msg)
  }

  if(!'sdimy' %in% colnames(slot(object, 'coordinates'))) {
    msg = 'Column "sdimy" for y spatial location was not found'
    errors = c(errors, msg)
  }

  if(!'cell_ID' %in% colnames(slot(object, 'coordinates'))) {
    msg = 'Column "cell_ID" for cell ID was not found'
    errors = c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}


## * definition ####
# spatLocsObj Class

#' @title S4 spatLocsObj Class
#' @description Framework to store spatial location information
#' @slot name name of spatLocsObj
#' @slot coordinates data.table of spatial coordinates/locations
#' @slot spat_unit spatial unit tag
#' @slot provenance origin of aggregated information (if applicable)
#' @export
setClass('spatLocsObj',
         contains = c('coordDataDT', 'spatData', 'provData', 'miscData'),
         slots = c(name = 'nullOrChar'),
         prototype = list(name = NULL),
         validity = check_spat_locs_obj)




# * show ####
# spatLocsObj Class

#' show method for spatLocsObj class
#' @param object spatial locations object
#' @aliases show,spatLocsObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('spatLocsObj'), function(object) {

    sdimx = sdimy = NULL

    cat("An object of class",  class(object), "\n")
    if(!is.null(slot(object, 'spat_unit'))) cat(paste0('for spatial unit: "', slot(object, 'spat_unit'), '"\n'))
    if(!is.null(slot(object, 'provenance'))) cat(paste0('provenance: ', slot(object, 'provenance'), '\n'))

    cat('   ------------------------\n\npreview:\n')
    if(!is.null(slot(object, 'coordinates'))) show(slot(object, 'coordinates'))

    cat('\nranges:\n')
    try(expr = print(sapply(slot(object, 'coordinates')[,.(sdimx,sdimy)], range)),
        silent = TRUE)

    cat('\n\n')

  })





## spatialNetworkObj Class ####

### * check ####
# spatialNetworkObj Class

#' @title Check spatialNetworkObj
#' @name check_spat_net_obj
#' @description Check function for S4 spatialNetworkObj
#' @param object S4 spatialNetworkObj to check
#' @keywords internal
check_spat_net_obj = function(object) {
  errors = character()
  method_slot = slot(object, 'method')
  length_method = length(method_slot)
  if(length_method > 1) {
    msg = paste0('method is length ', length_method, '. Should be 1')
    errors = c(errors, msg)
  }

  # if(is.null(method_slot)) {
  #   msg = 'A spatial network generation method must be given'
  #   errors = c(errors, msg)
  # }

  if(is.null(object@networkDT) & is.null(object@networkDT_before_filter)) {
    msg = 'No data in either networkDT or networkDT_before_filter slots.\nThis object contains no network information.\n'
    errors = c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}



### * definition ####
# spatialNetworkObj Class

#' @title S4 spatialNetworkObj Class
#' @description Framework to store spatial network information
#' @slot name name of spatialNetworkObj
#' @slot method method used to generate spatial network
#' @slot parameters additional method-specific parameters used during spatial network generation
#' @slot outputObj network geometry object
#' @slot networkDT data.table of network connections, distances, and weightings
#' @slot networkDT_before_filter unfiltered data.table  of network connections, distances, and weightings
#' @slot cellShapeObj network cell shape information
#' @slot crossSectionObjects crossSectionObjects (see \code{\link{create_crossSection_object}})
#' @slot spat_unit spatial unit tag
#' @slot provenance origin of aggregated information (if applicable)
#' @slot misc misc
#' @details The generic access operators work with the data within the \code{networkDT}
#' slot (filtered).
#' @export
setClass('spatialNetworkObj',
         contains = c('spatData', 'provData', 'miscData'),
         slots = c(name = 'nullOrChar',
                   method = 'nullOrChar',
                   parameters = 'nullOrList',
                   outputObj = 'ANY',
                   networkDT = 'nullOrDatatable',
                   networkDT_before_filter = 'nullOrDatatable',
                   cellShapeObj = 'ANY',
                   crossSectionObjects = 'ANY'),
         prototype = list(name = NULL,
                          method = NULL,
                          parameters = NULL,
                          outputObj = NULL,
                          networkDT = NULL,
                          networkDT_before_filter = NULL,
                          cellShapeObj = NULL,
                          crossSectionObjects = NULL),
         validity = check_spat_net_obj)


### * show ####
# spatialNetworkObj Class

#' show method for spatialNetworkObj class
#' @param object spatial network object
#' @aliases show,spatialNetworkObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('spatialNetworkObj'), function(object) {

    cat("An object of class",  class(object), "\n")
    if(!is.null(object@method)) cat('Contains spatial network generated with:', object@method, '\n')
    if(!is.na(object@spat_unit)) cat(paste0('for spatial unit: "', object@spat_unit, '"\n'))
    if(!is.na(object@provenance)) cat(paste0('provenance: "', object@provenance, '"\n'))

    if(!is.null(object@networkDT)) cat('  ', nrow(object@networkDT), 'connections (filtered)\n')
    if(!is.null(object@networkDT_before_filter)) cat('  ', nrow(object@networkDT_before_filter), 'connections (before filter)\n\n')

  })



### * Additional functions ####
# spatialNetworkObj Class

# S3 to S4 backwards compatibility

#' @title Spatial Networks
#' @name S3toS4spatNetObj
#' @description convert S3 spatialNetworkObj to S4
#' @param object S3 spatNetworkObj
#' @param spat_unit spatial unit metadata to append
#' @keywords internal
S3toS4spatNetObj = function(object,
                            spat_unit = NULL) {
  if(!isS4(object)) {
    object = new('spatialNetworkObj',
                 name = object$name,
                 method = object$method,
                 parameters = object$parameters,
                 outputObj = object$outputObj,
                 networkDT = object$networkDT,
                 networkDT_before_filter = object$networkDT_before_filter,
                 cellShapeObj = object$cellShapeObj,
                 crossSectionObjects = object$crossSectionObjects,
                 spat_unit = spat_unit,
                 misc = object$misc)
  }
  object
}




## crossSectionObj class ####
# See cross_section.R
# TODO







## spatialGridObj Class ####

### * check ####
# spatialGridObj Class

#' @title Check spatialGridObj
#' @name check_spat_grid_obj
#' @description Check function for S4 spatialGridObj
#' @param object S4 spatialGridObj to check
#' @keywords internal
check_spat_grid_obj = function(object) {
  errors = character()
  method_slot = slot(object, 'method')
  length_method = length(method_slot)
  if(length_method > 1) {
    msg = paste0('method is length ', length_method, '. Should be 1')
    errors = c(errors, msg)
  }

  # if(is.null(method_slot)) {
  #   msg = 'A grid generation method must be given'
  #   errors = c(errors, msg)
  # }

  if(is.null(object@gridDT)) {
    msg = 'No data in gridDT slot.\nThis object contains no spatial grid information\n'
    errors = c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}



### * definition ####
# spatialGridObj Class

#' @title S4 spatialGridObj Class
#' @description Framework to store spatial grid
#' @slot name name of spatialGridObj
#' @slot method method used to generate spatial grid
#' @slot parameters additional method-specific parameters used during spatial grid generation
#' @slot gridDT data.table holding the spatial grid information
#' @slot spat_unit spatial unit
#' @slot feat_type feature type
#' @slot provenance origin of aggregated information (if applicable)
#' @slot misc misc
#' @details
#' This is an S4 object that defines a spatial grid. The structure of the grid is stored as a
#' \code{data.table} within the \code{gridDT} slot and is defined by start and stop spatial
#' locations along the spatial axes. The \code{data.table} also includes names for each cell
#' of the grid and names for each of the spatial axis locations that make up the cell.
#' Grids can be annotated with both spatial and feature information
#' @export
setClass('spatialGridObj',
         slots = c(name = 'nullOrChar',
                   method = 'nullOrChar',
                   parameters = 'nullOrList',
                   gridDT = 'data.table',
                   spat_unit = 'nullOrChar',
                   feat_type = 'nullOrChar',
                   provenance = 'ANY',
                   misc = 'ANY'),
         prototype = list(name = NULL,
                          method = NULL,
                          parameters = NULL,
                          gridDT = NULL,
                          spat_unit = NULL,
                          feat_type = NULL,
                          provenance = NULL,
                          misc = NULL),
         validity = check_spat_grid_obj)



### * show ####
# spatialGridObj Class

#' show method for spatialGridObj class
#' @param object spatial grid object
#' @aliases show,spatialGridObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('spatialGridObj'), function(object) {

    # define for data.table
    x_start = x_end = y_start = y_end = z_start = z_end = NULL

    cat("An object of class",  class(object), "\n")
    if(!is.null(slot(object, 'spat_unit'))) cat('Contains annotations for spatial unit: "', slot(object, 'spat_unit'),'"', sep = '')
    if(!is.null(slot(object, 'feat_type'))) {
      cat(' and feature type: "', slot(object, 'feat_type'), '"\n', sep = '')
    } else cat('\n')

    # find grid spatial extent
    gridNames = colnames(slot(object, 'gridDT'))
    sdimx_max = slot(object, 'gridDT')[, max(x_start, x_end)]
    sdimx_min = slot(object, 'gridDT')[, min(x_start, x_end)]
    sdimy_max = slot(object, 'gridDT')[, max(y_start, y_end)]
    sdimy_min = slot(object, 'gridDT')[, min(y_start, y_end)]
    sdimx_uniques = slot(object, 'gridDT')[, length(unique(x_start))]
    sdimy_uniques = slot(object, 'gridDT')[, length(unique(y_start))]
    cat('Contains spatial grid defined for:\n  ', sdimx_uniques,'intervals from x range:',
        sdimx_min, 'to', sdimx_max, '\n  ' ,sdimy_uniques,'intervals from y range:',
        sdimy_min, 'to', sdimy_max)
    if('z_start' %in% gridNames & 'z_end' %in% gridNames) {
      sdimz_max = slot(object, 'gridDT')[, max(z_start, z_end)]
      sdimz_min = slot(object, 'gridDT')[, min(z_start, z_end)]
      sdimz_uniques = slot(object, 'gridDT')[, length(unique(z_start))]
      cat('\n  ', sdimz_uniques, 'intervals from z range:', sdimz_min, 'to', sdimz_max, '\n\n')
    } else {
      cat('\n\n')
    }

    if(!is.null(slot(object, 'method'))) cat('Contains spatial grid generated with:', slot(object, 'method'), '\n\n')

    if(!is.null(slot(object, 'parameters'))) {
      cat('Parameters used:\n')
      for(param in names(slot(object, 'parameters'))) {
        cat(paste0('  ',param, ': ', slot(object, 'parameters')[[param]], '\n'))
      }
      cat('\n')
    }

    if(!is.null(slot(object, 'misc'))) {
      cat('Additional included info:\n')
      print(names(slot(object, 'misc')))
      cat('\n')
    }
  })



## * Additional functions ####
# spatialGridObj Class

# S3 to S4 backwards compatibility

#' @title Spatially Binned Data
#' @name S3toS4spatGridObj
#' @description convert S3 spatialGridObj to S4
#' @param object S3 spatialGridObj
#' @keywords internal
S3toS4spatialGridObj = function(object) {
  if(!isS4(object)) {
    object = new('spatialGridObj',
                 name = object$name,
                 method = object$method,
                 parameters = object$parameters,
                 gridDT = object$gridDT,
                 misc = object$misc)
  }
  object
}



# SUBCELLULAR ####

## giottoPolygon class ####

# * definition ####
# giottoPolygon class

#' @title S4 giotto polygon Class
#' @description Giotto class to store and operate on polygon-like data
#' @concept giotto polygon class
#' @slot name name of polygon shapes
#' @slot spatVector terra spatVector to store polygon shapes
#' @slot spatVectorCentroids centroids of polygon shapes
#' @slot overlaps information about overlapping points and polygons
#' @details holds polygon data
#'
#' @export
giottoPolygon <- setClass(
  Class = "giottoPolygon",

  slots = c(
    name = "ANY",
    spatVector = "ANY",
    spatVectorCentroids = "ANY",
    overlaps = "ANY"
  ),

  prototype = list(
    name = NULL,
    spatVector = NULL,
    spatVectorCentroids = NULL,
    overlaps = NULL
  )
)



## giottoPoints class ####


## * definition ####
# giottoPoints class

#' @title S4 giotto points Class
#' @description Giotto class to store and operate on points data
#' @concept giotto points class
#' @slot feat_type name of feature type
#' @slot spatVector terra spatVector to store point shapes
#' @slot networks feature networks
#' @details Contains vector-type feature data
#'
#' @export
giottoPoints <- setClass(
  Class = "giottoPoints",

  slots = c(
    feat_type = "ANY",
    spatVector = "ANY",
    networks = "ANY"
  ),

  prototype = list(
    feat_type = NULL,
    spatVector = NULL,
    networks = NULL
  )
)


## featureNetwork class ####


## * definition ####
# featureNetwork class


#' @title S4 giotto feature network Class
#' @description Giotto class to store and operate on feature network
#' @concept giotto points network class
#' @slot name name of feature network
#' @slot network_datatable feature network in data.table format
#' @slot network_lookup_id table mapping numeric network ID to unique feature numerical IDs
#' @slot full fully connected network
#' @details contains feature network information
#'
#' @export
featureNetwork <- setClass(
  Class = "featureNetwork",

  slots = c(
    name = "ANY",
    network_datatable = "ANY",
    network_lookup_id = "ANY",
    full = "ANY"
  ),

  prototype = list(
    name = NULL,
    network_datatable = NULL,
    network_lookup_id = NULL,
    full = NULL
  )
)


# IMAGES ####

## giottoImage class ####

# * definition ####
# giottoImage class

#' @title S4 giottoImage Class
#' @description Framework of giotto object to store and work with spatial expression data
#' @concept giotto image object
#' @slot name name of Giotto image
#' @slot mg_object magick image object
#' @slot minmax minimum and maximum of associated spatial location coordinates
#' @slot boundaries x and y coordinate adjustments (default to 0)
#' @slot scale_factor image scaling relative to spatial locations
#' @slot resolution spatial location units covered per pixel
#' @slot file_path file path to the image if given
#' @slot OS_platform Operating System to run Giotto analysis on
#' @details
#' [\strong{mg_object}] Core object is any image that can be read by the magick package
#'
#' [\strong{boundaries}] Boundary adjustments can be used to manually or
#' automatically through a script adjust the image with the spatial data.
#'
#'
#' @export
giottoImage <- setClass(
  Class = "giottoImage",

  slots = c(
    name = "ANY",
    mg_object = "ANY",
    minmax = "ANY",
    boundaries = "ANY",
    scale_factor = "ANY",
    resolution = "ANY",
    file_path = "ANY",
    OS_platform = "ANY"
  ),

  prototype = list(
    name = NULL,
    mg_object = NULL,
    minmax = NULL,
    boundaries = NULL,
    scale_factor = NULL,
    resolution = NULL,
    file_path = NULL,
    OS_platform = NULL
  )
)



# * show ####
# giottoImage class

#' show method for giottoImage class
#' @param object giottoImage object
#' @aliases show,giottoImage-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giottoImage",
  definition = function(object) {

    cat("An object of class '",  class(object), "' with name ", object@name, "\n \n")

    cat("Min and max values are: \n",
        "Max on x-axis: ", object@minmax[['xmax_sloc']], "\n",
        "Min on x-axis: ", object@minmax[['xmin_sloc']], "\n",
        "Max on y-axis: ", object@minmax[['ymax_sloc']], "\n",
        "Min on y-axis: ", object@minmax[['ymin_sloc']], "\n",
        "\n")

    cat("Boundary adjustment are: \n",
        "Max adjustment on x-axis: ", object@boundaries[['xmax_adj']], "\n",
        "Min adjustment on x-axis: ", object@boundaries[['xmin_adj']], "\n",
        "Max adjustment on y-axis: ", object@boundaries[['ymax_adj']], "\n",
        "Min adjustment on y-axis: ", object@boundaries[['ymin_adj']], "\n",
        "\n")

    cat("Boundaries are: \n",
        "Image x-axis max boundary: ", object@minmax[['xmax_sloc']] + object@boundaries[['xmax_adj']], "\n",
        "Image x-axis min boundary: ", object@minmax[['xmin_sloc']] - object@boundaries[['xmin_adj']], "\n",
        "Image y-axis max boundary: ", object@minmax[['ymax_sloc']] + object@boundaries[['ymax_adj']], "\n",
        "Image y-axis min boundary: ", object@minmax[['ymin_sloc']] - object@boundaries[['ymin_adj']], "\n",
        "\n")

    cat("Scale factor: \n")
    print(object@scale_factor)

    cat("\n Resolution: \n")
    print(object@resolution)

    cat("\n File Path: \n")
    print(object@file_path)

    # print(object@mg_object)

  }
)


## giottoLargeImage class ####


## * definition ####
# giottoLargeImage class

#' @title S4 giottoLargeImage Class
#' @description class to handle images too large to load in normally through magick
#' @concept giotto object image
#' @slot name name of large Giotto image
#' @slot raster_object terra raster object
#' @slot extent tracks the extent of the raster object. Note that most processes should rely on the extent of the raster object instead of this.
#' @slot overall_extent terra extent object covering the original extent of image
#' @slot scale_factor image scaling relative to spatial locations
#' @slot resolution spatial location units covered per pixel
#' @slot max_intensity value to set as maximum intensity in color scaling
#' @slot min_intensity minimum value found
#' @slot is_int values are integers
#' @slot file_path file path to the image if given
#' @slot OS_platform Operating System to run Giotto analysis on
#' @export
giottoLargeImage <- setClass(
  Class = "giottoLargeImage",

  slots = c(
    name = "ANY",
    raster_object = "ANY",
    extent = "ANY",
    overall_extent = "ANY",
    scale_factor = "ANY",
    resolution = "ANY",
    max_intensity = "ANY",
    min_intensity = "ANY",
    is_int = "ANY",
    file_path = "ANY",
    OS_platform = "ANY"
  ),

  prototype = list(
    name = NULL,
    raster_object = NULL,
    extent = NULL,
    overall_extent = NULL,
    scale_factor = NULL,
    resolution = NULL,
    max_intensity = NULL,
    min_intensity = NULL,
    is_int = NULL,
    file_path = NULL,
    OS_platform = NULL
  )
)


## * show ####
# giottoLargeImage class

#' show method for giottoLargeImage class
#' @param object giottoLargeImage object
#' @aliases show,giottoLargeImage-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giottoLargeImage",
  definition = function(object) {

    cat("An object of class '",  class(object), "' with name ", object@name, "\n \n")

    cat("Image boundaries are: \n")
    print(terra::ext(object@raster_object)[1:4])

    cat("Original image boundaries are: \n")
    print(object@overall_extent[1:4])

    cat("\n Scale factor: \n")
    print(object@scale_factor)

    cat("\n Resolution: \n")
    print(object@resolution)

    cat('\n Estimated maximum intensity is: ', object@max_intensity, ' \n',
        'Estimated minimum intensity is: ', object@min_intensity, ' \n')

    if(object@is_int == TRUE) cat('Values are integers')
    if(object@is_int == FALSE) cat('Values are floating point')

    cat('\n File path is: ', object@file_path, '\n')

  }
)





