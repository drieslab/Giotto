# MISC ####
## * Define class unions ####

#' @title NULL or char class union
#' @description class to allow either NULL or character
#' @keywords internal
#' @noRd
setClassUnion('nullOrChar', c('NULL', 'character'))

#' @title NULL or list class union
#' @description class to allow either NULL or list
#' @keywords internal
#' @noRd
setClassUnion('nullOrList', c('NULL', 'list'))

#' @title NULL or data.table class union
#' @description class to allow either NULL or data.table
#' @keywords internal
#' @noRd
setClassUnion('nullOrDatatable', c('NULL', 'data.table'))

## * Define external classes ####

#' data.table S4 class for method dispatch
#' @name data.table-class
#' @aliases data.table
#' @family data.table
#' @exportClass data.table
#' @noRd
setOldClass('data.table')






# VIRTUAL CLASSES ####

# ** nameData Class ####
#'
setClass('nameData',
         contains = 'VIRTUAL',
         slots = list(name = 'character'),
         prototype = prototype(name = NA_character_))

# ** exprData Class ####
#' Basic class for classes with expression information
#'
setClass('exprData',
         contains = 'VIRTUAL',
         slots = list(exprMat = 'ANY'),
         prototype = prototype(exprMat = NULL))



# ** coordData Class ####
#' Basic class for classes with coordinate information
#'
#' coordDataDT is the specific flavor that deals with objects where the coordinate
#' information is stored within data.table objects and should work similarly to
#' data.table when interacting with some basic generic operators for data
#' retreival and setting.
setClass('coordDataDT',
         contains = 'VIRTUAL',
         slots = list(coordinates = 'data.table'),
         prototype = prototype(coordinates = data.table::data.table()))


# * Initialize ####
setMethod('initialize', 'coordDataDT',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@coordinates)) {
              .Object@coordinates = data.table::setalloccol(.Object@coordinates)
            }
            .Object
          })


# setClass('coordDataMT',
#          slots = list(coordinates = 'matrix'),
#          prototype = prototype(coordinates = matrix()))


# ** metaData Class ####
#' Basic class for classes with metadata information
#'
#' Classes that inherit from this class will contain a metadata slot that stores
#' information in a data.table and should work similarly to data.table when interacting
#' with some basic generic operators for data retrieval and setting
setClass('metaData',
         contains = 'VIRTUAL',
         slots = list(metaDT = 'data.table',
                      col_desc = 'character'),
         prototype = prototype(metaDT = data.table::data.table(),
                               col_desc = NA_character_))


setMethod('initialize', 'metaData',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@metaDT)) {
              .Object@metaDT = data.table::setalloccol(.Object@metaDT)
            }
            .Object
          })


# ** enrData ####
#' enrData
setClass('enrData',
         contains = 'VIRTUAL',
         slots = list(method = 'character',
                      enrichDT = 'nullOrDatatable'),
         prototype = prototype(method = NA_character_,
                               enrichDT = NULL))

setMethod('initialize', 'enrData',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@enrichDT)) {
              .Object@enrichDT = data.table::setalloccol(.Object@enrichDT)
            }
            .Object
          })



# ** nnData ####
setClass('nnData',
         contains = 'VIRTUAL',
         slots = list(nn_type = 'character',
                      igraph = 'ANY'),
         prototype = prototype(nn_type = NA_character_,
                               igraph = NULL))


# ** spatNetData ####
setClass('spatNetData',
         contains = 'VIRTUAL',
         slots = list(method = 'character',
                      parameters = 'ANY',
                      outputObj = 'ANY',
                      networkDT = 'nullOrDatatable',
                      networkDT_before_filter = 'nullOrDatatable',
                      cellShapeObj = 'ANY'),
         prototype = prototype(method = NA_character_,
                               parameters = NULL,
                               outputObj = NULL,
                               networkDT = NULL,
                               networkDT_before_filter = NULL,
                               cellShapeObj = NULL))

setMethod('initialize', 'spatNetData',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@networkDT)) {
              .Object@networkDT = data.table::setalloccol(.Object@networkDT)
            }
            if(!is.null(.Object@networkDT_before_filter)) {
              .Object@networkDT_before_filter = data.table::setalloccol(.Object@networkDT_before_filter)
            }
            .Object
          })


# ** spatGridData ####
setClass('spatGridData',
         contains = 'VIRTUAL',
         slots = list(method = 'character',
                      parameters = 'ANY',
                      gridDT = 'nullOrDatatable'),
         prototype = prototype(method = NA_character_,
                               parameters = NULL,
                               gridDT = NULL))


setMethod('initialize', 'spatGridData',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@gridDT)) {
              .Object@gridDT = data.table::setalloccol(.Object@gridDT)
            }
            .Object
          })


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
         contains = 'VIRTUAL',
         slots = list(provenance = 'ANY'),
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
         contains = c('provData', 'VIRTUAL'),
         slots = list(spat_unit = 'character'), # not allowed to be NULL
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
         contains = 'VIRTUAL',
         slots = list(feat_type = 'character'), # not allowed to be NULL
         prototype = prototype(feat_type = NA_character_))



# ** miscData Class ####
#' @title Basic class for additional miscellaneous information
#'
#' @description
#' Classes (such as dimObj) that can hold information from multiple types of methods
#' use the misc slot to hold additional information specific to each method.
#' Information may be stored within as S3 structures.
setClass('miscData',
         contains = 'VIRTUAL',
         slots = list(misc = 'ANY'),
         prototype = prototype(misc = NULL))








# SUBCLASSES ####

# ** spatFeatData ####
#' Superclass for classes that contain both spatial and feature data
setClass('spatFeatData',
         contains = c('spatData', 'featData', 'VIRTUAL'))









# CORE ####

## Giotto class ####


## * Check ###
# Giotto class


# #' @title Check Giotto Object
# #' @name check_giottoObj
# #' @description check function for S4 giotto object
# #' @param object giotto object to check
# #' @keywords internal
# #' @noRd
# check_giotto_obj = function(object) {
#
#
#   errors = character()
#   ch = box_chars()
#
#
#   ## Validate Nesting ##
#   ## ---------------- ##
#
#   slot_depths = giotto_slot_depths()
#   for(slot_i in seq(nrow(slot_depths))) {
#
#     slot_data = slot(object, slot_depths$slot[[slot_i]])
#     if(!is.null(slot_data)) {
#       if(depth(slot_data, method = 'min') != slot_depths$depth[[slot_i]]) {
#         msg = wrap_txt('Invalid nesting discovered for',
#                        slot_depths$slot[[slot_i]],
#                        'slot')
#         errors = c(errors, msg)
#       }
#     }
#   }




  # ## Match spatial units and feature types ##
  #
  # # Find existing spat_units in source data
  # uniqueSrcSU = list_expression(object)$spat_unit
  # uniqueSrcSU = unique(uniqueSrcSU, list_spatial_info_names(object))
  #
  # # Find existing feat_types in source data
  # uniqueSrcFT = list_expression(object)$feat_type
  # uniqueSrcFT = unique(uniqueSrcFT, list_feature_info_names(object))
  #
  # # Find existing spat_units and feat_types within other data slots
  # uniqueDatSU = lapply(gDataSlots, function(x) list_giotto_data(gobject = object, slot = x)$spat_unit)
  # uniqueDatFT = lapply(gDataSlots, function(x) list_giotto_data(gobject = object, slot = x)$feat_type)

  # # If any spat_units exist, ensure that associated objects have identical cell_IDs
  # if(length(uniqueSU) > 0) {
  #   for(SU in seq_along(uniqueSU)) {
  #
  #     all(lapply(gDataSlots, function(x) {
  #       list_giotto_data(gobject = object,
  #                        slot = x)
  #     }))
  #   }
  # }
#
#
#
#
#
#
#
#   if(length(errors) == 0) TRUE else errors
# }







#' @title Update giotto object
#' @name updateGiottoObject
#' @description Updates the giotto object for changes in structure for backwards
#' compatibility with earlier versions
#' @param gobject giotto object to update
#' @details
#' Supported updates:
#' \itemize{
#'   \item{3.2.0 update adding multiomics slot}
#'   \item{master branch to suite - TODO}
#' }
#' @examples
#' \dontrun{
#' gobject = updateGiottoObject(gobject)
#' }
#' @export
updateGiottoObject = function(gobject) {

  if(!inherits(gobject, 'giotto')) {
    stop(wrap_txt('This function is intended for updating giotto objects'))
  }

  # 3.2.0 release adds multiomics slot
  if(is.null(attr(gobject, 'multiomics'))) {
    attr(gobject, 'multiomics') = NA
    gobject@multiomics = NULL
  }

  return(gobject)
}



##### * Definition ####
# Giotto class
# ! Any slot modifications should also be reflected in packedGiotto class !

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
#' @slot multiomics multiomics integration results
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
    join_info = "ANY",
    multiomics = "ANY"

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
    parameters = list(),
    instructions = NULL,
    offset_file = NULL,
    OS_platform = NULL,
    join_info = NULL,
    multiomics = NULL
  )

  # validity = check_giotto_obj
)








##### * Initialize ####
#' @noRd
#' @keywords internal
setMethod('initialize', signature('giotto'), function(.Object, ...) {

  .Object = callNextMethod()

  # a = list(...)


# TODO
  ## set slots ##
  ## --------- ##

  # if('spatial_info' %in% names(a)) {
  #   .Object = setPolygonInfo(.Object, gpolygon = a$spatial_info)
  # }
  # if('expression' %in% names(a)) {
  #   .Object = setExpression(.Object, values = a$expression)
  # }




  # message('initialize Giotto run\n\n')  # debug


  ## set instructions ##
  ## ---------------- ##

  # set default instructions (no recursive initialize)
  if(is.null(instructions(.Object))) {
    instructions(.Object, initialize = FALSE) = createGiottoInstructions()
  }

  ## test python module availability ##
  python_modules = c('pandas', 'igraph', 'leidenalg', 'community', 'networkx', 'sklearn')
  my_python_path = instructions(.Object, 'python_path')
  for(module in python_modules) {
    if(reticulate::py_module_available(module) == FALSE) {
      warning('module: ', module, ' was not found with python path: ', my_python_path, '\n')
    }
  }



  ## Slot Detection ##
  ## -------------- ##

  # detect expression and subcellular data
  avail_expr = list_expression(.Object)
  avail_si = list_spatial_info(.Object)
  avail_fi = list_feature_info(.Object)

  used_spat_units = unique(c(avail_expr$spat_unit, avail_si$spat_info))
  used_feat_types = unique(c(avail_expr$feat_type, avail_fi$feat_info))

  # detect ID slots
  avail_cid = list_cell_id_names(.Object)
  avail_fid = list_cell_id_names(.Object)

  # detect metadata slots
  avail_cm = list_cell_metadata(.Object)
  avail_fm = list_feat_metadata(.Object)

  # detect spatial location slot
  avail_sl = list_spatial_locations(.Object)

  # detect nearest network slot
  avail_nn = list_nearest_networks(.Object)

  # detect dimension reduction slot
  avail_dr = list_dim_reductions(.Object)

  # detect spatial network slot
  avail_sn = list_spatial_networks(.Object)

  # detect spatial enrichment slot
  avail_se = list_spatial_enrichments(.Object)


  ## Perform any subobject updates ##
  ## ----------------------------- ##

  # Feature Info #
  if(!is.null(avail_fi)) {
    info_list = get_feature_info_list(.Object)
    # update S4 object if needed
    info_list = lapply(info_list, function(info) {
      try_val = try(validObject(info), silent = TRUE)
      if(inherits(try_val, 'try-error')) {
        info = updateGiottoPointsObject(info)
      }
      return(info)
    })
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    .Object = setFeatureInfo(gobject = .Object,
                             x = info_list,
                             verbose = FALSE,
                             initialize = FALSE)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  }


  # Spatial Info #
  if(!is.null(avail_si)) {
    info_list = get_polygon_info_list(.Object)

    # update S4 object if needed
    info_list = lapply(info_list, function(info) {
      try_val = try(validObject(info), silent = TRUE)
      if(inherits(try_val, 'try-error')) {
        info = updateGiottoPolygonObject(info)
      }
      return(info)
    })
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    .Object = setPolygonInfo(gobject = .Object,
                             x = info_list,
                             verbose = FALSE,
                             centroids_to_spatlocs = FALSE,
                             initialize = FALSE)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  }






  ## Set active/default spat_unit and feat_type ##
  ## ------------------------------------------ ##

  # detect if actives are set in giotto instructions
  active_su = try(instructions(.Object, 'active_spat_unit'), silent = TRUE)
  active_ft = try(instructions(.Object, 'active_feat_type'), silent = TRUE)

  # determine actives using defaults if data exists then set
  if(inherits(active_su, 'try-error')) {
    if(!is.null(avail_expr) | !is.null(avail_si)) {
      active_su = set_default_spat_unit(gobject = .Object)
      instructions(.Object, 'active_spat_unit', initialize = FALSE) = active_su
    }
  }
  if(inherits(active_ft, 'try-error')) {
    if(!is.null(avail_expr) | !is.null(avail_fi)) {
      active_ft = set_default_feat_type(gobject = .Object,
                                        spat_unit = active_su)
      instructions(.Object, 'active_feat_type', initialize = FALSE) = active_ft
    }
  }






  ## Set expression_feat ##
  ## ------------------- ##
  e_feat = used_feat_types
  if('rna' %in% e_feat) {
    rna_idx = which(e_feat == 'rna')
    e_feat = c(e_feat[rna_idx], e_feat[-rna_idx])
  }
  .Object@expression_feat = e_feat





  ## Ensure Consistent IDs ##
  ## --------------------- ##

  # cell IDs can be expected to be constant across a spatial unit

  # expression
  if(!is.null(avail_expr)) {
    unique_expr_sets = unique(avail_expr[, .(spat_unit, feat_type)])

    for(set_i in nrow(unique_expr_sets)) {
      exp_list = get_expression_values_list(
        gobject = .Object,
        spat_unit = unique_expr_sets$spat_unit[[set_i]],
        feat_type = unique_expr_sets$feat_type[[set_i]],
        output = 'exprObj',
        set_defaults = FALSE
      )

      exp_list_names = lapply(exp_list, spatIDs)
      list_match = sapply(exp_list_names, setequal, exp_list_names[[1L]])
      if(!all(list_match)) {
        print(list_match)
        warning(wrap_text(
          'spat_unit:', unique_expr_sets$spat_unit[[set_i]], '/',
          'feat_type:', unique_expr_sets$feat_type[[set_i]],
          '\nNot all expression matrices share the same cell_IDs'
        ))
      }
    }
  }




  # MIGHT BE CHANGED IN THE FUTURE
  # feat_IDs cannot be expected to be constant across spat units.







  ## ID initialization ##
  ## ----------------- ##

  # Must be after default spat_unit/feat_type are set.
  # feat_ID initialization depends on active spat_unit


  # Initialization of cell_ID and feat_ID slots. These slots hold their     #
  # respective IDs for each spatial unit and feature type respectively.     #
  #                                                                         #
  # cell_metadata and feat_metadata slots are initialized off these slots.  #
  #                                                                         #
  # expression information is PREFERRED for ID initialization.              #
  # subcellular information, being raw data may also be used.               #

  .Object = init_cell_and_feat_IDs(gobject = .Object)





  ## Metadata initialization ##
  ## ----------------------- ##

  # Initialization of all spat_unit/feat_type combinations if the metadata  #
  # does not currently exist.                                               #

  # provenance is always updated from matched expression info if existing



  for(spatial_unit in used_spat_units) {
    for(feature_type in used_feat_types) {

      provenance = NULL
      # get expression for provenance info
      if(!is.null(avail_expr)) {
        if(nrow(avail_expr[spat_unit == spatial_unit &
                           feat_type == feature_type]) != 0L) {
          provenance = prov(get_expression_values(
            gobject = .Object,
            spat_unit = spatial_unit,
            feat_type = feature_type,
            output = 'exprObj',
            set_defaults = FALSE
          ))
        }
      }

      # initialize if no metadata exists OR none for this spat/feat

      # cell metadata
      if(is.null(avail_cm)) {
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        .Object = set_cell_metadata(
            gobject = .Object,
            metadata = 'initialize',
            spat_unit = spatial_unit,
            feat_type = feature_type,
            verbose = FALSE,
            set_defaults = FALSE
        )
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      } else if(nrow(avail_cm[spat_unit == spatial_unit &
                              feat_type == feature_type]) == 0L) {
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        .Object = set_cell_metadata(
          gobject = .Object,
          metadata = 'initialize',
          spat_unit = spatial_unit,
          feat_type = feature_type,
          verbose = FALSE,
          set_defaults = FALSE
        )
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      }

      # feature metadata
      if(is.null(avail_fm)) {
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        .Object = set_feature_metadata(
          gobject = .Object,
          metadata = 'initialize',
          spat_unit = spatial_unit,
          feat_type = feature_type,
          verbose = FALSE,
          set_defaults = FALSE
        )
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      } else if(nrow(avail_fm[spat_unit == spatial_unit &
                              feat_type == feature_type]) == 0L) {
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        .Object = set_feature_metadata(
          gobject = .Object,
          metadata = 'initialize',
          spat_unit = spatial_unit,
          feat_type = feature_type,
          verbose = FALSE,
          set_defaults = FALSE
        )
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      }


      # update provenance (always happens for all metadata objects)
      if(is.null(provenance)) next() # skip if no provenance info

      cm = get_cell_metadata(gobject = .Object,
                             spat_unit = spatial_unit,
                             feat_type = feature_type,
                             output = 'cellMetaObj',
                             copy_obj = FALSE,
                             set_defaults = FALSE)
      fm = get_feature_metadata(gobject = .Object,
                                spat_unit = spatial_unit,
                                feat_type = feature_type,
                                output = 'featMetaObj',
                                copy_obj = FALSE,
                                set_defaults = FALSE)
      prov(cm) = provenance
      prov(fm) = provenance
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      .Object = set_cell_metadata(gobject = .Object, metadata = cm, verbose = FALSE)
      .Object = set_feature_metadata(gobject = .Object, metadata = fm, verbose = FALSE)
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    }
  }



  ## Metadata ##
  ## ------------- ##

  if(!is.null(avail_cm)) {
    check_cell_metadata(gobject = .Object) # modifies by reference
  }

  if(!is.null(avail_fm)) {
    check_feat_metadata(gobject = .Object) # modifies by reference
  }


  ## Spatial locations ##
  ## ----------------- ##

  if(!is.null(avail_expr) & !is.null(avail_sl)) {
    # 1. ensure spatial locations and expression matrices have the same cell IDs
    # 2. give cell IDs if not provided
    check_spatial_location_data(gobject = .Object) # modifies by reference
  }



  ## Spatial network ##
  ## --------------- ##

  if(!is.null(avail_sl) & !is.null(avail_sn)) {
    # 1. ensure vertices have same IDs as seen in spat_unit for gobject
    # 2. ensure spatial locations of same spat_unit exists
    check_spatial_networks(gobject = .Object)
  }



  ## Spatial enrichment ##
  ## ------------------ ##

  if(!is.null(avail_sl) & !is.null(avail_se)) {
    # 1. ensure IDs in enrichment match gobject for same spat_unit
    # 2. ensure spatial locations exist for same spat_unit
    check_spatial_enrichment(gobject = .Object)
  }



  ## Nearest networks ##
  ## ---------------- ##

  if(!is.null(avail_expr) & !is.null(avail_nn)) {
    check_nearest_networks(gobject = .Object)
  }



  ## Dimension reduction ##
  ## ------------------- ##

  if(!is.null(avail_dr)) {
    .Object = check_dimension_reduction(gobject = .Object)
  }



  ## Spatial info ##
  ## ------------ ##

  if(!is.null(avail_si) & !is.null(avail_sl)) {
    check_spatial_info(gobject = .Object)
  }





  ## validity check ##
  ## -------------- ##
  obj_valid = try(validObject(.Object), silent = TRUE)
  if(inherits(obj_valid, 'try-error')) {
    .Object = updateGiottoObject(.Object)
    validObject(.Object)
  }



  .Object

})










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
    spat_unit = feat_type = prints = name = img_type = name = NULL

    cat("An object of class",  class(object), "\n")


    # active spat_unit and feat_type
    active_su = try(instructions(object, 'active_spat_unit'), silent = TRUE)
    active_ft = try(instructions(object, 'active_feat_type'), silent = TRUE)
    if(!inherits(active_su, 'try-error')) {
      cat('>Active spat_unit: ', active_su, '\n')
    }
    if(!inherits(active_ft, 'try-error')) {
      cat('>Active feat_type: ', active_ft, '\n')
    }


    cat('[SUBCELLULAR INFO]\n')
    if(!is.null(object@spatial_info)) cat('polygons      :', wrap_txt(list_spatial_info_names(object)), '\n')
    if(!is.null(object@feat_info)) cat('features      :', wrap_txt(list_feature_info_names(object)), '\n')


    mini_avail_print = function(avail_dt) {
      if(!'spat_unit' %in% colnames(avail_dt)) {
        avail_dt[, spat_unit := '']
      } else avail_dt[, spat_unit := paste0('[', spat_unit, ']')]
      if(!'feat_type' %in% colnames(avail_dt)) {
        avail_dt[, feat_type := '']
      } else avail_dt[, feat_type := paste0('[', feat_type, ']')]
      avail_dt[, prints := paste0(spat_unit, feat_type)]

      unique_entry = avail_dt[, unique(prints)]
      for(entry in unique_entry) {
        cat('  ', entry, paste0(' ', wrap_txt(avail_dt[prints == entry, name])), '\n', sep = '')
      }
    }


    cat('[AGGREGATE INFO]\n')
    avail_expr = list_expression(object)
    if(!is.null(avail_expr)) {
      cat('expression -----------------------\n')
      mini_avail_print(avail_expr)
    }

    avail_sl = list_spatial_locations(object)
    if(!is.null(avail_sl)) {
      cat('spatial locations ----------------\n')
      mini_avail_print(avail_sl)
    }

    avail_sn = list_spatial_networks(object)
    if(!is.null(avail_sn)) {
      cat('spatial networks -----------------\n')
      mini_avail_print(avail_sn)
    }

    avail_dim = list_dim_reductions(object)
    if(!is.null(avail_dim)) {
      cat('dim reduction --------------------\n')
      mini_avail_print(avail_dim)
    }

    avail_nn = list_nearest_networks(object)
    if(!is.null(avail_nn)) {
      cat('nearest neighbor networks --------\n')
      mini_avail_print(avail_nn)
    }

    avail_im = list_images(object)
    if(!is.null(avail_im)) {
      cat('attached images ------------------\n')
      if('image' %in% avail_im$img_type) {
        if(sum(avail_im$img_type == 'image') > 3) {
          cat(  'giottoLargeImage :', sum(avail_im$img_type == 'image'), 'items...\n')
        } else {
          cat('giottoImage      :', wrap_txt(avail_im[img_type == 'image', name]), '\n')
        }
      }
      if('largeImage' %in% avail_im$img_type) {
        if(sum(avail_im$img_type == 'largeImage') > 3) {
          cat(  'giottoLargeImage :', sum(avail_im$img_type == 'largeImage'), 'items...\n')
        } else {
          cat(  'giottoLargeImage :', wrap_txt(avail_im[img_type == 'largeImage', name]), '\n')
        }
      }
    }


    cat(wrap_txt('\n\nUse objHistory() to see steps and params used\n\n'))
    invisible(x = NULL)
  }
)



# for use with wrap() generic
# not intended to be used until after unwrapped to giotto class
# does not inherit giotto to avoid any method inheritance
setClass(
  "packedGiotto",
  slots = c(
    packed_spatial_info = "ANY",
    packed_feat_info = "ANY",

    expression = "nullOrList",
    expression_feat = "ANY",
    spatial_locs = "ANY",
    cell_metadata = "ANY",
    feat_metadata = "ANY",
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
    join_info = "ANY",
    multiomics = "ANY"

  ),

  prototype = list(
    packed_spatial_info = NULL,
    packed_feat_info = NULL,

    expression = NULL,
    expression_feat = NULL,
    spatial_locs = NULL,
    cell_metadata = NULL,
    feat_metadata = NULL,
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
    join_info = NULL,
    multiomics = NULL
  )
)

setMethod("show", signature(object='packedGiotto'),
          function(object) {
            print(paste("This is a", class(object), "object. Use 'Giotto::unwrap()' to unpack it"))
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

  if(length(errors) == 0) TRUE else errors
}



## * Definition ####
# exprObj Class

#' @title S4 exprObj
#' @description Framework to store aggregated expression information
#' @slot name name of exprObj
#' @slot exprMat matrix of expression information
#' @slot spat_unit spatial unit of expression (e.g. 'cell')
#' @slot feat_type feature type of expression (e.g. 'rna', 'protein')
#' @slot provenance origin data of expression information (if applicable)
#' @slot misc misc
#' @export
setClass('exprObj',
         contains = c('nameData', 'exprData', 'spatFeatData', 'miscData'),
         validity = check_expr_obj)


## * Initialize
# setMethod('initialize', 'exprObj',
#           function(.Object, ...) {
#
#             # expand args
#             a = list(.Object = .Object, ...)
#
#             # evaluate data
#             if('exprMat' %in% names(a)) {
#               exprMat = a$exprMat
#               if(is.null(exprMat)) exprMat = matrix()
#               else {
#                 # Convert matrix input to preferred format
#                 exprMat = evaluate_expr_matrix(exprMat)
#               }
#
#               # return to arg list
#               a$exprMat = exprMat
#             }
#
#             .Object = do.call('methods'::'callNextMethod', a)
#             .Object
#           })


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

    cat("An object of class",  class(object), "\n")

    # print spat/feat and provenance info
    cat(paste0('for spatial unit: "', slot(object, 'spat_unit'), '" and feature type: "', slot(object, 'feat_type'),'" \n'))
    if(!is.null(slot(object, 'provenance'))) cat('  Provenance: ', unlist(slot(object, 'provenance')),'\n')

    cat('\ncontains:\n')
    # preview matrix

    # * Matrix sparseMatrix specific *
    if(inherits(slot(object, 'exprMat'), 'sparseMatrix')) {
      print_cap = capture.output(Matrix::printSpMatrix2(x = slot(object, 'exprMat'),
                                                        zero.print = ".",
                                                        col.names = FALSE,
                                                        note.dropping.colnames = FALSE,
                                                        suppRows = TRUE,
                                                        suppCols = TRUE,
                                                        width = 40,
                                                        maxp = 80))
      print_cap = print_cap[-which(print_cap == ' ..............................')]
      writeLines(gsub(pattern = "in show(.*?))'", replacement = '', x = print_cap))
      cat('\n First four colnames:')
      cat('\n', wrap_txt(head(colnames(slot(object, 'exprMat')), 4), strWidth = 40), '\n')
    } else if(inherits(slot(object, 'exprMat'), 'denseMatrix')) {
      abb_mat(object, nrows = 10, ncols = 10, header = FALSE)
    } else {
      # * other matrices *
      print(slot(object, 'exprMat'))
      cat('\n')
    }

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
      errors = c(errors, msg)
    }

    if(colnames(object@metaDT)[[1]] != 'cell_ID') {
      msg = '"cell_ID" column should be the first column.'
      errors = c(errors, msg)
    }

  }
  if(length(errors) == 0) TRUE else errors
}

# * Definition ####
#' @title S4 cellMetaObj
#' @description Framework to store cell metadata
#' @slot metaDT metadata info
#' @slot col_desc (optional) character vector describing columns of the metadata
#' @slot spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @slot feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @slot provenance origin data of aggregated expression information (if applicable)
#' @export
setClass('cellMetaObj',
         contains = c('metaData', 'spatFeatData'),
         validity = check_cell_meta_obj)


setMethod('show', signature('cellMetaObj'), function(object) {

  cat('An object of class', class(object), '\n')
  cat('Provenance:', slot(object, 'provenance'), '\n')
  # TODO print spat/feat
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
      errors = c(errors, msg)
    }

    if(colnames(object@metaDT)[[1]] != 'feat_ID') {
      msg = '"feat_ID" column should be the first column.'
      errors = c(errors, msg)
    }

  }
  if(length(errors) == 0) TRUE else errors
}

# * Definition ####
#' @title S4 featMetaObj
#' @description Framework to store feature metadata
#' @slot metaDT metadata info
#' @slot col_desc (optional) character vector describing columns of the metadata
#' @slot spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @slot feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @slot provenance origin data of aggregated expression information (if applicable)
#' @export
setClass('featMetaObj',
         contains = c('metaData', 'spatFeatData'),
         validity = check_feat_meta_obj)


setMethod('show', signature('featMetaObj'), function(object) {

  cat('An object of class', class(object), '\n')
  cat('Provenance:', slot(object, 'provenance'), '\n')
  # TODO print spat/feat
  if(!is.null(slot(object, 'metaDT'))) print(slot(object, 'metaDT'))

})


# DIMENSION REDUCTION ####

## dimObj Class ####



##### * Check #####
# dimObj Class

#' @title Check dimObj
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

  # This check applied using check_dimension_reduction()
  # if(!inherits(rownames(object@coordinates, 'character'))) {
  #   msg = 'Dim reduction coordinate rownames must be character'
  #   errors = c(errors, msg)
  # }

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
#' @slot reduction whether reduction was performed on 'feats' or 'cells'
#' @slot reduction_method method used to generate dimension reduction
#' @slot coordinates embedding coordinates
#' @slot misc method-specific additional outputs
#' @export
setClass('dimObj',
         contains = c('nameData', 'spatFeatData'),
         slots = c(reduction = 'character',
                   reduction_method = 'character',
                   coordinates = 'ANY',
                   misc = 'ANY'),
         prototype = list(reduction = NA_character_,
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

## * Definition ####
# nnNetObj Class

#' @title S4 nnNetObj
#' @description Framework to store nearest neighbor network information
#' @slot name name of nnNetObj
#' @slot nn_type type of nearest neighbor network
#' @slot igraph igraph object containing network information
#' @slot feat_type feature type of data
#' @slot spat_unit spatial unit of data
#' @slot provenance origin of aggregated information (if applicable)
#' @slot misc misc
#' @export
setClass('nnNetObj',
         contains = c('nameData', 'nnData', 'spatFeatData', 'miscData'))





## * Initialize
# setMethod('initialize', 'nnNetObj',
#           function(.Object, ...) {
#
#             # expand args
#             a = list(.Object = .Object, ...)
#
#             # evaluate data
#             if('igraph' %in% names(a)) {
#               igraph = a$igraph
#               if(is.null(igraph)) igraph = NULL
#               else {
#                 # Convert igraph input to preferred format
#                 igraph = evaluate_nearest_networks(igraph)
#               }
#
#               # return to arg list
#               a$igraph = igraph
#             }
#
#             .Object = do.call('methods'::'callNextMethod', a)
#             .Object
#           })





## * Show ####
#' show method for nnNetObj class
#' @param object nearest neigbor network object
#' @aliases show,nnNetObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('nnNetObj'), function(object) {

    cat("An object of class",  class(object), ':', object@name, '\n')
    if(!is.null(object@nn_type)) cat('--| Contains nearest neighbor network generated with:', object@nn_type, '\n')
    if(!is.null(object@feat_type) & !is.null(object@spat_unit)) {
      cat('----| for feat_type:', object@feat_type, '\n')
      cat('----|     spat_unit:', object@spat_unit, '\n')
    }
    if(!is.null(object@provenance)) cat('----|     provenance:', object@provenance, '\n\n')

    if(!is.null(object@igraph)) {
      print(object@igraph)
      cat('\n\n')
    }

    if(!is.null(object@misc)) {
      cat('Additional included info:\n')
      print(names(object@misc))
      cat('\n')
    }

  })









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

  # Allow check_spatial_location_data() to compensate for missing cell_ID
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
         contains = c('nameData', 'coordDataDT', 'spatData', 'miscData'),
         validity = check_spat_locs_obj)




## * Initialize
# setMethod('initialize', 'spatLocsObj',
#           function(.Object, ...) {
#
#             # expand args
#             a = list(.Object = .Object, ...)
#
#             # evaluate data
#             if('coordinates' %in% names(a)) {
#               coordinates = a$coordinates
#               if(is.null(coordinates)) {
#                 coordinates = data.table::data.table(
#                   sdimx = NA_real_,
#                   sdimy = NA_real_,
#                   cell_ID = NA_character_
#                 )
#               } else {
#                 coordinates = evaluate_spatial_locations(coordinates)
#               }
#
#               # return to arg list
#               a$coordinates = coordinates
#             }
#
#             .Object = do.call('methods'::'callNextMethod', a)
#             .Object
#           })




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
    cat(paste0('for spatial unit: "', slot(object, 'spat_unit'), '"\n'))
    if(!is.null(slot(object, 'provenance'))) cat('provenance:', slot(object, 'provenance'), '\n')

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
         contains = c('nameData', 'spatNetData' ,'spatData', 'miscData'),
         slots = c(crossSectionObjects = 'ANY'),
         prototype = list(crossSectionObjects = NULL),
         validity = check_spat_net_obj)


### * show ####
# spatialNetworkObj Class

#' show method for spatialNetworkObj class
#' @param object spatial network object
#' @aliases show,spatialNetworkObj-method
#' @docType methods
#' @importFrom methods show
#' @importFrom graphics segments

#' @rdname show-methods
setMethod(
  f = "show", signature('spatialNetworkObj'), function(object) {

    cat("An object of class",  class(object), "\n")
    if(!is.na(object@method)) cat('Contains spatial network generated with:', object@method, '\n')
    if(!is.na(object@spat_unit)) cat('for spatial unit: "', object@spat_unit, '"\n', sep = '')
    if(!is.null(object@provenance)) cat('provenance:', object@provenance, '\n')

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
         contains = c('nameData', 'spatGridData', 'spatFeatData', 'miscData'),
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
    cat('Contains annotations for spatial unit: "', slot(object, 'spat_unit'),'"', sep = '')
    if(!is.na(slot(object, 'feat_type'))) {
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



## spatEnrObj class ####

# * definition ####
# spatEnrObj class

#' @title S4 spatEnrObj Class
#' @description Framework to store spatial enrichment results
#' @slot name name of enrichment object
#' @slot method method used to perform spatial enrichment
#' @slot enrichDT spatial enrichment data.table
#' @slot spat_unit spatial unit
#' @slot feat_type feature type
#' @slot provenance provenance information
#' @slot misc misc
#' @export
setClass('spatEnrObj',
         contains = c('nameData', 'enrData', 'spatFeatData', 'miscData'))


# * show ####
# spatEnrObj Class

#' show method for spatEnrObj class
#' @param object spatial locations object
#' @aliases show,spatEnrObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('spatEnrObj'), function(object) {

    cat("An object of class",  class(object), "\n")
    cat(paste0('for spatial unit: "', slot(object, 'spat_unit'), '" and feature type: "',
               slot(object, 'feat_type'), '"\n'))
    if(!is.null(slot(object, 'provenance'))) cat('provenance:', slot(object, 'provenance'), '\n')

    cat('   ------------------------\n\npreview:\n')
    if(!is.null(slot(object, 'enrichDT'))) {
      enr_cols = ncol(slot(object, 'enrichDT'))
      if(enr_cols > 10L) {
        show(slot(object, 'enrichDT')[1:4, 1:10])
        cat(rep(' ', times = getOption('width')/2.5 - 10L), rep('.', 20L),'\n', sep = '')
        show(slot(object, 'enrichDT')[1:4, 'cell_ID'])
        cat('...', enr_cols - 11L, ' cols omitted\n', sep = '')
      } else {
        show(slot(object, 'enrichDT')[1:4])
      }
    }

    cat('\n...first 20 remaining colnames:\n')
    cat('\n', wrap_txt(head(colnames(slot(object, 'enrichDT'))[-c(1:10)], 20L), strWidth = 40L), '\n')

    cat('\n\n')

  })


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
#' @slot unique_ID_cache cached unique spatial IDs that should match the spatVector slot
#' @details holds polygon data
#'
#' @export
giottoPolygon = setClass(
  Class = "giottoPolygon",
  contains = c('nameData'),

  slots = c(
    spatVector = "ANY",
    spatVectorCentroids = "ANY",
    overlaps = "ANY",
    unique_ID_cache = 'character'
  ),

  prototype = list(
    spatVector = NULL,
    spatVectorCentroids = NULL,
    overlaps = NULL,
    unique_ID_cache = NA_character_
  )
)




#' @title Update giotto polygon object
#' @name updateGiottoPolygonObject
#' @param gpoly giotto polygon object
#' @export
updateGiottoPolygonObject = function(gpoly) {
  if(!inherits(gpoly, 'giottoPolygon')) {
    stop('This function is only for giottoPoints')
  }

  # 3.2.X adds cacheing of IDs
  if(is.null(attr(gpoly, 'unique_ID_cache'))) {
    attr(gpoly, 'unique_ID_cache') = unique(as.list(gpoly@spatVector)$poly_ID)
  }

  gpoly
}


# * show ####
setMethod('show', signature = 'giottoPolygon', function(object) {

  cat('An object of class giottoPolygon with name "', object@name, '"\n', sep = '')
  cat('Spatial Information:\n')
  print(object@spatVector)

  if(!is.null(object@spatVectorCentroids)) {
    cat(' centroids   : calculated\n')
  } else {
    cat(' centroids   : NULL\n')
  }

  if(!is.null(object@overlaps)) {
    cat(' overlaps    : calculated')
  } else {
    cat(' overlaps    : NULL')
  }

})



# for use with wrap() generic
setClass('packedGiottoPolygon',
         contains = c('nameData'),

         slots = c(
           packed_spatVector = 'ANY',
           packed_spatVectorCentroids = 'ANY',
           packed_overlaps = 'ANY',
           unique_ID_cache = 'character'
         ),
         prototype = list(
           packed_spatVector = NULL,
           packed_spatVectorCentroids = NULL,
           packed_overlaps = NULL,
           unique_ID_cache = NA_character_
         ))


setMethod("show", signature(object='packedGiottoPolygon'),
          function(object) {
            print(paste("This is a", class(object), "object. Use 'Giotto::unwrap()' to unpack it"))
          }
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
#' @slot unique_ID_cache cached unique feature IDs that should match the spatVector slot
#' @details Contains vector-type feature data
#'
#' @export
giottoPoints <- setClass(
  Class = "giottoPoints",
  contains = c('featData'),

  slots = c(
    spatVector = "ANY",
    networks = "ANY",
    unique_ID_cache = 'character'
  ),

  prototype = list(
    spatVector = NULL,
    networks = NULL,
    unique_ID_cache = NA_character_
  )
)




#' @title Update giotto points object
#' @name updateGiottoPointsObject
#' @param gpoints giotto points object
#' @export
updateGiottoPointsObject = function(gpoints) {
  if(!inherits(gpoints, 'giottoPoints')) {
    stop('This function is only for giottoPoints')
  }

  # 3.2.X adds cacheing of IDs
  if(is.null(attr(gpoints, 'unique_ID_cache'))) {
    attr(gpoints, 'unique_ID_cache') = unique(as.list(gpoints@spatVector)$feat_ID)
  }

  gpoints
}




# * show ####
setMethod('show', signature = 'giottoPoints', function(object) {

  cat('An object of class giottoPoints with feature type "', object@feat_type, '"\n', sep = '')
  cat('Feature Information:\n')
  print(object@spatVector)

  if(!is.null(object@networks)) {
    cat(' feat. net.  :')
    print(object@networks)
  }

})






# for use with wrap() generic
setClass(
  'packedGiottoPoints',

  slots = c(
    feat_type = 'character',
    packed_spatVector = 'ANY',
    networks = 'ANY',
    unique_ID_cache = 'character'
  ),
  prototype = list(
    feat_type = NA_character_,
    packed_spatVector = NULL,
    networks = NULL,
    unique_ID_cache = NA_character_
  )
)





setMethod("show", signature(object='packedGiottoPoints'),
          function(object) {
            print(paste("This is a", class(object), "object. Use 'Giotto::unwrap()' to unpack it"))
          }
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
  contains = 'nameData',

  slots = c(
    network_datatable = "ANY",
    network_lookup_id = "ANY",
    full = "ANY"
  ),

  prototype = list(
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



# constructor functions for S4 subobjects ####

#' @title Create S4 exprObj
#' @name createExprObj
#' @description Create an S4 exprObj
#' @param expression_data expression data
#' @param name name of exprObj
#' @param spat_unit spatial unit of expression (e.g. 'cell')
#' @param feat_type feature type of expression (e.g. 'rna', 'protein')
#' @param provenance origin data of expression information (if applicable)
#' @param misc misc
createExprObj = function(expression_data,
                         name = 'test',
                         spat_unit = 'cell',
                         feat_type = 'rna',
                         provenance = NULL,
                         misc = NULL) {

  exprMat = evaluate_expr_matrix(expression_data)

  create_expr_obj(name = name,
                  exprMat = exprMat,
                  spat_unit = spat_unit,
                  feat_type = feat_type,
                  provenance = provenance,
                  misc = misc)
}


#' @param exprMat matrix of expression information
#' @keywords internal
#' @noRd
create_expr_obj = function(name = 'test',
                           exprMat = NULL,
                           spat_unit = 'cell',
                           feat_type = 'rna',
                           provenance = NULL,
                           misc = NULL) {

  if(is.null(exprMat)) exprMat = matrix()

  return(new('exprObj',
             name = name,
             exprMat = exprMat,
             spat_unit = spat_unit,
             feat_type = feat_type,
             provenance = provenance,
             misc = misc))
}





#' @title Create S4 cellMetaObj
#' @name createCellMetaObj
#' @description Create an S4 cellMetaObj
#' @param metadata metadata info
#' @param col_desc (optional) character vector describing columns of the metadata
#' @param spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @param feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @param provenance origin data of aggregated expression information (if applicable)
#' @param verbose be verbose
#' @export
createCellMetaObj = function(metadata,
                             spat_unit = 'cell',
                             feat_type = 'rna',
                             provenance = NULL,
                             col_desc = NULL,
                             verbose = TRUE) {

  metadata = evaluate_cell_metadata(metadata = metadata,
                                    verbose = verbose)

  create_cell_meta_obj(metaDT = metadata,
                       col_desc = col_desc,
                       spat_unit = spat_unit,
                       feat_type = feat_type,
                       provenance = provenance)
}


#' @keywords internal
#' @noRd
create_cell_meta_obj = function(metaDT = NULL,
                                col_desc = NA_character_,
                                spat_unit = 'cell',
                                feat_type = 'rna',
                                provenance = NULL) {

  if(is.null(col_desc)) col_desc = NA_character_

  if(is.null(metaDT)) metaDT = data.table::data.table(cell_ID = NA_character_)

  return(new('cellMetaObj',
             metaDT = metaDT,
             col_desc = col_desc,
             spat_unit = spat_unit,
             provenance = provenance,
             feat_type = feat_type))
}





#' @title Create S4 featMetaObj
#' @name createFeatMetaObj
#' @description Create an S4 featMetaObj
#' @param metadata metadata info
#' @param col_desc (optional) character vector describing columns of the metadata
#' @param spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @param feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @param provenance origin data of aggregated expression information (if applicable)
#' @param verbose be verbose
#' @export
createFeatMetaObj = function(metadata,
                             spat_unit = 'cell',
                             feat_type = 'rna',
                             provenance = NULL,
                             col_desc = NULL,
                             verbose = TRUE) {

  metadata = evaluate_feat_metadata(metadata = metadata,
                                    verbose = verbose)

  create_feat_meta_obj(metaDT = metadata,
                       col_desc = col_desc,
                       spat_unit = spat_unit,
                       feat_type = feat_type,
                       provenance = provenance)
}



#' @keywords internal
#' @noRd
create_feat_meta_obj = function(metaDT = NULL,
                                col_desc = NA_character_,
                                spat_unit = 'cell',
                                feat_type = 'rna',
                                provenance = NULL) {

  if(is.null(col_desc)) col_desc = NA_character_

  if(is.null(metaDT)) metaDT = data.table::data.table(feat_ID = NA_character_)

  return(new('featMetaObj',
             metaDT = metaDT,
             col_desc = col_desc,
             spat_unit = spat_unit,
             provenance = provenance,
             feat_type = feat_type))
}







#' @title Create S4 dimObj
#' @name createDimObj
#' @description Create an S4 dimObj
#' @param coordinates embedding coordinates
#' @param name name of dimObj
#' @param reduction reduction on columns (e.g. cells) or rows (e.g. features)
#' @param reduction_method method used to generate dimension reduction
#' @param spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @param feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @param provenance origin data of aggregated expression information (if applicable)
#' @param misc misc
#' @param my_rownames (optional) if needed, set coordinates rowname values here
#' @export
createDimObj = function(coordinates,
                        name = 'test',
                        spat_unit = 'cell',
                        feat_type = 'rna',
                        method = NULL,
                        reduction = 'cells',
                        provenance = NULL,
                        misc = NULL,
                        my_rownames = NULL) {

  coordinates = evaluate_dimension_reduction(coordinates)

  create_dim_obj(name = name,
                 reduction = reduction,
                 reduction_method = method,
                 coordinates = coordinates,
                 spat_unit = spat_unit,
                 feat_type = feat_type,
                 provenance = provenance,
                 misc = misc,
                 my_rownames = my_rownames)
}


#' @keywords internal
#' @noRd
create_dim_obj = function(name = 'test',
                          reduction = 'cells',
                          reduction_method = NA_character_,
                          coordinates = NULL,
                          spat_unit = 'cell',
                          feat_type = 'rna',
                          provenance = NULL,
                          misc = NULL,
                          my_rownames = NULL) {

  if(is.null(reduction_method)) reduction_method = NA_character_

  number_of_dimensions = ncol(coordinates)
  colnames(coordinates) = paste0('Dim.', seq(number_of_dimensions))

  if(!is.null(my_rownames)) {
    rownames(coordinates) = as.character(my_rownames)
  }

  new('dimObj',
      name = name,
      reduction = reduction,
      reduction_method = reduction_method,
      coordinates = coordinates,
      spat_unit = spat_unit,
      feat_type = feat_type,
      provenance = if(is.null(provenance)) spat_unit else provenance, # assumed
      misc = misc)
}








#' @title Create S4 nnNetObj
#' @name createNearestNetObj
#' @description Create an S4 nnNetObj
#' @param name name of nnNetObj
#' @param nn_type type of nearest neighbor network
#' @param network igraph object or data.frame containing nearest neighbor
#' information (see details)
#' @slot spat_unit spatial unit of data
#' @slot feat_type feature type of data
#' @slot provenance origin of aggregated information (if applicable)
#' @param misc misc
#' @details igraph and dataframe-like inputs must include certain information.
#' For igraph, it must have, at minimum vertex 'name' attributes and 'distance'
#' edge attribute.
#' dataframe-like inputs must have 'from', 'to', and 'distance' columns
#' @export
createNearestNetObj = function(name = 'test',
                               network,
                               nn_type = NULL,
                               spat_unit = 'cell',
                               feat_type = 'rna',
                               provenance = NULL,
                               misc = NULL) {

  if(is.null(network)) igraph = NULL
  else {
    # convert igraph input to preferred format
    igraph = evaluate_nearest_networks(network)
  }

  create_nn_net_obj(name = name,
                    igraph = igraph,
                    nn_type = nn_type,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    provenance = provenance,
                    misc = misc)
}


#' @keywords internal
#' @noRd
create_nn_net_obj = function(name = 'test',
                             nn_type = NA_character_,
                             igraph = NULL,
                             spat_unit = 'cell',
                             feat_type = 'rna',
                             provenance = NULL,
                             misc = NULL) {

  if(is.null(nn_type)) nn_type = NA_character_

  new('nnNetObj',
      name = name,
      nn_type = nn_type,
      igraph = igraph,
      spat_unit = spat_unit,
      feat_type = feat_type,
      provenance = provenance,
      misc = misc)
}








#' @title Create S4 spatLocsObj
#' @name create_spat_locs_obj
#' @description Create an S4 spatLocsObj
#' @param coordinates spatial coordinates
#' @param name name of spatLocsObj
#' @param spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @param provenance origin data of aggregated expression information (if applicable)
#' @param misc misc
#' @export
createSpatLocsObj = function(coordinates,
                             name = 'test',
                             spat_unit = 'cell',
                             provenance = NULL,
                             misc = NULL,
                             verbose = TRUE) {

  # convert coordinates input to preferred format
  coordinates = evaluate_spatial_locations(spatial_locs = coordinates,
                                           verbose = verbose)

  create_spat_locs_obj(name = name,
                       coordinates = coordinates,
                       spat_unit = spat_unit,
                       provenance = provenance,
                       misc = misc)
}



#' @keywords internal
#' @noRd
create_spat_locs_obj = function(name = 'test',
                                coordinates = NULL,
                                spat_unit = 'cell',
                                provenance = NULL,
                                misc = NULL) {

  if(is.null(coordinates)) {
    coordinates = data.table::data.table(
      sdimx = NA_real_,
      sdimy = NA_real_,
      cell_ID = NA_character_
    )
  }

  # set cell_ID col if missing to conform to spatialLocationsObj validity
  # should already never be the case after evaluation
  if(!'cell_ID' %in% colnames(coordinates)) coordinates[, cell_ID := NA_character_]

  new('spatLocsObj',
      name = name,
      coordinates = coordinates,
      spat_unit = spat_unit,
      provenance = provenance,
      misc = misc)
}








#' @title Create S4 spatialNetworkObj
#' @name createSpatNetObj
#' @param network network data with connections, distances, and weightings
#' @param name name of spatialNetworkObj
#' @param networkDT_before_filter (optional) unfiltered data.table  of network connections, distances, and weightings
#' @param spat_unit spatial unit tag
#' @param method method used to generate spatial network
#' @param parameters (optional) additional method-specific parameters used during spatial network generation
#' @param outputObj (optional) network geometry object
#' @param cellShapeObj (optional) network cell shape information
#' @param crossSectionObjects (optional) crossSectionObjects (see \code{\link{create_crossSection_object}})
#' @param provenance (optional) origin of aggregated information (if applicable)
#' @param misc misc
#' @export
createSpatNetObj = function(network,
                            name = 'test',
                            networkDT_before_filter = NULL,
                            method = NULL,
                            spat_unit = 'cell',
                            provenance = NULL,
                            parameters = NULL,
                            outputObj = NULL,
                            cellShapeObj = NULL,
                            crossSectionObjects = NULL,
                            misc = NULL) {

  networkDT = evaluate_spatial_network(network)

  create_spat_net_obj(name = name,
                      method = method,
                      parameters = parameters,
                      outputObj = outputObj,
                      networkDT = networkDT,
                      networkDT_before_filter = networkDT_before_filter,
                      cellShapeObj = cellShapeObj,
                      crossSectionObjects = crossSectionObjects,
                      spat_unit = spat_unit,
                      provenance = provenance,
                      misc = misc)
}


#' @keywords internal
#' @noRd
create_spat_net_obj = function(name = 'test',
                               method = NA_character_,
                               parameters = NULL,
                               outputObj = NULL,
                               networkDT = NULL,
                               networkDT_before_filter = NULL,
                               cellShapeObj = NULL,
                               crossSectionObjects = NULL,
                               spat_unit = 'cell',
                               provenance = NULL,
                               misc = NULL ) {

  if(is.null(method)) method = NA_character_

  new('spatialNetworkObj',
      name = name,
      method = method,
      parameters = parameters,
      outputObj = outputObj,
      networkDT = networkDT,
      networkDT_before_filter = networkDT_before_filter,
      cellShapeObj = cellShapeObj,
      crossSectionObjects = crossSectionObjects,
      spat_unit = spat_unit,
      provenance = provenance,
      misc = misc)
}






#' @title Create S4 spatEnrObj
#' @name create_spat_enr_obj
#' @param enrichment_data spatial enrichment results, provided a dataframe-like object
#' @param name name of S4 spatEnrObj
#' @param method method used to generate spatial enrichment information
#' @param spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @param feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @param provenance origin data of aggregated expression information (if applicable)
#' @param misc misc additional information about the spatial enrichment or how it
#' was generated
#' @export
createSpatEnrObj = function(enrichment_data,
                            name = 'test',
                            spat_unit = 'cell',
                            feat_type = 'rna',
                            method = NULL,
                            provenance = NULL,
                            misc = NULL,
                            verbose = TRUE) {

  enrichDT = evaluate_spatial_enrichment(enrichment_data, verbose = verbose)

  create_spat_enr_obj(name = name,
                      method = method,
                      enrichDT = enrichment_data,
                      spat_unit = spat_unit,
                      feat_type = feat_type,
                      provenance = provenance,
                      misc = misc)
}


#' @keywords internal
#' @noRd
create_spat_enr_obj = function(name = 'test',
                               method = NA_character_,
                               enrichDT = NULL,
                               spat_unit = 'cell',
                               feat_type = 'rna',
                               provenance = NULL,
                               misc = NULL) {

  if(is.null(method)) method = NA_character_

  new('spatEnrObj',
      name = name,
      method = method,
      enrichDT = enrichDT,
      spat_unit = spat_unit,
      feat_type = feat_type,
      provenance = provenance,
      misc = misc)
}













#' @title Create S4 spatialGridObj
#' @name create_spat_grid_obj
#' @param name name of spatialGridObj
#' @param method method used to generate spatial grid
#' @param parameters additional method-specific parameters used during spatial grid generation
#' @param gridDT data.table holding the spatial grid information
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param provenance origin of aggregated information (if applicable)
#' @param misc misc
#' @keywords internal
create_spat_grid_obj = function(name = 'test',
                                method = NA_character_,
                                parameters = NULL,
                                gridDT = NULL,
                                spat_unit = 'cell',
                                feat_type = 'rna',
                                provenance = NULL,
                                misc = NULL) {

  return(new('spatGridObj',
             name = name,
             method = method,
             parameters = parameters,
             gridDT = gridDT,
             spat_unit = spat_unit,
             feat_type = feat_type,
             provenance = provenance,
             misc = misc))
}







#' @title Create feature network object
#' @name create_featureNetwork_object
#' @param name name to assign the created feature network object
#' @param network_datatable network data.table object
#' @param network_lookup_id network lookup id
#' @param full fully connected status
#' @keywords internal
create_featureNetwork_object = function(name = 'feat_network',
                                        network_datatable = NULL,
                                        network_lookup_id = NULL,
                                        full = NULL) {


  # create minimum giotto points object
  f_network = featureNetwork(name = name,
                             network_datatable = NULL,
                             network_lookup_id = NULL,
                             full = NULL)

  ## 1. check network data.table object
  if(!methods::is(network_datatable, 'data.table')) {
    stop("network_datatable needs to be a network data.table object")
  }
  f_network@network_datatable = network_datatable

  ## 2. provide network fully connected status
  f_network@full = full

  ## 3. provide feature network name
  f_network@name = name

  ## 4. network lookup id
  f_network@network_lookup_id = network_lookup_id

  # giotoPoints object
  return(f_network)

}


# create Giotto points from data.frame or spatVector

#' @title Create giotto points object
#' @name createGiottoPoints
#' @description Creates Giotto point object from a structured dataframe-like object
#' @param x spatVector or data.frame-like object with points coordinate information (x, y, feat_ID)
#' @param feat_type feature type
#' @param verbose be verbose
#' @param unique_IDs (optional) character vector of unique IDs present within
#' the spatVector data. Provided for cacheing purposes
#' @return giottoPoints
#' @concept polygon
#' @export
createGiottoPoints = function(x,
                              feat_type = 'rna',
                              verbose = TRUE,
                              unique_IDs = NULL) {

  if(inherits(x, 'data.frame')) {

    spatvec = create_spatvector_object_from_dfr(x = x,
                                                verbose = verbose)
    g_points = create_giotto_points_object(feat_type = feat_type,
                                           spatVector = spatvec,
                                           unique_IDs = unique_IDs)

  } else if(inherits(x, 'spatVector')) {

    g_points = create_giotto_points_object(feat_type = feat_type,
                                           spatVector = x,
                                           unique_IDs = unique_IDs)

  } else {

    stop('Class ', class(x), ' is not supported')

  }

  return(g_points)

}



#' @title Create giotto points object
#' @name create_giotto_points_object
#' @param feat_type feature type
#' @param spatVector terra spatVector object containing point data
#' @param networks (optional) feature network object
#' @param unique_IDs (optional) unique IDs in spatVector for cacheing
#' @keywords internal
create_giotto_points_object = function(feat_type = 'rna',
                                       spatVector = NULL,
                                       networks = NULL,
                                       unique_IDs = NULL) {

  if(is.null(feat_type)) feat_type = NA # compliance with featData class

  # create minimum giotto points object
  g_points = giottoPoints(feat_type = feat_type,
                          spatVector = NULL,
                          networks = NULL)

  ## 1. check terra spatVector object
  if(!inherits(spatVector, 'SpatVector')) {
    stop("spatVector needs to be a spatVector object from the terra package")
  }

  g_points@spatVector = spatVector

  ## 2. provide feature id
  g_points@feat_type = feat_type

  ## 3. feature_network object
  g_points@networks = networks

  ## 4. feat_ID cacheing
  if(is.null(unique_IDs)) {
    g_points@unique_ID_cache = featIDs(g_points)
  } else {
    g_points@unique_ID_cache = unique_IDs
  }

  # giottoPoints object
  return(g_points)

}








## extension of spatVector object
## name should match the cellular structure

#' @title Create a giotto polygon object
#' @name create_giotto_polygon_object
#' @param name name of polygon object
#' @param spatVector SpatVector of polygons
#' @param spatVectorCentroids (optional) SpatVector of polygon centroids
#' @param overlaps (optional) feature overlaps of polygons
#' @param unique_IDs unique polygon IDs for cacheing
#' @keywords internal
create_giotto_polygon_object = function(name = 'cell',
                                        spatVector = NULL,
                                        spatVectorCentroids = NULL,
                                        overlaps = NULL,
                                        unique_IDs = NULL) {


  # create minimum giotto
  g_polygon = giottoPolygon(name = name,
                            spatVector = NULL,
                            spatVectorCentroids = NULL,
                            overlaps = NULL)

  ## 1. check spatVector object
  if(!methods::is(spatVector, 'SpatVector')) {
    stop("spatVector needs to be a SpatVector object from the terra package")
  }

  g_polygon@spatVector = spatVector


  ## 2. centroids need to be of similar length as polygons
  if(!is.null(spatVectorCentroids)) {
    if(!methods::is(spatVectorCentroids, 'SpatVector')) {
      stop("spatVectorCentroids needs to be a spatVector object from the terra package")
    }

    l_centroids = nrow(terra::values(spatVectorCentroids))
    l_polygons = nrow(terra::values(spatVector))

    if(l_centroids == l_polygons) {
      g_polygon@spatVectorCentroids = spatVectorCentroids
    } else {
      stop('number of centroids does not equal number of polygons')
    }

  }

  ## 3. overlaps info
  g_polygon@overlaps = overlaps

  ## 4. spat_ID cacheing
  if(is.null(unique_IDs)) {
    g_polygon@unique_ID_cache = spatIDs(g_polygon)
  } else {
    g_polygon@unique_ID_cache = unique_IDs
  }

  # provide name
  g_polygon@name = name

  # giotto polygon object
  return(g_polygon)
}







# Possibly to be implemented ####
# icfObject


