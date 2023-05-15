

### ---------------------------------------------------------------- ###
# Get and set functions to get and set values the giotto class slots   #
#                                                                      #
# Accessor functions are paired by internal and externals.             #
# External functions call the internals to perform their work          #
#                                                                      #
# > externals responsibilities:                                        #
#     - provide default spat_unit/feat_type                            #
#                                                                      #
#     [getters]                                                        #
#        - provide subsetting by cell/feat ID                          #
#     [setters]                                                        #
#        - Parse input for nesting information (read fxns)             #
#        - Ensure input and gobject compatibility                      #
#                                                                      #
# > internals responsibilities:                                        #
#     - understand giotto nesting structure                            #
#     - provide default spat_unit/feat_type                            #
#        - (never when called by external)                             #
#                                                                      #
#     [getters]                                                        #
#        - perform copy of objects if needed                           #
#     [setters]                                                        #
#        - read S4 metadata for nesting info                           #
#        - provide 'initialize' options                                #
#        - provide method to REMOVE info using NULL                    #
#        - call initialize(giotto) if sensitive slot                   #
#                                                                      #
# > initialize generic responsibilities (see classes.R)                #
#     [setters]                                                        #
#        - Convert to correct type and formatting (eval fxns)          #
### ---------------------------------------------------------------- ###



#%%%%% NOTE: python and instructions accessors are currently in giotto.R %%%%%#




## common in internal functions ####


## Slot Depth Information ####
# Function to provide correct slot nesting depth definitions for easy testing
# Values provided are for when the slots are populated
giotto_slot_depths = function() {
  data.table::data.table(
    slot = c(
      'expression',
      'expression_feat',
      'spatial_locs',
      'spatial_info',
      'feat_info',
      'cell_metadata',
      'feat_metadata',
      'cell_ID',
      'feat_ID',
      'spatial_network',
      'spatial_enrichment',
      'dimension_reduction',
      'nn_network',
      'images',
      'largeImages'
      # 'spatial_grid',
      # 'parameters',
      # 'instructions',
      # 'offset_file',
      # 'OS_platform',
      # 'join_info'
      # 'multiomics'
    ),
    depth = c(
      3L, 0L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 3L, 5L, 4L, 1L, 1L
    )
  )
}











## Read S4 Nesting Tags ####
#' @noRd
#' @description SHOULD ONLY BE CALLED FROM ACCESSORS. RELIES ON SPECIFIC NAMES FOR
#' NESTING ELEMENTS IN PARENT FRAME. \cr
#' Reads the nesting information attached to giotto S4 subobjects and
#' compares against the already-existing default values in parent frame. Final
#' nesting values will be sent to the parent frame.
#' @param x a giotto S4 subobject
#' @param nest_elements named character vector of nesting elements for how the
#' information should be nested. (see details)
#' @param specified named logical vector for whether user specified input for a
#' nesting element
#' @details Nesting elements define the nesting structure within giotto slots.
#' Common examples are 'spat_unit', 'feat_type', and 'name' \cr
#' This function compares the nest_elements that are currently available vs those
#' that are suggested by the S4 subobject's appended information on the basis of
#' whether the existing nest_elements were directly specified by the user or if
#' they were simply default values \cr
#' If values were directly specified, then the external nest_elements values will
#' be used downstream and assigned into the relevant S4 subobject slots.
#' If the values were NOT specified then the subobject values will be used downstream.
#' Values will be directly pulled from and set to the parent frame, with the
#' exception of the S4 object itself.
#' @return modified S4 object
#' @keywords internal
read_s4_nesting = function(x) {

  p = parent.frame()

  s_names = methods::slotNames(x)


  # Determine nesting element to use. If parent frame variables are edited, it
  # will happen within the if statements.

  # if S4 objects will also be edited within if statements, but values will only
  # be sent back to parent frame at end of function

  if('spat_unit' %in% s_names) {
    if(isTRUE(p$nospec_unit)) {
      if(!is.na(spatUnit(x))) p$spat_unit = spatUnit(x)
    } else {
      spatUnit(x) = p$spat_unit
    }
  }

  if('feat_type' %in% s_names) {
    if(isTRUE(p$nospec_feat)) {
      if(!is.na(featType(x))) p$feat_type = featType(x)
    } else {
      featType(x) = p$feat_type
    }
  }

  if('name' %in% s_names) {
    if(isTRUE(p$nospec_name)) {
      if(!is.na(objName(x))) p$name = objName(x)
    } else {
      objName(x) = p$name
    }
  }

  if('provenance' %in% s_names) {
    if(is.null(p$provenance)) {
      if(!is.null(prov(x))) p$provenance = prov(x)
    } else {
      prov(x) = p$provenance
    }
  }

  if('reduction' %in% s_names) {
    if(isTRUE(p$nospec_red)) {
      if(!is.na(slot(x, 'reduction'))) p$reduction = slot(x, 'reduction')
    } else {
      slot(x, 'reduction') = p$reduction
    }
  }

  if('reduction_method' %in% s_names) {
    if(isTRUE(p$nospec_red_method)) {
      if(!is.na(slot(x, 'reduction_method'))) p$reduction_method = slot(x, 'reduction_method')
    } else {
      slot(x, 'reduction_method') = p$reduction_method
    }
  }

  if('nn_type' %in% s_names) {
    if(isTRUE(p$nospec_net)) {
      if(!is.na(slot(x, 'nn_type'))) p$nn_type = slot(x, 'nn_type')
    } else {
      slot(x, 'nn_type') = p$nn_type
    }
  }




  return(x)

}









## cell_ID slot ####

#' @title Get cell IDs for a given spatial unit
#' @name get_cell_id
#' @inheritParams data_access_params
#' @description Data for each spatial unit is expected to agree on a single set of cell_IDs
#' that are shared across any feature types. These cell_IDs are stored within the
#' giotto object's \code{cell_ID} slot. Getters and setters for this slot directly
#' retrieve (get) or replace (set) this slot.
#' @return character vector of cell_IDs
#' @seealso set_cell_id
#' @family functions to set data in giotto object
#' @keywords internal
get_cell_id = function(gobject,
                       spat_unit = NULL,
                       set_defaults = TRUE) {

  guard_against_notgiotto(gobject)
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
  }

  cell_IDs = slot(gobject, 'cell_ID')[[spat_unit]]
  cell_IDs = as.character(cell_IDs)

  return(cell_IDs)

}





#' @title Set cell IDs for a given spatial unit
#' @name set_cell_id
#' @inheritParams data_access_params
#' @param cell_IDs character vector of cell IDs to set. (See details)
#' @param verbose be verbose
#' @description Setter function for the cell_ID slot. Directly replaces (sets)
#' this slot
#' @details
#' Data for each spatial unit is expected to agree on a single set of cell_IDs
#' that are shared across any feature types. These cell_IDs are stored within the
#' giotto object's \code{cell_ID} slot. \cr
#'
#' Pass \code{NULL} to \code{cell_IDs} param in order to delete the entry. \cr
#' Pass \code{'initialize'} to \code{cell_IDs} param in order to initialize the
#' specified entry. \cr
#'
#' \strong{NOTE:} The main purpose of the setter is to initialize, as cell_ID values
#' are AUTOMATICALLY updated every time \code{initialize()} is called on the
#' giotto object.
#' @seealso get_cell_id
#' @return giotto object with set cell_ID slot
#' @family functions to set data in giotto object
#' @keywords internal
set_cell_id = function(gobject,
                       spat_unit = NULL,
                       cell_IDs,
                       set_defaults = TRUE,
                       verbose = TRUE) {

  guard_against_notgiotto(gobject)

  # set default spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
  }

  if(!is.null(cell_IDs)) {
    if(!inherits(cell_IDs, 'character')) stop('cell_IDs must be a character vector.')
  }

  # if input is 'initialize', RESET/reinitialize object
  if(identical(cell_IDs, 'initialize')) {
    if(isTRUE(verbose)) wrap_msg('Initializing', spat_unit, 'cell_IDs.')
    expr_avail = list_expression(gobject = gobject, spat_unit = spat_unit)
    si_avail = list_spatial_info(gobject = gobject)


    # get cell ID values
    if(spat_unit %in% expr_avail$spat_unit) { # preferred from expression

      cell_IDs = spatIDs(get_expression_values(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = expr_avail$feat_type[[1L]],
        values = expr_avail$name[[1L]],
        output = 'exprObj',
        set_defaults = TRUE
      ))

      # IDs = lapply(seq(nrow(expr_avail)), function(expr_i) {
      #   ex_ID = spatIDs(
      #     get_expression_values(
      #       gobject = gobject,
      #       spat_unit = spat_unit,
      #       feat_type = expr_avail$feat_type[[expr_i]],
      #       values = expr_avail$name[[expr_i]],
      #       output = 'exprObj'
      #     )
      #   )
      # })
      # cell_IDs = unique(unlist(IDs))
    } else if(spat_unit %in% si_avail$spat_info) { # fallback to spat_info

      cell_IDs = spatIDs(get_polygon_info(
        gobject = gobject,
        polygon_name = spat_unit,
        return_giottoPolygon = TRUE
      ))
    } else {
      # catch
      stop(wrap_txt('No data found to initialize cell_ID with', errWidth = TRUE))
    }
  }


  # set values
  cell_IDs = as.character(cell_IDs)
  slot(gobject, 'cell_ID')[[spat_unit]] = cell_IDs

  return(gobject)
}






## feat_ID slot ####

#' @title Get feat IDs for a given feature type
#' @name get_feat_id
#' @inheritParams data_access_params
#' @description Across a single modality/feature type, all feature information is
#' expected to share a single set of feat_IDs. These feat_IDs are stored within
#' the giotto object's \code{feat_ID} slot. Getters and setters for this slot
#' directly (get) or replace (set) this slot.
#' @seealso set_feat_id
#' @family functions to set data in giotto object
#' @keywords internal
get_feat_id = function(gobject,
                       feat_type = NULL,
                       set_defaults = TRUE) {
  guard_against_notgiotto(gobject)
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = NULL)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  feat_IDs = slot(gobject, 'feat_ID')[[feat_type]]
  feat_IDs = as.character(feat_IDs)
  return(feat_IDs)

}





#' @title Set feat IDs for a given feature type
#' @name set_feat_id
#' @inheritParams data_access_params
#' @param feat_IDs character vector of feature IDs to set.
#' @param verbose be verbose
#' @description Setter function for the feat_ID slot. Directly replaces (sets)
#' this slot
#' @details
#' Across a single modality/feature type, and within a spatial unit, all feature
#' information is expected to share a single set of feat_IDs. These feat_IDs
#' are stored within the giotto object's \code{feat_ID} slot separated by
#' feat_type. \cr
#'
#' Pass \code{NULL} to \code{feat_IDs} param in order to delete the entry. \cr
#' Pass \code{'initialize'} to \code{feat_IDs} param in order to initialize the
#' specified entry. \cr
#'
#' \strong{NOTE:} The main purpose of the setter is to initialize, as feat_ID values
#' are AUTOMATICALLY updated every time \code{initialize()} is called on the
#' giotto object.
#' @seealso get_feat_id
#' @return giotto object with set cell_ID slot
#' @family functions to set data in giotto object
#' @keywords internal
set_feat_id = function(gobject,
                       feat_type = NULL,
                       feat_IDs,
                       set_defaults = TRUE,
                       verbose = TRUE) {

  guard_against_notgiotto(gobject)

  if(isTRUE(set_defaults)) {
    if(identical(feat_IDs, 'initialize')) {
      spat_unit = suppressWarnings( # expected to be missing sometimes with init
        set_default_spat_unit(gobject = gobject,
                              spat_unit = NULL)
      )
    } else {
      spat_unit = set_default_spat_unit(gobject = gobject,
                                        spat_unit = NULL)
    }
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  if(!is.null(feat_IDs)) {
    if(!inherits(feat_IDs, 'character')) stop('feat_IDs must be a character vector.\n')
  }

  # initialize feat_ID
  if(identical(feat_IDs, 'initialize')) {
    expr_avail = list_expression(gobject = gobject, feat_type = feat_type)
    fi_avail = list_feature_info(gobject = gobject)

    if(feat_type %in% expr_avail$feat_type) { # preferred from expression
      # IDs = lapply(seq(nrow(expr_avail)), function(expr_i) {
      #   ex_ID = featIDs(
      #     get_expression_values(gobject = gobject,
      #                           spat_unit = expr_avail$spat_unit[[expr_i]],
      #                           feat_type = feat_type,
      #                           values = expr_avail$name[[expr_i]],
      #                           output = 'exprObj')
      #   )
      # })
      # feat_IDs = unique(unlist(IDs))

      feat_IDs = featIDs(get_expression_values(
        gobject = gobject,
        spat_unit = expr_avail$spat_unit[[1L]],
        feat_type = feat_type,
        values = expr_avail$name[[1L]],
        set_defaults = FALSE,
        output = 'exprObj'
      ))

    } else if(feat_type %in% fi_avail$feat_info) { # fallback to feature info

      feat_IDs = unique(featIDs(get_feature_info(
        gobject = gobject,
        feat_type = feat_type,
        return_giottoPoints = TRUE,
        set_defaults = FALSE
      )))
    } else {
      # catch
      stop(wrap_txt('No data found to intitialize feat_ID with', errWidth = TRUE))
    }
  }


  feat_IDs = as.character(feat_IDs)
  slot(gobject, 'feat_ID')[[feat_type]] = feat_IDs

  return(gobject)

}



## old metadata setters (originally from Giotto.R)

# #' @title Set cell metadata
# #' @name set_cell_metadata
# #' @description Set cell metadata
# #' @param gobject giotto object
# #' @param cell_metadata cell_metadata supplied as a nested list with spat_unit$feat_type$name
# #' @keywords internal
# set_cell_metadata = function(gobject,
#                              cell_metadata) {
#
#   # if metadata is not provided, then:
#   # create metadata for each spatial unit and feature type combination
#
#   # define for data.table :=
#   cell_ID = NULL
#
#   if(is.null(cell_metadata)) {
#
#     for(spat_unit in names(gobject@expression)) {
#
#       for(feat_type in names(gobject@expression[[spat_unit]])) {
#
#         if(is.null(gobject@spatial_info)) {
#           gobject@cell_metadata[[spat_unit]][[feat_type]] = data.table::data.table(cell_ID = gobject@cell_ID[[spat_unit]])
#         } else {
#           for(poly in names(gobject@spatial_info)) {
#             gobject@cell_metadata[[poly]][[feat_type]] = data.table::data.table(cell_ID = gobject@cell_ID[[poly]])
#           }
#         }
#       }
#     }
#
#   } else {
#
#     # extract all metadata information
#     # need to be nested list (feature type and spatial unit)
#     for(spat_unit in names(cell_metadata)) {
#
#       for(feat_type in names(cell_metadata[[spat_unit]])) {
#
#
#         gobject@cell_metadata[[spat_unit]][[feat_type]] = data.table::as.data.table(cell_metadata[[spat_unit]][[feat_type]])
#         gobject@cell_metadata[[spat_unit]][[feat_type]][, cell_ID := gobject@cell_ID[[spat_unit]]]
#
#         # put cell_ID first
#         all_colnames = colnames(gobject@cell_metadata[[spat_unit]][[feat_type]])
#         other_colnames = grep('cell_ID', all_colnames, invert = T, value = T)
#         gobject@cell_metadata[[spat_unit]][[feat_type]] = gobject@cell_metadata[[spat_unit]][[feat_type]][, c('cell_ID', other_colnames), with = FALSE]
#
#       }
#     }
#   }
#
#   return(gobject)
#
# }


# #' @title Set feature metadata
# #' @name set_feature_metadata
# #' @description Set feature metadata
# #' @keywords internal
# set_feature_metadata = function(gobject,
#                                 feat_metadata) {
#
#   # define for data.table :=
#   feat_ID = NULL
#
#   if(is.null(feat_metadata)) {
#
#     for(spat_unit in names(gobject@expression)) {
#       for(feat_type in names(gobject@expression[[spat_unit]])) {
#         gobject@feat_metadata[[spat_unit]][[feat_type]] = data.table::data.table(feat_ID = gobject@feat_ID[[feat_type]])
#       }
#     }
#
#   } else {
#
#     for(spat_unit in names(gobject@expression)) {
#       for(feat_type in names(gobject@expression[[spat_unit]])) {
#         gobject@feat_metadata[[spat_unit]][[feat_type]] = data.table::as.data.table(feat_metadata[[spat_unit]][[feat_type]])
#         gobject@feat_metadata[[spat_unit]][[feat_type]][, feat_ID := gobject@feat_ID[[feat_type]]]
#
#         # put feat_ID first
#         all_colnames = colnames(gobject@feat_metadata[[spat_unit]][[feat_type]])
#         other_colnames = grep('feat_ID', all_colnames, invert = T, value = T)
#         gobject@feat_metadata[[spat_unit]][[feat_type]] = gobject@feat_metadata[[spat_unit]][[feat_type]][, c('feat_ID', other_colnames), with = FALSE]
#
#       }
#     }
#   }
#
#   return(gobject)
#
# }


## cell metadata slot ####

#' @title Get cell metadata
#' @name get_cell_metadata
#' @inheritParams data_access_params
#' @param output return as either 'data.table' or 'cellMetaObj'
#' @description Get cell metadata from giotto object
#' @seealso pDataDT
#' @keywords internal
get_cell_metadata = function(gobject,
                             spat_unit = NULL,
                             feat_type = NULL,
                             output = c('cellMetaObj', 'data.table'),
                             copy_obj = TRUE,
                             set_defaults = TRUE) {

  output = match.arg(output, choices = c('cellMetaObj', 'data.table'))

  # 1. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 2. Find object - note that metadata objects do not have names
  cellMeta = gobject@cell_metadata[[spat_unit]][[feat_type]]

  if(inherits(cellMeta, 'list') | is.null(cellMeta)) stop('metadata referenced does not exist.')
  if(!inherits(cellMeta, 'cellMetaObj')) stop('metadata referenced is not cellMetaObj')

  if(isTRUE(copy_obj)) cellMeta[] = data.table::copy(cellMeta[])

  # 3. Return as desired object type
  if(output == 'cellMetaObj') {
    # return cellMetaObj
    return(cellMeta)

  } else if(output == 'data.table') {

    cellMeta = slot(cellMeta, 'metaDT')

    # return data.table
    return(cellMeta)
  }

}


#' @title getCellMetadata
#' @name getCellMetadata
#' @inheritParams data_access_params
#' @param output return as either 'data.table' or 'cellMetaObj'
#' @description Get cell metadata from giotto object
#' @seealso pDataDT
#' @export
getCellMetadata = function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL,
                           output = c('cellMetaObj', 'data.table'),
                           copy_obj = TRUE,
                           set_defaults = TRUE) {

  get_cell_metadata(gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    output = output,
                    copy_obj = copy_obj,
                    set_defaults = set_defaults)

}






#' @title Set cell metadata
#' @name setCellMetadata
#' @description Function to set cell metadata into giotto object
#' @inheritParams data_access_params
#' @param x cellMetaObj or list of cellMetaObj to set. Passing NULL will
#' reset a specified set of cell metadata in the giotto object.
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @return giotto object
#' @family functions to set data in giotto object
setCellMetadata = function(gobject,
                           x,
                           spat_unit = NULL,
                           feat_type = NULL,
                           provenance = NULL,
                           verbose = TRUE,
                           initialize = TRUE) {
  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x param (data to set) must be given',
                                        errWidth = TRUE))

  # check hierarchical slots
  used_su = list_cell_id_names(gobject)
  if(is.null(used_su))
    stop(wrap_txt('Add expression or spatial (polygon) information first'))

  # 1. Determine user inputs
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  .external_accessor_cellmeta = list(nospec_unit = nospec_unit,
                                     nospec_feat = nospec_feat)
  # checked by internal setter to determine if called by external

  # SINGLE INPUT
  # 2. if input is cellMetaObj or NULL, pass to internal
  if(is.null(x) | inherits(x, 'cellMetaObj')) {

    # pass to internal
    gobject = set_cell_metadata(
      gobject = gobject,
      metadata = x,
      spat_unit = spat_unit,
      feat_type = feat_type,
      provenance = provenance,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)
  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'cellMetaObj')) {

      # MULTIPLE INPUT
      # 3. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_cell_metadata(
          gobject = gobject,
          metadata = x[[obj_i]],
          spat_unit = spat_unit,
          feat_type = feat_type,
          provenance = provenance,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)
    }
  }

  # catch
  stop(wrap_txt('only cellMetaObj or lists of cellMetaObj accepted.
                For raw or external data, please first use readCellMetadata()'))

}






#' @title Set cell metadata
#' @name set_cell_metadata
#' @description Function to set cell metadata information into giotto object
#' @inheritParams data_access_params
#' @param provenance provenance information to set
#' @param metadata cellMetaObj or data.table containing cell metadata. Setting NULL will remove
#' the object. Passing 'initialize' will reset the object.
#' @param verbose be verbose
#' @return giotto object
#' @family functions to set data in giotto object
#' @keywords internal
set_cell_metadata = function(gobject,
                             metadata,
                             spat_unit = NULL,
                             feat_type = NULL,
                             provenance = NULL,
                             verbose = TRUE,
                             set_defaults = TRUE,
                             initialize = FALSE) {

  # data.table vars
  cell_ID = NULL

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(metadata)) stop(wrap_txt('metadata param must be given'))

  if(!inherits(metadata, c('cellMetaObj', 'NULL')) &
     !identical(metadata, 'initialize')) {
    stop(wrap_txt(deparse(substitute(metadata)), 'is not a cellMetaObj (set)
                  or NULL (remove)'))
  }


  # 1. determine if user input was supplied
  p = parent.frame() # get values if called from external
  call_from_external = exists('.external_accessor_cellmeta', where = p)

  if(call_from_external) {
    nospec_unit = p$.external_accessor_cellmeta$nospec_unit
    nospec_feat = p$.external_accessor_cellmeta$nospec_feat
  } else {
    nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
    nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  }


  # 2. set spat unit/ feat type if needed
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 3.1 if input is NULL, remove object
  if(is.null(metadata)) {
    if(isTRUE(verbose)) message('NULL passed to metadata.\n Removing specified metadata.')
    gobject@cell_metadata[[spat_unit]][[feat_type]] = NULL
    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }

  # 3.2 if input is 'initialize', RESET/reinitialize object
  if(inherits(metadata, 'character')) {
    if(identical(metadata, 'initialize')) {
      if(isTRUE(verbose)) message('Initializing specified metadata.')
      gobject@cell_metadata[[spat_unit]][[feat_type]] = create_cell_meta_obj(
        metaDT = data.table::data.table(cell_ID = get_cell_id(gobject,
                                                              spat_unit = spat_unit)),
        col_desc = c(cell_ID = 'cell-specific unique ID value'),
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = if(is.null(provenance)) spat_unit else provenance
      )
      if(isTRUE(initialize)) return(initialize(gobject))
      else return(gobject)
    }
  }


  # 4.1 import info if S4 object, else generate S4
  if(inherits(metadata, 'cellMetaObj')) {

    if(isTRUE(nospec_unit)) {
      if(!is.na(slot(metadata, 'spat_unit'))) spat_unit = slot(metadata, 'spat_unit')
    } else {
      slot(metadata, 'spat_unit') = spat_unit
    }
    if(isTRUE(nospec_feat)) {
      if(!is.na(slot(metadata, 'feat_type'))) feat_type = slot(metadata, 'feat_type')
    } else {
      slot(metadata, 'feat_type') = feat_type
    }
    if(!is.null(provenance)) {
      slot(metadata, 'provenance') = provenance
    }

  } else {

    # 4.2 if nested list structure, extract spat_unit/feat_type
    if(inherits(metadata, 'list')) {
      cellMetaObj_list = read_cell_metadata(gobject,
                                            metadata = metadata,
                                            provenance = if(is.null(provenance)) spat_unit else provenance)
      # recursively run
      for(obj_i in seq_along(cellMetaObj_list)) {
        # (provenance info set during prev. step)
        gobject = set_cell_metadata(gobject,
                                    metadata = cellMetaObj_list[[obj_i]])
      }
      if(isTRUE(initialize)) initialize(gobject)
      else return(gobject)
    }

    # 4.3 otherwise assume data.frame type object

    if(is.null(spat_unit) | is.null(feat_type)) stop(wrap_txt(
      'Add expression or polygon info first
        Alternatively, specify expected spat_unit and feat_type using activeSpatUnit() and activeFeatType()',
      errWidth = TRUE
    ))


    metadata = data.table::as.data.table(metadata)

    # if cell ID col is missing, try to automatically set
    if(is.null(metadata[['cell_ID']])) {
      id_error = try(metadata[, cell_ID := get_cell_id(gobject, spat_unit = spat_unit)], silent = TRUE)
      if(inherits(id_error, 'try-error')) stop('cannot automatically set metadata cell_ID based on gobject cell_ID slot.')
    } else if(spat_unit %in% list_cell_id_names(gobject)) {

      # if cell ID col is present in both, try to match
      if(!identical(metadata[, cell_ID], get_cell_id(gobject, spat_unit = spat_unit))) {
        stop('metadata cell_ID does not match that in gobject cell_ID slot for spat_unit "', spat_unit, '".\n')
      }

    }

    # put cell_ID first
    all_colnames = colnames(metadata)
    other_colnames = grep('cell_ID', all_colnames, invert = TRUE, value = TRUE)
    metadata = metadata[, c('cell_ID', other_colnames), with = FALSE]

    metadata = new('cellMetaObj',
                   metaDT = metadata,
                   col_desc = NA_character_, # unknown
                   spat_unit = spat_unit,
                   feat_type = feat_type,
                   provenance = if(is.null(provenance)) spat_unit else provenance)
  }

  # 5. check if nesting address is already used - just feat_type for metadata
  potential_names = list_cell_metadata(gobject, spat_unit = spat_unit)[, feat_type]
  if(feat_type %in% potential_names) {
    if(isTRUE(verbose)) wrap_msg('> Cell metadata for spat_unit "', spat_unit, '" and feat_type "',
                                 feat_type, '" already exists and will be replaced with new metadata.')
  }


  # old_meta = data.table::copy(gobject@cell_metadata[[spat_unit]][[feat_type]])
  #
  # new_cols = NULL
  # if(!is.null(old_meta) && inherits(old_meta, "data.table")){
  #   old_cols = colnames(old_meta)
  #   new_cols = colnames(meta_dt)
  #   cat("\nWarning: The following columns will be overwritten:",
  #   old_cols, "\nThey will be replaced with:", new_cols,"\n")
  #
  #   old_meta[, (old_cols) := NULL]
  #
  # } else {
  #   new_cols = colnames(meta_dt)
  #   cat("\nWriting data within columns:", new_cols,
  #   "\nto cell metadata\n")
  # }

  # suppressWarnings({
  #   gobject <- removeCellAnnotation(gobject = gobject,
  #                                   spat_unit = spat_unit,
  #                                   feat_type = feat_type,
  #                                   columns = old_cols,
  #                                   return_gobject = TRUE)
  #   gobject@cell_metadata[[spat_unit]][[feat_type]] = meta_dt
  # })
  # message("\nCell Metadata slot '",spat_unit, feat_type, "' set.\n")
  #

  # 6. set and return gobject
  gobject@cell_metadata[[spat_unit]][[feat_type]] = metadata

  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting cell metadata [', spatUnit(metadata), '][', featType(metadata), '] ',
    sep = ''
  )

  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}


## feature metadata slot ####

#' @title Get feature metadata
#' @name get_feature_metadata
#' @inheritParams data_access_params
#' @param output return as either 'data.table' or 'featMetaObj'
#' @param copy_obj whether to perform a deepcopy of the data.table information
#' @description Get feature metadata from giotto object
#' @seealso fDataDT
#' @keywords internal
get_feature_metadata = function(gobject,
                                spat_unit = NULL,
                                feat_type = NULL,
                                output = c('featMetaObj', 'data.table'),
                                copy_obj = TRUE,
                                set_defaults = TRUE) {

  output = match.arg(output, choices = c('featMetaObj', 'data.table'))

  # 1. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 2. Find object - note that metadata objects do not have names
  featMeta = gobject@feat_metadata[[spat_unit]][[feat_type]]
  if(is.null(featMeta)) stop('metadata referenced does not exist')
  if(!inherits(featMeta, 'featMetaObj')) stop('metadata referenced is not featMetaObj')

  # 3. Return as desired object type

  if(isTRUE(copy_obj)) featMeta[] = data.table::copy(featMeta[])

  if(output == 'featMetaObj') return(featMeta)
  if(output == 'data.table') return(featMeta[])

}


#' @title getFeatureMetadata
#' @name getFeatureMetadata
#' @inheritParams data_access_params
#' @param output return as either 'data.table' or 'featMetaObj'
#' @param copy_obj whether to perform a deepcopy of the data.table information
#' @description Get feature metadata from giotto object
#' @seealso fDataDT
#' @export
getFeatureMetadata = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              output = c('featMetaObj', 'data.table'),
                              copy_obj = TRUE,
                              set_defaults = TRUE) {

  get_feature_metadata(gobject = gobject,
                       spat_unit = spat_unit,
                       feat_type = feat_type,
                       output = output,
                       copy_obj = copy_obj,
                       set_defaults = set_defaults)

}






#' @title Set feature metadata
#' @name setFeatureMetadata
#' @description Function to set feature metadata into giotto object
#' @inheritParams data_access_params
#' @param x featMetaObj or list of featMetaObj to set. Passing NULL will
#' reset a specified set of feature metadata in the giotto object.
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @return giotto object
#' @family functions to set data in giotto object
setFeatureMetadata = function(gobject,
                              x,
                              spat_unit = NULL,
                              feat_type = NULL,
                              provenance = NULL,
                              verbose = TRUE,
                              initialize = TRUE) {
  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x param (data to set) must be given',
                                        errWidth = TRUE))

  # 1. Determine user inputs
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  .external_accessor_featmeta = list(nospec_unit = nospec_unit,
                                     nospec_feat = nospec_feat)
  # checked by internal setter to determine if called by external

  # SINGLE INPUT
  # 2. if input is featMetaObj or NULL, pass to internal
  if(is.null(x) | inherits(x, 'featMetaObj')) {

    # pass to internal
    gobject = set_feature_metadata(
      gobject = gobject,
      metadata = x,
      spat_unit = spat_unit,
      feat_type = feat_type,
      provenance = provenance,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)
  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'featMetaObj')) {

      # MULTIPLE INPUT
      # 3. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_feature_metadata(
          gobject = gobject,
          metadata = x[[obj_i]],
          spat_unit = spat_unit,
          feat_type = feat_type,
          provenance = provenance,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)
    }
  }

  # catch
  stop(wrap_txt('only featMetaObj or lists of featMetaObj accepted.
                For raw or external data, please first use readFeatMetadata()'))

}







#' @keywords internal
#' @noRd
set_feature_metadata = function(gobject,
                                metadata,
                                spat_unit = NULL,
                                feat_type = NULL,
                                provenance = NULL,
                                verbose = TRUE,
                                set_defaults = TRUE,
                                initialize = FALSE) {

  # data.table vars
  feat_ID = NULL

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(metadata)) stop(wrap_txt('metadata param must be given'))

  if(!inherits(metadata, c('featMetaObj', 'NULL')) &
     !identical(metadata, 'initialize')) {
    stop(wrap_txt(deparse(substitute(metadata)), 'is not featMetaObj (set)
                  or NULL (remove)'))
  }


  # 1. determine if user input was supplied
  p = parent.frame() # get values if called from external
  call_from_external = exists('.external_accessor_featmeta', where = p)

  if(call_from_external) {
    nospec_unit = p$.external_accessor_featmeta$nospec_unit
    nospec_feat = p$.external_accessor_featmeta$nospec_feat
  } else {
    nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
    nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  }


  # 2. set spat unit/ feat type if needed
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 3.1 if input is NULL, remove object
  if(is.null(metadata)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to metadata.\n Removing specified metadata.')
    gobject@feat_metadata[[spat_unit]][[feat_type]] = NULL
    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }

  # 3.2 if input is 'initialize', RESET/reinitialize object
  if(inherits(metadata, 'character')) {
    if(metadata == 'initialize') {
      if(isTRUE(verbose)) message('Initializing specified metadata.')
      gobject@feat_metadata[[spat_unit]][[feat_type]] = create_feat_meta_obj(
        metaDT = data.table::data.table(feat_ID = get_feat_id(gobject,
                                                              feat_type = feat_type)),
        col_desc = c(feat_ID = 'feature-specific unique ID value'),
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = if(is.null(provenance)) spat_unit else provenance
      )
      if(isTRUE(initialize)) return(initialize(gobject))
      else return(gobject)
    }
  }


  # 4.1 import info if S4 object, else generate S4
  if(inherits(metadata, 'featMetaObj')) {

    if(isTRUE(nospec_unit)) {
      if(!is.na(slot(metadata, 'spat_unit'))) spat_unit = slot(metadata, 'spat_unit')
    } else {
      slot(metadata, 'spat_unit') = spat_unit
    }
    if(isTRUE(nospec_feat)) {
      if(!is.na(slot(metadata, 'feat_type'))) feat_type = slot(metadata, 'feat_type')
    } else {
      slot(metadata, 'feat_type') = feat_type
    }
    if(!is.null(provenance)) {
      slot(metadata, 'provenance') = provenance
    }

  } else {

    # 4.2 if nested list structure, extract spat_unit/feat_type
    if(inherits(metadata, 'list')) {
      featMetaObj_list = readFeatMetadata(data_list = metadata,
                                          provenance = if(is.null(provenance)) spat_unit else provenance)
      # recursively run
      for(obj_i in seq_along(featMetaObj_list)) {
        # (provenance info set during prev. step)
        gobject = set_feature_metadata(gobject,
                                       metadata = featMetaObj_list[[obj_i]])
      }
      return(gobject)
    }

    # 4.3 otherwise assume data.frame type object

    if(is.null(spat_unit) | is.null(feat_type)) stop(wrap_txt(
      'Add expression or polygon info first
        Alternatively, specify expected spat_unit and feat_type using activeSpatUnit() and activeFeatType()',
      errWidth = TRUE
    ))


    metadata = data.table::as.data.table(metadata)

    # if feat ID col is missing, try to automatically set
    if(is.null(metadata[['feat_ID']])) {
      id_error = try(metadata[, feat_ID := get_feat_id(gobject, feat_type = feat_type)], silent = TRUE)
      if(inherits(id_error, 'try-error')) stop('cannot automatically set metadata feat_ID based on gobject feat_ID slot.')
    } else if(feat_type %in% list_feat_id_names(gobject)) {

      # if feat ID col is present in both, try to match
      if(!identical(metadata[, feat_ID], get_feat_id(gobject, feat_type = feat_type))) {
        stop('metadata feat_ID does not match that in gobject feat_ID slot for feat_type "', feat_type, '".\n')
      }

    }

    # put feat_ID first
    all_colnames = colnames(metadata)
    other_colnames = grep('feat_ID', all_colnames, invert = TRUE, value = TRUE)
    metadata = metadata[, c('feat_ID', other_colnames), with = FALSE]

    metadata = new('featMetaObj',
                   metaDT = metadata,
                   col_desc = NA_character_, # unknown
                   spat_unit = spat_unit,
                   feat_type = feat_type,
                   provenance = if(is.null(provenance)) spat_unit else provenance)
  }

  # 5. check if nesting address is already used - just feat_type for metadata
  potential_names = list_feat_metadata(gobject, spat_unit = spat_unit)[, feat_type]
  if(feat_type %in% potential_names) {
    if(isTRUE(verbose)) wrap_msg('> Feat metadata for spat_unit "', spat_unit, '" and feat_type "',
                                 feat_type, '" already exists and will be replaced with new metadata.')
  }


  # 6. return object
  gobject@feat_metadata[[spat_unit]][[feat_type]] = metadata

  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting feature metadata [', spatUnit(metadata), '][', featType(metadata), '] ',
    sep = ''
  )

  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}


## expression values slot ####




#' @title Get expression values
#' @name getExpression
#' @description Function to get expression values from giotto object
#' @inheritParams data_access_params
#' @param values expression values to extract (e.g. "raw", "normalized", "scaled")
#' @param output what object type to retrieve the expression as. Currently either
#' 'matrix' for the matrix object contained in the exprObj or 'exprObj' (default) for
#' the exprObj itself are allowed.
#' @return exprObj or matrix depending on output param
#' @family expression accessor functions
#' @family functions to get data from giotto object
#' @export
getExpression = function(gobject,
                         values = NULL,
                         spat_unit = NULL,
                         feat_type = NULL,
                         output = c('exprObj', 'matrix'),
                         set_defaults = TRUE) {

  # 0. Check input
  guard_against_notgiotto(gobject)
  output = match.arg(output, choices = c('exprObj', 'matrix'))



  # 1. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }



  # 2. Find object
  potential_values = list_expression_names(gobject = gobject,
                                           spat_unit = spat_unit,
                                           feat_type = feat_type)

  if(is.null(values)) values = potential_values[[1]]


  # 3. Get object in desired format
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values,
                                      output = output,
                                      set_defaults = FALSE)

  return(expr_values)


}






# Internal function to get expression values from giotto object.
# **Controls expression slot nesting and structure**
get_expression_values = function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 values = NULL,
                                 output = c('exprObj', 'matrix'),
                                 set_defaults = TRUE) {

  guard_against_notgiotto(gobject)

  output = match.arg(output, choices = c('exprObj', 'matrix'))

  # 1. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 2. Find object
  potential_values = list_expression_names(gobject = gobject,
                                           spat_unit = spat_unit,
                                           feat_type = feat_type)

  if(is.null(values)) values = potential_values[[1]]
  if(is.null(values)) {
    stop(wrap_txt('No expression values discovered by getter:',
                  '\nspat_unit:', spat_unit,
                  '\nfeat_type:', feat_type))
  }


  ## special checks/cases for giotto standard pipeline
  if(values == 'scaled' & !'scaled' %in% potential_values) {
    stop(wrap_txt('Scaled expression not found.
                  First run scaling (& normalization) step(s)', errWidth = TRUE))
  } else if(values == 'normalized' & !'normalized' %in% potential_values) {
    stop(wrap_txt('Normalized expression not found.
                  First run normalization step', errWidth = TRUE))
  } else if(values == 'custom' & !'custom' %in% potential_values) {
    stop(wrap_txt('Custom expression not found.
                  First add custom expression matrix', errWidth = TRUE))
  }

  if(!values %in% potential_values) stop(wrap_txt(
    'Requested expression info not found
    [spat_unit:', spat_unit, '] [feat_type:', feat_type, '] [values:', values, ']',
    sep = '',
    errWidth = TRUE
  ))


  # Get info from slot nesting
  expr_vals = gobject@expression[[spat_unit]][[feat_type]][[values]]

  if(output == 'exprObj') return(expr_vals)
  else if(output == 'matrix') return(expr_vals[])

}





#' @description Get all expression values for a specified spatial unit
#' and feature type
#' @keywords internal
#' @return list of exprObj
#' @noRd
get_expression_values_list = function(gobject,
                                      spat_unit = NULL,
                                      feat_type = NULL,
                                      output = c('exprObj', 'matrix'),
                                      set_defaults = TRUE) {

  guard_against_notgiotto(gobject)

  output = match.arg(output, choices = c('exprObj', 'matrix'))

  if(isTRUE(set_defaults)) {
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)
  }

  data_list = slot(gobject, 'expression')[[spat_unit]][[feat_type]]

  if(output == 'exprObj') {
    return(data_list)
  }
  if(output == 'matrix') {
    return(lapply(data_list, `[`))
  }
}







#' @title Set expression data
#' @name setExpression
#' @description Function to set expression values for giotto object.
#' @inheritParams data_access_params
#' @param x exprObj or list of exprObj to set. Passing NULL will remove a
#' specified set of expression data from the giotto object
#' @param name name for the expression information
#' @param provenance provenance information (optional)
#' information for the giotto object. Pass NULL to remove an expression object
#' @param verbose be verbose
#' @return giotto object
#' @family expression accessor functions
#' @family functions to set data in giotto object
#' @export
setExpression = function(gobject,
                         x,
                         spat_unit = NULL,
                         feat_type = NULL,
                         name = 'raw',
                         provenance = NULL,
                         verbose = TRUE,
                         initialize = TRUE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x param (data to set) must be given'))

  # 1. Determine user inputs
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
  .external_accessor_expr = list(nospec_unit = nospec_unit,
                                 nospec_feat = nospec_feat,
                                 nospec_name = nospec_name)
  # checked by internal setter to determine if called by external


  # SINGLE INPUT
  # 3. if input is exprObj or NULL, pass to internal
  if(is.null(x) | inherits(x, 'exprObj')) {

    # pass to internal
    gobject = set_expression_values(
      gobject = gobject,
      values = x,
      spat_unit = spat_unit,
      feat_type = feat_type,
      name = name,
      provenance = provenance,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)

  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'exprObj')) {

      # MULTIPLE INPUT
      # 4. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_expression_values(
          gobject = gobject,
          values = x[[obj_i]],
          spat_unit = spat_unit,
          feat_type = feat_type,
          name = name,
          provenance = provenance,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)

    }
  }

  # catch
  stop(wrap_txt('Only exprObj or lists of exprObj accepted.
                    For raw or external data, please first use readExprData()'))

}









#' @title Set expression values
#' @name set_expression_values
#' @description Function to set expression values for giotto object
#' @inheritParams data_access_params
#' @param name name for the expression slot
#' @param provenance provenance information (optional)
#' @param values exprObj If NULL, then the object will be removed.
#' @param verbose be verbose
#' @param initialize (default = FALSE) whether to initialize the gobject before
#' returning. Will be set to TRUE when called by the external
#' @return giotto object
#' @family expression accessor functions
#' @family functions to set data in giotto object
#' @export
set_expression_values = function(gobject,
                                 values,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 name = 'test',
                                 provenance = NULL,
                                 verbose = TRUE,
                                 set_defaults = TRUE,
                                 initialize = FALSE) {

  guard_against_notgiotto(gobject)

  if(!inherits(values, c('exprObj', 'NULL'))) {
    stop(wrap_txt(deparse(substitute(values)), 'is not exprObj (set)
                  or NULL (remove)'))
  }

  # 1. Determine user inputs
  p = parent.frame() # Get values if called from external
  call_from_external = exists('.external_accessor_expr', where = p)

  if(call_from_external) {
    nospec_unit = p$.external_accessor_expr$nospec_unit
    nospec_feat = p$.external_accessor_expr$nospec_feat
    nospec_name = p$.external_accessor_expr$nospec_name
  } else {
    nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
    nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
    nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
  }

  if(inherits(values, 'exprObj')) {
    if(!is.na(spatUnit(values)) &
       !is.na(featType(values)) &
       isTRUE(nospec_unit) &
       isTRUE(nospec_feat)) {
      set_defaults = FALSE
    }
  }


  # 2. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }


  # 3. if input is NULL, remove object (no initialize option)
  if(is.null(values)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to values param.
                                 Removing specified expression')
    gobject@expression[[spat_unit]][[feat_type]][[name]] = NULL

    # prune if empty
    if(length(gobject@expression[[spat_unit]][[feat_type]]) == 0) {
      gobject@expression[[spat_unit]] = NULL
      if(length(gobject@expression[[spat_unit]] == 0)) {
        gobject@expression = list()
      }
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }


  # 4 import data from S4 if available
  # NOTE: modifies spat_unit/feat_type/name/provenance/data slots
  values = read_s4_nesting(values)

  # 5. check if specified name has already been used
  potential_names = list_expression_names(gobject, spat_unit = spat_unit, feat_type = feat_type)
  if(name %in% potential_names) {
    if(isTRUE(verbose)) wrap_msg('> ', name, ' already exists and will be replaced with new values \n')
  }

  ## 6. update and return giotto object
  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting expression [', spatUnit(values), '][', featType(values), '] ',
    objName(values), sep = ''
  )

  gobject@expression[[spat_unit]][[feat_type]][[name]] = values
  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}







## spatial locations slot ####





#' @title Get spatial locations
#' @name getSpatialLocations
#' @description Function to get a spatial location data.table
#' @inheritParams data_access_params
#' @param spat_loc_name name of spatial locations (defaults to first name in spatial_locs slot, e.g. "raw")
#' @param output what object type to get the spatial locations as. Default is as
#' a 'spatLocsObj'. Returning as 'data.table' is also possible.
#' @param copy_obj whether to copy/duplicate when getting the object (default = TRUE)
#' @param verbose be verbose
#' @return data.table with coordinates or spatLocsObj depending on \code{output}
#' @family spatial location data accessor functions
#' @family functions to get data from giotto object
#' @export
getSpatialLocations = function(gobject,
                               spat_unit = NULL,
                               name = NULL,
                               output = c('spatLocsObj', 'data.table'),
                               copy_obj = TRUE,
                               verbose = TRUE,
                               set_defaults = TRUE) {

  # Pass to internal function
  spatloc = get_spatial_locations(gobject = gobject,
                                  spat_unit = spat_unit,
                                  spat_loc_name = name,
                                  output = output,
                                  copy_obj = copy_obj,
                                  verbose = verbose,
                                  set_defaults = set_defaults)

  return(spatloc)
}








#' @title Get spatial locations
#' @name get_spatial_locations
#' @description Function to get a spatial location data.table
#' @inheritParams data_access_params
#' @param spat_loc_name name of spatial locations (defaults to first name in spatial_locs slot, e.g. "raw")
#' @param output what object type to get the spatial locations as. Default is as
#' a 'spatLocsObj'. Returning as 'data.table' is also possible.
#' @param copy_obj whether to copy/duplicate when getting the object (default = TRUE)
#' @param verbose be verbose
#' @return data.table with coordinates or spatLocsObj depending on \code{output}
#' @family spatial location data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatial_locations = function(gobject,
                                 spat_unit = NULL,
                                 spat_loc_name = NULL,
                                 output = c('spatLocsObj', 'data.table'),
                                 copy_obj = TRUE,
                                 verbose = TRUE,
                                 set_defaults = TRUE) {


  output = match.arg(output, choices = c('spatLocsObj', 'data.table'))

  if(!is.null(spat_unit)) {
    potential_spat_unit = names(slot(gobject, 'spatial_locs'))
    if(!spat_unit %in% potential_spat_unit) stop(wrap_txt('No spatial locations for spatial unit "', spat_unit, '" exist in spatial_locs slot.'))
  }

  # spatial locations
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
  }


  # if NULL (not given) and spatial locations have been added, then use first one
  # if NULL (not given) and spatial locations have NOT been added, then keep NULL
  if(is.null(spat_loc_name)) {
    if(!is.null(slot(gobject, 'spatial_locs'))) {
      spat_loc_name = names(slot(gobject, 'spatial_locs')[[spat_unit]])[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      spat_loc_name = NULL
      cat('No spatial locations have been found \n')
      return(NULL)
    }
  }

  potential_names = names(slot(gobject, 'spatial_locs')[[spat_unit]])

  # get object to return in desired format
  if(spat_loc_name %in% potential_names) { # case if found
    if(output == 'data.table') { # case if return as data.table
      if(inherits(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]], 'spatLocsObj')) { # if is spatLocsObj
        if(isTRUE(copy_obj)) spatloc = data.table::copy(slot(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]], 'coordinates')) # copy
        else spatloc = slot(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]], 'coordinates') # no copy
      } else { # to be deprecated | if is data.table
        if(isTRUE(copy_obj)) spatloc = data.table::copy(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]]) # copy
        else spatloc = slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]] # no copy
      }
    } else if(output == 'spatLocsObj') { # case if return as spatLocsObj
      if(inherits(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]], 'spatLocsObj')) { # if is spatLocsObj
        if(isTRUE(copy_obj)) spatloc = copy(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]]) # copy
        else spatloc = slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]] # no copy
      } else { # to be deprecated | if is data.table
        if(isTRUE(copy_obj)) { # copy
          spatloc = new('spatLocsObj',
                        name = spat_loc_name,
                        spat_unit = spat_unit,
                        coordinates = data.table::copy(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]]),
                        provenance = spat_unit) # assume
        } else { # no copy
          spatloc = new('spatLocsObj',
                        name = spat_loc_name,
                        spat_unit = spat_unit,
                        coordinates = slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]],
                        provenance = spat_unit) # assume
        }
      }
    } else {
      if(isTRUE(verbose)) wrap_msg('Other spatial locations outputs not yet supported')
    }
    return(spatloc)
  } else { # case if not found
    stop(wrap_txt("The spatial locations with name ","'", spat_loc_name, "'"," can not be found \n"))
  }
}






#' @description Get all spatial locations for a specified spatial unit
#' @keywords internal
#' @return list of spatLocsObj or data.tables depending on output param
#' @noRd
get_spatial_locations_list = function(gobject,
                                      spat_unit = NULL,
                                      output = c('spatLocsObj', 'data.table'),
                                      copy_obj = TRUE,
                                      verbose = TRUE,
                                      set_defaults = TRUE) {

  guard_against_notgiotto(gobject)

  output = match.arg(output, choices = c('spatLocsObj', 'data.table'))

  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
  }

  data_list = slot(gobject, 'spatial_locs')[[spat_unit]]

  if(isTRUE(copy_obj)) {
    data_list = lapply(data_list, copy)
  }

  if(output == 'spatLocsObj') {
    return(data_list)
  }
  if(output == 'data.table') {
    return(lapply(data_list, `[`))
  }
}








#' @title Set spatial locations
#' @name setSpatialLocations
#' @description Function to set a spatial location slot
#' @inheritParams data_access_params
#' @param x spatLocsObj or list of spatLocsObj. Passing NULL will remove a
#' specified set of spatial locations data.
#' @param name name of spatial locations, default "raw"
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @details Spatial information will be set to the nested location described
#' by their tagged spat_unit and name information. An alternative location can also
#' be specified through the respective params in this function.
#' @return giotto object
#' @family spatial location data accessor functions
#' @family functions to set data in giotto object
#' @export
setSpatialLocations = function(gobject,
                               x,
                               spat_unit = NULL,
                               name = 'raw',
                               provenance = NULL,
                               verbose = TRUE,
                               initialize = TRUE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x (data to set) param must be given'))


  # check hierarchical slots
  avail_ex = list_expression(gobject)
  avail_si = list_spatial_info(gobject)
  if(is.null(avail_ex) & is.null(avail_si))
    stop(wrap_txt('Add expression or spatial (polygon) information first'))


  # 1. Determine user inputs
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
  .external_accessor_spatloc = list(nospec_unit = nospec_unit,
                                    nospec_name = nospec_name)
  # checked by internal setter to determine if called by external


  # NATIVE INPUT TYPES
  # 3. If input is spatLocsObj or NULL, pass to internal
  if(is.null(x) | inherits(x, 'spatLocsObj')) {

    # pass to internal
    gobject = set_spatial_locations(
      gobject = gobject,
      spatlocs = x,
      spat_unit = spat_unit,
      spat_loc_name = name,
      provenance = provenance,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)

  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'spatLocsObj')) {

      # MULTIPLE INPUT
      # 4. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_spatial_locations(
          gobject = gobject,
          spatlocs = x[[obj_i]],
          spat_unit = spat_unit,
          spat_loc_name = name,
          provenance = provenance,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)

    }
  }

  # catch
  stop(wrap_txt('Only spatLocsObj or lists of spatLocsObj accepted.
                For raw or external data, please first use readSpatLocsData()'))

}









#' @title Set spatial locations
#' @name set_spatial_locations
#' @description Function to set a spatial location slot
#' @inheritParams data_access_params
#' @param spatlocs spatial locations (accepts either \code{data.table} or
#' \code{spatLocsObj})
#' @param spat_loc_name name of spatial locations, default "raw"
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @details If a \code{spatLocsObj} is provided to \code{spatlocs} param then any
#' attached name and spat_unit info will be used for input to this function's
#' \code{spat_loc_name} and \code{spat_unit}params, BUT will be overridden by any
#' alternative specific inputs to those params. \cr
#' ie: a \code{spatLocsObj} with spat_unit slot == 'cell' will be automatically
#' nested by spat_unit 'cell' when using \code{set_spatial_locations} as long as
#' param \code{spat_unit = NULL}. BUT if param \code{spat_unit = 'nucleus'} then
#' the \code{spatLocsObj} will be nested by spat_unit 'nucleus' instead and
#' its spat_unit slot will be changed to 'nucleus'
#' @return giotto object
#' @family spatial location data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatial_locations = function(gobject,
                                 spatlocs,
                                 spat_unit = NULL,
                                 spat_loc_name = 'raw',
                                 provenance = NULL,
                                 verbose = TRUE,
                                 set_defaults = TRUE,
                                 initialize = FALSE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(spatlocs)) stop(wrap_txt('spatlocs param must be given'))

  # 0. pass to external if not native formats
  if(!inherits(spatlocs, c('spatLocsObj', 'NULL'))) {
    stop(wrap_txt(deparse(substitute(spatlocs)), 'is not spatLocsObj (set)
                  or NULL (remove)'))
    gobject = setSpatialLocations(gobject = gobject,
                                  x = spatlocs,
                                  spat_unit = spat_unit,
                                  name = spat_loc_name,
                                  provenance = provenance,
                                  verbose = verbose)
    return(gobject) # initialize done already
  }

  # 1. determine if input was supplied to spat_unit and spat_loc_name
  p = parent.frame() # Get values if called from external
  call_from_external = exists('.external_accessor_spatloc', where = p)

  if(isTRUE(call_from_external)) {
    nospec_unit = p$.external_accessor_spatloc$nospec_unit
    nospec_name = p$.external_accessor_spatloc$nospec_name
  } else {
    nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
    nospec_name = methods::hasArg(spat_loc_name)
  }

  # use name from now on for compatiblity with S4 reading
  name = spat_loc_name

  # 2. set spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
  }

  # 3. remove if input is NULL
  if(is.null(spatlocs)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to spatlocs.
                                 Removing specified spatial locations.')
    gobject@spatial_locs[[spat_unit]][[name]] = NULL

    # prune if empty
    if(length(gobject@spatial_locs[[spat_unit]]) == 0) {
      gobject@spatial_locs = NULL
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }

  # 4. import data from S4 if available
  # NOTE: modifies spat_unit/name/provenance
  spatlocs = read_s4_nesting(spatlocs)

  # 5. check if specified name has already been used
  potential_names = list_spatial_locations_names(gobject, spat_unit = spat_unit)
  if(name %in% potential_names) {
    if(isTRUE(verbose)) {
      wrap_msg('> ', name, ' already exists and will be replaced with new spatial locations \n')
    }
  }

  # 6. update and return giotto object
  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting spatial locations [', spatUnit(spatlocs), '] ',
    objName(spatlocs), sep = ''
  )

  gobject@spatial_locs[[spat_unit]][[name]] = spatlocs
  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}







## dimension reduction slot ####

#' @title Get dimension reduction
#' @name get_dimReduction
#' @inheritParams data_access_params
#' @param reduction reduction on cells or features (e.g. "cells", "feats")
#' @param reduction_method reduction method (e.g. "pca", "umap", "tsne")
#' @param name name of reduction results
#' @param output object type to return as. Either 'dimObj' (default) or 'matrix'
#' of the embedding coordinates.
#' @description Function to get a dimension reduction object
#' @return dim reduction object (default) or dim reduction coordinates
#' @family dimensional reduction data accessor functions
#' @family functions to get data from giotto object
#' @export
get_dimReduction = function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            reduction = c('cells', 'feats'),
                            reduction_method = c('pca', 'umap', 'tsne'),
                            name = 'pca',
                            output = c('dimObj', 'matrix'),
                            set_defaults = TRUE) {

  guard_against_notgiotto(gobject)

  # to be deprecated ('data.table' -> 'matrix')
  if(!identical(output, c('dimObj', 'matrix'))) {
    if(output == 'data.table') output = 'matrix'
  }
  output = match.arg(output, choices = c('dimObj', 'matrix'))

  # Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  ## check parameters
  reduction = match.arg(arg = reduction, choices = c('cells', 'feats'))
  reduction_method = match.arg(arg = reduction_method, choices = unique(c('pca', 'umap', 'tsne', reduction_method)))

  ## check reduction
  reduction_res = slot(gobject, 'dimension_reduction')[[reduction]][[spat_unit]][[feat_type]]
  if(is.null(reduction_res)) {
    stop('No dimension reduction for ', reduction, ' has been applied \n')
  }

  ## check method
  reduction_res = reduction_res[[reduction_method]]
  if(is.null(reduction_res)) {
    stop(reduction_method, ' has not been performed on this dataset \n')
  }

  ## check name for method
  reduction_res = reduction_res[[name]]
  if(is.null(reduction_res)) {
    stop(name, ': this name is not available for method: ', reduction_method, '\n')
  }

  ## S3 backwards compatibility
  if(!isS4(reduction_res)) reduction_res = S3toS4dimObj(reduction_res)
  silent = validObject(reduction_res) # variable hides TRUE print

  ## return object or coordinates
  if(output == 'dimObj') {
    return(reduction_res)
  } else if(output == 'matrix') {
    return(slot(reduction_res, 'coordinates'))
  } else {
    stop('other outputs not supported')
  }

}



#' @title Get dimension reduction
#' @name getDimReduction
#' @inheritParams data_access_params
#' @param reduction reduction on cells or features (e.g. "cells", "feats")
#' @param reduction_method reduction method (e.g. "pca", "umap", "tsne")
#' @param name name of reduction results
#' @param output object type to return as. Either 'dimObj' (default) or 'matrix'
#' of the embedding coordinates.
#' @description Function to get a dimension reduction object
#' @return dim reduction object (default) or dim reduction coordinates
#' @family dimensional reduction data accessor functions
#' @family functions to get data from giotto object
#' @export
getDimReduction = function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL,
                           reduction = c('cells', 'feats'),
                           reduction_method = c('pca', 'umap', 'tsne'),
                           name = 'pca',
                           output = c('dimObj', 'matrix'),
                           set_defaults = TRUE) {

  # pass to internal
  dimRed = get_dimReduction(gobject = gobject,
                            spat_unit = spat_unit,
                            feat_type = feat_type,
                            reduction = reduction,
                            reduction_method = reduction_method,
                            name = name,
                            output = output,
                            set_defaults = set_defaults)

  return(dimRed)
}








#' @description Get all dimension reductions for a specified spatial unit,
#' feature type, and reduction (either on 'cells' or 'feats')
#' @keywords internal
#' @return list of dimObj or matrix depending on output param
#' @noRd
get_dim_reduction_list = function(gobject,
                                  spat_unit = NULL,
                                  feat_type = NULL,
                                  reduction = c('cells', 'feats'),
                                  output = c('dimObj', 'matrix'),
                                  set_defaults = TRUE) {

  guard_against_notgiotto(gobject)

  reduction = match.arg(reduction, choices = c('cells', 'feats'))

  # to be deprecated ('data.table' -> 'matrix')
  if(!identical(output, c('dimObj', 'matrix'))) {
    if(output == 'data.table') output = 'matrix'
  }
  output = match.arg(output, choices = c('dimObj', 'matrix'))

  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  data_list = slot(gobject, 'dimension_reduction')[[reduction]][[spat_unit]][[feat_type]]

  data_list = unlist(data_list, recursive = TRUE, use.names = FALSE)
  data_list = assign_objnames_2_list(data_list)

  if(output == 'dimObj') {
    return(data_list)
  }
  if(output == 'matrix') {
    return(lapply(data_list, `[`))
  }
}







#' @title Set dimension reduction data
#' @name setDimReduction
#' @description Function to dimension reduction information into the Giotto
#' object.
#' @inheritParams data_access_params
#' @param x dimObj or list of dimObj to set. Passing NULL will remove a specified
#' set of dimension reduction information from the gobject
#' @param name name of reduction results
#' @param reduction reduction on cells or features
#' @param reduction_method reduction method (e.g. "pca")
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @return giotto object
#' @keywords autocomplete
#' @family dimensional reduction data accessor functions
#' @family functions to set data in giotto object
#' @export
setDimReduction = function(gobject,
                           x,
                           spat_unit = NULL,
                           feat_type = NULL,
                           name = 'pca',
                           reduction = c('cells', 'feats'),
                           reduction_method = c('pca', 'umap', 'tsne'),
                           provenance = NULL,
                           verbose = TRUE,
                           initialize = TRUE) {

  guard_against_notgiotto(gobject)
  reduction = match.arg(reduction, choices = c('cells', 'feats'))
  # reduction_method = match.arg(reduction_method, choices = c('pca', 'umap', 'tsne'))
  if(!methods::hasArg(x)) stop(wrap_txt('x (data to set) param must be given'))


  # check hierarchical slots
  avail_ex = list_expression(gobject)
  if(is.null(avail_ex)) stop(wrap_txt('Add expression information first'))



  # 1. Determine user inputs
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
  nospec_red = ifelse(is.null(match.call()$reduction), yes = TRUE, no = FALSE)
  nospec_red_method = ifelse(is.null(match.call()$reduction_method), yes = TRUE, no = FALSE)
  .external_accessor_dimred = list(nospec_unit = nospec_unit,
                                   nospec_feat = nospec_feat,
                                   nospec_name = nospec_name,
                                   nospec_red = nospec_red,
                                   nospec_red_method = nospec_red_method)
  # checked by internal setter to determine if called by external


  # NATIVE INPUT TYPES
  # 3. If input is dimObj or NULL, pass to internal
  if(is.null(x) | inherits(x, 'dimObj')) {

    # pass to internal
    gobject = set_dimReduction(
      gobject = gobject,
      dimObject = x,
      reduction = reduction,
      spat_unit = spat_unit,
      feat_type = feat_type,
      reduction_method = reduction_method,
      name = name,
      provenance = provenance,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)

  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'dimObj')) {

      # MULTIPLE INPUT
      # 4. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_dimReduction(
          gobject = gobject,
          dimObject = x[[obj_i]],
          reduction = reduction,
          spat_unit = spat_unit,
          feat_type = feat_type,
          reduction_method = reduction_method,
          name = name,
          provenance = provenance,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)
    }
  }

  # catch
  stop(wrap_txt('Only dimObj or lists of dimObj accepted.
                  For raw or external data, please first use readDimReducData()'))

}







#' @title Set dimension reduction
#' @name set_dimReduction
#' @description Function to set a dimension reduction slot
#' @inheritParams data_access_params
#' @param reduction reduction on cells or features
#' @param reduction_method reduction method (e.g. "pca")
#' @param name name of reduction results
#' @param dimObject dimension object result to set
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @return giotto object
#' @family dimensional reduction data accessor functions
#' @family functions to set data in giotto object
#' @export
set_dimReduction = function(gobject,
                            dimObject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            reduction = c('cells', 'feats'),
                            reduction_method = c('pca', 'umap', 'tsne'),
                            name = 'pca',
                            provenance = NULL,
                            verbose = TRUE,
                            set_defaults = TRUE,
                            initialize = FALSE) {

  guard_against_notgiotto(gobject)
  reduction = match.arg(reduction, choices = c('cells', 'feats'))
  # reduction_method = match.arg(reduction_method, choices = c('pca', 'umap', 'tsne'))
  if(!methods::hasArg(dimObject)) stop(wrap_txt('dimObject param must be given'))

  # 0. pass to external if not native format
  if(!inherits(dimObject, c('dimObj', 'NULL'))) {
    stop(wrap_txt(deparse(substitute(dimObject)), 'is not dimObj (set)
                  or NULL (remove)'))
  }

  # 1. Determine user inputs
  p = parent.frame() # Get values if called from external
  call_from_external = exists('.external_accessor_dimred', where = p)

  if(call_from_external) {
    nospec_unit = p$.external_accessor_dimred$nospec_unit
    nospec_feat = p$.external_accessor_dimred$nospec_feat
    nospec_name = p$.external_accessor_dimred$nospec_name
    nospec_red = p$.external_accessor_dimred$nospec_red
    nospec_red_method = p$.external_accessor_dimred$nospec_red_method
  } else {
    nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
    nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
    nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
    nospec_red = ifelse(is.null(match.call()$reduction), yes = TRUE, no = FALSE)
    nospec_red_method = ifelse(is.null(match.call()$reduction_method), yes = TRUE, no = FALSE)
  }


  # 2. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 3. if input is NULL, remove object
  if(is.null(dimObject)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to dimObject param.
                                 Removing specified expression')
    slot(gobject, 'dimension_reduction')[[reduction]][[spat_unit]][[feat_type]][[reduction_method]][[name]] = NULL

    # prune if empty
    if(length(gobject@dimension_reduction[[reduction]][[spat_unit]][[feat_type]][[reduction_method]]) == 0L) {
      gobject@dimension_reduction = NULL
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }

  # 4. import data from S4 if available
  # NOTE: modifies spat_unit/feat_type/name/provenance/reduction/reduction_method/data slots
  dimObject = read_s4_nesting(dimObject)


  ## 5. check if specified name has already been used
  potential_names = list_dim_reductions_names(gobject,
                                              spat_unit = spat_unit,
                                              feat_type = feat_type,
                                              data_type = reduction,
                                              dim_type = reduction_method)
  if(name %in% potential_names) {
    if(isTRUE(verbose)) wrap_msg('> ', name, ' already exists and will be replaced with new dimension reduction object \n')
  }



  ## 6. update and return giotto object
  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting dimension reduction [', spatUnit(dimObject), '][', featType(dimObject), '] ',
    objName(dimObject), sep = ''
  )

  slot(gobject, 'dimension_reduction')[[reduction]][[spat_unit]][[feat_type]][[reduction_method]][[name]] = dimObject
  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}






## nearest neighbor network slot ####

#' @title Get nearest network
#' @name get_NearestNetwork
#' @description Get a NN-network from a Giotto object
#' @inheritParams data_access_params
#' @param nn_network_to_use "kNN" or "sNN"
#' @param network_name name of NN network to be used
#' @param output return a igraph or data.table object. Default 'igraph'
#' @return igraph or data.table object
#' @family expression space nearest network accessor functions
#' @family functions to get data from giotto object
#' @export
get_NearestNetwork = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              nn_network_to_use = NULL,
                              network_name = NULL,
                              output = c('nnNetObj', 'igraph', 'data.table'),
                              set_defaults = TRUE) {


  output = match.arg(arg = output, choices = c('nnNetObj', 'igraph', 'data.table'))

  # 1.  Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 2 Find the object
  if(is.null(nn_network_to_use)) {
    nn_network_to_use = names(slot(gobject, 'nn_network')[[spat_unit]][[feat_type]])[[1]]
    if(is.null(nn_network_to_use)) {

      stop(wrap_txt('There is currently no nearest-neighbor network created for
                     spatial unit: "', spat_unit, '" and feature type "', feat_type, '".
                     First run createNearestNetwork()\n', sep = ''))
    } else {
      wrap_msg('The NN network type was not specified, default to the first: "',
               nn_network_to_use,'"', sep = '')
    }
  }

  if(is.null(network_name)) {
    network_name = names(slot(gobject, 'nn_network')[[spat_unit]][[feat_type]][[nn_network_to_use]])[[1]]
    if(is.null(network_name)) {
      stop(wrap_txt('There is currently no nearest-neighbor network built for spatial unit: "', spat_unit,
                    '" feature type: "', feat_type, '" and network type: "', nn_network_to_use,'"\n',
                    sep = ''))
    }else {
      wrap_msg('The NN network name was not specified, default to the first: "',
               network_name,'"', sep = '')
    }
  }

  # 3. get object in desired format

  nnNet = slot(gobject, 'nn_network')[[spat_unit]][[feat_type]][[nn_network_to_use]][[network_name]]
  if(is.null(nnNet)) {
    stop(wrap_txt('nn_network_to_use: "', nn_network_to_use, '" or network_name: "', network_name, '" does not exist.
                  Create a nearest-neighbor network first', sep = ''))
  }

  if(output == 'nnNetObj') {

    return(nnNet) # return nnNetObj

  } else if(output == 'igraph' | output == 'data.table') {
    nnNet = slot(nnNet, 'igraph')

    if(output == 'igraph') return(nnNet) # return igraph
    if(output == 'data.table') {
      nnNet = data.table::setDT(igraph::get.data.frame(x = nnNet))
      return(nnNet) # return data.table
    }
  }

}




#' @title Get nearest neighbor network
#' @name getNearestNetwork
#' @description Get a NN-network from a Giotto object
#' @inheritParams data_access_params
#' @param nn_type "kNN" or "sNN"
#' @param name name of NN network to be used
#' @param output return a igraph or data.table object. Default 'igraph'
#' @return igraph or data.table object
#' @family expression space nearest network accessor functions
#' @family functions to get data from giotto object
#' @export
getNearestNetwork = function(gobject,
                             spat_unit = NULL,
                             feat_type = NULL,
                             nn_type = NULL,
                             name = NULL,
                             output = c('nnNetObj', 'igraph', 'data.table'),
                             set_defaults = TRUE) {

  # pass to internal
  nn = get_NearestNetwork(gobject = gobject,
                          spat_unit = spat_unit,
                          feat_type = feat_type,
                          nn_network_to_use = nn_type,
                          network_name = name,
                          output = output,
                          set_defaults = set_defaults)

  return(nn)
}






#' @description Get all nearest neighbor networks for a specified spatial unit
#' and feature type
#' @keywords internal
#' @return list of nnNetObj, igraph, or data.table depending on output param
#' @noRd
get_nearest_network_list = function(gobject,
                                    spat_unit = NULL,
                                    feat_type = NULL,
                                    output = c('nnNetObj', 'igraph', 'data.table'),
                                    set_defaults = TRUE) {

  guard_against_notgiotto(gobject)

  output = match.arg(output, choices = c('nnNetObj', 'igraph', 'data.table'))

  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  data_list = slot(gobject, 'nn_network')[[spat_unit]][[feat_type]]

  data_list = unlist(data_list, recursive = TRUE, use.names = FALSE)
  data_list = assign_objnames_2_list(data_list)

  if(output == 'nnNetObj') {
    return(data_list)
  }
  if(output == 'igraph') {
    return(lapply(data_list, `[`))
  }
  if(output == 'data.table') {
    return(lapply(data_list, function(obj_i) {
      data.table::setDT(igraph::get.data.frame(x = obj_i[]))
    }))
  }
}








#' @title Set nearest neighbor network
#' @name setNearestNetwork
#' @description Set a NN-network for a Giotto object
#' @inheritParams data_access_params
#' @param x nnNetObj or list of nnNetObj. Passing NULL will remove a specified
#' set of nearest neighbor network information from the gobject
#' @param nn_type "kNN" or "sNN"
#' @param name name of NN network to be used
#' yet supported.
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @return giotto object
#' @family expression space nearest network accessor functions
#' @family functions to set data in giotto object
#' @export
setNearestNetwork = function(gobject,
                             x,
                             spat_unit = NULL,
                             feat_type = NULL,
                             nn_type = 'sNN',
                             name = 'sNN.pca',
                             provenance = NULL,
                             verbose = TRUE,
                             initialize = TRUE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x (data to set) param must be given'))


  # check hierarchical slots
  avail_dr = list_dim_reductions(gobject)
  if(is.null(avail_dr))
    stop(wrap_txt('Add dimension reduction information first'))


  # 1. Determine user inputs
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  nospec_net = ifelse(is.null(match.call()$nn_type), yes = TRUE, no = FALSE)
  nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
  .external_accessor_nn = list(nospec_unit = nospec_unit,
                               nospec_feat = nospec_feat,
                               nospec_net = nospec_net,
                               nospec_name = nospec_name)
  # checked by internal setter to determine if called by external


  # NATIVE INPUT TYPES
  # 3. If input is nnNetObj or NULL, pass to internal
  if(is.null(x) | inherits(x, 'nnNetObj')) {

    # pass to internal
    gobject = set_NearestNetwork(
      gobject = gobject,
      nn_network = x,
      spat_unit = spat_unit,
      feat_type = feat_type,
      nn_network_to_use = nn_type,
      network_name = name,
      provenance = provenance,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)

  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'nnNetObj')) {

      # MULTIPLE INPUT
      # 4. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_NearestNetwork(
          gobject = gobject,
          nn_network = x[[obj_i]],
          spat_unit = spat_unit,
          feat_type = feat_type,
          nn_network_to_use = nn_type,
          network_name = name,
          provenance = provenance,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)

    }
  }
  # catch
  stop(wrap_txt('Only nnNetObj or lists of nnNetObj accepted.
                For raw or external data, please first use readNearestNetData()'))

}








#' @title Set nearest network
#' @name set_NearestNetwork
#' @description Set a NN-network for a Giotto object
#' @inheritParams data_access_params
#' @param nn_network_to_use "kNN" or "sNN"
#' @param network_name name of NN network to be used
#' @param nn_network nnNetObj or igraph nearest network object. Data.table not
#' yet supported.
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @return giotto object
#' @family expression space nearest network accessor functions
#' @family functions to set data in giotto object
#' @export
set_NearestNetwork = function(gobject,
                              nn_network,
                              spat_unit = NULL,
                              feat_type = NULL,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              provenance = NULL,
                              verbose = TRUE,
                              set_defaults = TRUE,
                              initialize = FALSE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(nn_network)) stop(wrap_txt('nn_network param must be given'))

  # 0. stop if not native formats
  if(!inherits(nn_network, c('nnNetObj', 'NULL'))) {
    stop(wrap_txt(deparse(substitute(nn_network)), 'is not nnNetObj (set).
                  or NULL (remove)'))
  }


  # 1. determine user input
  p = parent.frame() # Get values if called from external
  call_from_external = exists('.external_accessor_nn', where = p)

  if(isTRUE(call_from_external)) {
    nospec_unit = p$.external_accessor_nn$nospec_unit
    nospec_feat = p$.external_accessor_nn$nospec_feat
    nospec_net = p$.external_accessor_nn$nospec_net
    nospec_name = p$.external_accessor_nn$nospec_name
  } else {
    nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
    nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
    nospec_net = ifelse(is.null(match.call()$nn_network_to_use), yes = TRUE, no = FALSE)
    nospec_name = ifelse(is.null(match.call()$network_name), yes = TRUE, no = FALSE)
  }

  # change var name for use with read_s4_nesting()
  name = network_name
  nn_type = nn_network_to_use


  # 2. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 3. If input is null, remove object
  if(is.null(nn_network)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to nn_network.
                                 Removing specified nearest neighbor network.')
    gobject@nn_network[[spat_unit]][[feat_type]][[nn_type]][[name]] = NULL

    # prune if empty
    if(length(gobject@nn_network[[spat_unit]][[feat_type]][[nn_type]]) == 0L) {
      gobject@nn_network = NULL
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    return(gobject)
  }


  # 4. import data from S4 if available
  # NOTE: modifies spat_unit/feat_type/name/nn_type/nn_network/provenance
  nn_network = read_s4_nesting(nn_network)
  # Use updated values in nn_type and name instead of network_name and nn_network_to_use


  ## 5. check if specified name has already been used
  potential_names = list_nearest_networks_names(gobject,
                                                spat_unit = spat_unit,
                                                feat_type = feat_type,
                                                nn_type = nn_type)
  if(name %in% potential_names) {
    if(isTRUE(verbose)) wrap_msg('> "', name, '" already exists and will be replaced with new nearest neighbor network')
  }

  ## 6. update and return giotto object
  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting nearest neighbor network [', spatUnit(nn_network), '][', featType(nn_network), '] ',
    objName(nn_network), sep = ''
  )

  gobject@nn_network[[spat_unit]][[feat_type]][[nn_type]][[name]] = nn_network
  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}








## spatial network slot ####

#' @title Get spatial network
#' @name get_spatialNetwork
#' @description Function to get a spatial network
#' @inheritParams data_access_params
#' @param name name of spatial network
#' @param output object type to return as. Options: 'spatialNetworkObj' (default),
#' 'networkDT' and 'networkDT_before_filter' for data.table outputs.
#' @param copy_obj whether to copy/duplicate when getting the object (default = TRUE)
#' @param verbose be verbose
#' @family spatial network data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatialNetwork = function(gobject,
                              spat_unit = NULL,
                              name = NULL,
                              output = c('spatialNetworkObj',
                                         'networkDT',
                                         'networkDT_before_filter',
                                         'outputObj'),
                              set_defaults = TRUE,
                              copy_obj = TRUE,
                              verbose = TRUE) {


  output = match.arg(output, choices = c('spatialNetworkObj',
                                         'networkDT',
                                         'networkDT_before_filter',
                                         'outputObj'))


  # Set spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
  }

  # set name to get if missing
  if(is.null(name)) name = names(slot(gobject, 'spatial_network')[[spat_unit]])[[1]]

  # check if given name is present
  if (!is.element(name, names(slot(gobject, 'spatial_network')[[spat_unit]]))){
    if(isTRUE(verbose)) msg = wrap_txt('spatial network', name, 'has not been created. Returning NULL.
                                       check which spatial networks exist with showGiottoSpatNetworks()')
    warning(msg)
    return(NULL)
  }else{
    networkObj = slot(gobject, 'spatial_network')[[spat_unit]][[name]]

    # S3 backwards compatibility
    if(!isS4(networkObj)) networkObj = S3toS4spatNetObj(networkObj)
    silent = validObject(networkObj) # Variable used to hide TRUE print

  }

  if(copy_obj) {
    networkObj@networkDT = data.table::copy(networkObj@networkDT)
    if(!is.null(networkObj@networkDT_before_filter)) {
      networkObj@networkDT_before_filter = data.table::copy(networkObj@networkDT_before_filter)
    }
  }

  if (output == 'spatialNetworkObj'){
    return(networkObj)
  } else if(output == 'networkDT'){
    return(slot(networkObj, 'networkDT'))
  } else if(output == 'networkDT_before_filter') {
    return(slot(networkObj, 'networkDT_before_filter'))
  } else if(output == 'outputObj') {
    return(slot(networkObj, 'outputObj'))
  }
}




#' @title Get spatial network
#' @name getSpatialNetwork
#' @description Function to get a spatial network
#' @inheritParams data_access_params
#' @param name name of spatial network
#' @param output object type to return as. Options: 'spatialNetworkObj' (default),
#' 'networkDT' and 'networkDT_before_filter' for data.table outputs.
#' @param copy_obj whether to copy/duplicate when getting the object (default = TRUE)
#' @param verbose be verbose
#' @family spatial network data accessor functions
#' @family functions to get data from giotto object
#' @export
getSpatialNetwork = function(gobject,
                             spat_unit = NULL,
                             name = NULL,
                             output = c('spatialNetworkObj',
                                        'networkDT',
                                        'networkDT_before_filter',
                                        'outputObj'),
                             set_defaults = TRUE,
                             copy_obj = TRUE,
                             verbose = TRUE) {

  # Pass to internal function
  network = get_spatialNetwork(gobject = gobject,
                               spat_unit = spat_unit,
                               name = name,
                               output = output,
                               set_defaults = set_defaults,
                               copy_obj = copy_obj,
                               verbose = verbose)

  return(network)
}






#' @description Get all spatial networks for a specified spatial unit
#' @keywords internal
#' @return list of dimObj or data.table depending on output param
#' @noRd
get_spatial_network_list = function(gobject,
                                    spat_unit = NULL,
                                    output = c('spatialNetworkObj',
                                               'networkDT',
                                               'networkDT_before_filter',
                                               'outputObj'),
                                    set_defaults = TRUE,
                                    copy_obj = TRUE) {

  guard_against_notgiotto(gobject)

  output = match.arg(output, choices = c('spatialNetworkObj',
                                         'networkDT',
                                         'networkDT_before_filter',
                                         'outputObj'))

  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
  }

  data_list = slot(gobject, 'spatial_network')[[spat_unit]]

  data_list = unlist(data_list, recursive = TRUE, use.names = FALSE)
  data_list = assign_objnames_2_list(data_list)

  # copy object
  if(isTRUE(copy_obj)) data_list = lapply(data_list, copy)


  # return object list
  if(output == 'spatialNetworkObj') {
    return(data_list)
  }
  if(output == 'networkDT') {
    return(lapply(data_list, `[`))
  }
  if(output == 'networkDT_before_filter') {
    return(lapply(data_list, slot, 'networkDT_before_filter'))
  }
  if(output == 'outputObj') {
    return(lapply(data_list, slot, 'outputObj'))
  }
}





#' @title Set spatial network
#' @name setSpatialNetwork
#' @description Function to set a spatial network
#' @inheritParams data_access_params
#' @param x spatialNetworkObj or list of spatialNetworkObj to set. Passing NULL
#' removes a specified set of spatial network information from the gobject.
#' @param name name of spatial network
#' @param provenance provenance name
#' @param verbose be verbose
#' @return giotto object
#' @family spatial network data accessor functions
#' @family functions to set data in giotto object
#' @export
setSpatialNetwork = function(gobject,
                             x,
                             spat_unit = NULL,
                             name = NULL,
                             provenance = NULL,
                             verbose = TRUE,
                             initialize = TRUE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x param (data to set) must be given'))


  # check hierarchical slots
  avail_ex = list_expression(gobject)
  avail_sl = list_spatial_locations(gobject)
  if(is.null(avail_ex)) stop(wrap_txt('Add expression and spatial location information first'))
  if(is.null(avail_sl)) stop(wrap_txt('Add spatial location information first'))


  # 1. Determine user inputs
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
  .external_accessor_sn = list(nospec_unit = nospec_unit,
                               nospec_name = nospec_name)
  # checked by internal setter to determine if called by external


  # NATIVE INPUT TYPES
  # 3. If input is spatialNetworkObj or NULL, pass to internal
  if(is.null(x) | inherits(x, 'spatialNetworkObj')) {

    # pass to internal
    gobject = set_spatialNetwork(
      gobject = gobject,
      spatial_network = x,
      spat_unit = spat_unit,
      name = name,
      provenance = provenance,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)

  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'spatialNetworkObj')) {

      # MULTIPLE INPUT
      # 4. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_spatialNetwork(
          gobject = gobject,
          spatial_network = x[[obj_i]],
          spat_unit = spat_unit,
          name = name,
          provenance = provenance,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)
    }
  }

  # catch
  stop(wrap_txt('Only spatialNetworkObj or lists of spatialNetworkObj accepted.
                For raw or external data, please first use readSpatNetData()'))

}







#' @title Set spatial network
#' @name set_spatialNetwork
#' @description Function to set a spatial network
#' @inheritParams data_access_params
#' @param name name of spatial network
#' @param provenance provenance name
#' @param spatial_network spatial network
#' @param verbose be verbose
#' @return giotto object
#' @family spatial network data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatialNetwork = function(gobject,
                              spatial_network,
                              spat_unit = NULL,
                              name = NULL,
                              provenance = NULL,
                              verbose = TRUE,
                              set_defaults = TRUE,
                              initialize = FALSE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(spatial_network)) stop(wrap_txt('spatial_network param must be given'))

  # 0. stop if not native formats
  if(!inherits(spatial_network, c('spatialNetworkObj', 'NULL'))) {
    stop(wrap_txt(deparse(substitute(spatial_network)), 'is not spatialNetworkObj (set).
                  or NULL (remove)'))
  }

  # 1. determine if input was supplied to spat_unit and name
  p = parent.frame() # Get values if called from external
  call_from_external = exists('.external_accessor_sn', where = p)

  if(isTRUE(call_from_external)) {
    nospec_unit = p$.external_accessor_sn$nospec_unit
    nospec_name = p$.external_accessor_sn$nospec_name
  } else {
    nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
    nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
  }


  # 2. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
  }

  # 3. If input is NULL, remove object
  if(is.null(spatial_network)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to spatial_network.
                                 Removing specified spatial network.')
    gobject@spatial_network[[spat_unit]][[name]] = NULL

    # prune if empty
    if(length(gobject@spatial_network[[spat_unit]]) == 0L) {
      gobject@spatial_network = NULL
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }


  # 4. import data from S4 if available
  # NOTE: modifies spat_unit/name/provenance/spatial_network
  spatial_network = read_s4_nesting(spatial_network)


  # 5. check if specified name has already been used
  if(isTRUE(verbose)) {
    potential_names = list_spatial_networks_names(gobject = gobject,
                                                  spat_unit = spat_unit)
    if(name %in% potential_names) {
      wrap_msg('> "', name, '" already exists and will be replaced with new spatial network')
    }
  }


  # 6. update and return giotto object
  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting spatial network [', spatUnit(spatial_network), '] ',
    objName(spatial_network), sep = ''
  )

  slot(gobject, 'spatial_network')[[spat_unit]][[name]] = spatial_network
  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}







## spatial grid slot ####

#' @title Get spatial grid
#' @name get_spatialGrid
#' @description Function to get spatial grid
#' @inheritParams data_access_params
#' @param name name of spatial grid
#' @param return_grid_Obj return grid object (default = FALSE)
#' @family spatial grid data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatialGrid = function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL,
                           name = NULL,
                           return_grid_Obj = FALSE,
                           set_defaults = TRUE) {


  # Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # **To be deprecated** - check for old nesting
  if(is.null(names(gobject@spatial_grid[[spat_unit]][[feat_type]])[[1]])) { # If gobject has nothing for this feat_type
    available = list_spatial_grids(gobject,
                                   spat_unit = spat_unit)
    if(nrow(available > 0 & is.null(available$feat_type))) { # If ANY old nesting objects are discovered (only reports old nestings if any detected)
      if(is.null(name)) gridObj = gobject@spatial_grid[[spat_unit]][[available$name[[1]]]]
      else gridObj = gobject@spatial_grid[[spat_unit]][[name]]
      if(inherits(gridObj, 'spatialGridObj')) {
        if(is.null(name)) message('The grid name was not specified, default to the first: "', available$name[[1]],'"')
        # S3 backwards compatibility
        if(!isS4(gridObj)) gridObj = S3toS4spatialGridObj(gridObj)
        silent = validObject(gridObj) # variable used to hide TRUE print

        gridDT = slot(gridObj, 'gridDT')
        if (return_grid_Obj == TRUE){
          return(gridObj)
        }else{
          return(gridDT)
        }
      } else {
        stop('There is currently no spatial grid created for
             spatial unit: "', spat_unit,'" and feature type "', feat_type,'".
             First run createSpatialGrid()')
      }
    }
  } # ** deprecation end**

  # Automatically select first grid for given spat_unit/feat_type combination
  if(is.null(name)) {
    name = names(slot(gobject, 'spatial_grid')[[spat_unit]][[feat_type]])[[1]]
    message('The grid name was not specified, default to the first: "', name, '"')
  } else if (!is.element(name, names(slot(gobject, 'spatial_grid')[[spat_unit]][[feat_type]]))){

    message = sprintf("spatial grid %s has not been created. Returning NULL.
                      check which spatial grids exist with showGiottoSpatGrids()\n", name)
    warning(message)
    return(NULL)
  }

  # Get spatialGridObj
  gridObj = slot(gobject, 'spatial_grid')[[spat_unit]][[feat_type]][[name]]

  # S3 backwards compatibility
  if(!isS4(gridObj)) gridObj = S3toS4spatialGridObj(gridObj)
  silent = validObject(gridObj) # variable used to hide TRUE print

  gridDT = slot(gridObj, 'gridDT')

  if (return_grid_Obj == TRUE){
    return(gridObj)
  }else{
    return(gridDT)
  }
}



#' @title Get spatial grid
#' @name getSpatialGrid
#' @description Function to get spatial grid
#' @inheritParams data_access_params
#' @param name name of spatial grid
#' @param return_grid_Obj return grid object (default = FALSE)
#' @family spatial grid data accessor functions
#' @family functions to get data from giotto object
#' @export
getSpatialGrid = function(gobject,
                          spat_unit = NULL,
                          feat_type = NULL,
                          name = NULL,
                          return_grid_Obj = FALSE,
                          set_defaults = TRUE) {

  # Pass to internal function
  grid = get_spatialGrid(gobject = gobject,
                         spat_unit = spat_unit,
                         feat_type = feat_type,
                         name = name,
                         return_grid_Obj = return_grid_Obj,
                         set_defaults = set_defaults)

  return(grid)
}

#' @title Set spatial grid
#' @name set_spatialGrid
#' @description Function to set a spatial grid
#' @inheritParams data_access_params
#' @param spatial_grid spatial grid object
#' @param name name of spatial grid
#' @param verbose be verbose
#' @return giotto object
#' @family spatial grid data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatialGrid = function(gobject,
                           spatial_grid,
                           spat_unit = NULL,
                           feat_type = NULL,
                           name = NULL,
                           verbose = TRUE,
                           set_defaults = TRUE) {


  # 1. check input
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  nospec_name = ifelse(is.null(name), yes = TRUE, no = FALSE)

  # 2. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 3. if input is null, remove object
  if(is.null(spatial_grid)) {
    if(isTRUE(verbose)) message('NULL passed to metadata.\n Removing specified metadata.')
    gobject@spatial_grid[[spat_unit]][[feat_type]][[name]] = NULL
  }

  # 4. import information from S4 if possible
  if(inherits(spatial_grid, 'spatialGridObj')) {

    spatial_grid = read_s4_nesting(spatial_grid)

    # if(isTRUE(nospec_unit)) {
    #   if(!is.na(slot(spatial_grid, 'spat_unit'))) spat_unit = slot(spatial_grid, 'spat_unit')
    #   else slot(spatial_grid, 'spat_unit') = spat_unit
    # } else {
    #   slot(spatial_grid, 'spat_unit') = spat_unit
    # }
    # if(isTRUE(nospec_feat)) {
    #   if(!is.na(slot(spatial_grid, 'feat_type'))) feat_type = slot(spatial_grid, 'feat_type')
    #   else slot(spatial_grid, 'feat_type') = feat_type
    # } else {
    #   slot(spatial_grid, 'feat_type') = feat_type
    # }
    # if(isTRUE(nospec_name)) {
    #   if(!is.na(slot(spatial_grid, 'name'))) name = slot(spatial_grid, 'name')
    #   else slot(spatial_grid, 'name') = name
    # } else {
    #   slot(spatial_grid, 'name') = name
    # }

  } else {
    stop('spatial_grid must be a spatialGridObj')
  }

  ## 5. check if specified name has already been used
  if(isTRUE(verbose)) {
    potential_names = names(slot(gobject, 'spatial_grid')[[spat_unit]][[feat_type]])
    if(name %in% potential_names) {
      wrap_msg('> "', name, '" already exists and will be replaced with new spatial grid \n')
    }
  }


  ## TODO: 2. check input for spatial grid
  if(!inherits(spatial_grid, 'spatialGridObj')) stop('spatial_grid to set must be S4 "spatialGridObj"\n')
  silent = validObject(spatial_grid) # Variable only used to hide TRUE prints

  ## 6. update and return giotto object
  slot(gobject, 'spatial_grid')[[spat_unit]][[feat_type]][[name]] = spatial_grid

  return(gobject)

}

#' @title Set spatial grid
#' @name setSpatialGrid
#' @description Function to set a spatial grid
#' @inheritParams data_access_params
#' @param spatial_grid spatial grid object
#' @param name name of spatial grid
#' @param verbose be verbose
#' @return giotto object
#' @family spatial grid data accessor functions
#' @family functions to set data in giotto object
#' @export
setSpatialGrid = function(gobject,
                           spatial_grid,
                           spat_unit = NULL,
                           feat_type = NULL,
                           name = NULL,
                           verbose = TRUE,
                           set_defaults = TRUE) {

  # Pass to internal function
  gobject = set_spatialGrid(gobject = gobject,
                            spatial_grid = spatial_grid,
                            spat_unit = spat_unit,
                            feat_type = feat_type,
                            name = name,
                            verbose = verbose,
                            set_defaults = set_defaults)

  return(gobject)
}

## polygon cell info ####

#' @title Get polygon info
#' @name get_polygon_info
#' @description Get giotto polygon spatVector
#' @param gobject giotto object
#' @param polygon_name name of polygons. Default "cell"
#' @param polygon_overlap include polygon overlap information
#' @param return_giottoPolygon (Defaults to FALSE) Return as giottoPolygon S4 object
#' @family polygon info data accessor functions
#' @family functions to get data from giotto object
#' @export
get_polygon_info = function(gobject,
                            polygon_name = NULL,
                            polygon_overlap = NULL,
                            return_giottoPolygon = FALSE) {


  potential_names = names(slot(gobject, 'spatial_info'))
  if(is.null(potential_names)) stop('Giotto object contains no polygon information')

  # If polygon_name is not given...
  if(is.null(polygon_name)) {
    if('cell' %in% potential_names) {
      polygon_name = 'cell' # Default to 'cell' as polygon_name if available
    } else {
      polygon_name = potential_names[1] # Select 1st available name if 'cell' is missing
      message('No polygon information named "cell" discovered.\n selecting first available ("',polygon_name,'")')
    }
  }

  if(isTRUE(return_giottoPolygon)) return(slot(gobject, 'spatial_info')[[polygon_name]])

  if(!polygon_name %in% potential_names) {
    stop('There is no polygon information with name ', polygon_name, '\n')
  } else {

    if(is.null(polygon_overlap)) {
      poly_info = gobject@spatial_info[[polygon_name]]@spatVector
    } else {
      potential_overlaps = names(gobject@spatial_info[[polygon_name]]@overlaps)

      if(!polygon_overlap %in% potential_overlaps) {
        stop('There is no polygon overlap information with name ', polygon_overlap, '\n')
      } else {
        poly_info = gobject@spatial_info[[polygon_name]]@overlaps[[polygon_overlap]]
      }
    }
    return(poly_info)
  }
}






#' @title Get polygon info
#' @name getPolygonInfo
#' @description Get giotto polygon spatVector
#' @param gobject giotto object
#' @param polygon_name name of polygons. Default is "cell"
#' @param polygon_overlap include polygon overlap information
#' @param return_giottoPolygon (Defaults to FALSE) Return as giottoPolygon S4 object
#' @family polygon info data accessor functions
#' @family functions to get data from giotto object
#' @export
getPolygonInfo = function(gobject = NULL,
                          polygon_name = NULL,
                          polygon_overlap = NULL,
                          return_giottoPolygon = FALSE) {
  if (!inherits(gobject, 'giotto')){
    wrap_msg("Unable to get polygon spatVector from non-Giotto object.")
    stop(wrap_txt("Please provide a Giotto object to the gobject argument.",
                  errWidth = TRUE))
  }

  poly_info = get_polygon_info(gobject = gobject,
                               polygon_name = polygon_name,
                               polygon_overlap = polygon_overlap,
                               return_giottoPolygon = return_giottoPolygon)

  return (poly_info)
}





#' @description Get list of all polygon info
#' @keywords internal
#' @return list of giottoPolygon or SpatVector depending on return_giottoPolygon
#' param
#' @noRd
get_polygon_info_list = function(gobject,
                                 return_giottoPolygon = TRUE) {

  guard_against_notgiotto(gobject)

  data_list = slot(gobject, 'spatial_info')

  # return objects
  if(isTRUE(return_giottoPolygon)) {
    return(data_list)
  } else {
    return(lapply(data_list, `[`))
  }

}






#' @title Set polygon info
#' @name setPolygonInfo
#' @description Set polygon information into Giotto object
#' @inheritParams data_access_params
#' @param x single object or named list of objects to set as polygon
#' information (see details)
#' @param name (optional, character) name to assign to polygon and spatial unit
#' that polygon might define. Only used for single giottoPolygon objects. Names
#' are taken from a named list for multiple polygons.
#' @param centroids_to_spatlocs if centroid information is discovered, whether
#' to additionally set them as a set of spatial locations (default = FALSE)
#' @return giotto object
#' @details Inputs can be provided as either single objects or named lists of
#' objects. If the list is not named, then a generic name of the template
#' 'cell_i' will be applied. \cr
#' If an input is a character string, then it is assumed that it is a
#' filepath. \cr
#' For required formatting when reading tabular data or objects, see
#' \code{\link{createGiottoPolygonsFromDfr}} details.
#' @family polygon info data accessor functions
#' @family functions to set data in giotto object
#' @export
setPolygonInfo = function(gobject,
                          x,
                          name = 'cell',
                          centroids_to_spatlocs = FALSE,
                          verbose = TRUE,
                          initialize = TRUE) {

  # data.table vars
  poly_ID = y = NULL

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x param (data to be set) must be given'))


  # 1. determine user inputs
  nospec_name = !methods::hasArg(name)
  .external_accessor_poly = list(nospec_name = nospec_name)
  # checked by internal setter to determine if called by external



  # NATIVE INPUT TYPES
  # 2. If input is giottoPolygon or NULL, pass to internal
  if(inherits(x, c('giottoPolygon', 'NULL'))) {
    # pass to internal
    gobject = set_polygon_info(
      gobject = gobject,
      polygon_name = name,
      gpolygon = x,
      verbose = verbose,
      initialize = !isTRUE(centroids_to_spatlocs) & initialize # delay so centroids can be added
    )

    # Attach centroids if found
    if(inherits(x, 'giottoPolygon') & isTRUE(centroids_to_spatlocs)) {
      if(!is.null(x@spatVectorCentroids)) {

        centroids = x@spatVectorCentroids
        centroidsDT = spatVector_to_dt(centroids)
        centroidsDT_loc = centroidsDT[, .(poly_ID, x, y)]
        colnames(centroidsDT_loc) = c('cell_ID', 'sdimx', 'sdimy')

        locsObj = create_spat_locs_obj(name = 'raw',
                                       coordinates = centroidsDT_loc,
                                       spat_unit = x@name, # tag same spat_unit as poly
                                       provenance = x@name,
                                       misc = NULL)

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        .external_accessor_spatloc = list(
          # set spatlocs 'spat_unit' using 'name' if it was EXPLICITLY
          # supplied to setPolygonInfo
          # otherwise, set spatlocs 'spat_unit' as x@name
          nospec_unit = nospec_name,
          # set spatlocs name based on locsObj@name
          nospec_name = TRUE
        )
        gobject = set_spatial_locations(gobject,
                                        spatlocs = locsObj,
                                        spat_unit = name, # useif explicit here
                                        verbose = verbose,
                                        set_defaults = FALSE,
                                        initialize = initialize)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      }
    }
    return(gobject)

  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'giottoPolygon')) {

      # 3. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_polygon_info(
          gobject = gobject,
          gpolygon = x[[obj_i]],
          polygon_name = name,
          verbose = verbose,
          initialize = !isTRUE(centroids_to_spatlocs) & initialize
        )

        # Attach centroids if found
        if(inherits(x[[obj_i]], 'giottoPolygon') & isTRUE(centroids_to_spatlocs)) {
          if(!is.null(x[[obj_i]]@spatVectorCentroids)) {

            centroids = x[[obj_i]]@spatVectorCentroids
            centroidsDT = spatVector_to_dt(centroids)
            centroidsDT_loc = centroidsDT[, .(poly_ID, x, y)]
            colnames(centroidsDT_loc) = c('cell_ID', 'sdimx', 'sdimy')

            locsObj = create_spat_locs_obj(name = 'raw',
                                           coordinates = centroidsDT_loc,
                                           spat_unit = x[[obj_i]]@name, # tag same spat_unit as poly
                                           provenance = x[[obj_i]]@name, # TODO change this if polygons get prov
                                           misc = initialize)

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            .external_accessor_spatloc = list(
              # set spatlocs 'spat_unit' using 'name' if it was EXPLICITLY
              # supplied to setPolygonInfo
              # otherwise, set spatlocs 'spat_unit' as x@name
              nospec_unit = nospec_name,
              # set spatlocs name based on locsObj@name
              nospec_name = TRUE
            )
            gobject = set_spatial_locations(gobject,
                                            spatlocs = locsObj,
                                            spat_unit = name, # useif explicit here
                                            verbose = verbose,
                                            set_defaults = FALSE,
                                            initialize = initialize)
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
          }
        }
      }
      return(gobject)

    }
  }

  # catch
  stop(wrap_txt('Only giottoPolygon or lists of giottoPolygon accepted.
                For raw or external data, please first use readPolygonData()',
                errWidth = TRUE))
}






#' @title Set polygon info
#' @name set_polygon_info
#' @description Set giotto polygon spatVector
#' @inheritParams data_access_params
#' @param polygon_name name of polygons. Default "cell" (only used when gpolygon
#' is length of 1)
#' @param gpolygon giottoPolygon object
#' @return giotto object
#' @family polygon info data accessor functions
#' @family functions to set data in giotto object
#' @export
set_polygon_info = function(gobject,
                            gpolygon,
                            polygon_name = 'cell',
                            verbose = TRUE,
                            initialize = FALSE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(gpolygon)) stop(wrap_txt('gpolygon param must be given'))

  # 0. stop if not native formats
  if(inherits(gpolygon, 'list')) {
    if(!all(sapply(gpolygon, class) == 'giottoPolygon')) {
      stop(wrap_txt('If providing a list to internal setter, only lists of',
                    'giottoPolygon objects are permitted',
                    errWidth = TRUE))
    }
  }
  if(!inherits(gpolygon, c('giottoPolygon', 'NULL', 'list'))) {
    stop(wrap_txt(deparse(substitute(gpolygon)), 'is not a giottoPolygon (set),
                  list of giottoPolygons (set),
                  or NULL (remove)'))
  }

  # 1. determine user input
  p = parent.frame() # get values if called from external
  call_from_external = exists('.external_accessor_poly', where = p)

  if(isTRUE(call_from_external)) {
    nospec_name = p$.external_accessor_poly$nospec_name
  } else {
    nospec_name = !methods::hasArg(polygon_name)
  }


  # use name instead of polygon_name for compatibility with S4 reading
  name = polygon_name

  # 2. set default spat_unit
  # not needed when general default exists ('cell')

  # 3.1 if input is NULL, remove object
  if(is.null(gpolygon)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to gpolygon.
                                 Removing specified polygon information')
    gobject@spatial_info[[name]] = NULL

    # prune if empty
    if(length(gobject@spatial_info[[name]]) == 0L) {
      gobject@spatial_info = NULL
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }

  # 3.2 if input is list, set list
  if(inherits(gpolygon, 'list')) {
    # ensure list names are accurate
    gp_names = names(gpolygon)
    if(is.null(gp_names)) stop(wrap_txt('if "gpolygon" is a list, then it must be a named list',
                                        errWidth = TRUE))
    if(any(is.na(gp_names))) stop(wrap_txt('No NA values allowed in "gpolygon" list names'))
    dup_bool = duplicated(gp_names)
    if(any(dup_bool)) stop(wrap_txt('Duplicated list names:', gp_names[dup_bool],
                                      '\nAll gpolygon list names must be unique'))

    # iterate through list
    for(gp_name in gp_names) {
      gpolygon[[gp_name]]@name = gp_name

      ## check if specified name has already been used
      potential_names = names(gobject@spatial_info)
      if(gp_name %in% potential_names) {
        if(verbose) wrap_msg('> "', gp_name, '" already exists and will be replaced with new giotto polygon \n')
      }

      # set items
      gobject@spatial_info[[gp_name]] = gpolygon[[gp_name]]
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }



  # 4. import data from S4 if available (for single objects)
  # NOTE: modifies name/gpolygon
  gpolygon = read_s4_nesting(gpolygon)


  # 5. check if specified name has already been used
  potential_names = names(gobject@spatial_info)
  if(name %in% potential_names) {
    if(verbose) wrap_msg('> "', name, '" already exists and will be replaced with new giotto polygon \n')
  }


  ## 6. update and return giotto object
  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting polygon info [', objName(gpolygon), '] ', sep = ''
  )

  gobject@spatial_info[[name]] = gpolygon
  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}





## feature info ####

#' @title Get feature info
#' @name getFeatureInfo
#' @description Get giotto points spatVector
#' @inheritParams data_access_params
#' @param return_giottoPoints return as a giottoPoints object
#' @family feature info data accessor functions
#' @family functions to get data from giotto object
#' @export
getFeatureInfo = function(gobject = gobject,
                          feat_type = NULL,
                          return_giottoPoints = FALSE,
                          set_defaults = TRUE) {
  if (!inherits(gobject, 'giotto')){
    wrap_msg("Unable to get giotto points spatVector feature info from non-Giotto object.")
    stop(wrap_txt("Please provide a Giotto object to the gobject argument.",
                  errWidth = TRUE))
  }
  feat_info = get_feature_info(gobject = gobject,
                               feat_type = feat_type,
                               return_giottoPoints = return_giottoPoints,
                               set_defaults = set_defaults)
  return(feat_info)
}

#' @title Get feature info
#' @name get_feature_info
#' @param return_giottoPoints return as a giottoPoints object
#' @description Get giotto points spatVector
#' @inheritParams data_access_params
#' @family feature info data accessor functions
#' @family functions to get data from giotto object
#' @return a SpatVector (default) or giottoPoints object depending on value of
#' return_giottoPoints
#' @export
get_feature_info = function(gobject,
                            feat_type = NULL,
                            set_defaults = TRUE,
                            return_giottoPoints = FALSE) {

  guard_against_notgiotto(gobject)

  # specify feat_type
  if(isTRUE(set_defaults)) {
    feat_type = set_default_feat_type(gobject = gobject,
                                      feat_type = feat_type)
  }

  potential_names = names(gobject@feat_info)

  if(!feat_type %in% potential_names) {
    stop('There is no feature information with name ', feat_type, '\n')
  } else {
    feat_info = gobject@feat_info[[feat_type]]
    if(return_giottoPoints) {
      return(feat_info)
    } else {
      return(feat_info@spatVector)
    }
  }
}






#' @description Get list of all feature information
#' @keywords internal
#' @return list of giottoPoints or SpatVector depending on return_giottoPoints
#' param
#' @noRd
get_feature_info_list = function(gobject,
                                 return_giottoPoints = TRUE) {

  guard_against_notgiotto(gobject)

  data_list = slot(gobject, 'feat_info')

  # return objects
  if(isTRUE(return_giottoPoints)) {
    return(data_list)
  } else {
    return(lapply(data_list, `[`))
  }
}






#' @title Set feature info
#' @name setFeatureInfo
#' @description Set giotto polygon spatVector for features
#' @inheritParams data_access_params
#' @param x giottoPoints object or list of giottoPoints to set. Passing NULL
#' will remove the specified giottoPoints object from the giotto object
#' @param verbose be verbose
#' @return giotto object
#' @family feature info data accessor functions
#' @family functions to set data in giotto object
#' @export
setFeatureInfo = function(gobject,
                          x,
                          feat_type = NULL,
                          verbose = TRUE,
                          initialize = TRUE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x param (data to set) must be given'))

  # 1. Determine user inputs
  nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  .external_accessor_point = list(nospec_feat = nospec_feat)
  # checked by internal setter to determine if called by external


  # Extra defaults: expression, feature_info, spat_info specific
  default_feat = if(is.null(gobject@expression_feat)) 'rna' else gobject@expression_feat[[1L]]


  # NATIVE INPUT TYPES
  # 2. if input is giottoPoints or NULL, pass to internal
  if(is.null(x) | inherits(x, 'giottoPoints')) {

    # pass to internal
    gobject = set_feature_info(
      gobject = gobject,
      gpoints = x,
      feat_type = feat_type,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)

  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'giottoPoints')) {

      # MULTIPLE INPUT
      # 3. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_feature_info(
          gobject = gobject,
          gpoints = x[[obj_i]],
          feat_type = feat_type,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)
    }
  }

  # catch
  stop(wrap_txt('Only giottoPoints or lists of giottoPoints accepted.
                For raw or external data, please first use readFeatureInfo()'))

}








#' @title Set feature info
#' @name set_feature_info
#' @description Set giotto polygon spatVector for features
#' @inheritParams data_access_params
#' @param gpoints giotto points object
#' @param gpolygon typo do not use
#' @param verbose be verbose
#' @return giotto object
#' @family feature info data accessor functions
#' @family functions to set data in giotto object
#' @export
set_feature_info = function(gobject,
                            gpoints,
                            feat_type = NULL,
                            verbose = TRUE,
                            set_defaults = TRUE,
                            initialize = FALSE,
                            gpolygon = NULL) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(gpoints) & !methods::hasArg(gpolygon)) stop(wrap_txt('gpoints param must be given'))

  if(!is.null(gpolygon)) { # deprecation
    warning(wrap_txt('do not use gpolygon param. Use gpoints instead'))
    if(is.null(gpoints)) gpoints = gpolygon
  }


  # 0. stop if not native formats
  if(inherits(gpoints, 'list')) {
    if(!all(sapply(gpoints, class) == 'giottoPoints')) {
      stop(wrap_txt('If providing a list to internal setter, only lists of',
                    'giottoPoints objects are permitted',
                    errWidth = TRUE))
    }
  }
  if(!inherits(gpoints, c('giottoPoints', 'NULL', 'list'))) {
    stop(wrap_txt(deparse(substitute(gpoints)), 'is not giottoPoints (set),
                  list of giottoPoints (set),
                  or NULL (remove)'))
  }


  # 1. determine user input
  p = parent.frame() # get values if called from external
  call_from_external = exists('.external_accessor_point', where = p)

  if(isTRUE(call_from_external)) {
    nospec_feat = p$.external_accessor_point$nospec_feat
  } else {
    nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  }


  # 2. set default feat_type
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 3.1 if input is NULL, remove object
  if(is.null(gpoints)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to gpoints.
                                 Removing specified feature information.')
    gobject@feat_info[[feat_type]] = NULL

    # prune if empty
    if(length(gobject@feat_info[[feat_type]]) == 0L) {
      gobject@feat_info = NULL
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }

  # 3.2 if input is list, set list
  if(inherits(gpoints, 'list')) {
    # ensure list names are accurate
    gp_names = names(gpoints)
    if(is.null(gp_names)) stop(wrap_txt('If "gpoints" is a list, then it must be a named list',
                                        errWidth = TRUE))
    if(any(is.na(gp_names))) stop(wrap_txt('No NA values allowed in "gpoints" list names"'))
    dup_bool = duplicated(gp_names)
    if(any(dup_bool)) stop(wrap_txt('Duplicated list names:', gp_names[dup_bool],
                                    '\nAll gpoints list names must be unique'))
    for(gp_name in gp_names) {
      featType(gpoints[[gp_name]]) = gp_name
    }

    # replacements warning already given during extract points list (external setter)
    # remove items to replace
    for(gp_name in gp_names) {
      gobject@feat_info[[gp_name]] = gpoints[[gp_name]]
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }


  # 4. import data from S4 if available
  # NOTE: modifies feat_type/gpoints
  gpoints = read_s4_nesting(gpoints)


  ## 5. check if specified name has already been used
  potential_names = names(gobject@feat_info)
  if(feat_type %in% potential_names) {
    if(isTRUE(verbose)) wrap_msg('> "', feat_type, '" already exists and will be replaced with new giotto points \n')
  }

  ## 6. update and return giotto object
  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting feature info [', featType(gpoints), '] ', sep = ''
  )

  gobject@feat_info[[feat_type]] = gpoints
  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}









## spatial enrichment slot ####


#' @title Get spatial enrichment
#' @name get_spatial_enrichment
#' @description Function to get a spatial enrichment data.table
#' @inheritParams data_access_params
#' @param enrichm_name name of spatial enrichment results. Default "DWLS"
#' @return data.table with fractions
#' @family spatial enrichment data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatial_enrichment = function(gobject,
                                  spat_unit = NULL,
                                  feat_type = NULL,
                                  enrichm_name = 'DWLS',
                                  output = c('spatEnrObj', 'data.table'),
                                  copy_obj = TRUE,
                                  set_defaults = TRUE) {


  output = match.arg(output, choices = c('spatEnrObj', 'data.table'))

  # Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # spatial locations
  # if NULL (not given) and spatial locations have been added, then use first one
  # if NULL (not given) and spatial locations have NOT been added, then keep NULL
  if(is.null(enrichm_name)) {
    if(!is.null(gobject@spatial_enrichment)) {
      enrichm_name = list_spatial_enrichments_names(gobject,
                                                    spat_unit = spat_unit,
                                                    feat_type = feat_type)[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      enrichm_name = NULL
      cat('No spatial enrichment results have been found \n')
      return(NULL)
    }
  }


  potential_names = list_spatial_enrichments_names(gobject,
                                                   spat_unit = spat_unit,
                                                   feat_type = feat_type)

  if(enrichm_name %in% potential_names) {
    enr_res = gobject@spatial_enrichment[[spat_unit]][[feat_type]][[enrichm_name]]

    if(isTRUE(copy_obj)) enr_res[] = data.table::copy(enr_res[])

    if(output == 'spatEnrObj') {
      return(enr_res)
    } else if(output == 'data.table') {
      return(enr_res[])
    }

  } else {
    stop("The spatial enrichment result with name ","'", enrichm_name, "'"," can not be found \n")
  }
}

#' @title Get spatial enrichment
#' @name getSpatialEnrichment
#' @description Function to get a spatial enrichment data.table
#' @inheritParams data_access_params
#' @param name name of spatial enrichment results. Default "DWLS"
#' @return data.table with fractions
#' @family spatial enrichment data accessor functions
#' @family functions to get data from giotto object
#' @export
getSpatialEnrichment = function(gobject,
                                spat_unit = NULL,
                                feat_type = NULL,
                                name = 'DWLS',
                                output = c('spatEnrObj', 'data.table'),
                                copy_obj = TRUE,
                                set_defaults = TRUE) {

  # Pass to internal function
  enr_res = get_spatial_enrichment(gobject = gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   enrichm_name = name,
                                   output = output,
                                   copy_obj = copy_obj,
                                   set_defaults = set_defaults)

  return(enr_res)
}





#' @description Get all spatial enrichments for a specified spatial unit and
#' feature type
#' @keywords internal
#' @return list of spatEnrObj or data.table depending on output param
#' @noRd
get_spatial_enrichment_list = function(gobject,
                                       spat_unit = NULL,
                                       feat_type = NULL,
                                       output = c('spatEnrObj', 'data.table'),
                                       copy_obj = TRUE,
                                       set_defaults = TRUE) {

  guard_against_notgiotto(gobject)

  output = match.arg(output, choices = c('spatEnrObj', 'data.table'))

  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  data_list = slot(gobject, 'spatial_enrichment')[[spat_unit]][[feat_type]]

  # copy object
  if(isTRUE(copy_obj)) data_list = lapply(data_list, copy)


  # return object
  if(output == 'spatEnrObj') {
    return(data_list)
  }
  if(output == 'data.table') {
    return(lapply(data_list, `[`))
  }
}






#' @title Set spatial enrichment
#' @name setSpatialEnrichment
#' @description Function to set a spatial enrichment slot
#' @inheritParams data_access_params
#' @param name name of spatial enrichment results. Default "DWLS"
#' @param x spatEnrObj or list of spatEnrObj to set. Passing NULL will remove
#' a specified set of spatial enrichment information from the gobject.
#' @param verbose be verbose
#' @return giotto object
#' @family spatial enrichment data accessor functions
#' @family functions to set data in giotto object
#' @export
setSpatialEnrichment = function(gobject,
                                x,
                                spat_unit = NULL,
                                feat_type = NULL,
                                name = 'enrichment',
                                provenance = NULL,
                                verbose = TRUE,
                                initialize = TRUE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(x)) stop(wrap_txt('x param (data to set) must be given'))

  # check hierarchical slots
  avail_ex = list_expression(gobject)
  avail_sl = list_spatial_locations(gobject)
  if(is.null(avail_ex)) stop(wrap_txt('Add expression and spatial information first'))
  if(is.null(avail_sl)) stop(wrap_txt('Add spatial location information first'))

  # 1. determine user inputs
  nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
  nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
  nospec_name = ifelse(is.null(match.call()$name), yes = TRUE, no = FALSE)
  .external_accessor_spatenr = list(nospec_unit = nospec_unit,
                                    nospec_feat = nospec_feat,
                                    nospec_name = nospec_name)
  # checked by internal setter to determine if called by external


  # NATIVE INPUT TYPES
  # 3. if input is spatEnrObj or NULL, pass to internal
  if(is.null(x) | inherits(x, 'spatEnrObj')) {

    # pass to internal
    gobject = set_spatial_enrichment(
      gobject = gobject,
      spatenrichment = x,
      spat_unit = spat_unit,
      feat_type = feat_type,
      enrichm_name = name,
      verbose = verbose,
      set_defaults = FALSE,
      initialize = initialize
    )
    return(gobject)

  } else if(inherits(x, 'list')) {

    # check list items are native
    if(all(sapply(x, class) == 'spatEnrObj')) {

      # MULTIPLE INPUT
      # 4. iteratively set
      for(obj_i in seq_along(x)) {

        # if(isTRUE(verbose)) message('[', obj_i, ']')

        gobject = set_spatial_enrichment(
          gobject = gobject,
          spatenrichment = x[[obj_i]],
          spat_unit = spat_unit,
          feat_type = feat_type,
          enrichm_name = name,
          verbose = verbose,
          set_defaults = FALSE,
          initialize = initialize
        )
      }
      return(gobject)

    }
  }
  # catch
  stop(wrap_txt('Only spatEnrObj or lists of spatEnrObj accepted.
                For raw or external data, please first use readSpatEnrichData()'))

}








#' @title Set spatial enrichment
#' @name set_spatial_enrichment
#' @description Function to set a spatial enrichment slot
#' @inheritParams data_access_params
#' @param enrichm_name name of spatial enrichment results. Default "DWLS"
#' @param spatenrichment spatial enrichment results
#' @param verbose be verbose
#' @return giotto object
#' @family spatial enrichment data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatial_enrichment = function(gobject,
                                  spatenrichment,
                                  spat_unit = NULL,
                                  feat_type = NULL,
                                  enrichm_name = 'enrichment',
                                  provenance = NULL,
                                  verbose = TRUE,
                                  set_defaults = TRUE,
                                  initialize = TRUE) {

  guard_against_notgiotto(gobject)
  if(!methods::hasArg(spatenrichment)) stop(wrap_txt('spatenrichment param must be given'))

  # 0. stop if not native formats
  if(!inherits(spatenrichment, c('spatEnrObj', 'NULL'))) {
    stop(wrap_txt(deparse(substitute(spatenrichment)), 'is not spatEnrObj (set).
                  or NULL (remove)'))
  }

  # 1. Check user input
  p = parent.frame() # Get values if called from external
  call_from_external = exists('.external_accessor_spatenr', where = p)

  if(isTRUE(call_from_external)) {
    nospec_unit = p$.external_accessor_spatenr$nospec_unit
    nospec_feat = p$.external_accessor_spatenr$nospec_feat
    nospec_name = p$.external_accessor_spatenr$nospec_name
  } else {
    nospec_unit = ifelse(is.null(spat_unit), yes = TRUE, no = FALSE)
    nospec_feat = ifelse(is.null(feat_type), yes = TRUE, no = FALSE)
    nospec_name = ifelse(is.null(match.call()$enrichm_name), yes = TRUE, no = FALSE)
  }

  # change var name to be compatible with read_S4_nesting()
  name = enrichm_name


  # 2. Set feat_type and spat_unit
  if(isTRUE(set_defaults)) {
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)
  }

  # 3. Remove object if input is NULL
  if(is.null(spatenrichment)) {
    if(isTRUE(verbose)) wrap_msg('NULL passed to spatenrichment.
                                 Removing specified spatial enrichment.')
    gobject@spatial_enrichment[[spat_unit]][[feat_type]][[name]] = NULL

    # prune if empty
    if(length(gobject@spatial_enrichment[[spat_unit]][[feat_type]]) == 0L) {
      gobject@spatial_enrichment = NULL
    }

    if(isTRUE(initialize)) return(initialize(gobject))
    else return(gobject)
  }

  # 4. Import info from S4 if given
  # NOTE: modifies spat_unit/feat_type/name/provenance/spatenrichment
  spatenrichment = read_s4_nesting(spatenrichment)
  # use updated values in name instead of enrichm_name


  # 5. check if specified name has already been used
  if(isTRUE(verbose)) {
    potential_names = list_spatial_enrichments_names(gobject = gobject,
                                                     spat_unit = spat_unit,
                                                     feat_type = feat_type)
    if(name %in% potential_names) {
      wrap_msg('> "', name, '" already exists and will be replaced with new spatial enrichment results')
    }
  }

  # 6. update and return giotto object
  if(isTRUE(verbose) & isTRUE(call_from_external)) wrap_msg(
    'Setting spatial enrichment [', spatUnit(spatenrichment), '][', featType(spatenrichment), '] ',
    objName(spatenrichment), sep = ''
  )

  gobject@spatial_enrichment[[spat_unit]][[feat_type]][[name]] = spatenrichment
  if(isTRUE(initialize)) return(initialize(gobject))
  else return(gobject)

}






## MG image slot ####

#' @title Get \emph{magick}-based giotto \code{image}
#' @name get_giottoImage_MG
#' @description Get a giottoImage from a giotto object
#' @param gobject giotto object
#' @param name name of giottoImage \code{\link{showGiottoImageNames}}
#' @return a giottoImage
#' @keywords internal
get_giottoImage_MG = function(gobject = NULL,
                              name = NULL) {

  if(is.null(gobject)) stop('The giotto object holding the giottoImage needs to be provided \n')
  g_image_names = names(gobject@images)
  if(is.null(g_image_names)) stop('No giottoImages have been found \n')

  if(is.null(name)) {
    name = g_image_names[1]
  }

  if(!name %in% g_image_names) stop(name, ' was not found among the image names, see showGiottoImageNames()')

  g_image = gobject@images[[name]]

  return(g_image)
}



#' @title Set \emph{magick}-based giotto \code{image}
#' @name set_giottoImage_MG
#' @description Set a giottoImage for a giotto object with no additional modifications
#' @param gobject giotto object
#' @param image_object a giottoImage object
#' @param name name to assign giottoImage
#' @param verbose be verbose
#' @return giotto object
#' @keywords internal
set_giottoImage_MG = function(gobject,
                              image_object,
                              name = NULL,
                              verbose = TRUE) {

  # Check params
  if(is.null(gobject)) stop('gobject must be given \n')
  if(is.null(image_object)) stop('image_object to be attached must be given \n')

  # Default to name present in image object name slot
  if(is.null(name)) name = image_object@name

  # Find existing names
  potential_names = list_images_names(gobject = gobject, img_type = 'image')

  if(verbose == TRUE) {
    if(name %in% potential_names) wrap_msg('> "' ,name, '" already exists and will be replaced with new image object \n')
  }

  gobject@images[[name]] = image_object

  return(gobject)

}



## large image slot ####



#' @title Get \emph{terra}-based giotto \code{largeImage}
#' @name get_giottoLargeImage
#' @description Set a giottoLargeImage from a giottoObject
#' @param gobject giotto object
#' @param name name of giottoLargeImage \code{\link{showGiottoImageNames}}
#' @return a giottoLargeImage
#' @keywords internal
get_giottoLargeImage = function(gobject = NULL,
                                name = NULL) {

  if(is.null(gobject)) stop('The giotto object holding the giottoLargeImage needs to be provided \n')
  g_image_names = names(gobject@largeImages)
  if(is.null(g_image_names)) stop('No giottoLargeImages have been found \n')

  if(is.null(name)) {
    name = g_image_names[1]
  }

  if(!name %in% g_image_names) stop(name,' was not found among the largeImage names. See showGiottoImageNames() \n')

  g_imageL = gobject@largeImages[[name]]

  return(g_imageL)
}




#' @title Set \emph{terra}-based giotto \code{largeImage}
#' @name set_giottoLargeImage
#' @description Set a giottoLargeImage for a giotto object with no additional modifications
#' @param gobject giotto object
#' @param largeImage_object a giottoLargeImage object
#' @param name name to assign giottoLargeImage
#' @param verbose be verbose
#' @return giotto object
#' @keywords internal
set_giottoLargeImage = function(gobject,
                                largeImage_object,
                                name = NULL,
                                verbose = TRUE) {

  # Check params
  if(is.null(gobject)) stop('gobject must be given \n')
  if(is.null(largeImage_object)) stop('largeImage_object to be attached must be given \n')

  # Default to name present in image object name slot
  if(is.null(name)) name = largeImage_object@name

  # Find existing names
  potential_names = list_images_names(gobject = gobject, img_type = 'largeImage')

  if(verbose == TRUE) {
    if(name %in% potential_names) wrap_msg('> "' ,name, '" already exists and will be replaced with new image object \n')
  }

  gobject@largeImages[[name]] = largeImage_object

  return(gobject)

}



## all image slots ####



#' @title Get giotto image object
#' @name get_giottoImage
#' @description Get giotto image object from gobject
#' @param gobject giotto object
#' @param image_type type of giotto image object. Either "image" or "largeImage"
#' @param name name of a giotto image object \code{\link{showGiottoImageNames}}
#' @return a giotto image object
#' @family image data accessor functions
#' @family functions to get data from giotto object
#' @export
get_giottoImage = function(gobject = NULL,
                           image_type = c('image','largeImage'),
                           name = NULL) {


  # Check image type
  image_type = match.arg(image_type, choices = c('image','largeImage'))

  # Select get function
  if(image_type == 'image') {
    g_img = get_giottoImage_MG(gobject = gobject,
                               name = name)
  }
  if(image_type == 'largeImage') {
    g_img = get_giottoLargeImage(gobject = gobject,
                                 name = name)
  }
  return(g_img)
}

#' @title Get giotto image object
#' @name getGiottoImage
#' @description Get giotto image object from gobject
#' @param gobject giotto object
#' @param image_type type of giotto image object. Either "image" or "largeImage"
#' @param name name of a giotto image object \code{\link{showGiottoImageNames}}
#' @return a giotto image object
#' @family image data accessor functions
#' @family functions to get data from giotto object
#' @export
getGiottoImage = function(gobject = NULL,
                          image_type = c('image','largeImage'),
                          name = NULL) {
  if (!inherits(gobject, 'giotto')){
    wrap_msg("Unable to get Giotto Image from non-Giotto object.")
    stop(wrap_txt("Please provide a Giotto object to the gobject argument.",
                  errWidth = TRUE))
  }

  g_img = get_giottoImage(gobject = gobject,
                          image_type = image_type,
                          name = name)

  return (g_img)

}






#' @description Get list of all giottoImages
#' @keywords internal
#' @noRd
get_giotto_image_list = function(gobject,
                                 image_type = c('image', 'largeImage')) {

  guard_against_notgiotto(gobject)

  image_type = match.arg(image_type, choices = c('image', 'largeImage'))

  # return object
  if(image_type == 'image') {
    return(
      slot(gobject, 'images')
    )
  }
  if(image_type == 'largeImage') {
    return(
      slot(gobject, 'largeImages')
    )
  }
}








#' @title Set giotto image object
#' @name set_giottoImage
#' @description Directly attach a giotto image to giotto object
#' @details \emph{\strong{Use with care!}} This function directly attaches giotto image
#'   objects to the gobject without further modifications of spatial positioning values
#'   within the image object that are generally needed in order for them to
#'   plot in the correct location relative to the other modalities of spatial data. \cr
#'   For the more general-purpose method of attaching image objects, see \code{\link{addGiottoImage}}
#' @param gobject giotto object
#' @param image giotto image object to be attached without modification to the
#'   giotto object
#' @param image_type type of giotto image object. Either "image" or "largeImage"
#' @param name name of giotto image object
#' @param verbose be verbose
#' @return giotto object
#' @family image data accessor functions
#' @family functions to set data in giotto object
#' @seealso \code{\link{addGiottoImage}}
#' @export
set_giottoImage = function(gobject = NULL,
                           image = NULL,
                           image_type = NULL,
                           name = NULL,
                           verbose = TRUE) {


  # Check image type
  image_type = match.arg(image_type, choices = c('image','largeImage'))

  # Select set function
  if(image_type == 'image') {
    gobject = set_giottoImage_MG(gobject = gobject,
                                 image_object = image,
                                 name = name,
                                 verbose = verbose)
  }
  if(image_type == 'largeImage') {
    gobject = set_giottoLargeImage(gobject = gobject,
                                   largeImage_object = image,
                                   name = name,
                                   verbose = verbose)
  }
  return(gobject)
}

#' @title Set giotto image object
#' @name setGiottoImage
#' @description Directly attach a giotto image to giotto object
#' @details \emph{\strong{Use with care!}} This function directly attaches giotto image
#'   objects to the gobject without further modifications of spatial positioning values
#'   within the image object that are generally needed in order for them to
#'   plot in the correct location relative to the other modalities of spatial data. \cr
#'   For the more general-purpose method of attaching image objects, see \code{\link{addGiottoImage}}
#' @param gobject giotto object
#' @param image giotto image object to be attached without modification to the
#'   giotto object
#' @param image_type type of giotto image object. Either "image" or "largeImage"
#' @param name name of giotto image object
#' @param verbose be verbose
#' @return giotto object
#' @family image data accessor functions
#' @family functions to set data in giotto object
#' @seealso \code{\link{addGiottoImage}}
#' @export
setGiottoImage = function(gobject = NULL,
                          image = NULL,
                          image_type = NULL,
                          name = NULL,
                          verbose = TRUE){

  if (!inherits(gobject, 'giotto')){
    wrap_msg("Unable to set Giotto Image to non-Giotto object.")
    stop(wrap_txt("Please provide a Giotto object to the gobject argument.",
                  errWidth = TRUE))
  } else if (is.null(image)) {
    wrap_msg("Warning: image argument set to NULL. Replacing current image slot with NULL will remove the image.")
  } else if( !"giottoImage" %in% image || !"giottoLargeImage" %in% image) {
    wrap_msg("Unable to set non-giottoImage objects. Please ensure a giottoImage or giottoLargeImage is provided to this function.")
    wrap_msg("See createGiottoImage or createGiottoLargeImage for more details.")
    stop(wrap_txt("Unable to set non-giottoImage object.",
                  errWidth = TRUE))
  }

  gobject = set_giottoImage(gobject = gobject,
                          image = image,
                          image_type = image_type,
                          name = name,
                          verbose = verbose)

  return (gobject)

}


## Show functions ####

#' @title showGiottoExpression
#' @name showGiottoExpression
#' @description shows the available matrices
#' @param gobject giotto object
#' @param nrows number of rows to print for each matrix (ignored for sparse matrices)
#' @param ncols number of columns to print for each matrix (ignored for sparse matrices)
#' @return prints the name and small subset of available matrices
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoExpression = function(gobject, nrows = 4, ncols = 4) {

  # import print styling
  ch = box_chars()
  ct = color_tag()

  # 1. check inputs
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. get availability matrix
  available_data = list_expression(gobject = gobject)
  if(is.null(available_data)) {
    cat('No expression data available \n')
  } else {

    # 3.1 set up object printouts
    objPrints = list()
    for(obj_i in seq(nrow(available_data))) {

      # get object
      dataObj = get_expression_values(gobject = gobject,
                                      values = available_data$name[[obj_i]],
                                      spat_unit = available_data$spat_unit[[obj_i]],
                                      feat_type = available_data$feat_type[[obj_i]],
                                      output = 'exprObj')

      # collect object prints
      if(inherits(dataObj[], 'sparseMatrix')) {
        objPrints[[obj_i]] = capture.output(show(dataObj))
      } else if(inherits(dataObj[], c('denseMatrix', 'matrix', 'data.frame'))) {
        objPrints[[obj_i]] = capture.output(abb_mat(dataObj, nrows, ncols))
      } else {
        # directly print slot (catch)
        objPrints[[obj_i]] = capture.output(show(slot(dataObj, 'exprMat')))
      }

    }

    # object printblock edits
    objPrints = lapply(objPrints, function(x) paste0(ch$s, x)) # add indent
    objPrints = lapply(objPrints, function(x) paste(x, collapse = ('\n'))) # linearize print

    # append to availability table
    available_data$values = unlist(objPrints)

    # 3.2 setup general prints
    if(isTRUE(use_color_text())) {
      available_data$spat_unit = paste0('Spatial unit "', ct$b, available_data$spat_unit, ct$x, '"')
      available_data$feat_type = paste0('Feature type "', ct$r, available_data$feat_type, ct$x, '"')
      available_data$name = paste0('Expression data "', ct$t, available_data$name, ct$x, '" values:')
    } else {
      available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
      available_data$feat_type = paste0('Feature type "', available_data$feat_type, '"')
      available_data$name = paste0('Expression data "', available_data$name, '" values:')
    }


    # 4. print information
    print_leaf(level_index = 1,
               availableDT = available_data,
               inherit_last = TRUE,
               indent = '')

  }

}



#' @title showGiottoCellMetadata
#' @name showGiottoCellMetadata
#' @description shows the available cell metadata
#' @param gobject giotto object
#' @param nrows number of rows to print for each metadata
#' @return prints the name and small subset of available metadata
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoCellMetadata = function(gobject,
                                  nrows = 4) {

  # import print styling
  ch = box_chars()
  ct = color_tag()

  # 1. check inputs
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. get availability matrix
  available_data = list_cell_metadata(gobject = gobject)
  if(is.null(available_data)) {
    cat('No cell metadata available \n')
  } else {

    # 3.1 set up object printouts
    objPrints = objRows = list()
    for(obj_i in seq(nrow(available_data))) {

      # get object
      dataObj = get_cell_metadata(gobject = gobject,
                                  spat_unit = available_data$spat_unit[[obj_i]],
                                  feat_type = available_data$feat_type[[obj_i]],
                                  output = 'cellMetaObj',
                                  copy_obj = TRUE)

      # collect object prints
      objRows[[obj_i]] = nrow(dataObj[])

      objPrints[[obj_i]] =
        dataObj[1:if(nrows <= objRows[[obj_i]]) nrows else objRows[[obj_i]],] %>%
        print %>%
        capture.output

    }

    # object printblock edits
    objPrints = lapply(objPrints, function(x) paste0(ch$s, x)) # add indent
    objPrints = lapply(objPrints, function(x) paste(x, collapse = ('\n'))) # linearize print

    # append to availability table
    available_data$values = unlist(objPrints)

    # 3.2 setup general prints
    if(isTRUE(use_color_text())) {
      available_data$spat_unit = paste0('Spatial unit "', ct$b, available_data$spat_unit, ct$x, '"')
      available_data$feat_type = paste0('Feature type "', ct$r, available_data$feat_type, ct$x, '"')
    } else {
      available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
      available_data$feat_type = paste0('Feature type "', available_data$feat_type, '"')
    }


    # 4. print information
    print_leaf(level_index = 1,
               availableDT = available_data,
               inherit_last = TRUE,
               indent = '')

  }

}


#' @title showGiottoFeatMetadata
#' @name showGiottoFeatMetadata
#' @description shows the available feature metadata
#' @param gobject giotto object
#' @param nrows number of rows to print for each metadata
#' @return prints the name and small subset of available metadata
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoFeatMetadata = function(gobject,
                                  nrows = 4) {

  # import print styling
  ch = box_chars()
  ct = color_tag()

  # 1. check inputs
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. get availability matrix
  available_data = list_feat_metadata(gobject = gobject)
  if(is.null(available_data)) {
    cat('No feature metadata available \n')
  } else {

    # 3.1 set up object printouts
    objPrints = objRows = list()
    for(obj_i in seq(nrow(available_data))) {

      # get object
      dataObj = get_feature_metadata(gobject = gobject,
                                     spat_unit = available_data$spat_unit[[obj_i]],
                                     feat_type = available_data$feat_type[[obj_i]],
                                     output = 'featMetaObj',
                                     copy_obj = TRUE)

      # collect object prints
      objRows[[obj_i]] = nrow(dataObj[])

      objPrints[[obj_i]] =
        dataObj[1:if(nrows <= objRows[[obj_i]]) nrows else objRows[[obj_i]],] %>%
        print %>%
        capture.output

    }

    # object printblock edits
    objPrints = lapply(objPrints, function(x) paste0(ch$s, x)) # add indent
    objPrints = lapply(objPrints, function(x) paste(x, collapse = ('\n'))) # linearize print

    # append to availability table
    available_data$values = unlist(objPrints)

    # 3.2 setup general prints
    if(isTRUE(use_color_text())) {
      available_data$spat_unit = paste0('Spatial unit "', ct$b, available_data$spat_unit, ct$x, '"')
      available_data$feat_type = paste0('Feature type "', ct$r, available_data$feat_type, ct$x, '"')
    } else {
      available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
      available_data$feat_type = paste0('Feature type "', available_data$feat_type, '"')
    }


    # 4. print information
    print_leaf(level_index = 1,
               availableDT = available_data,
               inherit_last = TRUE,
               indent = '')

  }
}



#' @title showGiottoSpatLocs
#' @name showGiottoSpatLocs
#' @description shows the available spatial locations
#' @param gobject giotto object
#' @param nrows number of rows to print for each spatial location data.table
#' @return prints the name and small subset of available data.table
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatLocs = function(gobject,
                              nrows = 4) {

  # import print styling
  ch = box_chars()
  ct = color_tag()

  # 1. check inputs
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. get availability matrix
  available_data = list_spatial_locations(gobject = gobject)
  if(is.null(available_data)) {
    cat('No spatial locations available \n')
  } else {

    # 3.1 set up object printouts
    objPrints = objRows = list()
    for(obj_i in seq(nrow(available_data))) {

      # get object
      dataObj = get_spatial_locations(gobject = gobject,
                                      spat_unit = available_data$spat_unit[[obj_i]],
                                      spat_loc_name = available_data$name[[obj_i]],
                                      output = 'spatLocsObj',
                                      copy_obj = TRUE)

      # collect object prints
      objRows[[obj_i]] = nrow(dataObj[])

      objPrints[[obj_i]] = capture.output(abb_spatlocs(dataObj, nrows))

    }

    # object printblock edits
    objPrints = lapply(objPrints, function(x) paste0(ch$s, x)) # add indent
    objPrints = lapply(objPrints, function(x) paste(x, collapse = ('\n'))) # linearize print

    # append to availability table
    available_data$values = unlist(objPrints)

    # 3.2 setup general prints
    if(isTRUE(use_color_text())) {
      available_data$spat_unit = paste0('Spatial unit "', ct$b, available_data$spat_unit, ct$x, '"')
      available_data$name = paste0('S4 spatLocsObj "', ct$t, available_data$name, ct$x, '" coordinates:')
    } else {
      available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
      available_data$name = paste0('S4 spatLocsObj "', available_data$name, '" coordinates:')
    }
    for(obj_i in seq(nrow(available_data))) {
      available_data$name[[obj_i]] = paste0(available_data$name[[obj_i]],
                                            ch$s, '(', objRows[[obj_i]], ' rows)')
    }

    # 4. print information
    print_leaf(level_index = 1,
               availableDT = available_data,
               inherit_last = TRUE,
               indent = '')

  }

  # if(is.null(gobject)) stop('A giotto object needs to be provided \n')
  #
  # available_data = list_spatial_locations(gobject = gobject)
  # if(is.null(available_data)) cat('No spatial locations available \n')
  #
  # for(spatial_unit in unique(available_data$spat_unit)) {
  #
  #   cat('Spatial unit: ', spatial_unit, ' \n\n')
  #
  #   for(spatlocname in available_data[available_data$spat_unit == spatial_unit,]$name) {
  #     if(inherits(gobject@spatial_locs[[spatial_unit]][[spatlocname]], 'data.frame')) {
  #       cat('--> Name: ', spatlocname, ' \n\n')
  #       print(gobject@spatial_locs[[spatial_unit]][[spatlocname]][1:nrows,])
  #       cat('\n')
  #     }
  #     if(inherits(gobject@spatial_locs[[spatial_unit]][[spatlocname]], 'spatLocsObj')) {
  #       cat('--> Name: ', spatlocname, ' \n\n')
  #       cat('An object of class spatLocsObj\n')
  #       cat('Provenance: ', slot(gobject@spatial_locs[[spatial_unit]][[spatlocname]], 'provenance'), ' \n')
  #       print(slot(gobject@spatial_locs[[spatial_unit]][[spatlocname]], 'coordinates')[1:nrows,])
  #       cat('\n')
  #     }
  #   }
  # }

}


#' @title showGiottoSpatEnrichments
#' @name showGiottoSpatEnrichments
#' @description shows the available spatial enrichment results
#' @param gobject giotto object
#' @param nrows number of rows to print for each spatial enrichment data.table
#' @return prints the name and small subset of available data.table
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatEnrichments = function(gobject,
                                     nrows = 4) {

  # define for data.table [] subsetting
  spat_unit = NULL
  feat_type = NULL

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_spatial_enrichments(gobject = gobject)

  if(is.null(available_data)) cat('No spatial enrichments available \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial unit: ', spatial_unit, ' \n\n')

    for(feature_type in available_data[spat_unit == spatial_unit][['feat_type']]) {

      cat('--> Feature type: ', feature_type, ' \n\n')

      for(spatenrichname in available_data[spat_unit == spatial_unit][feat_type == feature_type][['name']]) {

        cat('----> Name ', spatenrichname, ': \n\n')

        print(gobject@spatial_enrichment[[spatial_unit]][[feature_type]][[spatenrichname]][][1:nrows,])

      }

    }

  }


}



#' @title showGiottoDimRed
#' @name showGiottoDimRed
#' @description shows the available dimension reductions
#' @param gobject giotto object
#' @param nrows number of coordinates rows to print
#' @param ncols number of coordinates columns to print
#' @return prints the name and small subset of available dimension reduction coordinates
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoDimRed = function(gobject,
                            nrows = 3,
                            ncols = 2) {

  # Define for data.table
  data_type = NULL

  # import print styling
  ch = box_chars()
  ct = color_tag()

  # 1. Check inputs
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. Get availability matrix
  available_data = list_dim_reductions(gobject)
  if(is.null(available_data)) {
    cat('No dimensional reductions available \n')
  } else {

    # 3.1 Set up object printouts
    objPrints = objRows = objCols = list()
    for(obj_i in seq(nrow(available_data))) {

      # Get object
      dataObj = get_dimReduction(gobject = gobject,
                                 reduction = available_data$data_type[[obj_i]],
                                 spat_unit = available_data$spat_unit[[obj_i]],
                                 feat_type = available_data$feat_type[[obj_i]],
                                 reduction_method = available_data$dim_type[[obj_i]],
                                 name = available_data$name[[obj_i]],
                                 output = 'data.table')

      # Collect object prints
      objRows[[obj_i]] = nrow(dataObj)
      objCols[[obj_i]] = ncol(dataObj)

      objPrints[[obj_i]] =
        dataObj[1:if(nrows <= objRows[[obj_i]]) nrows else objRows[[obj_i]],
                1:if(ncols <= objCols[[obj_i]]) ncols else objCols[[obj_i]]] %>%
        print %>%
        capture.output

    }

    # object printblock edits
    objPrints = lapply(objPrints, function(x) paste0(ch$s, x)) # Add indent
    objPrints = lapply(objPrints, function(x) paste(x, collapse = ('\n'))) # Linearize print

    # Append to availability table
    available_data$values = unlist(objPrints)

    # 3.2 Setup general prints
    if(isTRUE(use_color_text())) {
      available_data$spat_unit = paste0('Spatial unit "', ct$b, available_data$spat_unit, ct$x, '"')
      available_data$feat_type = paste0('Feature type "', ct$r, available_data$feat_type, ct$x, '"')
      available_data$dim_type = paste0('Dim reduction type "', ct$p, available_data$dim_type, ct$x, '"')
      available_data$name = paste0('S4 dimObj "', ct$t, available_data$name, ct$x, '" coordinates:')
    } else {
      available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
      available_data$feat_type = paste0('Feature type "', available_data$feat_type, '"')
      available_data$dim_type = paste0('Dim reduction type "', available_data$dim_type, '"')
      available_data$name = paste0('S4 dimObj "', available_data$name, '" coordinates:')
    }
    for(obj_i in seq(nrow(available_data))) {
      available_data$name[[obj_i]] = paste0(available_data$name[[obj_i]],
                                            ch$s ,'(', objRows[[obj_i]], ' rows ', objCols[[obj_i]], ' cols)')
    }

    # 4. Print information
    for(data_type_red in unique(available_data$data_type)) {
      data_type_subset = available_data$data_type == data_type_red

      if(isTRUE(use_color_text())) {
        if(data_type_red == 'feats') cat(paste0('Dim reduction on ', ct$y, 'features', ct$x, ':'))
        if(data_type_red == 'cells') cat(paste0('Dim reduction on ', ct$y, 'cells' , ct$x, ':'))
      } else {
        if(data_type_red == 'feats') cat(paste0('Dim reduction on ', 'features', ':'))
        if(data_type_red == 'cells') cat(paste0('Dim reduction on ', 'cells' , ':'))
      }

      cat('\n',
          '-------------------------',
          '\n\n.\n')

      print_leaf(level_index = 2, # skip over dim reduction layer
                 availableDT = available_data[data_type == data_type_red],
                 inherit_last = TRUE,
                 indent = '')

    }

  }

}




#' @title showGiottoNearestNetworks
#' @name showGiottoNearestNetworks
#' @description shows the available nearest neighbor networks
#' @param gobject giotto object
#' @param nrows number of network rows to print
#' @return prints the name and small subset of available nearest neighbor network info
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoNearestNetworks = function(gobject,
                                     nrows = 3) {

  # import print styling
  ch = box_chars()
  ct = color_tag()

  # 1. check input
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. get availability matrix
  available_data = list_nearest_networks(gobject)
  if(is.null(available_data)) {
    cat('No nearest neighbor networks available \n')
  } else {

    # 3.1 Set up object printouts
    objPrints = objRows = list()
    for(obj_i in seq(nrow(available_data))) {

      # Get object
      dataObj = get_NearestNetwork(gobject = gobject,
                                   spat_unit = available_data$spat_unit[[obj_i]],
                                   feat_type = available_data$feat_type[[obj_i]],
                                   nn_network_to_use = available_data$nn_type[[obj_i]],
                                   network_name = available_data$name[[obj_i]],
                                   output = 'data.table')

      # Collect object prints
      objRows[[obj_i]] = nrow(dataObj)

      objPrints[[obj_i]] =
        dataObj[1:if(nrows <= objRows[[obj_i]]) nrows else objRows[[obj_i]],] %>%
        print %>%
        capture.output

    }

    # object printblock edits
    objPrints = lapply(objPrints, function(x) paste0(ch$s, x)) # Add indent
    objPrints = lapply(objPrints, function(x) paste(x, collapse = ('\n'))) # Linearize print

    # Append to availability table
    available_data$values = unlist(objPrints)

    # 3.2 Setup general prints
    if(isTRUE(use_color_text())) {
      available_data$spat_unit = paste0('Spatial unit "', ct$b, available_data$spat_unit, ct$x, '"')
      if(!is.null(available_data$feat_type)) {
        available_data$feat_type = paste0('Feature type "', ct$r, available_data$feat_type, ct$x, '"')  # Check to be deprecated
      } else warning('Only networks from the deprecated nesting will be shown')
      available_data$nn_type = paste0('NN network type "', ct$p, available_data$nn_type, ct$x, '"')
      available_data$name = paste0('S4 nnNetObj "', ct$t, available_data$name, ct$x, '"')
    } else {
      available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
      if(!is.null(available_data$feat_type)) {
        available_data$feat_type = paste0('Feature type "', available_data$feat_type, '"')  # Check to be deprecated
      } else warning('Only networks from the deprecated nesting will be shown')
      available_data$nn_type = paste0('NN network type "', available_data$nn_type, '"')
      available_data$name = paste0('S4 nnNetObj "', available_data$name, '"')
    }
    for(obj_i in seq(nrow(available_data))) {
      available_data$name[[obj_i]] = paste0(available_data$name[[obj_i]],
                                            ch$s ,'(', objRows[[obj_i]], ' rows)')
    }

    # 4. Print information
    print_leaf(level_index = 1,
               availableDT = available_data,
               inherit_last = TRUE,
               indent = '')

  }

}




#' @title showGiottoSpatialInfo
#' @name showGiottoSpatialInfo
#' @description show the available giotto spatial polygon information
#' @param gobject giotto object
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatialInfo = function(gobject) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_spatial_info(gobject = gobject)
  if(is.null(available_data)) cat('No spatial info available \n')

  for(info in available_data$spat_info) {

    cat("For Spatial info: ", info, "\n\n")
    print(gobject@spatial_info[[info]])
    cat("-----------------------------")
    cat("\n \n")
  }

}


#' @title showGiottoFeatInfo
#' @name showGiottoFeatInfo
#' @description show the available giotto spatial feature information
#' @param gobject giotto object
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoFeatInfo = function(gobject) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_feature_info(gobject = gobject)
  if(is.null(available_data)) cat('No feature info available \n')

  for(info in available_data$feat_info) {

    cat("For Feature info: ", info, "\n\n")
    print(gobject@feat_info[[info]])
    cat("-----------------------------")
    cat("\n \n")
  }

}




#' @title showGiottoSpatNetworks
#' @name showGiottoSpatNetworks
#' @description Prints the available spatial networks that are attached to the Giotto object
#' @param gobject a giotto object
#' @param nrows number of rows to print
#' @return prints names and small subset of available spatial network info
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatNetworks = function(gobject,
                                  nrows = 4) {

  # import print styling
  ch = box_chars()
  ct = color_tag()

  # 1. Check input
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. Get availability matrix
  available_data = list_spatial_networks(gobject = gobject)
  if(is.null(available_data)) {
    cat('No spatial networks are available \n')
  } else {

    # 3.1 Set up object printouts
    objPrints = objRows = list()
    for(obj_i in seq(nrow(available_data))) {

      # Get object
      dataObj = get_spatialNetwork(gobject = gobject,
                                   spat_unit = available_data$spat_unit[[obj_i]],
                                   name = available_data$name[[obj_i]],
                                   output = 'networkDT')

      # Collect object prints
      objRows[[obj_i]] = nrow(dataObj)

      objPrints[[obj_i]] =
        dataObj[1:if(nrows <= objRows[[obj_i]]) nrows else objRows[[obj_i]],] %>%
        print %>%
        capture.output

    }

    # object printblock edits
    objPrints = lapply(objPrints, function(x) paste0(ch$s, x)) # Add indent
    objPrints = lapply(objPrints, function(x) paste(x, collapse = ('\n'))) # Linearize print

    # Append to availability table
    available_data$values = unlist(objPrints)

    # 3.2 Setup general prints
    if(isTRUE(use_color_text())) {
      available_data$spat_unit = paste0('Spatial unit "', ct$b, available_data$spat_unit, ct$x, '"')
      available_data$name = paste0('S4 spatialNetworkObj "', ct$t, available_data$name, ct$x, '"')
    } else {
      available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
      available_data$name = paste0('S4 spatialNetworkObj "', available_data$name, '"')
    }
    for(obj_i in seq(nrow(available_data))) {
      available_data$name[[obj_i]] = paste0(available_data$name[[obj_i]],
                                            ch$s ,'(', objRows[[obj_i]], ' rows)')
    }

    # 4. Print information
    print_leaf(level_index = 1,
               availableDT = available_data,
               inherit_last = TRUE,
               indent = '')
  }

}


#' @title Show networks
#' @name showNetworks
#' @inheritDotParams showGiottoSpatNetworks
#' @seealso \code{\link{showGiottoSpatNetworks}}
#' @export
showNetworks = function(...) {

  .Deprecated(new = "showGiottoSpatNetworks")

  showGiottoSpatNetworks(...)

}


#' @title showGiottoSpatGrids
#' @name showGiottoSpatGrids
#' @description Prints the available spatial grids that are attached to the Giotto object
#' @param gobject giotto object
#' @param nrows number of rows to print
#' @return prints name of available spatial grids
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatGrids = function(gobject,
                               nrows = 4) {

  # import boxchars
  ch = box_chars()

  # 1. check input
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. get availability matrix
  available_data = list_spatial_grids(gobject = gobject)
  if(is.null(available_data)) {
    cat('No available spatial grids \n')
  } else {

    # 3.1 Set up object printouts
    objPrints = objRows = list()
    for(obj_i in seq(nrow(available_data))) {

      # Get object
      dataObj = get_spatialGrid(gobject = gobject,
                                spat_unit = available_data$spat_unit[[obj_i]],
                                name = available_data$name[[obj_i]],
                                return_grid_Obj = FALSE)

      # Collect object prints
      objRows[[obj_i]] = nrow(dataObj)

      objPrints[[obj_i]] =
        dataObj[1:if(nrows <= objRows[[obj_i]]) nrows else objRows[[obj_i]],] %>%
        print %>%
        capture.output

    }

    # object printblock edits
    objPrints = lapply(objPrints, function(x) paste0(ch$s, x)) # Add indent
    objPrints = lapply(objPrints, function(x) paste(x, collapse = ('\n'))) # Linearize print

    # Append to availability table
    available_data$values = unlist(objPrints)

    # 3.2 Setup general prints
    available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
    if(!is.null(available_data$feat_type)) {
      available_data$feat_type = paste0('Feature type "', available_data$feat_type, '"') # Check to be deprecated
    } else warning('Only networks from the deprecated nesting will be shown')
    available_data$name = paste0('S4 spatialGridObj "', available_data$name, '"')
    for(obj_i in seq(nrow(available_data))) {
      available_data$name[[obj_i]] = paste0(available_data$name[[obj_i]],
                                            ch$s ,'(', objRows[[obj_i]], ' rows)')
    }

    # 4. Print information
    print_leaf(level_index = 1,
               availableDT = available_data,
               inherit_last = TRUE,
               indent = '')

  }

}


#' @title Show Spatial Grids
#' @name showGrids
#' @inheritDotParams showGiottoSpatGrids
#' @seealso \code{\link{showGiottoSpatGrids}}
#' @export
showGrids = function(...) {

  .Deprecated(new = "showGiottoSpatGrids")

  showGiottoSpatGrids(...)

}


#' @title showGiottoImageNames
#' @name showGiottoImageNames
#' @description Prints the available giotto images that are attached to the Giotto object
#' @param gobject a giotto object
#' @return prints names of available giotto image objects
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoImageNames = function(gobject) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_images(gobject = gobject)
  if(is.null(available_data)) cat('No available images \n')

  for(image_type in unique(available_data$img_type)) {

    cat('Image type:', image_type, '\n\n')

    for(image_name in available_data[available_data$img_type == image_type,]$name) {

      cat('--> Name:', image_name, '\n')

    }
    cat('\n')
  }
}



# List functions ####

#' @title list_giotto_data
#' @name list_giotto_data
#' @description list the available data within specified giotto object slot
#' @param gobject giotto object
#' @param slot giotto object slot of interest (e.g. "expression", "spatial_locs", etc.)
#' @param ... additional params to pass
#' @return names and locations of data within giotto object slot
#' @keywords internal
list_giotto_data = function(gobject = NULL,
                            slot = NULL,
                            ...) {

  if(slot == 'expression') return(list_expression(gobject = gobject,...))
  if(slot == 'cell_metadata') return(list_cell_metadata(gobject = gobject,...))
  if(slot == 'feat_metadata') return(list_feat_metadata(gobject = gobject,...))
  if(slot == 'spatial_locs') return(list_spatial_locations(gobject = gobject,...))
  if(slot == 'spatial_enrichment') return(list_spatial_enrichments(gobject = gobject,...))
  if(slot == 'dimension_reduction') return(list_dim_reductions(gobject = gobject,...))
  if(slot == 'nn_network') return(list_nearest_networks(gobject = gobject,...))
  if(slot == 'spatial_info') return(list_spatial_info(gobject = gobject))
  if(slot == 'feat_info') return(list_feature_info(gobject = gobject))
  if(slot == 'spatial_network') return(list_spatial_networks(gobject = gobject,...))
  if(slot == 'spatial_grid') return(list_spatial_grids(gobject = gobject,...))
  if(slot == 'images') return(list_images_names(gobject = gobject, img_type = 'image'))
  if(slot == 'largeImages') return(list_images_names(gobject = gobject, img_type = 'largeImage'))

}


#' @title list_expression
#' @name list_expression
#' @description lists the available matrices
#' @inheritParams data_access_params
#' @return names and locations of available matrices as data.table. col order matters.
list_expression = function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL) {

  availableExpr = data.table()
  for(spatial_unit in names(gobject@expression)) {
    for(feature_type in names(gobject@expression[[spatial_unit]])) {
      for(mat_i in names(gobject@expression[[spatial_unit]][[feature_type]])) {
        availableExpr = rbind(availableExpr,
                              list(spat_unit = spatial_unit,
                                   feat_type = feature_type,
                                   name = mat_i))
      }
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableExpr$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableExpr$feat_type == feat_type else feat_type_subset = TRUE

  availableExpr = availableExpr[spat_unit_subset & feat_type_subset,]

  # return data.table (NULL if empty)
  if(nrow(availableExpr) == 0) return(NULL)
  else return(availableExpr)
}



#' @title list_expression_names
#' @name list_expression_names
#' @description lists the available matrices names for a given spatial unit and feature type
#' @inheritParams data_access_params
#' @return vector with names of available matrices
list_expression_names = function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL) {

  if(is.null(spat_unit)) stop('spat_unit must be given\n')
  if(is.null(feat_type)) stop('feat_type must be given\n')

  expression_names = names(gobject@expression[[spat_unit]][[feat_type]])

  return(expression_names)
}



#' @title List cell ID names
#' @name list_cell_id_names
#' @description lists the available cell id names. In effect, these names are the
#' spat_units and poly info in the gobject
#' @inheritParams data_access_params
#' @return vector with names of available sets of cell_IDs
list_cell_id_names = function(gobject) {
  return(names(gobject@cell_ID))
}


#' @title List feat ID names
#' @name list_feat_id_names
#' @description lists the available feat id names In effect, these names are the
#' feat_types and feature info in the gobject
#' @inheritParams data_access_params
#' @return vector with names of available sets of feat_IDs
list_feat_id_names = function(gobject) {
  return(names(gobject@feat_ID))
}


#' @title list_cell_metadata
#' @name list_cell_metadata
#' @description lists the available cell metadata.
#' @inheritParams data_access_params
#' @return names and locations of available cell metadata as data.table
list_cell_metadata = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              return_uniques = FALSE) {

  availableCMet = data.table()
  uniques = list()
  for(spatial_unit in names(gobject@cell_metadata)) {
    uniques$spat_unit = c(uniques$spat_unit, spatial_unit)
    for(feature_type in names(gobject@cell_metadata[[spatial_unit]])) {
      uniques$feat_type = c(uniques$feat_type, feature_type)
      availableCMet = rbind(availableCMet,
                            list(spat_unit = spatial_unit,
                                 feat_type = feature_type))
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableCMet$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableCMet$feat_type == feat_type else feat_type_subset = TRUE

  availableCMet = availableCMet[spat_unit_subset & feat_type_subset,]

  if(!isTRUE(return_uniques)) {
    # return data.table (NULL if empty)
    if(nrow(availableCMet) == 0) return(NULL)
    else return(availableCMet)
  } else {
    return(lapply(uniques, unique))
  }
}



#' @title list_feat_metadata
#' @name list_feat_metadata
#' @description lists the available feature metadata
#' @inheritParams data_access_params
#' @return names and locations of available feature metadata as data.table
list_feat_metadata = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              return_uniques = FALSE) {

  availableFMet = data.table()
  uniques = list()
  for(spatial_unit in names(gobject@feat_metadata)) {
    uniques$spat_unit = c(uniques$spat_unit, spatial_unit)
    for(feature_type in names(gobject@feat_metadata[[spatial_unit]])) {
      uniques$feat_type = c(uniques$feat_type, feature_type)
      availableFMet = rbind(availableFMet,
                            list(spat_unit = spatial_unit,
                                 feat_type = feature_type))
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableFMet$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableFMet$feat_type == feat_type else feat_type_subset = TRUE

  availableFMet = availableFMet[spat_unit_subset & feat_type_subset,]

  if(!isTRUE(return_uniques)) {
    # return data.table (NULL if empty)
    if(nrow(availableFMet) == 0) return(NULL)
    else return(availableFMet)
  } else {
    return(lapply(uniques, unique))
  }
}



#' @title list_spatial_locations
#' @name list_spatial_locations
#' @description shows the available spatial locations
#' @inheritParams data_access_params
#' @return names and locations of available data.table as data.table
list_spatial_locations = function(gobject,
                                  spat_unit = NULL,
                                  return_uniques = FALSE) {

  availableSpatLocs = data.table()
  uniques = list()
  for(spatial_unit in names(gobject@spatial_locs)) {
    uniques$spat_unit = c(uniques$spat_unit, spatial_unit)
    for(spatloc_name in names(gobject@spatial_locs[[spatial_unit]])) {
      uniques$name = c(uniques$name, spatloc_name)
      if(inherits(slot(gobject, 'spatial_locs')[[spatial_unit]][[spatloc_name]], c('data.table', 'spatLocsObj'))) {
        availableSpatLocs = rbind(availableSpatLocs,
                                  list(spat_unit = spatial_unit,
                                       name = spatloc_name))
      }
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableSpatLocs$spat_unit == spat_unit else spat_unit_subset = TRUE

  availableSpatLocs = availableSpatLocs[spat_unit_subset,]

  if(!isTRUE(return_uniques)) {
    if(nrow(availableSpatLocs) == 0) return(NULL)
    else return(availableSpatLocs)
  } else {
    return(lapply(uniques, unique))
  }
}




#' @title list_spatial_locations_names
#' @name list_spatial_locations_names
#' @description lists the available spatial location names for a given spatial unit
#' @inheritParams data_access_params
#' @return vector with names of available spatial locations
list_spatial_locations_names = function(gobject,
                                        spat_unit = NULL) {

  if(is.null(spat_unit)) stop('spat_unit must be given\n')

  spatlocs_names = names(gobject@spatial_locs[[spat_unit]])

  return(spatlocs_names)
}



#' @title list_spatial_enrichments
#' @name list_spatial_enrichments
#' @description return the available spatial enrichment results
#' @inheritParams data_access_params
#' @return names and locations of available data as data.table
list_spatial_enrichments = function(gobject,
                                    spat_unit = NULL,
                                    feat_type = NULL) {

  availableSpatEnr = data.table()

  for(spatial_unit in names(gobject@spatial_enrichment)) {

    for(feature_type in names(gobject@spatial_enrichment[[spatial_unit]])) {

      for(spatenr_name in names(gobject@spatial_enrichment[[spatial_unit]][[feature_type]])) {

        availableSpatEnr = rbind(availableSpatEnr,
                                 list(spat_unit = spatial_unit,
                                      feat_type = feature_type,
                                      name = spatenr_name))
      }

    }

  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableSpatEnr$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableSpatEnr$feat_type == feat_type else feat_type_subset = TRUE

  availableSpatEnr = availableSpatEnr[spat_unit_subset & feat_type_subset,]

  if(nrow(availableSpatEnr) == 0) return(NULL)
  else return(availableSpatEnr)
}





#' @title list_spatial_enrichments_names
#' @name list_spatial_enrichments_names
#' @description returns the available spatial enrichment names for a given spatial unit
#' @inheritParams data_access_params
#' @return vector of names for available spatial enrichments
list_spatial_enrichments_names = function(gobject,
                                          spat_unit = NULL,
                                          feat_type = NULL) {

  if(is.null(spat_unit)) stop('spat_unit must be given\n')
  if(is.null(feat_type)) stop('feat_type must be given\n')

  spatenr_names = names(gobject@spatial_enrichment[[spat_unit]][[feat_type]])

  return(spatenr_names)
}





#' @title list_dim_reductions
#' @name list_dim_reductions
#' @description return the available dimension reductions
#' @inheritParams data_access_params
#' @param data_type "cells" or "feats" data used in dim reduction
#' @param dim_type dimensional reduction method (e.g. "pca", "umap")
#' @return names and locations of dimension reduction as a data.table
list_dim_reductions = function(gobject,
                               data_type = NULL,
                               spat_unit = NULL,
                               feat_type = NULL,
                               dim_type = NULL) {

  availableDimRed = data.table()
  for(dataType in names(slot(gobject, 'dimension_reduction'))) {
    for(spatUnit in names(slot(gobject, 'dimension_reduction')[[dataType]])) {
      for(featType in names(slot(gobject, 'dimension_reduction')[[dataType]][[spatUnit]])) {
        for(dimType in names(slot(gobject, 'dimension_reduction')[[dataType]][[spatUnit]][[featType]])) {
          for(subType in names(slot(gobject, 'dimension_reduction')[[dataType]][[spatUnit]][[featType]][[dimType]])) {
            if(inherits(slot(gobject, 'dimension_reduction')[[dataType]][[spatUnit]][[featType]][[dimType]][[subType]], 'dimObj')) {
              availableDimRed = rbind(availableDimRed,
                                      list(data_type = dataType,
                                           spat_unit = spatUnit,
                                           feat_type = featType,
                                           dim_type = dimType,
                                           name = subType))
            }
          }
        }
      }
    }
  }

  # check if a specific category is desired
  if(!is.null(data_type)) data_type_subset = availableDimRed$data_type == data_type else data_type_subset = TRUE
  if(!is.null(spat_unit)) spat_unit_subset = availableDimRed$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableDimRed$feat_type == feat_type else feat_type_subset = TRUE
  if(!is.null(dim_type)) dimred_type_subset = availableDimRed$dim_type == dim_type else dimred_type_subset = TRUE

  availableDimRed = availableDimRed[data_type_subset & spat_unit_subset & feat_type_subset & dimred_type_subset,]

  # NULL if there is no data
  if(nrow(availableDimRed) == 0) return(NULL)
  else return(availableDimRed)
}



#' @title list_dim_reductions_names
#' @name list_dim_reductions_names
#' @description return the available dimension reductions object names
#' @inheritParams data_access_params
#' @param data_type cells or feats dim reduction
#' @param dim_type dimensional reduction type (method)
#' @return names of dimension reduction object
#' @details function that can be used to find which names have been used
list_dim_reductions_names = function(gobject,
                                     data_type = 'cells',
                                     spat_unit = NULL,
                                     feat_type = NULL,
                                     dim_type = NULL) {

  if(is.null(spat_unit)) stop('spat_unit must be given\n')
  if(is.null(feat_type)) stop('feat_type must be given\n')
  if(is.null(dim_type)) stop('dim_type must be given\n')

  dim_red_object_names = names(slot(gobject, 'dimension_reduction')[[data_type]][[spat_unit]][[feat_type]][[dim_type]])

  return(dim_red_object_names)
}



#' @title list_nearest_networks
#' @name list_nearest_networks
#' @description return the available nearest neighbor network information
#' @inheritParams data_access_params
#' @param nn_type nearest neighbor method (e.g. "sNN", "kNN")
#' @return names and locations of nearest neighbor networks as a data.table
list_nearest_networks = function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 nn_type = NULL,
                                 return_uniques = FALSE) {

  availableNN = data.table()
  uniques = list()
  for(spatUnit in names(slot(gobject, 'nn_network'))) {
    uniques$spat_unit = c(uniques$spat_unit, spatUnit)
    for(featType in names(slot(gobject, 'nn_network')[[spatUnit]])) {
      uniques$feat_type = c(uniques$feat_type, featType)
      for(nnType in names(slot(gobject, 'nn_network')[[spatUnit]][[featType]])) {
        uniques$nn_type = c(uniques$nn_type, nnType)
        for(nnNet in names(slot(gobject, 'nn_network')[[spatUnit]][[featType]][[nnType]])) {
          uniques$name = c(uniques$name, nnNet)
          if(inherits(slot(gobject, 'nn_network')[[spatUnit]][[featType]][[nnType]][[nnNet]], c('igraph', 'nnData')))
          availableNN = rbind(availableNN,
                              list(spat_unit = spatUnit,
                                   feat_type = featType,
                                   nn_type = nnType,
                                   name = nnNet))
        }
      }
    }
  }

  # **To be deprecated**
  # nn network has gained feat_type nesting. Check back one layer
  if(!all(uniques$nn_type %in% availableNN$nn_type)) {
    # Check for vaid igraph objects at lower nesting
    availableNN_old = data.table()
    for(spatUnit in names(slot(gobject, 'nn_network'))) {
      for(nnType in names(slot(gobject, 'nn_network')[[spatUnit]])) {
        for(nnNet in names(slot(gobject, 'nn_network')[[spatUnit]][[nnType]])) {
          if(inherits(slot(gobject, 'nn_network')[[spatUnit]][[nnType]][[nnNet]], 'igraph')) {
            availableNN_old = rbind(availableNN_old,
                                    list(spat_unit = spatUnit,
                                         nn_type = nnType,
                                         name = nnNet))
          }
        }
      }
    }
    if(nrow(availableNN_old > 0)) {
      message('Deprecated nesting found within nn_network slot:')
      print(availableNN_old)
      warning('Deprecated nesting discovered within Giotto nn_network slot. Consider remaking the object or changing the nesting to the suggested.')

      for(net in seq(nrow(availableNN_old))) {
        # Assign default feature type for each spat_unit
        featType = set_default_feat_type(gobject,
                                         spat_unit = availableNN_old$spat_unit[[net]])
        # Place object in new location
        gobject@nn_network[[availableNN_old$spat_unit[[net]]]][[featType]][[availableNN_old$nn_type[[net]]]][[availableNN_old$name[[net]]]] =
          gobject@nn_network[[availableNN_old$spat_unit[[net]]]][[availableNN_old$nn_type[[net]]]][[availableNN_old$name[[net]]]]
        # Remove old object so that it is not detected by this list function
        gobject@nn_network[[availableNN_old$spat_unit[[net]]]][[availableNN_old$nn_type[[net]]]][[availableNN_old$name[[net]]]] = NULL
      }
      # Recursive call on new nesting structure
      message('Suggested new nesting:')
      availableNN_suggest = list_nearest_networks(gobject)
      print(availableNN_suggest)
      cat('\n')
      availableNN = availableNN_old
    }
  } # **deprecation end**

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableNN$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableNN$feat_type == feat_type else feat_type_subset = TRUE
  if(!is.null(nn_type)) nn_type_subset = availableNN$nn_type == nn_type else nn_type_subset = TRUE

  availableNN = availableNN[spat_unit_subset & feat_type_subset & nn_type_subset,]

  if(!isTRUE(return_uniques)) {
    # NULL if there is no data
    if(nrow(availableNN) == 0) return(NULL)
    else return(availableNN)
  } else {
    return(lapply(uniques, unique))
  }
}



#' @title list_nearest_networks_names
#' @name list_nearest_networks_names
#' @description return the available nearest neighbor network object names
#' @inheritParams data_access_params
#' @param nn_type nearest neighbor method (e.g. "sNN", "kNN")
#' @return names of nearest neighbor network object
#' @details function that can be used to find which names have been used
list_nearest_networks_names = function(gobject,
                                       spat_unit = NULL,
                                       feat_type = NULL,
                                       nn_type = NULL) {

  if(is.null(spat_unit)) stop('spat_unit must be given\n')
  if(is.null(feat_type)) stop('feat_type must be given\n')
  if(is.null(nn_type)) stop('nn_type must be given\n')

  nn_object_names = names(slot(gobject, 'nn_network')[[spat_unit]][[feat_type]][[nn_type]])

  return(nn_object_names)
}



#' @title list_spatial_info
#' @name list_spatial_info
#' @description return the available giotto spatial polygon information
#' @param gobject giotto object
#' @return names of available spatial polygon information
list_spatial_info = function(gobject) {

  availableSpatInfo = data.table()
  for(info in names(gobject@spatial_info)) {
    availableSpatInfo = rbind(availableSpatInfo,
                              list(spat_info = info))
  }

  if(nrow(availableSpatInfo) == 0) return(NULL)
  else return(availableSpatInfo)
}




#' @title list_spatial_info_names
#' @name list_spatial_info_names
#' @description return the available names for giotto spatial polygon information
#' @param gobject giotto object
#' @return vector with names of available spatial polygon information
list_spatial_info_names = function(gobject) {

  spat_info_names = names(gobject@spatial_info)

  return(spat_info_names)
}



#' @title list_feature_info
#' @name list_feature_info
#' @description return the available giotto spatial feature information
#' @param gobject giotto object
#' @return names of available feature information
list_feature_info = function(gobject) {

  availableFeatInfo = data.table()
  for(info in names(gobject@feat_info)) {
    availableFeatInfo = rbind(availableFeatInfo,
                              list(feat_info = info))
  }

  if(nrow(availableFeatInfo) == 0) return(NULL)
  else return(availableFeatInfo)
}


#' @title list_feature_info_names
#' @name list_feature_info_names
#' @description return the available names for giotto feature information
#' @param gobject giotto object
#' @return vector with names of available feature information
list_feature_info_names = function(gobject) {

  feat_info_names = names(gobject@feat_info)

  return(feat_info_names)
}



#' @title list_spatial_networks
#' @name list_spatial_networks
#' @description return the available spatial networks that are attached to the Giotto object
#' @inheritParams data_access_params
#' @return data.table of names and locations of available spatial networks, col order matters or list of unique nestings
list_spatial_networks = function(gobject,
                                 spat_unit = NULL,
                                 return_uniques = FALSE) {

  availableSpatNetworks = data.table()
  uniques = list()
  for(spatial_unit in names(gobject@spatial_network)) {
    uniques$spat_unit = c(uniques$spat_unit, spatial_unit)
    for(spat_network_name in names(gobject@spatial_network[[spatial_unit]])) {
      uniques$name = c(uniques$name, spat_network_name)
      if(inherits(gobject@spatial_network[[spatial_unit]][[spat_network_name]], 'spatialNetworkObj'))
      availableSpatNetworks = rbind(availableSpatNetworks,
                                    list(spat_unit = spatial_unit,
                                         name = spat_network_name))
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableSpatNetworks$spat_unit == spat_unit else spat_unit_subset = TRUE

  availableSpatNetworks = availableSpatNetworks[spat_unit_subset,]

  if(!isTRUE(return_uniques)) {
    if(nrow(availableSpatNetworks) == 0) return(NULL)
    else return(availableSpatNetworks)
  } else {
    return(lapply(uniques, unique))
  }

}


#' @title list_spatial_networks_names
#' @name list_spatial_networks_names
#' @description return the available names for giotto feature information
#' @inheritParams data_access_params
#' @return vector with names of available feature information
list_spatial_networks_names = function(gobject,
                                       spat_unit = NULL) {

  if(is.null(spat_unit)) stop('spat_unit must be given\n')

  spat_network_names = names(gobject@spatial_network[[spat_unit]])

  return(spat_network_names)
}




#' @title list_spatial_grids
#' @name list_spatial_grids
#' @description return the available spatial grids that are attached to the Giotto object
#' @inheritParams data_access_params
#' @return data.table of names and locations of available spatial grids. col order matters
list_spatial_grids = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              return_uniques = FALSE) {

  availableSpatGrids = data.table()
  uniques = list()
  for(spatial_unit in names(gobject@spatial_grid)) {
    uniques$spat_unit = c(uniques$spat_unit, spatial_unit)
    for(feature_type in names(gobject@spatial_grid[[spatial_unit]])) {
      uniques$feat_type = c(uniques$feat_type, feature_type)
      for(grid_names in names(gobject@spatial_grid[[spatial_unit]][[feature_type]])) {
        uniques$name = c(uniques$name, grid_names)
        if(inherits(gobject@spatial_grid[[spatial_unit]][[feature_type]][[grid_names]], 'spatialGridObj')) {
          availableSpatGrids = rbind(availableSpatGrids,
                                     list(spat_unit = spatial_unit,
                                          feat_type = feature_type,
                                          name = grid_names))
        }
      }
    }
  }

  # **To be deprecated**
  # spatial_grid has gained feat_type nesting. Check back one layer
  if(!all(uniques$name %in% availableSpatGrids$name)) {
    # Check for valid spatialGridObj objects at lower nesting
    availableSpatGrids_old = data.table()
    for(spatial_unit in names(gobject@spatial_grid)) {
      for(grid_names in names(gobject@spatial_grid[[spatial_unit]])) {
        if(inherits(gobject@spatial_grid[[spatial_unit]][[grid_names]], 'spatialGridObj')) {
          availableSpatGrids_old = rbind(availableSpatGrids_old,
                                         list(spat_unit = spatial_unit,
                                              name = grid_names))
        }
      }
    }
    if(nrow(availableSpatGrids_old > 0)) {
      message('Deprecated nesting discovered within spatial_grid slot:')
      print(availableSpatGrids_old)
      warning('Deprecated nesting discovered within Giotto spatial_grid slot. Consider remaking the object or changing the nesting to the suggested.')
      for(grid in seq(nrow(availableSpatGrids_old))) {
        # Assign default feature type for each spat_unit
        feature_type = set_default_feat_type(gobject,
                                             spat_unit = availableSpatGrids_old$spat_unit[[grid]])
        # Place object in new location
        gobject@spatial_grid[[availableSpatGrids_old$spat_unit[[grid]]]][[feature_type]][[availableSpatGrids_old$name[[grid]]]] =
          gobject@spatial_grid[[availableSpatGrids_old$spat_unit[[grid]]]][[availableSpatGrids_old$name[[grid]]]]
        # Remove old object so that it is not detected by this list function
        gobject@spatial_grid[[availableSpatGrids_old$spat_unit[[grid]]]][[availableSpatGrids_old$name[[grid]]]] = NULL
      }
      # Recursive call on new nesting structure
      message('Suggested new nesting:')
      availableSpatGrids_suggest = list_spatial_grids(gobject)
      print(availableSpatGrids_suggest)
      cat('\n')
      availableSpatGrids = availableSpatGrids_old
    }
  } # **deprecation end**


  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableSpatGrids$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableSpatGrids$feat_type == feat_type else feat_type_subset = TRUE

  availableSpatGrids = availableSpatGrids[spat_unit_subset & feat_type_subset,]

  if(!isTRUE(return_uniques)) {
    if(nrow(availableSpatGrids) == 0) return(NULL)
    else return(availableSpatGrids)
  } else {
    return(lapply(uniques, unique))
  }
}




#' @title list_spatial_grids_names
#' @name list_spatial_grids_names
#' @description return the available spatial grids name for a given spatial unit that are attached to the Giotto object
#' @inheritParams data_access_params
#' @param return_uniques return unique nesting names (ignores if final object exists/is correct class)
#' @return vector with names of available spatial grids names
list_spatial_grids_names = function(gobject,
                                    spat_unit = NULL,
                                    feat_type = NULL,
                                    return_uniques = FALSE) {

  if(is.null(spat_unit)) stop('spat_unit must be given\n')
  if(is.null(feat_type)) stop('feat_type must be given\n')

  spat_grid_names = names(gobject@spatial_grid[[spat_unit]][[feat_type]])

  return(spat_grid_names)
}


#' @title list_images
#' @name list_images
#' @description Prints the available giotto images that are attached to the Giotto object
#' @param gobject giotto object
#' @param img_type "image" or "largeImage"
#' @return data.table of giotto image names attached to the giotto object
list_images = function(gobject,
                       img_type = NULL) {

  availableImages = data.table()

  g_image_names = names(slot(gobject, 'images'))
  g_limage_names = names(slot(gobject, 'largeImages'))

  for(image_type in c('image', 'largeImage')) {
    if(image_type == 'image') {
      for(name in g_image_names) {
        availableImages = rbind(availableImages,
                                list(img_type = image_type,
                                     name = name))
      }
    }
    if(image_type == 'largeImage') {
      for(name in g_limage_names) {
        availableImages = rbind(availableImages,
                                list(img_type = image_type,
                                     name = name))
      }
    }
  }

  # check if a specific category is desired
  if(!is.null(img_type)) img_type_subset = availableImages$img_type == img_type else img_type_subset = TRUE

  availableImages = availableImages[img_type_subset,]

  # NULL if there is no data
  if(nrow(availableImages) == 0) return(NULL)
  else return(availableImages)
}



#' @title list_images_names
#' @name list_images_names
#' @description return the available image names for a given image type that are attached to the Giotto object
#' @param gobject a giotto object
#' @param img_type "image" or "largeImage"
#' @return vector with names of available image names
list_images_names = function(gobject,
                             img_type) {

  if(!img_type %in% c('image', 'largeImage')) stop('img_type must be either "image" or "largeImage\n"')

  if(img_type == 'image') img_names = names(gobject@images)
  if(img_type == 'largeImage') img_names = names(gobject@largeImages)

  return(img_names)
}


