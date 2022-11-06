
# Commonly used params
# Use @inheritParams data_access when documenting these parameters
#' @title Giotto object data accessors
#' @description Access or examine slots within the giotto object
#' @name data_access
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param return_uniques return unique nesting names (ignores if final object exists/is correct class)
#' @param output what format in which to get information (e.g. "data.table")
#' @keywords internal
NULL

## Get and set functions to get and set values in one of the giotto class slots ##

#%%%%% NOTE: python and instructions accessors are currently in giotto.R %%%%%#








## cell_ID slot ####

#' @title Get cell IDs for a given spatial unit
#' @name get_cell_id
#' @inheritParams data_access
#' @description Data for each spatial unit is expected to agree on a single set of cell_IDs
#' that are shared across any feature types. These cell_IDs are stored within the
#' giotto object's \code{cell_ID} slot. Getters and setters for this slot directly
#' retrieve (get) or replace (set) this slot.
#' @return character vector of cell_IDs
#' @seealso set_cell_id
#' @family functions to set data in giotto object
#' @keywords internal
get_cell_id = function(gobject,
                       spat_unit = NULL) {

  if(!inherits(gobject, 'giotto')) stop('Only Giotto Objects are supported for this function.')
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  return(slot(gobject, 'cell_ID')[[spat_unit]])

}


#' @title Set cell IDs for a given spatial unit
#' @name set_cell_id
#' @inheritParams data_access
#' @param cell_IDs character vector of cell IDs to set
#' @description Data for each spatial unit is expected to agree on a single set of cell_IDs
#' that are shared across any feature types. These cell_IDs are stored within the
#' giotto object's \code{cell_ID} slot. Getters and setters for this slot directly
#' retrieve (get) or replace (set) this slot. Set \code{cell_IDs} \code{NULL} in order
#' to delete the entry.
#' @seealso get_cell_id
#' @return giotto object with set cell_ID slot
#' @family functions to set data in giotto object
#' @keywords internal
set_cell_id = function(gobject,
                       spat_unit = NULL,
                       cell_IDs) {
  if(!inherits(gobject, 'giotto')) stop('Only Giotto Objects are supported for this function.')
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  if(!is.null(cell_IDs)) {
    if(!inherits(cell_IDs, 'character')) stop('cell_IDs must be a character vector.')
  }

  slot(gobject, 'cell_ID')[[spat_unit]] = cell_IDs

  return(gobject)
}



## feat_ID slot ####

#' @title Get feat IDs for a given feature type
#' @name get_feat_id
#' @inheritParams data_access
#' @description Across a single modality/feature type, all feature information is
#' expected to share a single set of feat_IDs. These feat_IDs are stored within
#' the giotto object's \code{feat_ID} slot. Getters and setters for this slot
#' directly (get) or replace (set) this slot.
#' @seealso set_feat_id
#' @family functions to set data in giotto object
#' @keywords internal
get_feat_id = function(gobject,
                       feat_type = NULL) {
  if(!inherits(gobject, 'giotto')) stop('Only Giotto Objects are supported for this function.')
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = NULL)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  return(slot(gobject, 'feat_ID')[[feat_type]])

}


#' @title Set feat IDs for a given feature type
#' @name set_feat_id
#' @inheritParams data_access
#' @description Across a single modality/feature type, all feature information is
#' expected to share a single set of feat_IDs. These feat_IDs are stored within
#' the giotto object's \code{feat_ID} slot. Getters and setters for this slot
#' directly (get) or replace (set) this slot.
#' @seealso get_feat_id
#' @return giotto object with set cell_ID slot
#' @family functions to set data in giotto object
#' @keywords internal
set_feat_id = function(gobject,
                       feat_type = NULL,
                       feat_IDs) {
  if(!inherits(gobject, 'giotto')) stop('Only Giotto Objects are supported for this function.')
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = NULL)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  if(!is.null(feat_IDs)) {
    if(!inherits(feat_IDs, 'character')) stop('feat_IDs must be a character vector.')
  }

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
#' @inheritParams data_access
#' @param output return as either 'data.table' or 'cellMetaObj'
#' @description Get cell metadata from giotto object
#' @seealso pDataDT
#' @keywords internal
get_cell_metadata = function(gobject,
                             spat_unit = NULL,
                             feat_type = NULL,
                             output = c('data.table', 'cellMetaObj')) {

  output = match.arg(output, choices = c('data.table', 'cellMetaObj'))

  # 1. Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 2. Find object - note that metadata objects do not have names
  cellMeta = gobject@cell_metadata[[spat_unit]][[feat_type]]

  # 3. Return as desired object type
  if(output == 'cellMetaObj') {
    if(inherits(cellMeta, 'data.table')) { # ** TO BE DEPRECATED **
      cellMeta = new('cellMetaObj',
                     metaDT = cellMeta,
                     col_desc = NA_character_, # no info
                     spat_unit = spat_unit,
                     feat_type = feat_type,
                     provenance = spat_unit) # assumed
    }

    if(inherits(cellMeta, 'list') | is.null(cellMeta)) stop('metadata referenced does not exist.')
    if(!inherits(cellMeta, 'cellMetaObj')) stop('metadata referenced is not data.table or cellMetaObj')

    # return cellMetaObj
    return(cellMeta)
  } else if(output == 'data.table') {
    if(inherits(cellMeta, 'cellMetaObj')) {
      cellMeta = slot(cellMeta, 'metaDT')
    }

    if(inherits(cellMeta, 'list') | is.null(cellMeta)) stop('metadata referenced does not exist.')
    if(!inherits(cellMeta, 'data.table')) stop('metadata referenced is not data.table or cellMetaObj')

    # return data.table
    return(cellMeta)
  }

}





#' @title Set cell metadata
#' @name set_cell_metadata
#' @description Function to set cell metadata information into giotto object
#' @inheritParams data_access
#' @param provenance provenance information to set
#' @param metadata cellMetaObj or data.table containing cell metadata. Setting NULL will remove
#' the object. Passing 'initialize' will reset the object.
#' @param verbose be verbose
#' @return giotto object
#' @family functions to set data in giotto object
set_cell_metadata = function(gobject,
                             spat_unit = NULL,
                             feat_type = NULL,
                             provenance = NULL,
                             metadata,
                             verbose = TRUE) {

  if(!inherits(gobject, 'giotto')) stop("Only Giotto Objects are supported for this function.")

  # 1. determine if user input was supplied
  if(is.null(spat_unit)) nospec_unit = TRUE
  if(is.null(feat_type)) nospec_feat = TRUE

  # 2. set spat unit/ feat type if needed
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 3.1 if input is NULL, remove object
  if(is.null(metadata)) {
    if(isTRUE(verbose)) message('NULL passed to metadata.\n Removing specified metadata.')
    gobject@cell_metadata[[spat_unit]][[feat_type]] = NULL
    return(gobject)
  }

  # 3.2 if input is 'initialize', RESET/reinitialize object
  if(inherits(metadata, 'character')) {
    if(metadata == 'initialize') {
      if(isTRUE(verbose)) message('Initializing specified metadata.')
      gobject@cell_metadata[[spat_unit]][[feat_type]] = new('cellMetaObj',
                                                            metaDT = data.table::data.table(cell_ID = get_cell_id(gobject,
                                                                                                                  spat_unit = spat_unit)),
                                                            col_desc = c(cell_ID = 'cell-specific unique ID value'),
                                                            spat_unit = spat_unit,
                                                            feat_type = feat_type,
                                                            provenance = if(is.null(provenance)) spat_unit else provenance)
      return(gobject)
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
      return(gobject)
    }

    # 4.3 otherwise...
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
    if(isTRUE(verbose)) message('Cell metadata for spat_unit "', spat_unit, '" and feat_type "',
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
  gobject@cell_metadata[[spat_unit]][[feat_type]] = metadata
  return(gobject)

}


## feature metadata slot ####

#' @title Get feature metadata
#' @name get_feature_metadata
#' @inheritParams data_access
#' @param output return as either 'data.table' or 'featMetaObj'
#' @description Get feature metadata from giotto object
#' @seealso fDataDT
#' @keywords internal
get_feature_metadata = function(gobject,
                                spat_unit = NULL,
                                feat_type = NULL,
                                output = c('data.table', 'featMetaObj')) {

  output = match.arg(output, choices = c('data.table', 'featMetaObj'))

  # 1. Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 2. Find object - note that metadata objects do not have names
  featMeta = gobject@feat_metadata[[spat_unit]][[feat_type]]

  # 3. Return as desired object type
  if(output == 'featMetaObj') {
    if(inherits(featMeta, 'data.table')) { # ** TO BE DEPRECATED **
      featMeta = new('featMetaObj',
                     metaDT = featMeta,
                     col_desc = NA_character_, # no info
                     spat_unit = spat_unit,
                     feat_type = feat_type,
                     provenance = spat_unit) # assumed
    }

    if(inherits(featMeta, 'list') | is.null(featMeta)) stop('metadata referenced does not exist.')
    if(!inherits(featMeta, 'featMetaObj')) stop('metadata referenced is not data.table or featMetaObj')

    # return featMetaObj
    return(featMeta)
  } else if(output == 'data.table') {
    if(inherits(featMeta, 'featMetaObj')) {
      featMeta = slot(featMeta, 'metaDT')
    }

    if(inherits(featMeta, 'list') | is.null(featMeta)) stop('metadata referenced does not exist.')
    if(!inherits(featMeta, 'data.table')) stop('metadata referenced is not data.table or featMetaObj')

    # return data.table
    return(featMeta)
  }

}




#' @title Set feature metadata
#' @name set_feature_metadata
#' @description Function to set feature metadata into giotto object
#' @inheritParams data_access
#' @param metadata featMetaObj or data.table object containing feature metadata.
#' Setting NULL will remove the object.
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @return giotto object
#' @family functions to set data in giotto object
set_feature_metadata = function(gobject,
                                spat_unit = NULL,
                                feat_type = NULL,
                                provenance = NULL,
                                metadata,
                                verbose = TRUE) {

  if(!inherits(gobject, 'giotto')) stop("Only Giotto Objects are supported for this function.")

  # 1. determine if user input was supplied
  if(is.null(spat_unit)) nospec_unit = TRUE
  if(is.null(feat_type)) nospec_feat = TRUE

  # 2. set spat unit/ feat type if needed
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 3.1 if input is NULL, remove object
  if(is.null(metadata)) {
    if(isTRUE(verbose)) message('NULL passed to metadata.\n Removing specified metadata.')
    gobject@feat_metadata[[spat_unit]][[feat_type]] = NULL
    return(gobject)
  }

  # 3.2 if input is 'initialize', RESET/reinitialize object
  if(inherits(metadata, 'character')) {
    if(metadata == 'initialize') {
      if(isTRUE(verbose)) message('Initializing specified metadata.')
      gobject@feat_metadata[[spat_unit]][[feat_type]] = new('featMetaObj',
                                                            metaDT = data.table::data.table(feat_ID = get_feat_id(gobject,
                                                                                                                  feat_type = feat_type)),
                                                            col_desc = c(feat_ID = 'feature-specific unique ID value'),
                                                            spat_unit = spat_unit,
                                                            feat_type = feat_type,
                                                            provenance = if(is.null(provenance)) spat_unit else provenance)
      return(gobject)
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

  } else {

    # 4.2 if nested list structure, extract spat_unit/feat_type
    if(inherits(metadata, 'list')) {
      featMetaObj_list = read_feature_metadata(gobject,
                                               metadata = metadata,
                                               provenance = if(is.null(provenance)) spat_unit else provenance)
      # recursively run
      for(obj_i in seq_along(featMetaObj_list)) {
        # (provenance info set during prev. step)
        gobject = set_feature_metadata(gobject,
                                       metadata = featMetaObj_list[[obj_i]])
      }
      return(gobject)
    }

    # 4.3 otherwise...
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
    if(isTRUE(verbose)) message('Feat metadata for spat_unit "', spat_unit, '" and feat_type "',
                                feat_type, '" already exists and will be replaced with new metadata.')
  }

  gobject@feat_metadata[[spat_unit]][[feat_type]] = metadata
  return(gobject)


}


## expression values slot ####

#' @title  Get expression values
#' @name  get_expression_values
#' @description Function to get expression values from giotto object
#' @inheritParams data_access
#' @param values expression values to extract (e.g. "raw", "normalized", "scaled")
#' @param output what object type to retrieve the expression as. Currently either
#' 'Matrix' for dgeMatrix or dgCMatrix or 'exprObj' are allowed.
#' @return expression matrix
#' @family expression accessor functions
#' @family functions to get data from giotto object
#' @export
get_expression_values = function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 values,
                                 output = c('Matrix', 'exprObj')) {


  output = match.arg(output, choices = c('Matrix', 'exprObj'))

  # 1. Set feat_type and spat_unit

  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 2. Find object

  potential_values = names(gobject@expression[[spat_unit]][[feat_type]])

  ## special cases for giotto standard pipeline
  if(values == 'scaled' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('run first scaling (& normalization) step')
  } else if(values == 'normalized' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('run first normalization step')
  } else if(values == 'custom' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('first add custom expression matrix')
  }

  if(!values %in% potential_values) stop("The spatial unit ", spat_unit ," for expression matrix ", feat_type, " and with name ","'", values, "'"," can not be found \n")

  # 3. Get object in desired format

  expr_values = gobject@expression[[spat_unit]][[feat_type]][[values]]

  if(output == 'exprObj') {
    if(inherits(expr_values, c('dgeMatrix'))) { # ** TO BE DEPRECATED **
      expr_values = new('exprObj',
                        name = values,
                        exprMat = expr_values,
                        sparse = FALSE,
                        spat_unit = spat_unit,
                        feat_type = feat_type,
                        provenance = spat_unit, # assumed
                        misc = NULL)
    } else if(inherits(expr_values, c('dgCMatrix'))) { # ** TO BE DEPRECATED **
      expr_values = new('exprObj',
                        name = values,
                        exprMat = expr_values,
                        sparse = TRUE,
                        spat_unit = spat_unit,
                        feat_type = feat_type,
                        provenance = spat_unit, # assumed
                        misc = NULL)
    } else {
      expr_values = new('exprObj', # general catch
                        name = values,
                        exprMat = expr_values,
                        sparse = NA,            # unknown
                        spat_unit = spat_unit,
                        feat_type = feat_type,
                        provenance = spat_unit, # assumed
                        misc = NULL)
    }

    if(!inherits(expr_values, 'exprObj')) stop('Cannot convert to exprObj')

    # return exprObj
    return(expr_values)

  } else if(output == 'Matrix') {

    if(inherits(expr_values, 'exprObj')) expr_values = slot(expr_values, 'exprMat')

    # return 'Matrix'
    if(inherits(expr_values, c('dgeMatrix', 'dgCMatrix'))) {
      return(expr_values)
    } else {
      warning('matrix type conversions not yet implemented. Returning the current expression matrix')
      return(expr_values)
    }

  }



}


#' @title select_expression_values
#' @name select_expression_values
#' @inheritDotParams get_expression_values
#' @seealso \code{\link{get_expression_values}}
#' @keywords internal
select_expression_values = function(...) {

  .Deprecated(new = "get_expression_values")

  get_expression_values(...)

}


#' @title  Set expression values
#' @name  set_expression_values
#' @description Function to set expression values for giotto object
#' @inheritParams data_access
#' @param name name for the expression slot
#' @param values exprObj or matrix of expression values. If NULL, then the object
#' will be removed.
#' @param verbose be verbose
#' @return giotto object
#' @family expression accessor functions
#' @family functions to set data in giotto object
#' @export
set_expression_values = function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 name = 'test',
                                 values,
                                 verbose = TRUE) {

  # 1. Determine user inputs
  if(is.null(spat_unit)) nospec_unit = TRUE
  if(is.null(feat_type)) nospec_feat = TRUE
  if(is.null(match.call()$name)) nospec_name = TRUE

  # 2. Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 3. if input is NULL, remove object
  if(is.null(values)) {
    if(isTRUE(verbose)) message('NULL passed to values param.\n Removing specified expression')
    gobject@expression[[spat_unit]][[feat_type]][[name]] = values
    return(gobject)
  }

  # 4. import data from S4 if available
  if(inherits(values, 'exprObj')) {

    if(isTRUE(nospec_unit)) {
      if(!is.na(slot(values, 'spat_unit'))) spat_unit = slot(values, 'spat_unit')
    } else {
      slot(values, 'spat_unit') = spat_unit
    }
    if(isTRUE(nospec_feat)) {
      if(!is.na(slot(values, 'feat_type'))) feat_type = slot(values, 'feat_type')
    } else {
      slot(values, 'feat_type') = feat_type
    }
    if(isTRUE(nospec_name)) {
      if(!is.na(slot(values, 'name'))) name = slot(values, 'name')
    } else {
      slot(values, 'name') = name
    }

  } else {
    values = new('exprObj',
                 name = name,
                 exprMat = values,
                 sparse = NA,          # unknown
                 spat_unit = spat_unit,
                 feat_type = feat_type,
                 provenance = spat_unit, # assumed
                 misc = NULL)
  }

  ## 5. check if specified name has already been used
  potential_names = list_expression_names(gobject, spat_unit = spat_unit, feat_type = feat_type)
  if(name %in% potential_names) {
    if(isTRUE(verbose)) message(name, ' already exist and will be replaced with new values \n')
  }

  ## 6. update and return giotto object
  gobject@expression[[spat_unit]][[feat_type]][[name]] = values
  return(gobject)

}





## spatial locations slot ####

#' @title Get spatial locations
#' @name get_spatial_locations
#' @description Function to get a spatial location data.table
#' @inheritParams data_access
#' @param spat_loc_name name of spatial locations (defaults to first name in spatial_locs slot, e.g. "raw")
#' @param return_spatlocs_Obj return spatLocsObj (default = FALSE)
#' @param copy_obj whether to copy/duplicate when getting the object (default = TRUE)
#' @return data.table with coordinates or spatLocsObj depending on \code{return_spatlocs_Obj}
#' @family spatial location data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatial_locations = function(gobject,
                                 spat_unit = NULL,
                                 spat_loc_name = NULL,
                                 return_spatlocs_Obj = FALSE,
                                 copy_obj = TRUE) {

  if(!is.null(spat_unit)) {
    potential_spat_unit = names(slot(gobject, 'spatial_locs'))
    if(!spat_unit %in% potential_spat_unit) stop('No spatial locations for spatial unit "', spat_unit, '" exist in spatial_locs slot.')
  }

  # spatial locations
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

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
    if(!isTRUE(return_spatlocs_Obj)) { # case if return as data.table
      if(inherits(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]], 'spatLocsObj')) { # if is spatLocsObj
        if(isTRUE(copy_obj)) spatloc = data.table::copy(slot(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]], 'coordinates')) # copy
        else spatloc = slot(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]], 'coordinates') # no copy
      } else { # to be deprecated | if is data.table
        if(isTRUE(copy_obj)) spatloc = data.table::copy(slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]]) # copy
        else spatloc = slot(gobject, 'spatial_locs')[[spat_unit]][[spat_loc_name]] # no copy
      }
    } else { # case if return as spatLocsObj
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
    }
    return(spatloc)
  } else { # case if not found
    stop("The spatial locations with name ","'", spat_loc_name, "'"," can not be found \n")
  }
}



#' @title select_spatial_locations
#' @name select_spatial_locations
#' @inheritDotParams get_spatial_locations
#' @seealso \code{\link{get_spatial_locations}}
#' @keywords internal
select_spatial_locations = function(...) {

  .Deprecated(new = "get_spatial_locations")

  get_spatial_locations(...)

}




#' @title Set spatial locations
#' @name set_spatial_locations
#' @description Function to set a spatial location slot
#' @inheritParams data_access
#' @param spatlocs spatial locations (accepts either \code{data.table} or
#' \code{spatLocsObj})
#' @param spat_loc_name name of spatial locations, default "raw"
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
                                 verbose = TRUE) {

  # 1. determine if input was supplied to spat_unit and spat_loc_name
  if(is.null(spat_unit)) nospec_unit = TRUE
  if(is.null(match.call()$spat_loc_name)) nospec_name = TRUE

  # 2. spatial locations
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  # 3. remove if input is NULL
  if(is.null(spatlocs)) {
    if(isTRUE(verbose)) message('NULL passed to spatlocs.\n Removing specified spatial locations.')
    gobject@spatial_locs[[spat_unit]][[spat_loc_name]] = spatlocs
    return(gobject)
  }

  # 4. import data
  # if is spatLocsObj, decide if to use tagged spat_unit and name info
  if(inherits(spatlocs, 'spatLocsObj')) {
    if(isTRUE(nospec_unit)) { # case if no input - prioritize tagged info
      if(!is.na(slot(spatlocs, 'spat_unit'))) spat_unit = slot(spatlocs, 'spat_unit')
      else slot(spatlocs, 'spat_unit') = spat_unit
    } else { # case if input given - use input
      slot(spatlocs, 'spat_unit') = spat_unit
    }
    if(isTRUE(nospec_name)) { # case if no input - prioritize tagged info
      if(!is.na(slot(spatlocs, 'name'))) spat_loc_name = slot(spatlocs, 'name')
      else slot(spatlocs, 'name') = spat_loc_name
    } else { # case if input given - use input
      slot(spatlocs, 'name') = spat_loc_name
    }

  } else {
    spatlocs = new('spatLocsObj',
                   name = spat_loc_name,
                   spat_unit = spat_unit,
                   coordinates = spatlocs,
                   provenance = spat_unit) # assume
  }

  # 5. check if specified name has already been used
  potential_names = list_spatial_locations_names(gobject, spat_unit = spat_unit)
  if(spat_loc_name %in% potential_names) {
    if(isTRUE(verbose)) {
    message(spat_loc_name, ' already exists and will be replaced with new spatial locations \n')
    }
  }

  # 6. update and return giotto object
  gobject@spatial_locs[[spat_unit]][[spat_loc_name]] = spatlocs
  return(gobject)

}






## dimension reduction slot ####

#' @title Get dimension reduction
#' @name get_dimReduction
#' @inheritParams data_access
#' @param reduction reduction on cells or features (e.g. "cells", "feats")
#' @param reduction_method reduction method (e.g. "pca", "umap", "tsne")
#' @param name name of reduction results
#' @param return_dimObj return full dimension object result. Default = FALSE
#' @description Function to get a dimension reduction object
#' @return dim reduction coordinates (default) or dim reduction object
#' @family dimensional reduction data accessor functions
#' @family functions to get data from giotto object
#' @export
get_dimReduction = function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            reduction = c('cells', 'feats'),
                            reduction_method = c('pca', 'umap', 'tsne'),
                            name = 'pca',
                            return_dimObj = FALSE) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

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
  if(return_dimObj == TRUE) {
    return(reduction_res)
  } else {
    return(slot(reduction_res, 'coordinates'))
  }

}


#' @title select_dimReduction
#' @name select_dimReduction
#' @inheritDotParams get_dimReduction
#' @seealso \code{\link{get_dimReduction}}
#' @keywords internal
select_dimReduction = function(...) {

  .Deprecated(new = "get_dimReduction")

  get_dimReduction(...)

}


#' @title Set dimension reduction
#' @name set_dimReduction
#' @description Function to set a dimension reduction slot
#' @inheritParams data_access
#' @param reduction reduction on cells or features
#' @param reduction_method reduction method (e.g. "pca")
#' @param name name of reduction results
#' @param dimObject dimension object result to set
#' @return giotto object
#' @family dimensional reduction data accessor functions
#' @family functions to set data in giotto object
#' @export
set_dimReduction <- function(gobject,
                             spat_unit = NULL,
                             feat_type = NULL,
                             reduction = c('cells', 'genes'),
                             reduction_method = c('pca', 'umap', 'tsne'),
                             name = 'pca',
                             dimObject) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## 1. check if specified name has already been used
  potential_names = names(slot(gobject, 'dimension_reduction')[[reduction]][[spat_unit]][[feat_type]][[reduction_method]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new dimension reduction object \n')
  }

  ## TODO: 2. check input for dimension reduction object
  if(!inherits(dimObject, 'dimObj')) stop('Object to set must be a giotto S4 "dimObj"\n')
  silent = validObject(dimObject) # variable hides TRUE print

  ## 3. update and return giotto object
  slot(gobject, 'dimension_reduction')[[reduction]][[spat_unit]][[feat_type]][[reduction_method]][[name]] = dimObject

  return(gobject)

}







## nearest neighbor network slot ####

#' @title Get nearest network
#' @name get_NearestNetwork
#' @description Get a NN-network from a Giotto object
#' @inheritParams data_access
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
                              output = c('igraph', 'data.table', 'nnNetObj')) {

  output = match.arg(arg = output, choices = c('igraph', 'data.table', 'nnNetObj'))

  # 1.  Set feat_type and spat_unit

  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 2 Find the object

  # **to be deprecated - check for old nesting**
  if(is.null(names(gobject@nn_network[[spat_unit]][[feat_type]])[[1]])) { # If gobject has nothing for this feat_type
    available = list_nearest_networks(gobject,
                                      spat_unit = spat_unit,
                                      nn_type = nn_network_to_use)
    if(!is.null(available)) {
      if(nrow(available) > 0 & is.null(available$feat_type)) { # If ANY old nesting objects are discovered (only reports old nestings if detected)
        if(is.null(network_name)) igraph_object = gobject@nn_network[[spat_unit]][[available$nn_type[[1]]]][[available$name[[1]]]]
        else igraph_object = gobject@nn_network[[spat_unit]][[available$nn_type[[1]]]][[network_name]]
        if(inherits(igraph_object, 'igraph')) {
          if(is.null(nn_network_to_use)) message('The NN network type was not specified, default to the first: "', available$nn_type[[1]],'"')
          if(is.null(network_name)) message('The NN network name was not specified, default to the first: "', available$name[[1]],'"')
          ## convert igraph to data.table
          if(output == 'data.table') {
            igraph_object = data.table::as.data.table(igraph::get.data.frame(x = igraph_object))
            return(igraph_object)
          }
          return(igraph_object)
        } else {
          stop('There is currently no nearest-neighbor network created for
           spatial unit: "', spat_unit, '" and feature type "', feat_type, '".
           First run createNearestNetwork()\n')
        }
      }
    } else {
      stop('There is currently no nearest-neighbor network created for
           spatial unit: "', spat_unit, '" and feature type "', feat_type, '".
           First run createNearestNetwork()\n')
    }
  } # **deprecation end**


  # automatic nearest network selection
  if(is.null(nn_network_to_use)) {
    nn_network_to_use = names(gobject@nn_network[[spat_unit]][[feat_type]])[[1]]
    if(is.null(nn_network_to_use)) {

      stop('There is currently no nearest-neighbor network created for
           spatial unit: "', spat_unit, '" and feature type "', feat_type, '".
           First run createNearestNetwork()\n')
    } else {
      message('The NN network type was not specified, default to the first: "', nn_network_to_use,'"')
    }
  }

  if(is.null(network_name)) {
    network_name = names(gobject@nn_network[[spat_unit]][[feat_type]][[nn_network_to_use]])[[1]]
    if(is.null(network_name)) {
      stop('There is currently no nearest-neighbor network built for spatial unit: "', spat_unit,
           '" feature type: "', feat_type, '" and network type: "', nn_network_to_use,'"\n')
    }else {
      message('The NN network name was not specified, default to the first: "', network_name,'"')
    }
  }

  # 3. get object in desired format

  nnNet = gobject@nn_network[[spat_unit]][[feat_type]][[nn_network_to_use]][[network_name]]
  if(is.null(nnNet)) {
    stop('nn_network_to_use: "', nn_network_to_use, '" or network_name: "', network_name, '" does not exist.
          Create a nearest-neighbor network first')
  }

  if(output == 'nnNetObj') {
    if(inherits(nnNet, 'igraph')) { # ** TO BE DEPRECATED **
      nnNet = new('nnNetObj',
                  name = network_name,
                  nn_type = nn_network_to_use,
                  igraph = nnNet,
                  spat_unit = spat_unit,
                  feat_type = feat_type,
                  provenance = spat_unit, # assumed if nnNet is igraph
                  misc = NULL)
    }

    if(!inherits(nnNet, 'nnNetObj')) stop('Specified nnNet is neither igraph nor nnNetObj')

    # return nnNetObj
    return(nnNet)

  } else if(output == 'igraph' | output == 'data.table') {
    if(inherits(nnNet, 'nnNetObj')) { # ** TO BE DEPRECATED ** Will always be assumed to be the case moving forward
      nnNet = slot(nnNet, 'igraph')
    }

    if(output == 'igraph') return(nnNet) # return igraph

    if(output == 'data.table') {
      nnNet = data.table::as.data.table(igraph::get.data.frame(x = nnNet))
      return(nnNet) # return data.table
    }
  }

}


#' @title Extract nearest network
#' @name extractNearestNetwork
#' @inheritDotParams get_NearestNetwork
#' @seealso \code{\link{get_NearestNetwork}}
#' @keywords internal
extractNearestNetwork = function(...) {

  .Deprecated(new = "get_NearestNetwork")

  get_NearestNetwork(...)

}


#' @title Select nearest network
#' @name select_NearestNetwork
#' @inheritDotParams get_NearestNetwork
#' @seealso \code{\link{get_NearestNetwork}}
#' @keywords internal
select_NearestNetwork = function(...) {

  .Deprecated(new = "get_NearestNetwork")

  get_NearestNetwork(...)

}


#' @title Set nearest network
#' @name set_NearestNetwork
#' @description Set a NN-network for a Giotto object
#' @inheritParams data_access
#' @param nn_network_to_use "kNN" or "sNN"
#' @param network_name name of NN network to be used
#' @param nn_network nnNetObj or igraph nearest network object. Data.table not
#' yet supported.
#' @param verbose be verbose
#' @return giotto object
#' @family expression space nearest network accessor functions
#' @family functions to set data in giotto object
#' @export
set_NearestNetwork = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              nn_network,
                              verbose = TRUE) {

  # 1. determine user input
  if(is.null(spat_unit)) nospec_unit = TRUE
  if(is.null(feat_type)) nospec_feat = TRUE
  if(is.null(match.call()$nn_network_to_use)) nospec_net = TRUE
  if(is.null(match.call()$network_name)) nospec_name = TRUE

  # 2. Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # 3. If input is null, remove object
  if(is.null(nn_network)) {
    if(isTRUE(verbose)) message('NULL passed to nn_network.\n Removing specified nearest neighbor network.')
    gobject@nn_network[[spat_unit]][[feat_type]][[nn_network_to_use]][[network_name]] = nn_network
    return(gobject)
  }

  # 4. import data from S4 if available
  if(inherits(nn_network, 'nnNetObj')) {

    if(isTRUE(nospec_unit)) {
      if(!is.na(slot(nn_network, 'spat_unit'))) spat_unit = slot(nn_network, 'spat_unit')
    } else {
      slot(nn_network, 'spat_unit') = spat_unit
    }
    if(isTRUE(nospec_feat)) {
      if(!is.na(slot(nn_network, 'feat_type'))) feat_type = slot(nn_network, 'feat_type')
    } else {
      slot(nn_network, 'feat_type') = feat_type
    }
    if(isTRUE(nospec_net)) {
      if(!is.na(slot(nn_network, 'nn_type'))) nn_network_to_use = slot(nn_network, 'nn_type')
    } else {
      slot(nn_network, 'nn_type') = nn_network_to_use
    }
    if(isTRUE(nospec_name)) {
      if(!is.na(slot(nn_network, 'name'))) network_name = slot(nn_network, 'name')
    } else {
      slot(nn_network, 'name') = network_name
    }

  } else {
    # TODO convert to igraph if data.table class

    nn_network = new('nnNetObj',
                     name = network_name,
                     nn_type = nn_network_to_use,
                     igraph = nn_network,
                     spat_unit = spat_unit,
                     feat_type = feat_type,
                     provenance = spat_unit, # assumed
                     misc = NULL)
  }


  ## 5. check if specified name has already been used
  potential_names = list_nearest_networks_names(gobject,
                                                spat_unit = spat_unit,
                                                feat_type = feat_type,
                                                nn_type = nn_network_to_use)
  if(network_name %in% potential_names) {
    if(isTRUE(verbose)) message('"',network_name, '" already exists and will be replaced with new nearest neighbor network')
  }

  ## 6. update and return giotto object
  gobject@nn_network[[spat_unit]][[feat_type]][[nn_network_to_use]][[network_name]] = nn_network
  return(gobject)

}






## spatial network slot ####

#' @title Get spatial network
#' @name get_spatialNetwork
#' @description Function to get a spatial network
#' @inheritParams data_access
#' @param name name of spatial network
#' @param return_network_Obj return network object (default = FALSE)
#' @family spatial network data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatialNetwork <- function(gobject,
                               spat_unit = NULL,
                               name = NULL,
                               return_network_Obj = FALSE) {



  # Set spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  # set name to get if missing
  if(is.null(name)) name = names(slot(gobject, 'spatial_network')[[spat_unit]])[[1]]

  # check if given name is present
  if (!is.element(name, names(slot(gobject, 'spatial_network')[[spat_unit]]))){
    message = sprintf("spatial network %s has not been created. Returning NULL.
                      check which spatial networks exist with showGiottoSpatNetworks()\n", name)
    warning(message)
    return(NULL)
  }else{
    networkObj = slot(gobject, 'spatial_network')[[spat_unit]][[name]]

    # S3 backwards compatibility
    if(!isS4(networkObj)) networkObj = S3toS4spatNetObj(networkObj)
    silent = validObject(networkObj) # Variable used to hide TRUE print

    networkDT = slot(networkObj, 'networkDT')
  }

  if (return_network_Obj == TRUE){
    return(networkObj)
  }else{
    return(networkDT)
  }
}


#' @title Select spatial network
#' @name select_spatialNetwork
#' @inheritDotParams get_spatialNetwork
#' @seealso \code{\link{get_spatialNetwork}}
#' @keywords internal
select_spatialNetwork = function(...) {

  .Deprecated(new = "get_spatialNetwork")

  get_spatialNetwork(...)

}


#' @title Set spatial network
#' @name set_spatialNetwork
#' @description Function to set a spatial network
#' @inheritParams data_access
#' @param name name of spatial network
#' @param spatial_network spatial network
#' @return giotto object
#' @family spatial network data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatialNetwork <- function(gobject,
                               spat_unit = NULL,
                               name = NULL,
                               spatial_network) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_network[[spat_unit]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new spatial network \n')
  }

  ## TODO: 2. check input for spatial network
  if(!inherits(spatial_network, 'spatialNetworkObj')) stop('spatial_network to set must be S4 "spatialNetworkObj"\n')
  ## Variable used to hide TRUE prints
  silent = validObject(spatial_network)

  ## 3. update and return giotto object
  slot(gobject, 'spatial_network')[[spat_unit]][[name]] = spatial_network
  return(gobject)

}




## spatial grid slot ####

#' @title Get spatial grid
#' @name get_spatialGrid
#' @description Function to get spatial grid
#' @inheritParams data_access
#' @param name name of spatial grid
#' @param return_grid_Obj return grid object (default = FALSE)
#' @family spatial grid data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatialGrid <- function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            name = NULL,
                            return_grid_Obj = FALSE) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

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


#' @title Select spatial grid
#' @name select_spatialGrid
#' @inheritDotParams get_spatialGrid
#' @seealso \code{\link{get_spatialGrid}}
#' @keywords internal
select_spatialGrid = function(...) {

  .Deprecated(new = "get_spatialGrid")

  get_spatialGrid(...)

}


#' @title Set spatial grid
#' @name set_spatialGrid
#' @description Function to set a spatial grid
#' @inheritParams data_access
#' @param name name of spatial grid
#' @param spatial_grid spatial grid object
#' @return giotto object
#' @family spatial grid data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatialGrid <- function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            name = NULL,
                            spatial_grid) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## 1. check if specified name has already been used
  potential_names = names(slot(gobject, 'spatial_grid')[[spat_unit]][[feat_type]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new spatial grid \n')
  }

  ## TODO: 2. check input for spatial grid
  if(!inherits(spatial_grid, 'spatialGridObj')) stop('spatial_grid to set must be S4 "spatialGridObj"\n')
  silent = validObject(spatial_grid) # Variable only used to hide TRUE prints

  ## 3. update and return giotto object
  slot(gobject, 'spatial_grid')[[spat_unit]][[feat_type]][[name]] = spatial_grid

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


#' @title Select polygon info
#' @name select_polygon_info
#' @inheritDotParams get_polygon_info
#' @seealso \code{\link{get_polygon_info}}
#' @keywords internal
select_polygon_info = function(...) {

  .Deprecated(new = "get_polygon_info")

  get_polygon_info(...)

}



#' @title Set polygon info
#' @name set_polygon_info
#' @description Set giotto polygon spatVector
#' @param gobject giotto object
#' @param polygon_name name of polygons. Default "cell"
#' @param gpolygon giotto polygon
#' @param verbose verbosity
#' @return giotto object
#' @family polygon info data accessor functions
#' @family functions to set data in giotto object
#' @export
set_polygon_info = function(gobject,
                            polygon_name = 'cell',
                            gpolygon,
                            verbose = TRUE) {



  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_info)
  if(polygon_name %in% potential_names) {
     if(verbose) cat(polygon_name, ' already exist and will be replaced with new giotto polygon \n')
  }

  ## TODO: 2. check input for giotto polygon


  ## 3. update and return giotto object
  gobject@spatial_info[[polygon_name]] = gpolygon
  return(gobject)

}






## feature info ####

#' @title Get feature info
#' @name get_feature_info
#' @description Get giotto points spatVector
#' @inheritParams data_access
#' @family feature info data accessor functions
#' @family functions to get data from giotto object
#' @export
get_feature_info = function(gobject,
                            feat_type = NULL) {

  # specify feat_type
  feat_type = set_default_feat_type(gobject = gobject,
                                    feat_type = feat_type)

  potential_names = names(gobject@feat_info)

  if(!feat_type %in% potential_names) {
    stop('There is no feature information with name ', feat_type, '\n')
  } else {
    feat_info = gobject@feat_info[[feat_type]]@spatVector
    return(feat_info)
  }
}

#' @title Select feature info
#' @name select_feature_info
#' @inheritDotParams get_feature_info
#' @seealso \code{\link{get_feature_info}}
#' @keywords internal
select_feature_info = function(...) {

  .Deprecated(new = "get_feature_info")

  get_feature_info(...)

}


#' @title Set feature info
#' @name set_feature_info
#' @description Set giotto polygon spatVector for features
#' @inheritParams data_access
#' @param gpolygon giotto polygon
#' @return giotto object
#' @family feature info data accessor functions
#' @family functions to set data in giotto object
#' @export
set_feature_info = function(gobject,
                            feat_type = NULL,
                            gpolygon) {

  # specify feat_type
  if(is.null(feat_type)) {
    feat_type = gobject@expression_feat[[1]]
  }

  ## 1. check if specified name has already been used
  potential_names = names(gobject@feat_info)
  if(feat_type %in% potential_names) {
    cat(feat_type, ' already exist and will be replaced with new giotto polygon \n')
  }

  ## TODO: 2. check input for giotto polygon


  ## 3. update and return giotto object
  gobject@feat_info[[feat_type]] = gpolygon
  return(gobject)

}



## spatial enrichment slot ####


#' @title Get spatial enrichment
#' @name get_spatial_enrichment
#' @description Function to get a spatial enrichment data.table
#' @inheritParams data_access
#' @param enrichm_name name of spatial enrichment results. Default "DWLS"
#' @return data.table with fractions
#' @family spatial enrichment data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatial_enrichment <- function(gobject,
                                   spat_unit = NULL,
                                   feat_type = NULL,
                                   enrichm_name = 'DWLS') {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # spatial locations
  # if NULL (not given) and spatial locations have been added, then use first one
  # if NULL (not given) and spatial loactions have NOT been added, then keep NULL
  if(is.null(enrichm_name)) {
    if(!is.null(gobject@spatial_enrichment)) {
      enrichm_name = names(gobject@spatial_enrichment[[spat_unit]][[feat_type]])[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      enrichm_name = NULL
      cat('No spatial enrichment results have been found \n')
      return(NULL)
    }
  }



  potential_names = names(gobject@spatial_enrichment[[spat_unit]][[feat_type]])

  if(enrichm_name %in% potential_names) {
    enr_res = data.table::copy(gobject@spatial_enrichment[[spat_unit]][[feat_type]][[enrichm_name]])
    return(enr_res)
  } else {
    stop("The spatial enrichment result with name ","'", enrichm_name, "'"," can not be found \n")
  }
}


#' @title Set spatial enrichment
#' @name set_spatial_enrichment
#' @description Function to set a spatial enrichment slot
#' @inheritParams data_access
#' @param enrichm_name name of spatial enrichment results. Default "DWLS"
#' @param spatenrichment spatial enrichment results
#' @return giotto object
#' @family spatial enrichment data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatial_enrichment <- function(gobject,
                                   spat_unit = NULL,
                                   feat_type = NULL,
                                   enrichm_name = 'enrichment',
                                   spatenrichment) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_enrichment[[spat_unit]][[feat_type]][[enrichm_name]])
  if(enrichm_name %in% potential_names) {
    cat(enrichm_name, ' already exist and will be replaced with new spatial enrichment results \n')
  }

  ## TODO: 2. check input for spatial locations


  ## 3. update and return giotto object
  gobject@spatial_enrichment[[spat_unit]][[feat_type]][[enrichm_name]] = spatenrichment
  return(gobject)

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
    if(name %in% potential_names) cat(name, ' already exists and will be replaced with new image object \n')
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
    if(name %in% potential_names) cat(name, 'already exists and will be replaced with new image object \n')
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



## Show functions ####

#' @title showGiottoExpression
#' @name showGiottoExpression
#' @description shows the available matrices
#' @param gobject giotto object
#' @param nrows number of rows to print for each matrix
#' @param ncols number of columns to print for each matrix
#' @return prints the name and small subset of available matrices
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoExpression = function(gobject, nrows = 4, ncols = 4) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_expression(gobject = gobject)
  if(is.null(available_data)) cat('No expression data available \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial unit: ', spatial_unit, ' \n\n')

    for(feature_type in unique(available_data[available_data$spat_unit == spatial_unit,]$feat_type)) {

      cat('--> Feature: ', feature_type, ' \n\n')

      for(mat_i in available_data[available_data$spat_unit == spatial_unit & available_data$feat_type == feature_type,]$name) {

        cat('----> Name: ', mat_i, 'matrix: \n')

        dimensions = dim(gobject@expression[[spatial_unit]][[feature_type]][[mat_i]])
        nrows = min(nrows, dimensions[1])
        ncols = min(ncols, dimensions[2])

        print(gobject@expression[[spatial_unit]][[feature_type]][[mat_i]][1:nrows, 1:ncols])
        cat('\n')
      }
    }
  }
}



#' @title showGiottoCellMetadata
#' @name showGiottoCellMetadata
#' @description shows the available cell metadata
#' @param gobject giotto object
#' @param nrows number of rows to print for each metadata
#' @param ncols number of columns to print for each metadata
#' @return prints the name and small subset of available metadata
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoCellMetadata = function(gobject,
                                  nrows = 4,
                                  ncols = 4) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_cell_metadata(gobject = gobject)
  if(is.null(available_data)) cat('No cell metadata available \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial unit: ', spatial_unit, ' \n\n')

    for(feature_type in unique(available_data[available_data$spat_unit == spatial_unit,]$feat_type)) {

      cat('--> Feature: ', feature_type, ' \n\n')

      dimensions = dim(gobject@cell_metadata[[spatial_unit]][[feature_type]])
      nrows = min(nrows, dimensions[1])
      ncols = min(ncols, dimensions[2])

      print(gobject@cell_metadata[[spatial_unit]][[feature_type]][1:nrows, 1:ncols])
      cat('\n')
    }
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

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_spatial_locations(gobject = gobject)
  if(is.null(available_data)) cat('No spatial locations available \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial unit: ', spatial_unit, ' \n\n')

    for(spatlocname in available_data[available_data$spat_unit == spatial_unit,]$name) {
      if(inherits(gobject@spatial_locs[[spatial_unit]][[spatlocname]], 'data.frame')) {
        cat('--> Name: ', spatlocname, ' \n\n')
        print(gobject@spatial_locs[[spatial_unit]][[spatlocname]][1:nrows,])
        cat('\n')
      }
      if(inherits(gobject@spatial_locs[[spatial_unit]][[spatlocname]], 'spatLocsObj')) {
        cat('--> Name: ', spatlocname, ' \n\n')
        print(slot(gobject@spatial_locs[[spatial_unit]][[spatlocname]], 'coordinates')[1:nrows,])
        cat('\n')
      }
    }
  }

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
        print(gobject@spatial_enrichment[[spatial_unit]][[feature_type]][[spatenrichname]][1:nrows,])

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

  # import boxchars
  ch = box_chars()

  # 1. Check inputs
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. Get availability matrix
  available_data = list_dim_reductions(gobject)
  if(is.null(available_data)) cat('No dimensional reductions available \n')

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
                               return_dimObj = FALSE)

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
  available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
  available_data$feat_type = paste0('Feature type "', available_data$feat_type, '"')
  available_data$dim_type = paste0('Dim reduction type "', available_data$dim_type, '"')
  available_data$name = paste0('S4 dimObj "', available_data$name, '" coordinates:')
  for(obj_i in seq(nrow(available_data))) {
    available_data$name[[obj_i]] = paste0(available_data$name[[obj_i]],
                                          ch$s ,'(', objRows[[obj_i]], ' rows ', objCols[[obj_i]], ' cols)')
  }

  # 4. Print information
  for(data_type_red in unique(available_data$data_type)) {
    data_type_subset = available_data$data_type == data_type_red

    if(data_type_red == 'feats') cat('Dim reduction on features:')
    if(data_type_red == 'cells') cat('Dim reduction on cells:')

    cat('\n',
        '-------------------------',
        '\n\n.\n')

    print_leaf(level_index = 2, # skip over dim reduction layer
               availableDT = available_data[data_type == data_type_red],
               inherit_last = TRUE,
               indent = '')

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

  # import boxchars
  ch = box_chars()

  # 1. check input
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. get availability matrix
  available_data = list_nearest_networks(gobject)
  if(is.null(available_data)) cat('No nearest neighbor networks available \n')

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
  available_data$spat_unit = paste0('Spatial unit "', available_data$spat_unit, '"')
  if(!is.null(available_data$feat_type)) {
    available_data$feat_type = paste0('Feature type "', available_data$feat_type, '"')  # Check to be deprecated
  } else warning('Only networks from the deprecated nesting will be shown')
  available_data$nn_type = paste0('NN network type "', available_data$nn_type, '"')
  available_data$name = paste0('S3 igraph "', available_data$name, '"')
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

  # import boxchars
  ch = box_chars()

  # 1. Check input
  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  # 2. Get availability matrix
  available_data = list_spatial_networks(gobject = gobject)
  if(is.null(available_data)) cat('No spatial networks are available \n')

  # 3.1 Set up object printouts
  objPrints = objRows = list()
  for(obj_i in seq(nrow(available_data))) {

    # Get object
    dataObj = get_spatialNetwork(gobject = gobject,
                                 spat_unit = available_data$spat_unit[[obj_i]],
                                 name = available_data$name[[obj_i]],
                                 return_network_Obj = FALSE)

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
  available_data$name = paste0('S4 spatialNetworkObj "', available_data$name, '"')
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
  if(is.null(available_data)) cat('No available spatial grids \n')

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
#' @inheritParams data_access
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
#' @inheritParams data_access
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
#' @inheritParams data_access
#' @return vector with names of available sets of cell_IDs
list_cell_id_names = function(gobject) {
  return(names(gobject@cell_ID))
}


#' @title List feat ID names
#' @name list_feat_id_names
#' @description lists the available feat id names In effect, these names are the
#' feat_types and feature info in the gobject
#' @inheritParams data_access
#' @return vector with names of available sets of feat_IDs
list_feat_id_names = function(gobject) {
  return(names(gobject@feat_ID))
}


#' @title list_cell_metadata
#' @name list_cell_metadata
#' @description lists the available cell metadata.
#' @inheritParams data_access
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

  availableCMet = availableCMet[spat_unit_subset,]

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
#' @inheritParams data_access
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

  availableFMet = availableFMet[spat_unit_subset,]

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
#' @inheritParams data_access
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
  if(!is.null(spat_unit)) {
    availableSpatLocs = availableSpatLocs[availableSpatLocs$spat_unit == spat_unit,]
  }

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
#' @inheritParams data_access
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
#' @inheritParams data_access
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
  if(!is.null(spat_unit)) {
    availableSpatEnr = availableSpatEnr[availableSpatEnr$spat_unit == spat_unit,]
  }

  if(nrow(availableSpatEnr) == 0) return(NULL)
  else return(availableSpatEnr)
}





#' @title list_spatial_enrichments_names
#' @name list_spatial_enrichments_names
#' @description returns the available spatial enrichment names for a given spatial unit
#' @inheritParams data_access
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
#' @inheritParams data_access
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
#' @inheritParams data_access
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
#' @inheritParams data_access
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
          if(inherits(slot(gobject, 'nn_network')[[spatUnit]][[featType]][[nnType]][[nnNet]], 'igraph'))
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
#' @inheritParams data_access
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
#' @inheritParams data_access
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
#' @inheritParams data_access
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
#' @inheritParams data_access
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
#' @inheritParams data_access
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


