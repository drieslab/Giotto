
#' @include classes.R
#' @importFrom methods setMethod
NULL

# -------------------------------------------------- #




# giotto ####
#' @noRd
#' @keywords internal
setMethod('initialize', signature('giotto'), function(.Object, ...) {

  .Object = methods::callNextMethod()
  .Object = updateGiottoObject(.Object)


  # DT vars
  spat_unit = feat_type = NULL

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
        warning(wrap_txt(
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
  methods::validObject(.Object)



  .Object

})










# Virtual Classes ####






## metaData ####
setMethod('initialize', 'metaData',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@metaDT)) {
              .Object@metaDT = data.table::setalloccol(.Object@metaDT)
            }
            .Object
          })



## enrData ####
setMethod('initialize', 'enrData',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@enrichDT)) {
              .Object@enrichDT = data.table::setalloccol(.Object@enrichDT)
            }
            .Object
          })





## spatNetData ####
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





## coordDataDT ####
setMethod('initialize', 'coordDataDT',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@coordinates)) {
              .Object@coordinates = data.table::setalloccol(.Object@coordinates)
            }
            .Object
          })






## spatGridData ####
setMethod('initialize', 'spatGridData',
          function(.Object, ...) {
            .Object = methods::callNextMethod()
            # prepare DT for set by reference
            if(!is.null(.Object@gridDT)) {
              .Object@gridDT = data.table::setalloccol(.Object@gridDT)
            }
            .Object
          })


















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



