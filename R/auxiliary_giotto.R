

## Giotto auxiliary functions ####


#' @title set_default_spat_unit
#' @name set_default_spat_unit
#' @param gobject gobject
#' @param spat_unit spatial unit
#' @keywords internal
set_default_spat_unit = function(gobject,
                                 spat_unit = NULL) {


  # set spatial unit
  if(is.null(spat_unit)) {

    spat_unit = getOption('giotto.spat_unit')
    if(is.null(spat_unit)) {
      if(!is.null(gobject@expression) & length(gobject@expression) > 0) {
        spat_unit = names(gobject@expression)[[1]]
      } else if(!is.null(gobject@spatial_info)){
        spat_unit = names(gobject@spatial_info)[[1]]
      } else {
        warning('No default for spat_unit could be set \n')
      }
    }

  }

  return(spat_unit)

}



#' @title set_default_feat_type
#' @name set_default_feat_type
#' @param gobject gobject
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @keywords internal
set_default_feat_type = function(gobject,
                                 feat_type = NULL,
                                 spat_unit) {


  # set spatial unit
  if(is.null(feat_type)) {

    feat_type = getOption('giotto.feat_type')

    if(is.null(feat_type)) {
      if(!is.null(gobject@expression) & length(gobject@expression) > 0) {
        feat_type = names(gobject@expression[[spat_unit]])[[1]]
        if(is.null(feat_type)) stop('valid spat_unit input needed \n')
      } else if(!is.null(gobject@feat_info)){
        feat_type = names(gobject@feat_info)[[1]]
      } else {
        warning('No default for feat_type could be set \n')
      }
    }

  }

  return(feat_type)
}


#' @title mean_expr_det_test
#' @keywords internal
mean_expr_det_test = function(mymatrix, detection_threshold = 1) {
  mean_expr_detected = unlist(apply(X = mymatrix, MARGIN = 1, FUN = function(x) {
    detected_x = x[x > detection_threshold]
    mean(detected_x)
  }))
}

#' @title libNorm_giotto
#' @keywords internal
libNorm_giotto <- function(mymatrix, scalefactor){

  libsizes = colSums_flex(mymatrix)

  norm_expr = t_flex(t_flex(mymatrix)/ libsizes)*scalefactor

}

#' @title logNorm_giotto
#' @keywords internal
logNorm_giotto = function(mymatrix, base, offset) {

  if(methods::is(mymatrix, 'DelayedMatrix')) {
    mymatrix = log(mymatrix + offset)/log(base)
  } else if(methods::is(mymatrix, 'dgCMatrix')) {
    mymatrix@x = log(mymatrix@x + offset)/log(base) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    mymatrix@x = log(mymatrix@x + offset)/log(base)
  } else {
    mymatrix = log(as.matrix(mymatrix) + offset)/log(base)
  }

  return(mymatrix)
}



#' @title pDataDT
#' @name pDataDT
#' @description show cell metadata
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @return data.table with cell metadata
#' @export
pDataDT <- function(gobject,
                    spat_unit = NULL,
                    feat_type = NULL) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  if(!inherits(gobject, c('ExpressionSet', 'SCESet', 'seurat', 'giotto'))) {
    stop('only works with ExpressionSet (-like) objects')
  }

  if(inherits(gobject, c('ExpressionSet', 'SCESet'))) {
    return(data.table::as.data.table(Biobase::pData(gobject)))
  }
  else if(inherits(gobject, 'giotto')) {
    return(gobject@cell_metadata[[spat_unit]][[feat_type]])
  }
  else if(inherits(gobject, 'seurat')) {
    return(data.table::as.data.table(gobject@meta.data))
  }

}

#' @title fDataDT
#' @name fDataDT
#' @description show gene metadata
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @return data.table with gene metadata
#' @export
fDataDT <- function(gobject,
                    spat_unit = NULL,
                    feat_type = NULL) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  if(!inherits(gobject, c('ExpressionSet', 'SCESet', 'giotto'))) {
    stop('only works with ExpressionSet (-like) objects')
  }
  else if(inherits(gobject, 'giotto')) {
    return(gobject@feat_metadata[[spat_unit]][[feat_type]])
  }
  return(data.table::as.data.table(Biobase::fData(gobject)))

}


#' @title create_average_DT
#' @description calculates average gene expression for a cell metadata factor (e.g. cluster)
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param meta_data_name name of metadata column to use
#' @param expression_values which expression values to use
#' @return data.table with average gene epression values for each factor
#' @keywords internal
create_average_DT <- function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              meta_data_name,
                              expression_values = c('normalized', 'scaled', 'custom')) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # expression values to be used
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_data = get_expression_values(gobject = gobject,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit,
                                    values = values)


  # metadata
  cell_metadata = pDataDT(gobject,
                          feat_type = feat_type,
                          spat_unit = spat_unit)
  myrownames = rownames(expr_data)

  savelist <- list()
  for(group in unique(cell_metadata[[meta_data_name]])) {

    name = paste0('cluster_', group)

    temp = expr_data[, cell_metadata[[meta_data_name]] == group]
    temp_DT = rowMeans_flex(temp)

    savelist[[name]] <- temp_DT
  }

  finalDF = do.call('cbind', savelist)
  rownames(finalDF) = myrownames

  return(as.data.frame(finalDF))
}

#' @title create_average_detection_DT
#' @description calculates average gene detection for a cell metadata factor (e.g. cluster)
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param meta_data_name name of metadata column to use
#' @param expression_values which expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @return data.table with average gene epression values for each factor
#' @keywords internal
create_average_detection_DT <- function(gobject,
                                        feat_type = NULL,
                                        spat_unit = NULL,
                                        meta_data_name,
                                        expression_values = c('normalized', 'scaled', 'custom'),
                                        detection_threshold = 0) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # expression values to be used
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_data = get_expression_values(gobject = gobject,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit,
                                    values = values)

  # metadata
  cell_metadata <- pDataDT(gobject,
                           feat_type = feat_type,
                           spat_unit = spat_unit)
  myrownames <- rownames(expr_data)

  savelist <- list()
  for(group in unique(cell_metadata[[meta_data_name]])) {

    name = paste0('cluster_', group)

    temp = expr_data[, cell_metadata[[meta_data_name]] == group]
    temp = as.matrix(temp)

    if(is.matrix(temp)) {
      temp_DT = rowSums_flex(temp > detection_threshold)/ncol(temp)
    } else {
      temp_DT = as.numeric(temp > detection_threshold)
    }

    savelist[[name]] <- temp_DT
  }

  finalDF = do.call('cbind', savelist)
  rownames(finalDF) = myrownames

  return(as.data.frame(finalDF))
}




#' @title Subset expression data
#' @name subset_expression_data
#' @description Subset expression data from giotto object
#' @keywords internal
subset_expression_data = function(gobject,
                                  filter_bool_feats,
                                  filter_bool_cells,
                                  feat_type,
                                  spat_unit) {


  for(spat_unit_name in names(gobject@expression)) {

    for(feat_type_name in names(gobject@expression[[spat_unit_name]])) {

      if(feat_type_name == feat_type & spat_unit_name == spat_unit) {

        # filter features and cells
        expression_names = names(gobject@expression[[spat_unit_name]][[feat_type_name]])

        for(expr_name in expression_names) {

          # for HDF5Array
          if(methods::is(gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]], 'HDF5Array')) {
            gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]] = DelayedArray::realize(gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]][filter_bool_feats, filter_bool_cells], "HDF5Array")
          }

          gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]] = gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]][filter_bool_feats, filter_bool_cells]
        }

      } else if(feat_type_name == feat_type & spat_unit_name != spat_unit){

        # filter only features, but NOT cells
        expression_names = names(gobject@expression[[spat_unit_name]][[feat_type_name]])

        for(expr_name in expression_names) {

          # for HDF5Array
          if(methods::is(gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]], 'HDF5Array')) {
            gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]] = DelayedArray::realize(gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]][filter_bool_feats, ], "HDF5Array")
          }

          gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]] = gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]][filter_bool_feats, ]
        }

      } else if(feat_type_name != feat_type & spat_unit_name == spat_unit) {

        # filter only cells, but NOT features
        expression_names = names(gobject@expression[[spat_unit_name]][[feat_type_name]])

        for(expr_name in expression_names) {

          # for HDF5Array
          if(methods::is(gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]], 'HDF5Array')) {
            gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]] = DelayedArray::realize(gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]][, filter_bool_cells], "HDF5Array")
          }

          gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]] = gobject@expression[[spat_unit_name]][[feat_type_name]][[expr_name]][, filter_bool_cells]
        }

      } else {

        # don't filter
        expression_names = names(gobject@expression[[spat_unit_name]][[feat_type_name]])

      }

    }


  }

  return(gobject)


}

#' @title Subset spatial locations
#' @name subset_spatial_locations
#' @description Subset location data from giotto object
#' @keywords internal
subset_spatial_locations = function(gobject,
                                    filter_bool_cells,
                                    spat_unit) {

  # only subset cell_ID if the spatial unit is the same (e.g. cell)

  for(spat_unit_name in names(gobject@spatial_locs)) {
    if(spat_unit_name == spat_unit) {
      #print(spat_unit)
      #print(spat_unit_name)
      for(spatlocname in names(gobject@spatial_locs[[spat_unit_name]])) {
        gobject@spatial_locs[[spat_unit_name]][[spatlocname]] = gobject@spatial_locs[[spat_unit_name]][[spatlocname]][filter_bool_cells]
      }
    }
  }
  return(gobject)
}


#' @title Subset cell metadata
#' @name subset_cell_metadata
#' @description Subset cell metadata from giotto object
#' @keywords internal
subset_cell_metadata = function(gobject,
                                feat_type,
                                filter_bool_cells,
                                spat_unit) {

  if(!is.null(gobject@cell_metadata)) {
    # only subset cell_ID for selected spatial unit
    for(feat_type in names(gobject@cell_metadata[[spat_unit]])) {
      gobject@cell_metadata[[spat_unit]][[feat_type]] = gobject@cell_metadata[[spat_unit]][[feat_type]][filter_bool_cells,]
    }
  }

  return(gobject)
}


#' @title Subset feature metadata
#' @name subset_feature_metadata
#' @description Subset feature metadata from giotto object
#' @keywords internal
subset_feature_metadata = function(gobject,
                                   feat_type,
                                   spat_unit,
                                   filter_bool_feats) {


  if(!is.null(gobject@feat_metadata)) {

    for(spat_unit in names(gobject@feat_metadata)) {
      # only subset features of the feat type, but do it for all spaital regions
      gobject@feat_metadata[[spat_unit]][[feat_type]] = gobject@feat_metadata[[spat_unit]][[feat_type]][filter_bool_feats,]
    }

  }

  return(gobject)
}



#' @title Subset spatial network
#' @name subset_spatial_network
#' @description subset spatial networks from giotto object
#' @keywords internal
subset_spatial_network = function(gobject,
                                  spat_unit,
                                  cells_to_keep) {

  # define for data.table [] subset
  to = NULL
  from = NULL
  
  # cell spatial network
  if(!is.null(gobject@spatial_network)) {
    for(spat_unit_name in names(gobject@spatial_network)) {
      if(spat_unit_name == spat_unit) {
        for(network in names(gobject@spatial_network[[spat_unit_name]])) {
          gobject@spatial_network[[spat_unit_name]][[network]]$networkDT =   gobject@spatial_network[[spat_unit_name]][[network]]$networkDT[to %in% cells_to_keep & from %in% cells_to_keep]
        }
      }
    }
  }

  return(gobject)
}



#' @title Subset dimension reduction
#' @name subset_dimension_reduction
#' @description Subset dimension reduction results from giotto object
#' @keywords internal
subset_dimension_reduction = function(gobject,
                                      spat_unit,
                                      feat_type,
                                      cells_to_keep) {


  if(!is.null(gobject@dimension_reduction$cells)) {

    for(spat_unit_name in names(gobject@dimension_reduction[['cells']])) {


      for(feat_type_name in names(gobject@dimension_reduction[['cells']][[spat_unit_name]])) {


        for(dim_method in names(gobject@dimension_reduction[['cells']][[spat_unit_name]][[feat_type_name]])) {

          if(spat_unit_name == spat_unit & feat_type_name == feat_type) {

            dim_red_names = list_dim_reductions_names(gobject = gobject, data_type = 'cells',
                                                      spat_unit = spat_unit_name, feat_type = feat_type_name,
                                                      dim_type = dim_method)

            for(selected_name in dim_red_names) {

              old_coord = get_dimReduction(gobject = gobject,
                                           spat_unit = spat_unit_name,
                                           feat_type = feat_type_name, reduction = 'cells',
                                           reduction_method = dim_method, name = selected_name,
                                           return_dimObj = FALSE)

              new_coord = old_coord[rownames(old_coord) %in% cells_to_keep,]
              gobject@dimension_reduction[['cells']][[spat_unit_name]][[feat_type_name]][[dim_method]][[selected_name]][['coordinates']] = new_coord

            }

          }

        }

      }

    }

  }

  return(gobject)
}




#' @title Subset nearest network
#' @name subset_nearest_network
#' @description Subset nearest network results from giotto object
#' @keywords internal
subset_nearest_network = function(gobject,
                                  spat_unit,
                                  filter_bool_cells) {


  ## nn network ##
  if(!is.null(gobject@nn_network[['cells']])) {

    for(spat_unit_name in names(gobject@nn_network[['cells']])) {

      if(spat_unit_name == spat_unit) {

        for(knn_name in names(gobject@nn_network[['cells']][[spat_unit_name]][['kNN']])) {

          # layout
          old_layout = gobject@nn_network[['cells']][[spat_unit_name]][['kNN']][[knn_name]][['layout']]

          if(!is.null(old_layout)) {
            new_layout = old_layout[filter_bool_cells,]
            gobject@nn_network[['cells']][[spat_unit]][['kNN']][[knn_name]][['layout']] = new_layout
          }

          # igraph object
          old_graph = gobject@nn_network[['cells']][[spat_unit_name]][['kNN']][[knn_name]][['igraph']]
          vertices_to_keep = igraph::V(old_graph)[filter_bool_cells]
          new_subgraph = igraph::subgraph(graph = old_graph, v = vertices_to_keep)
          gobject@nn_network[['cells']][[spat_unit_name]][['kNN']][[knn_name]][['igraph']] = new_subgraph
        }

        for(snn_name in names(gobject@nn_network[['cells']][[spat_unit_name]][['sNN']])) {

          # layout
          old_layout = gobject@nn_network[['cells']][[spat_unit_name]][['sNN']][[snn_name]][['layout']]

          if(!is.null(old_layout)) {
            new_layout = old_layout[filter_bool_cells,]
            gobject@nn_network[['cells']][[spat_unit_name]][['sNN']][[snn_name]][['layout']] = new_layout
          }

          # igraph object
          old_graph = gobject@nn_network[['cells']][[spat_unit_name]][['sNN']][[snn_name]][['igraph']]
          vertices_to_keep = igraph::V(old_graph)[filter_bool_cells]
          new_subgraph = igraph::subgraph(graph = old_graph, v = vertices_to_keep)
          gobject@nn_network[['cells']][[spat_unit_name]][['sNN']][[snn_name]][['igraph']] = new_subgraph
        }

      }
    }

  }

  return(gobject)
}



#' @title Subset spatial enrichment
#' @name subset_spatial_enrichment
#' @description Subset spatial enrichment results from giotto object
#' @keywords internal
subset_spatial_enrichment = function(gobject,
                                     spat_unit,
                                     feat_type,
                                     filter_bool_cells) {

  if(!is.null(gobject@spatial_enrichment)) {
    for(spat_unit_name in names(gobject@spatial_enrichment)) {

      for(feat_type_name in names(gobject@spatial_enrichment[[spat_unit_name]])) {

        if(spat_unit_name == spat_unit & feat_type_name == feat_type) {
          for(spat_enrich_name in names(gobject@spatial_enrichment[[spat_unit_name]][[feat_type_name]])) {

            gobject@spatial_enrichment[[spat_unit_name]][[feat_type_name]][[spat_enrich_name]] = gobject@spatial_enrichment[[spat_unit_name]][[feat_type_name]][[spat_enrich_name]][filter_bool_cells]

          }
        }

      }


    }
  }

  return(gobject)
}






#' @title Subset giotto polygon object
#' @name subset_giotto_polygon_object
#' @description Subset a single giotto polygon object
#' @keywords internal
subset_giotto_polygon_object = function(gpolygon,
                                        cell_ids,
                                        feat_ids,
                                        feat_type) {


  if(!is.null(gpolygon@spatVector)) {
    cell_id_bool = gpolygon@spatVector$poly_ID %in% cell_ids
    gpolygon@spatVector = gpolygon@spatVector[cell_id_bool]
  }

  if(!is.null(gpolygon@spatVectorCentroids)) {
    cell_id_bool = gpolygon@spatVectorCentroids$poly_ID %in% cell_ids
    gpolygon@spatVectorCentroids = gpolygon@spatVectorCentroids[cell_id_bool]
  }

  if(!is.null(gpolygon@overlaps)) {

    for(feat in names(gpolygon@overlaps)) {
      cell_id_bool = gpolygon@overlaps[[feat]]$poly_ID %in% cell_ids
      gpolygon@overlaps[[feat]] = gpolygon@overlaps[[feat]][cell_id_bool]

      if(feat == feat_type) {
        feat_id_bool = gpolygon@overlaps[[feat]]$feat_ID %in% feat_ids
        gpolygon@overlaps[[feat]] = gpolygon@overlaps[[feat]][feat_id_bool]
      }

    }


  }

  return(gpolygon)

}



#' @title Subset spatial info data
#' @name subset_spatial_info_data
#' @description Subset  all spatial info (polygon) data
#' @keywords internal
subset_spatial_info_data = function(spatial_info,
                                    cell_ids,
                                    poly_info = 'cell',
                                    feat_ids,
                                    feat_type = NULL) {


  # set feat type
  if(is.null(feat_type)) {
    feat_type = 'rna'
  }

  res_list = list()
  for(spat_info in names(spatial_info)) {

    cat('for ', spat_info, '\n')

    if(spat_info %in% poly_info) {

      cat('--> ', spat_info, ' found back in polygon layer: ', poly_info, '\n')

      spat_subset = subset_giotto_polygon_object(spatial_info[[spat_info]],
                                                 cell_ids = cell_ids,
                                                 feat_ids = feat_ids,
                                                 feat_type = feat_type)
      #print(spat_subset)
      res_list[[spat_info]] = spat_subset

    } else {

      if(!is.null(spatial_info[[spat_info]]@overlaps)) {

        for(feat in names(spatial_info[[spat_info]]@overlaps)) {

          if(feat == feat_type) {

            feat_id_bool = spatial_info[[spat_info]]@overlaps[[feat]]$feat_ID %in% feat_ids

            spatial_info[[spat_info]]@overlaps[[feat]] = spatial_info[[spat_info]]@overlaps[[feat]][feat_id_bool]
          }
        }
      }

      res_list[[spat_info]] = spatial_info[[spat_info]]
    }


  }
  return(res_list)
}


# subset giotto points

#' @title Subset giotto points object
#' @name subset_giotto_points_object
#' @description Subset a single giotto points object
#' @details Subset on feature ids and on x,y coordinates
#' @keywords internal
subset_giotto_points_object = function(gpoints,
                                       feat_ids = NULL,
                                       x_min = NULL,
                                       x_max = NULL,
                                       y_min = NULL,
                                       y_max = NULL) {

  # define for data.table [] subset
  x = NULL
  y = NULL
  
  if(!is.null(gpoints@spatVector)) {

    if(!is.null(feat_ids)) {
      feat_id_bool = gpoints@spatVector$feat_ID %in% feat_ids
      gpoints@spatVector = gpoints@spatVector[feat_id_bool]
    }

    if(!any(is.null(c(x_min, x_max, y_min, y_max)))) {

      myspatvector = gpoints@spatVector
      spatDT = spatVector_to_dt(myspatvector)

      spatDT_subset = spatDT[x >= x_min & x <= x_max & y >= y_min & y <= y_max]
      myspatvector_subset = dt_to_spatVector_points(dt = spatDT_subset)

      gpoints@spatVector = myspatvector_subset
    }

  }

  return(gpoints)

}



#' @title Subset feature info data
#' @name subset_feature_info_data
#' @description Subset  all spatial feature (points) data
#' @keywords internal
subset_feature_info_data = function(feat_info,
                                    feat_ids,
                                    feat_type = 'rna',
                                    x_min = NULL,
                                    x_max = NULL,
                                    y_min = NULL,
                                    y_max = NULL) {

  res_list = list()
  for(feat in names(feat_info)) {

    if(feat == feat_type) {

      feat_subset = subset_giotto_points_object(feat_info[[feat]],
                                                feat_ids = feat_ids,
                                                x_min = x_min,
                                                x_max = x_max,
                                                y_min = y_min,
                                                y_max = y_max)
      res_list[[feat]] = feat_subset

    } else {
      res_list[[feat]] = feat_info[[feat]]
    }


  }
  return(res_list)
}











#' @title subsetGiotto
#' @description Subsets Giotto object including previous analyses.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type to use
#' @param cell_ids cell IDs to keep
#' @param feat_ids feature IDs to keep
#' @param gene_ids deprecated, use feat_ids
#' @param poly_info polygon information to use
#' @param x_max,x_min,y_max,y_min minimum and maximum x and y coordinates to keep for feature coordinates
#' @param verbose be verbose
#' @param toplevel_params parameters to extract
#' @return giotto object
#' @details Subsets a Giotto object for a specific spatial unit and feature type
#' @export
subsetGiotto <- function(gobject,
                         spat_unit = NULL,
                         feat_type = NULL,
                         cell_ids = NULL,
                         feat_ids = NULL,
                         gene_ids = NULL,
                         poly_info = NULL,
                         x_max = NULL,
                         x_min = NULL,
                         y_max = NULL,
                         y_min = NULL,
                         verbose = TRUE,
                         toplevel_params = 2) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # set poly_info
  if(is.null(poly_info)) {
    poly_info = spat_unit
  }

  ## deprecated arguments
  if(!is.null(gene_ids)) {
    feat_ids = gene_ids
    warning('gene_ids argument is deprecated, use feat_ids argument in the future \n')
  }

  g_cell_IDs = gobject@cell_ID[[spat_unit]]
  g_feat_IDs = gobject@feat_ID[[feat_type]]

  #print(g_cell_IDs[1:5])
  #print(g_feat_IDs[1:5])

  ## filter index
  if(!is.null(cell_ids)) {
    filter_bool_cells = g_cell_IDs %in% cell_ids
  } else filter_bool_cells = g_cell_IDs %in% g_cell_IDs
  if(!is.null(feat_ids)) {
    filter_bool_feats = g_feat_IDs %in% feat_ids
  } else filter_bool_feats = g_feat_IDs %in% g_feat_IDs

  cells_to_keep = g_cell_IDs[filter_bool_cells]
  feats_to_keep = g_feat_IDs[filter_bool_feats]

  #print(cells_to_keep[1:5])
  #print(length(filter_bool_cells))
  #print(feats_to_keep[1:5])
  #print(length(filter_bool_feats))

  if(verbose) cat('completed 1: preparation \n')


  #print(gobject@spatial_locs[[spat_unit]][['raw']])

  ## FILTER ##
  # filter expression data
  gobject = subset_expression_data(gobject = gobject,
                                   filter_bool_feats = filter_bool_feats,
                                   filter_bool_cells = filter_bool_cells,
                                   feat_type = feat_type,
                                   spat_unit = spat_unit)

  if(verbose) cat('completed 2: subset expression data \n')

  #print(gobject@spatial_locs[[spat_unit]][['raw']])

  # filter spatial locations
  #print(spat_unit)
  gobject = subset_spatial_locations(gobject = gobject,
                                     filter_bool_cells = filter_bool_cells,
                                     spat_unit = spat_unit)

  if(verbose) cat('completed 3: subset spatial locations \n')


  # filter cell_ID and gene_ID
  gobject@cell_ID[[spat_unit]] = gobject@cell_ID[[spat_unit]][filter_bool_cells]
  gobject@feat_ID[[feat_type]] = gobject@feat_ID[[feat_type]][filter_bool_feats]


  if(verbose) cat('completed 4: subset cell (spatial units) and feature IDs \n')


  ## cell & feature metadata ##
  # cell metadata
  gobject = subset_cell_metadata(gobject = gobject,
                                  feat_type = feat_type,
                                  filter_bool_cells = filter_bool_cells,
                                  spat_unit = spat_unit)

  if(verbose) cat('completed 5: subset cell metadata \n')

  # feature metadata
  gobject = subset_feature_metadata(gobject = gobject,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit,
                                    filter_bool_feats = filter_bool_feats)

  if(verbose) cat('completed 6: subset feature metadata \n')

  # data.table variables
  to = from = V = NULL

  ## spatial network & grid ##
  # cell spatial network
  gobject = subset_spatial_network(gobject = gobject,
                                   spat_unit = spat_unit,
                                   cells_to_keep = cells_to_keep)


  if(verbose) cat('completed 7: subset spatial network(s) \n')

  # spatial grid
  # need to be recomputed


  ## dimension reduction ##
  # cell dim reduction
  gobject = subset_dimension_reduction(gobject = gobject,
                                       spat_unit = spat_unit,
                                       feat_type = feat_type,
                                       cells_to_keep = cells_to_keep)

  if(verbose == TRUE) cat('completed 8: subsetted dimension reductions \n')


  ## nn network ##
  gobject = subset_nearest_network(gobject = gobject,
                                   spat_unit = spat_unit,
                                   filter_bool_cells = filter_bool_cells)

  if(verbose == TRUE) cat('completed 9: subsetted nearest network(s) \n')


  ## spatial enrichment ##
  gobject = subset_spatial_enrichment(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      filter_bool_cells = filter_bool_cells)

  if(verbose == TRUE) cat('completed 10: subsetted spatial enrichment results \n')

  ## spatial info
  if(!is.null(gobject@spatial_info)) {
    gobject@spatial_info = subset_spatial_info_data(spatial_info = gobject@spatial_info,
                                                    feat_type = feat_type,
                                                    cell_ids = cells_to_keep,
                                                    feat_ids = feats_to_keep,
                                                    poly_info = poly_info)

    if(verbose == TRUE) cat('completed 11: subsetted spatial information data \n')
  }


  ## feature info
  if(!is.null(gobject@feat_info)) {
    gobject@feat_info = subset_feature_info_data(feat_info = gobject@feat_info,
                                                 feat_ids = feats_to_keep,
                                                 feat_type = feat_type,
                                                 x_max = x_max,
                                                 x_min = x_min,
                                                 y_max = y_max,
                                                 y_min = y_min)

    if(verbose == TRUE) cat('completed 12: subsetted spatial feature data \n')
  }



  ## update parameters used ##
  nframes = sys.nframe()
  if(verbose == TRUE) cat('number of frames: ', nframes, '\n')

  parent = sys.parent()
  if(verbose == TRUE) cat('sys parent: ', parent, '\n')

  parameters_info = update_giotto_params(gobject,
                                         description = '_subset',
                                         return_gobject = FALSE,
                                         toplevel = toplevel_params)

  # extra parameters to include
  cells_removed = length(filter_bool_cells[filter_bool_cells==FALSE])
  feats_removed = length(filter_bool_feats[filter_bool_feats==FALSE])

  parameters_list = parameters_info[['plist']]
  update_name = parameters_info[['newname']]

  parameters_list[[update_name]] = c(parameters_list[[update_name]],
                                     'cells removed' = cells_removed,
                                     'feats removed' = feats_removed)
  gobject@parameters = parameters_list


  return(gobject)

}




#' @title Subset by spatial locations
#' @name subsetGiottoLocs
#' @description Subsets Giotto object based on spatial locations
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type to use
#' @param spat_loc_name name of spatial locations to use
#' @param x_max,x_min,y_max,y_min,z_max,z_min minimum and maximum x, y, and z coordinates
#'   to subset to
#' @param poly_info polygon information to use
#' @param return_gobject return Giotto object
#' @param verbose be verbose
#' @return giotto object
#' @details Subsets a Giotto based on spatial locations and for one provided spatial unit
#' if return_gobject = FALSE, then a filtered combined metadata data.table will be returned
#' @export
subsetGiottoLocs = function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            spat_loc_name = NULL,
                            x_max = NULL,
                            x_min = NULL,
                            y_max = NULL,
                            y_min = NULL,
                            z_max = NULL,
                            z_min = NULL,
                            poly_info = 'cell',
                            return_gobject = TRUE,
                            verbose = FALSE) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # Check spatial params
  spatError = NULL
  if(!is.null(x_min) && !is.null(x_max)) if(x_min > x_max) spatError = append(spatError, 'x_max must be larger than x_min \n')
  if(!is.null(y_min) && !is.null(y_max)) if(y_min > y_max) spatError = append(spatError, 'y_max must be larger than y_min \n')
  if(!is.null(z_min) && !is.null(z_max)) if(z_min > z_max) spatError = append(spatError, 'z_max must be larger than z_min \n')
  if(!is.null(spatError)) stop(spatError)


  # function requires spat_loc_name
  if(is.null(spat_loc_name)) {
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = names(gobject@spatial_locs[[spat_unit]])[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      spat_loc_name = NULL
      cat('No spatial locations have been found \n')
      return(NULL)
    }
  }

  comb_metadata = combineMetadata(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  spat_loc_name = spat_loc_name)
  comb_colnames =  colnames(comb_metadata)

  # x spatial dimension
  if('sdimx' %in% comb_colnames) {
    if(is.null(x_max)) x_max = max(comb_metadata[['sdimx']])
    if(is.null(x_min)) x_min = min(comb_metadata[['sdimx']])

    comb_metadata = comb_metadata[get('sdimx') < x_max & get('sdimx') > x_min]
  }

  # y spatial dimension
  if('sdimy' %in% comb_colnames) {
    if(is.null(y_max)) y_max = max(comb_metadata[['sdimy']])
    if(is.null(y_min)) y_min = min(comb_metadata[['sdimy']])

    comb_metadata = comb_metadata[get('sdimy') < y_max & get('sdimy') > y_min]
  }

  # z spatial dimension
  if('sdimz' %in% comb_colnames) {
    if(is.null(z_max)) z_max = max(comb_metadata[['sdimz']])
    if(is.null(z_min)) z_min = min(comb_metadata[['sdimz']])

    comb_metadata = comb_metadata[get('sdimz') < z_max & get('sdimz') > z_min]
  }

  if(return_gobject == TRUE) {

    filtered_cell_IDs = comb_metadata[['cell_ID']]

    subset_object = subsetGiotto(gobject = gobject,
                                 spat_unit = spat_unit,
                                 feat_type = feat_type,
                                 cell_ids = filtered_cell_IDs,
                                 poly_info = poly_info,
                                 x_max = x_max,
                                 x_min = x_min,
                                 y_max = y_max,
                                 y_min = y_min,
                                 verbose = verbose)

    return(subset_object)

  } else {
    return(comb_metadata)
  }

}




#' @title Subset by spatial locations -- multi
#' @name subsetGiottoLocsMulti
#' @description Subsets Giotto object based on spatial locations
#' @inheritParams subsetGiottoLocs
#' @return giotto object
#' @details Subsets a Giotto based on spatial locations for multiple spatial units
#' if return_gobject = FALSE, then a filtered combined metadata data.table will be returned
#' @export
subsetGiottoLocsMulti = function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 spat_loc_name = NULL,
                                 x_max = NULL,
                                 x_min = NULL,
                                 y_max = NULL,
                                 y_min = NULL,
                                 z_max = NULL,
                                 z_min = NULL,
                                 poly_info = NULL,
                                 return_gobject = TRUE,
                                 verbose = TRUE) {




  res_list = list()

  for(spat_unit_selected in spat_unit) {

    poly_info_selected = poly_info[[spat_unit_selected]]

    cat('\n \n')

    if(verbose) cat('Start subset on location for spatial unit: ', spat_unit_selected,
                    'and polygon information layers: ', poly_info_selected, '\n')


    if(return_gobject == TRUE) {
      gobject = subsetGiottoLocs(gobject = gobject,
                                 spat_unit = spat_unit_selected,
                                 feat_type = feat_type,
                                 spat_loc_name = spat_loc_name,
                                 x_max = x_max,
                                 x_min = x_min,
                                 y_max = y_max,
                                 y_min = y_min,
                                 z_max = z_max,
                                 z_min = z_min,
                                 poly_info = poly_info_selected,
                                 return_gobject = return_gobject,
                                 verbose = verbose)
    } else {

      res_list[[spat_unit_selected]] = subsetGiottoLocs(gobject = gobject,
                                                        spat_unit = spat_unit_selected,
                                                        feat_type = feat_type,
                                                        spat_loc_name = spat_loc_name,
                                                        x_max = x_max,
                                                        x_min = x_min,
                                                        y_max = y_max,
                                                        y_min = y_min,
                                                        z_max = z_max,
                                                        z_min = z_min,
                                                        poly_info = poly_info_selected,
                                                        return_gobject = return_gobject,
                                                        verbose = verbose)

    }
  }

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    return(res_list)
  }

}




#' @title filterDistributions
#' @name filterDistributions
#' @description show gene or cell distribution after filtering on expression threshold
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param detection consider features (e.g. genes) or cells
#' @param plot_type type of plot
#' @param nr_bins number of bins for histogram plot
#' @param fill_color fill color for plots
#' @param scale_axis ggplot transformation for axis (e.g. log2)
#' @param axis_offset offset to be used together with the scaling transformation
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot object
#' @export
filterDistributions <- function(gobject,
                                feat_type = NULL,
                                spat_unit = NULL,
                                expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                                expression_threshold = 1,
                                detection = c('feats', 'cells'),
                                plot_type = c('histogram', 'violin'),
                                nr_bins = 30,
                                fill_color = 'lightblue',
                                scale_axis = 'identity',
                                axis_offset = 0,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'filterDistributions') {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # expression values to be used
  values = match.arg(expression_values, unique(c('raw', 'normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      feat_type = feat_type,
                                      spat_unit = spat_unit,
                                      values = values)

  # plot distribution for genes or cells
  detection = match.arg(detection, c('feats', 'cells'))

  # plot type
  plot_type = match.arg(plot_type, c('histogram', 'violin'))

  # variables
  V1 = NULL

  # for genes
  if(detection == 'feats') {

    feat_detection_levels = data.table::as.data.table(rowSums_flex(expr_values >= expression_threshold))

    if(plot_type == 'violin') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_violin(data = feat_detection_levels, ggplot2::aes(x = 'feats', y = V1+axis_offset),
                                      fill = fill_color)
      pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(y = 'feat detected in # of cells', x = '')

    } else if(plot_type == 'histogram') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_histogram(data = feat_detection_levels, ggplot2::aes(x = V1+axis_offset),
                                         color = 'white', bins = nr_bins, fill = fill_color)
      pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(x = 'feat detected in # of cells')

    }

    # for cells
  } else if(detection == 'cells') {

    cell_detection_levels = data.table::as.data.table(colSums_flex(expr_values >= expression_threshold))

    if(plot_type == 'violin') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_violin(data = cell_detection_levels, ggplot2::aes(x = 'cells', y = V1+axis_offset),
                                      fill = fill_color)
      pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(y = 'feats detected per cell', x = '')

    } else if(plot_type == 'histogram') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_histogram(data = cell_detection_levels, ggplot2::aes(x = V1+axis_offset),
                                         color = 'white', bins = nr_bins, fill = fill_color)
      pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(x = 'feats detected per cell')

    }
  }

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }

}



#' @title filterCombinations
#' @name filterCombinations
#' @description Shows how many genes and cells are lost with combinations of thresholds.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param expression_thresholds all thresholds to consider a gene expressed
#' @param feat_det_in_min_cells minimum # of cells that need to express a feature
#' @param gene_det_in_min_cells deprecated, use feat_det_in_min_cells
#' @param min_det_feats_per_cell minimum # of features that need to be detected in a cell
#' @param min_det_genes_per_cell deprecated, use min_det_feats_per_cell
#' @param scale_x_axis ggplot transformation for x-axis (e.g. log2)
#' @param x_axis_offset x-axis offset to be used together with the scaling transformation
#' @param scale_y_axis ggplot transformation for y-axis (e.g. log2)
#' @param y_axis_offset y-axis offset to be used together with the scaling transformation
#' @param show_plot show plot
#' @param return_plot return only ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return list of data.table and ggplot object
#' @details Creates a scatterplot that visualizes the number of genes and cells that are
#' lost with a specific combination of a gene and cell threshold given an arbitrary cutoff
#' to call a gene expressed. This function can be used to make an informed decision at the
#' filtering step with filterGiotto.
#' @export
filterCombinations <- function(gobject,
                               feat_type = NULL,
                               spat_unit = NULL,
                               expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                               expression_thresholds = c(1, 2),
                               feat_det_in_min_cells = c(5, 50),
                               gene_det_in_min_cells = NULL,
                               min_det_feats_per_cell = c(200, 400),
                               min_det_genes_per_cell = NULL,
                               scale_x_axis = 'identity',
                               x_axis_offset = 0,
                               scale_y_axis = 'identity',
                               y_axis_offset = 0,
                               show_plot = TRUE,
                               return_plot = FALSE,
                               save_plot = NA,
                               save_param =  list(),
                               default_save_name = 'filterCombinations') {


  ## deprecated arguments
  if(!is.null(gene_det_in_min_cells)) {
    feat_det_in_min_cells = gene_det_in_min_cells
    warning('gene_det_in_min_cells is deprecated, use feat_det_in_min_cells in the future \n')
  }
  if(!is.null(min_det_genes_per_cell)) {
    min_det_feats_per_cell = min_det_genes_per_cell
    warning('min_det_genes_per_cell is deprecated, use min_det_feats_per_cell in the future \n')
  }

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  # expression values to be used
  values = match.arg(expression_values, unique(c('raw', 'normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      feat_type = feat_type,
                                      spat_unit = spat_unit,
                                      values = values)

  # feat and cell minimums need to have the same length
  if(length(feat_det_in_min_cells) != length(min_det_feats_per_cell)) {
    stop('\n feat_det_in_min_cells and min_det_feats_per_cell need to be the same size \n')
  }

  # compute the number of removed feats and cells
  result_list = list()
  for(thresh_i in 1:length(expression_thresholds)) {

    threshold = expression_thresholds[thresh_i]

    det_feats_res = list()
    det_cells_res = list()
    for(combn_i in 1:length(feat_det_in_min_cells)) {

      min_cells_for_feat = feat_det_in_min_cells[combn_i]
      min_feats_per_cell = min_det_feats_per_cell[combn_i]


      # first remove feats
      filter_index_feats = rowSums_flex(expr_values >= threshold) >= min_cells_for_feat
      removed_feats = length(filter_index_feats[filter_index_feats == FALSE])
      det_cells_res[[combn_i]] = removed_feats

      # then remove cells
      filter_index_cells = colSums_flex(expr_values[filter_index_feats, ] >= threshold) >= min_feats_per_cell
      removed_cells = length(filter_index_cells[filter_index_cells == FALSE])
      det_feats_res[[combn_i]] = removed_cells
    }

    temp_dt = data.table::data.table('threshold' = threshold,
                                     removed_feats = unlist(det_cells_res),
                                     removed_cells = unlist(det_feats_res))

    result_list[[thresh_i]] = temp_dt

  }

  result_DT = do.call('rbind', result_list)

  # data.table variables
  # feat_detected_in_min_cells = min_detected_feats_per_cell = combination = NULL

  # data.table variables
  feat_detected_in_min_cells = min_detected_feats_per_cell = combination = NULL

  result_DT[['feat_detected_in_min_cells']] = feat_det_in_min_cells
  result_DT[['min_detected_feats_per_cell']] = min_det_feats_per_cell
  result_DT[['combination']] = paste0(result_DT$feat_detected_in_min_cells,'-',result_DT$min_detected_feats_per_cell)

  result_DT = result_DT[,.(threshold,
                           feat_detected_in_min_cells,
                           min_detected_feats_per_cell,
                           combination,
                           removed_feats,
                           removed_cells)]

  maximum_x_value = max(result_DT[['removed_cells']], na.rm = T)
  maximum_y_value = max(result_DT[['removed_feats']], na.rm = T)

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_line(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                      y = removed_feats+y_axis_offset,
                                                      group = as.factor(threshold)), linetype = 2)
  pl <- pl + ggplot2::geom_point(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                       y = removed_feats+y_axis_offset,
                                                       color = as.factor(threshold)))
  pl <- pl + scale_color_discrete(guide = guide_legend(title = 'threshold(s)'))
  pl <- pl + ggrepel::geom_text_repel(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                            y = removed_feats+y_axis_offset,
                                                            label = combination))
  pl <- pl + ggplot2::scale_x_continuous(trans = scale_x_axis, limits = c(0, maximum_x_value))
  pl <- pl + ggplot2::scale_y_continuous(trans = scale_y_axis, limits = c(0, maximum_y_value))
  pl <- pl + ggplot2::labs(x = 'number of removed cells', y = 'number of removed feats')


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  } else {
    return(list(results = result_DT, ggplot = pl))
  }

}


#' @title filterGiotto
#' @name filterGiotto
#' @description filter Giotto object based on expression threshold
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param feat_det_in_min_cells minimum # of cells that need to express a feature
#' @param gene_det_in_min_cells deprecated, use feat_det_in_min_cells
#' @param min_det_feats_per_cell minimum # of features that need to be detected in a cell
#' @param min_det_genes_per_cell deprecated, use min_det_feats_per_cell
#' @param poly_info polygon information to use
#' @param verbose verbose
#' @return giotto object
#' @details The function \code{\link{filterCombinations}} can be used to explore the effect of different parameter values.
#' @export
filterGiotto <- function(gobject,
                         spat_unit = NULL,
                         feat_type = NULL,
                         expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                         expression_threshold = 1,
                         feat_det_in_min_cells = 100,
                         gene_det_in_min_cells = NULL,
                         min_det_feats_per_cell = 100,
                         min_det_genes_per_cell = NULL,
                         poly_info = 'cell',
                         verbose = TRUE) {

  ## deprecated arguments
  if(!is.null(gene_det_in_min_cells)) {
    feat_det_in_min_cells = gene_det_in_min_cells
    warning('gene_det_in_min_cells is deprecated, use feat_det_in_min_cells in the future \n')
  }
  if(!is.null(min_det_genes_per_cell)) {
    min_det_feats_per_cell = min_det_genes_per_cell
    warning('min_det_genes_per_cell is deprecated, use min_det_feats_per_cell in the future \n')
  }


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  # expression values to be used
  values = match.arg(expression_values, unique(c('raw', 'normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      feat_type = feat_type,
                                      spat_unit = spat_unit,
                                      values = values)


  # approach:
  # 1. first remove genes that are not frequently detected
  # 2. then remove cells that do not have sufficient detected genes

  ## filter features
  filter_index_feats = rowSums_flex(expr_values >= expression_threshold) >= feat_det_in_min_cells
  selected_feat_ids = gobject@feat_ID[[feat_type]][filter_index_feats]


  ## filter cells
  filter_index_cells = colSums_flex(expr_values[filter_index_feats, ] >= expression_threshold) >= min_det_feats_per_cell
  selected_cell_ids = gobject@cell_ID[[spat_unit]][filter_index_cells]


  print(selected_cell_ids[1:4])
  print(selected_feat_ids[1:4])

  newGiottoObject = subsetGiotto(gobject = gobject,
                                 feat_type = feat_type,
                                 spat_unit = spat_unit,
                                 cell_ids = selected_cell_ids,
                                 feat_ids = selected_feat_ids,
                                 poly_info = poly_info,
                                 verbose = verbose)

  ## print output ##
  removed_feats = length(filter_index_feats[filter_index_feats == FALSE])
  total_feats   = length(filter_index_feats)

  removed_cells = length(filter_index_cells[filter_index_cells == FALSE])
  total_cells   = length(filter_index_cells)

  if(verbose == TRUE) {
    cat('Feature type: ', feat_type, '\n')
    cat('Number of cells removed: ', removed_cells, ' out of ', total_cells, '\n')
    cat('Number of feats removed: ', removed_feats, ' out of ', total_feats, '\n')
  }


  ## update parameters used ##
  
  # Do not update downstream of processGiotto
  # Parameters will be updated within processGiotto
  try({
    upstream_func = sys.call(-2)
    fname = as.character(upstream_func[[1]])
    if (fname == 'processGiotto') return(newGiottoObject)
  },
  silent = TRUE)
  
  
  # If this function call is not downstream of processGiotto, update normally
  newGiottoObject = update_giotto_params(newGiottoObject, description = '_filter')
  
  return(newGiottoObject)


}




#' @title RNA standard normalization
#' @name rna_standard_normalization
#' @description standard function for RNA normalization
#' @keywords internal
rna_standard_normalization = function(gobject,
                                      raw_expr,
                                      feat_type,
                                      spat_unit,
                                      library_size_norm = TRUE,
                                      scalefactor = 6e3,
                                      log_norm = TRUE,
                                      log_offset = 1,
                                      logbase = 2,
                                      scale_feats = TRUE,
                                      scale_cells = TRUE,
                                      scale_order = c('first_feats', 'first_cells'),
                                      verbose = TRUE) {

  # check feature type compatibility
  if(!feat_type %in% c('rna', 'RNA')) {
    warning('Caution: Standard normalization was developed for RNA data \n')
  }

  feat_names = rownames(raw_expr)
  col_names = colnames(raw_expr)

  ## 1. library size normalize
  if(library_size_norm == TRUE) {
    norm_expr = libNorm_giotto(mymatrix = raw_expr, scalefactor = scalefactor)
  } else {
    norm_expr = raw_expr
  }

  ## 2. lognormalize
  if(log_norm == TRUE) {
    norm_expr = logNorm_giotto(mymatrix = norm_expr,  base = logbase, offset = log_offset)
  } else {
    norm_expr = norm_expr
  }

  ## 3. scale
  if(scale_feats == TRUE & scale_cells == TRUE) {

    scale_order = match.arg(arg = scale_order, choices = c('first_feats', 'first_cells'))

    if(scale_order == 'first_feats') {
      if(verbose == TRUE) cat('\n first scale feats and then cells \n')

      norm_scaled_expr = t_flex(standardise_flex(x = t_flex(norm_expr), center = TRUE, scale = TRUE))
      norm_scaled_expr = standardise_flex(x = norm_scaled_expr, center = TRUE, scale = TRUE)

      #if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
      #norm_scaled_expr = t(Rfast::standardise(x = t(norm_expr), center = TRUE, scale = TRUE))
      #norm_scaled_expr = Rfast::standardise(x = norm_scaled_expr, center = TRUE, scale = TRUE)

    } else if(scale_order == 'first_cells') {
      if(verbose == TRUE) cat('\n first scale cells and then feats \n')

      norm_scaled_expr = standardise_flex(x = norm_expr, center = TRUE, scale = TRUE)
      norm_scaled_expr = t_flex(standardise_flex(x = t_flex(norm_scaled_expr), center = TRUE, scale = TRUE))

      #if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
      #norm_scaled_expr = Rfast::standardise(x = norm_expr, center = TRUE, scale = TRUE)
      #norm_scaled_expr = t(Rfast::standardise(x = t(norm_scaled_expr), center = TRUE, scale = TRUE))

    } else {
      stop('\n scale order must be given \n')
    }

  } else if(scale_feats == TRUE) {

    norm_scaled_expr = t_flex(standardise_flex(x = t_flex(norm_expr), center = TRUE, scale = TRUE))

    #if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
    #norm_scaled_expr = t(Rfast::standardise(x = t(norm_expr), center = TRUE, scale = TRUE))

  } else if(scale_cells == TRUE) {

    norm_scaled_expr = standardise_flex(x = norm_expr, center = TRUE, scale = TRUE)

    #if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
    #norm_scaled_expr = Rfast::standardise(x = norm_expr, center = TRUE, scale = TRUE)

  } else {
    norm_scaled_expr = NULL
  }


  ## 4. add cell and gene names back
  if(!is.null(norm_expr)) {
    rownames(norm_expr) = feat_names
    colnames(norm_expr) = col_names
  }
  if(!is.null(norm_scaled_expr)) {
    rownames(norm_scaled_expr) = feat_names
    colnames(norm_scaled_expr) = col_names
  }

  # return Giotto object

  gobject = set_expression_values(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  name = 'normalized',
                                  values = norm_expr)

  gobject = set_expression_values(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  name = 'scaled',
                                  values = norm_scaled_expr)

  #gobject@expression[[spat_unit]][[feat_type]][['normalized']] = norm_expr
  #gobject@expression[[spat_unit]][[feat_type]][['scaled']] = norm_scaled_expr

  return(gobject)
}



#' @title RNA osmfish normalization
#' @name rna_osmfish_normalization
#' @description function for RNA normalization according to osmFISH paper
#' @keywords internal
rna_osmfish_normalization = function(gobject,
                                     raw_expr,
                                     feat_type,
                                     spat_unit,
                                     name = 'custom',
                                     verbose = TRUE) {

  # check feature type compatibility
  if(!feat_type %in% c('rna', 'RNA')) {
    warning('Caution: osmFISH normalization was developed for RNA in situ data \n')
  }

  # 1. normalize per gene with scale-factor equal to number of genes
  norm_feats = (raw_expr/rowSums_flex(raw_expr)) * nrow(raw_expr)
  # 2. normalize per cells with scale-factor equal to number of cells
  norm_feats_cells = t_flex((t_flex(norm_feats)/colSums_flex(norm_feats)) * ncol(raw_expr))

  # return results to Giotto object
  if(verbose == TRUE) message('\n osmFISH-like normalized data will be returned to the', name, 'Giotto slot \n')

  gobject = set_expression_values(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  name = name,
                                  values = norm_feats_cells)

  #gobject@expression[[spat_unit]][[feat_type]][[name]] = norm_feats_cells

  return(gobject)
}


#' @title RNA pearson residuals normalization
#' @name rna_pears_resid_normalization
#' @description function for RNA normalization according to Lause/Kobak et al paper
#' Adapted from https://gist.github.com/hypercompetent/51a3c428745e1c06d826d76c3671797c#file-pearson_residuals-r
#' @keywords internal
rna_pears_resid_normalization = function(gobject,
                                         raw_expr,
                                         feat_type,
                                         spat_unit,
                                         theta = 100,
                                         name = 'scaled',
                                         verbose = TRUE) {


  # print message with information #
  if(verbose) message("using 'Lause/Kobak' method to normalize count matrix If used in published research, please cite:
  Jan Lause, Philipp Berens, Dmitry Kobak (2020).
                      'Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data' ")


  # check feature type compatibility
  if(!feat_type %in% c('rna', 'RNA')) {
    warning('Caution: pearson residual normalization was developed for RNA count normalization \n')
  }

  if(methods::is(raw_expr, 'HDF5Matrix')) {

    counts_sum0 = methods::as(matrix(MatrixGenerics::colSums2(raw_expr),nrow=1),"HDF5Matrix")
    counts_sum1 = methods::as(matrix(MatrixGenerics::rowSums2(raw_expr),ncol=1),"HDF5Matrix")
    counts_sum  = sum(raw_expr)

    #get residuals
    mu = (counts_sum1 %*% counts_sum0) / counts_sum
    z  = (raw_expr - mu) / sqrt(mu + mu^2/theta)

    #clip to sqrt(n)
    n = ncol(raw_expr)
    z[z > sqrt(n)]  = sqrt(n)
    z[z < -sqrt(n)] = -sqrt(n)

  } else {


    counts_sum0 = methods::as(matrix(Matrix::colSums(raw_expr),nrow=1),"dgCMatrix")
    counts_sum1 = methods::as(matrix(Matrix::rowSums(raw_expr),ncol=1),"dgCMatrix")
    counts_sum  = sum(raw_expr)

    #get residuals
    mu = (counts_sum1 %*% counts_sum0) / counts_sum
    z  = (raw_expr - mu) / sqrt(mu + mu^2/theta)

    #clip to sqrt(n)
    n = ncol(raw_expr)
    z[z > sqrt(n)]  = sqrt(n)
    z[z < -sqrt(n)] = -sqrt(n)

  }

  # return results to Giotto object
  if(verbose == TRUE) message('\n Pearson residual normalized data will be returned to the ', name, ' Giotto slot \n')

  gobject = set_expression_values(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  name = name,
                                  values = z)

  #gobject@expression[[spat_unit]][[feat_type]][[name]] = z

  return(gobject)

}




#' @title normalizeGiotto
#' @name normalizeGiotto
#' @description fast normalize and/or scale expresion values of Giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param norm_methods normalization method to use
#' @param library_size_norm normalize cells by library size
#' @param scalefactor scale factor to use after library size normalization
#' @param log_norm transform values to log-scale
#' @param log_offset offset value to add to expression matrix, default = 1
#' @param logbase log base to use to log normalize expression values
#' @param scale_feats z-score genes over all cells
#' @param scale_genes deprecated, use scale_feats
#' @param scale_cells z-score cells over all genes
#' @param scale_order order to scale feats and cells
#' @param theta theta parameter for the pearson residual normalization step
#' @param update_slot slot or name to use for the results from osmFISH and pearson residual normalization
#' @param verbose be verbose
#' @return giotto object
#' @details Currently there are two 'methods' to normalize your raw counts data.
#'
#' A. The standard method follows the standard protocol which can be adjusted using
#' the provided parameters and follows the following order: \cr
#' \itemize{
#'   \item{1. Data normalization for total library size and scaling by a custom scale-factor.}
#'   \item{2. Log transformation of data.}
#'   \item{3. Z-scoring of data by genes and/or cells.}
#' }
#' B. The normalization method as provided by the osmFISH paper is also implemented: \cr
#' \itemize{
#'   \item{1. First normalize genes, for each gene divide the counts by the total gene count and
#' multiply by the total number of genes.}
#'   \item{2. Next normalize cells, for each cell divide the normalized gene counts by the total
#' counts per cell and multiply by the total number of cells.}
#' }
#' C. The normalization method as provided by Lause/Kobak et al is also implemented: \cr
#' \itemize{
#'   \item{1. First calculate expected values based on Pearson correlations.}
#'   \item{2. Next calculate z-scores based on observed and expected values.}
#' }
#' By default the latter two results will be saved in the Giotto slot for scaled expression,
#'  this can be changed by changing the update_slot parameters
#' @export
normalizeGiotto <- function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            expression_values = 'raw',
                            norm_methods = c('standard', 'pearson_resid', 'osmFISH'),
                            library_size_norm = TRUE,
                            scalefactor = 6e3,
                            log_norm = TRUE,
                            log_offset = 1,
                            logbase = 2,
                            scale_feats = TRUE,
                            scale_genes = NULL,
                            scale_cells = TRUE,
                            scale_order = c('first_feats', 'first_cells'),
                            theta = 100,
                            update_slot = 'scaled',
                            verbose = TRUE) {



  ## deprecated arguments
  if(!is.null(scale_genes)) {
    scale_feats = scale_genes
    warning('scale_genes is deprecated, use scale_feats in the future \n')
  }

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## default is to start from raw data
  values = match.arg(expression_values, unique(c('raw', expression_values)))
  raw_expr = get_expression_values(gobject = gobject,
                                   feat_type = feat_type,
                                   spat_unit = spat_unit,
                                   values = values)

  norm_methods = match.arg(arg = norm_methods, choices = c('standard', 'pearson_resid', 'osmFISH'))

  # normalization according to standard methods
  if(norm_methods == 'standard') {

    gobject = rna_standard_normalization(gobject = gobject,
                                         raw_expr = raw_expr,
                                         feat_type = feat_type,
                                         spat_unit = spat_unit,
                                         library_size_norm = library_size_norm,
                                         scalefactor = scalefactor,
                                         log_norm = log_norm,
                                         log_offset = log_offset,
                                         logbase = logbase,
                                         scale_feats = scale_feats,
                                         scale_cells = scale_cells,
                                         scale_order = scale_order,
                                         verbose = verbose)


  }

  else if(norm_methods == 'osmFISH') {

    gobject = rna_osmfish_normalization(gobject = gobject,
                                        raw_expr = raw_expr,
                                        feat_type = feat_type,
                                        spat_unit = spat_unit,
                                        name = update_slot,
                                        verbose = verbose)

  }

  else if(norm_methods == 'pearson_resid') {

    gobject = rna_pears_resid_normalization(gobject = gobject,
                                            raw_expr = raw_expr,
                                            feat_type = feat_type,
                                            spat_unit = spat_unit,
                                            theta = theta,
                                            name = update_slot,
                                            verbose = verbose)

  }

  ## update parameters used ##
  
  # Do not update downstream of processGiotto
  # Parameters will be updated within processGiotto
  try({
    upstream_func = sys.call(-2)
    fname = as.character(upstream_func[[1]])
    if (fname == 'processGiotto') return(gobject)
  },
  silent = TRUE)
  
  
  # If this function call is not downstream of processGiotto, update normally
  gobject = update_giotto_params(gobject, description = '_normalize')

  return(gobject)

}



#' @title Adjust expression values
#' @name adjustGiottoMatrix
#' @description Adjust expression values to account for known batch effects or technological covariates.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param batch_columns metadata columns that represent different batch (max = 2)
#' @param covariate_columns metadata columns that represent covariates to regress out
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param update_slot expression slot that will be updated (default = custom)
#' @return giotto object
#' @details This function implements the \code{\link[limma]{removeBatchEffect}} function to
#' remove known batch effects and to adjust expression values according to provided covariates.
#' @export
#'
adjustGiottoMatrix <- function(gobject,
                               spat_unit = NULL,
                               feat_type = NULL,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               batch_columns = NULL,
                               covariate_columns = NULL,
                               return_gobject = TRUE,
                               update_slot = c('custom')) {
  
  # Catch for both batch and covariate being null
  if (is.null(batch_columns) & is.null(covariate_columns)){
    stop('\nMetadata for either different batches or covariates must be provided.')
  }
  
  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # metadata
  cell_metadata = pDataDT(gobject,
                          feat_type = feat_type,
                          spat_unit = spat_unit)

  if(!is.null(batch_columns)) {
    if(!all(batch_columns %in% colnames(cell_metadata))) {
      stop('\n batch column name(s) were not found in the cell metadata \n')
    }
  }

  if(!is.null(covariate_columns)) {
    if(!all(covariate_columns %in% colnames(cell_metadata))) {
      stop('\n covariate column name(s) were not found in the cell metadata \n')
    }
  }

  update_slot = match.arg(update_slot, c('normalized', 'scaled', 'custom', update_slot))

  # expression values to be used
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_data = get_expression_values(gobject = gobject,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit,
                                    values = values)


  # batch columns
  if(!is.null(batch_columns)) {
    batch_column_1 = cell_metadata[[ batch_columns[1] ]]
    if(length(batch_columns) > 1) {
      batch_column_2 = cell_metadata[[ batch_columns[2] ]]
    } else {
      batch_column_2 = NULL
    }
  } else {
    batch_column_1 = NULL
    batch_column_2 = NULL
  }

  # covariate columns
  if(!is.null(covariate_columns)) {
    covariates = as.matrix(cell_metadata[, covariate_columns, with = F])
  } else {
    covariates = NULL
  }



  # TODO: implement ResidualMatrix to work with a delayed matrix
  adjusted_matrix = limma::removeBatchEffect(x = expr_data,
                                             batch = batch_column_1,
                                             batch2 =  batch_column_2,
                                             covariates = covariates)

  if(return_gobject == TRUE) {

    gobject = set_expression_values(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type,
                                    name = update_slot,
                                    values = adjusted_matrix)
    
    ## update parameters used ##
    
    # Do not update downstream of processGiotto
    # Parameters will be updated within processGiotto
    try({
      test = sys.call(-2)
      fname = as.character(test[[1]])
      if (fname == 'processGiotto') return(gobject)
    },
    silent = TRUE)
    
    
    # If this function call is not downstream of processGiotto, update normally
    gobject = update_giotto_params(gobject, description = '_adj_matrix')
    
    return(gobject)

  } else {

    return(adjusted_matrix)

  }

}


#' @title processGiotto
#' @name processGiotto
#' @description Wrapper for the different Giotto object processing functions
#' @param gobject giotto object
#' @param filter_params additional parameters to filterGiotto
#' @param norm_params additional parameters to normalizeGiotto
#' @param stat_params additional parameters to addStatistics
#' @param adjust_params additional parameters to adjustGiottoMatrix; set to NULL if not required
#' @param verbose be verbose (default is TRUE)
#' @return giotto object
#' @details See \code{\link{filterGiotto}}, \code{\link{normalizeGiotto}},
#' \code{\link{addStatistics}}, and \code{\link{adjustGiottoMatrix}}. For more
#' information about the different parameters in each step. If you do not provide
#' them it will use the default values. If no adjustment is required, adjust_params must be set to NULL
#' @export
#'
processGiotto = function(gobject,
                         filter_params = list(),
                         norm_params = list(),
                         stat_params = list(),
                         adjust_params = list(),
                         verbose = TRUE) {

  # filter Giotto
  if(verbose == TRUE) cat('1. start filter step \n')
  if(!inherits(filter_params, 'list')) stop('filter_params need to be a list of parameters for filterGiotto \n')
  gobject = do.call('filterGiotto', c(gobject = gobject, filter_params))

  # normalize Giotto
  if(verbose == TRUE) cat('2. start normalization step \n')
  if(!inherits(norm_params, 'list')) stop('norm_params need to be a list of parameters for normalizeGiotto \n')
  gobject = do.call('normalizeGiotto', c(gobject = gobject, norm_params))

  # add Statistics
  if(verbose == TRUE) cat('3. start cell and gene statistics step \n')
  if(!inherits(stat_params, 'list')) stop('stat_params need to be a list of parameters for addStatistics \n')
  stat_params[['return_gobject']] = TRUE # force this to be true
  gobject = do.call('addStatistics', c(gobject = gobject, stat_params))

  # adjust Giotto, if applicable
  if(!is.null(adjust_params)){
    if(verbose == TRUE) cat('4. start adjusted matrix step \n')
    if(!inherits(adjust_params, 'list')) stop('adjust_params need to be a list of parameters for adjustGiottoMatrix \n')
    adjust_params[['return_gobject']] = TRUE # force this to be true
    gobject = do.call('adjustGiottoMatrix', c(gobject = gobject, adjust_params))
  }

  gobject = update_giotto_params(gobject, description = '_process')
  
  return(gobject)

}







## * ####
## Gene & Cell metadata functions ####


#' @title Annotate giotto clustering
#' @name annotateGiotto
#' @description Converts cluster results into a user provided annotation.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param annotation_vector named annotation vector (names = cluster ids)
#' @param cluster_column cluster column to convert to annotation names
#' @param name new name for annotation column
#' @return giotto object
#' @details You need to specifify which (cluster) column you want to annotate
#' and you need to provide an annotation vector like this:
#' \itemize{
#'   \item{1. identify the cell type of each cluster}
#'   \item{2. create a vector of these cell types, e.g. cell_types =  c('T-cell', 'B-cell', 'Stromal')}
#'   \item{3. provide original cluster names to previous vector, e.g. names(cell_types) = c(2, 1, 3)}
#' }
#' @export
annotateGiotto <- function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL,
                           annotation_vector = NULL,
                           cluster_column = NULL,
                           name = 'cell_types') {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # data.table: set global variable
  temp_cluster_name = NULL

  if(is.null(annotation_vector) | is.null(cluster_column)) {
    stop('\n You need to provide both a named annotation vector and the corresponding cluster column  \n')
  }

  cell_metadata = pDataDT(gobject,
                          feat_type = feat_type,
                          spat_unit = spat_unit)

  # 1. verify if cluster column exist
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n Cluster column is not found in cell metadata \n')
  }

  # 2. verify if each cluster has an annotation
  uniq_names = names(annotation_vector)
  uniq_clusters = unique(cell_metadata[[cluster_column]])
  missing_annotations = uniq_clusters[!uniq_clusters %in% uniq_names]
  no_matching_annotations = uniq_names[!uniq_names %in% uniq_clusters]

  if(length(missing_annotations) > 0) {
    cat('Not all clusters have an accompanying annotation in the annotation_vector: \n',
        'These names are missing: ', as.character(missing_annotations), '\n',
        'These annotations have no match: ', as.character(no_matching_annotations), '\n')
    stop('Annotation interrupted \n')
  }


  # 3. remove previous annotation name if it's the same
  # but only if new name is not the same as cluster to be used
  if(name %in% colnames(cell_metadata)) {
    cat('\n annotation name ', name,' was already used \n',
        'and will be overwritten \n')

    cell_metadata[, temp_cluster_name := annotation_vector[[as.character(get(cluster_column))]], by = 1:nrow(cell_metadata)]
    cell_metadata[, (name) := NULL]

  } else {

    cell_metadata[, temp_cluster_name := annotation_vector[[as.character(get(cluster_column))]], by = 1:nrow(cell_metadata)]
  }

  data.table::setnames(cell_metadata, old = 'temp_cluster_name', new = name)
  gobject@cell_metadata[[spat_unit]][[feat_type]] = cell_metadata

  return(gobject)


}



#' @title Remove cell annotation
#' @name removeCellAnnotation
#' @description Removes cell annotation from a Giotto object for a specific feature modality (default = 'rna')
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param columns names of columns to remove
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object
#' @details if \code{return_gobject = FALSE}, it will return the cell metadata
#' @export
removeCellAnnotation <- function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 columns = NULL,
                                 return_gobject = TRUE) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  if(is.null(columns)) {
    stop('\t You need to provide a vector of metadata column names to remove \t')
  }

  gobject@cell_metadata[[spat_unit]][[feat_type]][, (columns) := NULL]

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    gobject@cell_metadata[[spat_unit]][[feat_type]]
  }

}


#' @title Remove feature annotation
#' @name removeFeatAnnotation
#' @description Removes feature annotation from a Giotto object for a specific feature modality
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param columns names of columns to remove
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object
#' @details if \code{return_gobject = FALSE}, it will return the gene metadata
#' @export
removeFeatAnnotation <- function(gobject,
                                 spat_unit = NULL,
                                 feat_type = NULL,
                                 columns = NULL,
                                 return_gobject = TRUE) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  if(is.null(columns)) {
    stop('\t You need to provide a vector of metadata column names to remove \t')
  }

  gobject@feat_metadata[[spat_unit]][[feat_type]][, (columns) := NULL]

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    gobject@gene_metadata[[spat_unit]][[feat_type]]
  }

}



#' @title Add cell metadata
#' @name addCellMetadata
#' @description Adds cell metadata to the giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param new_metadata new cell metadata to use (data.table, data.frame, ...)
#' @param vector_name (optional) custom name if you provide a single vector
#' @param by_column merge metadata based on \emph{cell_ID} column in \code{\link{pDataDT}} (default = FALSE)
#' @param column_cell_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details You can add additional cell metadata in two manners:
#' \itemize{
#'   \item{1. Provide a data.table or data.frame with cell annotations in the same order as the \emph{cell_ID} column in pDataDT(gobject) }
#'   \item{2. Provide a data.table or data.frame with cell annotations and specify which column contains the cell IDs, these cell IDs need to match with the \emph{cell_ID} column in pDataDT(gobject)}
#' }
#' @export
addCellMetadata <- function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            new_metadata,
                            vector_name = NULL,
                            by_column = FALSE,
                            column_cell_ID = NULL) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  cell_metadata = data.table::copy(gobject@cell_metadata[[spat_unit]][[feat_type]])
  ordered_cell_IDs = gobject@cell_ID[[spat_unit]]

  if(is.vector(new_metadata) | is.factor(new_metadata)) {
    original_name = deparse(substitute(new_metadata))
    new_metadata = data.table::as.data.table(new_metadata)

    if(!is.null(vector_name) & is.character(vector_name)) {
      colnames(new_metadata) = vector_name
    } else {
      colnames(new_metadata) = original_name
    }

  } else {
    new_metadata = data.table::as.data.table(new_metadata)
  }

  if(is.null(column_cell_ID)) {
    column_cell_ID = 'cell_ID'
  }

  # overwrite columns with same name
  new_col_names = colnames(new_metadata)
  new_col_names = new_col_names[new_col_names != column_cell_ID]
  old_col_names = colnames(cell_metadata)
  old_col_names = old_col_names[old_col_names != 'cell_ID']
  same_col_names = new_col_names[new_col_names %in% old_col_names]


  if(length(same_col_names) >= 1) {
    cat('\n these column names were already used: ', same_col_names, '\n',
        'and will be overwritten \n')
    cell_metadata[, (same_col_names) := NULL]
  }



  if(by_column == FALSE) {
    cell_metadata = cbind(cell_metadata, new_metadata)
  } else {
    if(is.null(column_cell_ID)) stop('You need to provide cell_ID column')
    cell_metadata = data.table::merge.data.table(cell_metadata,
                                                 by.x = 'cell_ID',
                                                 new_metadata,
                                                 by.y = column_cell_ID,
                                                 all.x = TRUE)
  }
  
  # data.table variables
  cell_ID = NULL

  # reorder
  cell_metadata = cell_metadata[match(ordered_cell_IDs, cell_ID)]

  gobject@cell_metadata[[spat_unit]][[feat_type]] = cell_metadata
  return(gobject)
}


#' @title Add feature metadata
#' @name addFeatMetadata
#' @description Adds feature metadata to the giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param new_metadata new metadata to use
#' @param by_column merge metadata based on \emph{feat_ID} column in \code{\link{fDataDT}}
#' @param column_feat_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details You can add additional feature metadata in two manners: \cr
#' 1. Provide a data.table or data.frame with feature annotations in the same order as the \emph{feat_ID} column in fDataDT(gobject) \cr
#' 2. Provide a data.table or data.frame with feature annotations and specify which column contains the feature IDs,
#' these feature IDs need to match with the \emph{feat_ID} column in fDataDT(gobject)
#' @export
addFeatMetadata <- function(gobject,
                            feat_type = NULL,
                            spat_unit = NULL,
                            new_metadata,
                            by_column = F,
                            column_feat_ID = NULL) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # data.table variables
  feat_ID = NULL

  feat_metadata = gobject@feat_metadata[[spat_unit]][[feat_type]]
  ordered_feat_IDs = gobject@feat_ID[[feat_type]]

  if(by_column == FALSE) {
    feat_metadata = cbind(feat_metadata, new_metadata)
  } else {
    if(is.null(column_feat_ID)) stop('You need to provide feat ID column')
    feat_metadata <- data.table::merge.data.table(feat_metadata,
                                                  by.x = 'feat_ID',
                                                  new_metadata,
                                                  by.y = column_feat_ID,
                                                  all.x = T)
  }

  # reorder
  feat_metadata = feat_metadata[match(ordered_feat_IDs, feat_ID)]

  gobject@feat_metadata[[spat_unit]][[feat_type]] = feat_metadata
  return(gobject)
}




#' @title Add gene metadata
#' @name addGeneMetadata
#' @description adds gene metadata to the giotto object
#' @param gobject giotto object
#' @param new_metadata new metadata to use
#' @param by_column merge metadata based on gene_ID column in \code{\link{fDataDT}}
#' @param column_gene_ID column name of new metadata to use if \code{by_column = TRUE}
#' @return giotto object
#' @details You can add additional gene metadata in two manners:
#' 1. Provide a data.table or data.frame with gene annotations in the same order as the gene_ID column in fDataDT(gobject)
#' 2. Provide a data.table or data.frame with gene annotations and specify which column contains the gene IDs,
#' these gene IDs need to match with the gene_ID column in fDataDT(gobject)
#' @export
addGeneMetadata <- function(gobject,
                            new_metadata,
                            by_column = F,
                            column_gene_ID = NULL) {

  warning("Deprecated and replaced by addFeatMetadata")

  result = addFeatMetadata(gobject = gobject,
                           feat_type = NULL,
                           new_metadata = new_metadata,
                           by_column = by_column,
                           column_feat_ID = column_gene_ID)

  return(result)
}



#' @title Add feature statistics
#' @name addFeatStatistics
#' @description Adds feature statistics to the giotto object
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to feature metadata:
#' \itemize{
#'   \item{nr_cells: }{Denotes in how many cells the feature is detected}
#'   \item{per_cells: }{Denotes in what percentage of cells the feature is detected}
#'   \item{total_expr: }{Shows the total sum of feature expression in all cells}
#'   \item{mean_expr: }{Average feature expression in all cells}
#'   \item{mean_expr_det: }{Average feature expression in cells with detectable levels of the gene}
#' }
#' @export
addFeatStatistics <- function(gobject,
                              feat_type = NULL,
                              spat_unit = NULL,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # expression values to be used
  expression_values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_data = get_expression_values(gobject = gobject,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit,
                                    values = expression_values)

  # calculate stats
  feat_stats = data.table::data.table(feats = rownames(expr_data),
                                      nr_cells = rowSums_flex(expr_data > detection_threshold),
                                      perc_cells = (rowSums_flex(expr_data > detection_threshold)/ncol(expr_data))*100,
                                      total_expr = rowSums_flex(expr_data),
                                      mean_expr = rowMeans_flex(expr_data))

  # data.table variables
  mean_expr_det = NULL

  mean_expr_detected = mean_expr_det_test(expr_data, detection_threshold = detection_threshold)
  feat_stats[, mean_expr_det := mean_expr_detected]


  if(return_gobject == TRUE) {

    # remove previous statistics
    feat_metadata = fDataDT(gobject,
                            feat_type = feat_type,
                            spat_unit = spat_unit)
    metadata_names = colnames(feat_metadata)

    if('nr_cells' %in% metadata_names) {
      cat('\n feat statistics has already been applied once, will be overwritten \n')
      feat_metadata[, c('nr_cells', 'perc_cells', 'total_expr', 'mean_expr', 'mean_expr_det') := NULL]
      gobject@feat_metadata[[feat_type]][[spat_unit]] = feat_metadata
    }

    gobject = addFeatMetadata(gobject = gobject,
                              feat_type = feat_type,
                              spat_unit = spat_unit,
                              new_metadata = feat_stats,
                              by_column = T,
                              column_feat_ID = 'feats')

    ## update parameters used ##

    # parent function name
    cl = sys.call(-1)

    # Do not update downstream of processGiotto
    # Parameters will be updated within processGiotto
    try({
      upstream_func = sys.call(-3)
      fname = as.character(upstream_func[[1]])
      if (fname == 'processGiotto') return(gobject)
    },
    silent = TRUE)
    
    
    # If this function call is not downstream of processGiotto, update normally    
    if(is.null(cl)) {
      gobject = update_giotto_params(gobject, description = '_feat_stats')
    } else {
      fname = as.character(cl[[1]])
      if(fname == 'addStatistics') {
        gobject = update_giotto_params(gobject, description = '_feat_stats', toplevel = 3)
      } else {
        gobject = update_giotto_params(gobject, description = '_feat_stats')
      }
    }


    return(gobject)


  } else {
    return(feat_stats)
  }

}



#' @title Add gene statistics
#' @name addGeneStatistics
#' @description adds gene statistics to the giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to gene metadata:
#' \itemize{
#'   \item{nr_cells: }{Denotes in how many cells the gene is detected}
#'   \item{per_cells: }{Denotes in what percentage of cells the gene is detected}
#'   \item{total_expr: }{Shows the total sum of gene expression in all cells}
#'   \item{mean_expr: }{Average gene expression in all cells}
#'   \item{mean_expr_det: }{Average gene expression in cells with detectable levels of the gene}
#' }
#' @export
addGeneStatistics <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  warning("Deprecated and replaced by addFeatStatistics")

  result = addFeatStatistics(gobject = gobject,
                             feat_type = NULL,
                             expression_values = expression_values,
                             detection_threshold = detection_threshold,
                             return_gobject = return_gobject)

}






#' @title addCellStatistics
#' @name addCellStatistics
#' @description adds cells statistics to the giotto object
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to cell metadata:
#' \itemize{
#'   \item{nr_feats: }{Denotes in how many features are detected per cell}
#'   \item{perc_feats: }{Denotes what percentage of features is detected per cell}
#'   \item{total_expr: }{Shows the total sum of feature expression per cell}
#' }
#' @export
addCellStatistics <- function(gobject,
                              feat_type = NULL,
                              spat_unit = NULL,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # expression values to be used
  expression_values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_data = get_expression_values(gobject = gobject,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit,
                                    values = expression_values)

  # calculate stats
  #print('ok 1')
  cell_stats = data.table::data.table(cells = colnames(expr_data),
                                      nr_feats = colSums_flex(expr_data > detection_threshold),
                                      perc_feats = (colSums_flex(expr_data > detection_threshold)/nrow(expr_data))*100,
                                      total_expr = colSums_flex(expr_data))

  if(return_gobject == TRUE) {

    # remove previous statistics
    cell_metadata = pDataDT(gobject,
                            feat_type = feat_type,
                            spat_unit = spat_unit)
    metadata_names = colnames(cell_metadata)
    if('nr_feats' %in% metadata_names) {
      cat('\n cells statistics has already been applied once, will be overwritten \n')
      cell_metadata[, c('nr_feats', 'perc_feats', 'total_expr') := NULL]
      gobject@cell_metadata[[feat_type]][[spat_unit]] = cell_metadata
    }


    #print('ok 2')
    #print(cell_stats)
    gobject = addCellMetadata(gobject = gobject,
                              feat_type = feat_type,
                              spat_unit = spat_unit,
                              new_metadata = cell_stats,
                              by_column = TRUE,
                              column_cell_ID = 'cells')

    ## update parameters used ##
    
    # parent function name
    cl = sys.call(-1)
    
    # Do not update downstream of processGiotto
    # Parameters will be updated within processGiotto
    try({
      upstream_func = sys.call(-3)
      fname = as.character(upstream_func[[1]])
      if (fname == 'processGiotto') return(gobject)
    },
    silent = TRUE)
    
    # If this function call is not downstream of processGiotto, update normally
    if(is.null(cl)) {
      gobject = update_giotto_params(gobject, description = '_cell_stats')
    } else {

      fname = as.character(cl[[1]])
      if(fname == 'addStatistics') {
        gobject = update_giotto_params(gobject, description = '_cell_stats', toplevel = 3)
      } else {
        gobject = update_giotto_params(gobject, description = '_cell_stats')
      }

    }


    return(gobject)


  } else {
    return(cell_stats)
  }

}


#' @title addStatistics
#' @name addStatistics
#' @description Adds feature and cell statistics to the giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a feature detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE, else a list with results
#' @details See \code{\link{addFeatStatistics}} and \code{\link{addCellStatistics}}
#' @export
addStatistics <- function(gobject,
                          feat_type = NULL,
                          spat_unit = NULL,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          detection_threshold = 0,
                          return_gobject = TRUE) {



  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # get feats statistics
  feat_stats = addFeatStatistics(gobject = gobject,
                                 feat_type = feat_type,
                                 spat_unit = spat_unit,
                                 expression_values = expression_values,
                                 detection_threshold = detection_threshold,
                                 return_gobject = return_gobject)

  if(return_gobject == TRUE) {
    gobject = feat_stats
  }

  # get cell statistics
  cell_stats = addCellStatistics(gobject = gobject,
                                 feat_type = feat_type,
                                 spat_unit = spat_unit,
                                 expression_values = expression_values,
                                 detection_threshold = detection_threshold,
                                 return_gobject = return_gobject)

  if(return_gobject == TRUE) {
    gobject = cell_stats
    return(gobject)
  } else {
    return(feat_stats = feat_stats, cell_stats = cell_stats)
  }

}




#' @title addFeatsPerc
#' @name addFeatsPerc
#' @description Calculates the total percentage of (normalized) counts for a subset of selected genes
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param feats vector of selected features
#' @param vector_name column name as seen in \code{\link{pDataDT}}
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if \code{return_gobject = TRUE}, else a vector with % results
#' @export
addFeatsPerc = function(gobject,
                        spat_unit = NULL,
                        feat_type = NULL,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        feats = NULL,
                        vector_name = 'feat_perc',
                        return_gobject = TRUE) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # tests
  if(is.null(feats)) {
    stop('You need to provide a vector of feat names \n')
  }

  if(!methods::is(gobject, 'giotto')) {
    stop('You need to provide a giotto object \n')
  }


  # expression values to be used
  expression_values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_data = get_expression_values(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type,
                                    values = expression_values)


  totalsum = colSums_flex(expr_data)
  feat_sum = colSums_flex(expr_data[rownames(expr_data) %in% feats,])
  perc_feats = round((feat_sum/totalsum)*100, 2)

  if(return_gobject == TRUE) {
    temp_gobj = addCellMetadata(gobject = gobject,
                                spat_unit = spat_unit,
                                feat_type = feat_type,
                                new_metadata = perc_feats,
                                vector_name = vector_name,
                                by_column = F)

    ## update parameters used ##
    temp_gobj = update_giotto_params(temp_gobj, description = '_feats_perc')

    return(temp_gobj)
  } else {
    return(perc_feats)
  }

}



#' @title addGenesPerc
#' @name addGenesPerc
#' @description calculates the total percentage of (normalized) counts for a subset of selected genes
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param genes vector of selected genes
#' @param vector_name column name as seen in \code{\link{pDataDT}}
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if \code{return_gobject = TRUE}, else a vector with % results
#' @export
addGenesPerc = function(gobject,
                        spat_unit = NULL,
                        feat_type = NULL,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes = NULL,
                        vector_name = 'gene_perc',
                        return_gobject = TRUE) {

  warning("Deprecated and replaced by addFeatsPerc")


  result = addFeatsPerc(gobject = gobject,
                        feat_type = NULL,
                        expression_values = expression_values,
                        feats = genes,
                        vector_name = vector_name,
                        return_gobject = return_gobject)

  return(result)
}





## * ####
## Giotto auxiliary functions ####


#' @title showProcessingSteps
#' @name showProcessingSteps
#' @description shows the sequential processing steps that were performed
#' on a Giotto object in a summarized format
#' @param gobject giotto object
#' @return list of processing steps and names
#' @export
showProcessingSteps <- function(gobject) {

  parameters = gobject@parameters

  cat('Processing steps: \n \n')

  for(step in names(parameters)) {
    cat('\n', step, '\n')

    sub_step = parameters[[step]]

    if(any(grepl('name', names(sub_step)) == TRUE)) {

      selected_names = grep('name', names(sub_step), value = T)
      cat('\t name info: ', sub_step[selected_names], '\n')

    }

  }


}




#' @title create_cluster_matrix
#' @name create_cluster_matrix
#' @description creates aggregated matrix for a given clustering column
#' @keywords internal
create_cluster_matrix <- function(gobject,
                                  spat_unit = NULL,
                                  feat_type = NULL,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  cluster_column,
                                  feat_subset = NULL,
                                  gene_subset = NULL) {


  # data.table variables
  feats = NULL

  ## deprecated arguments
  if(!is.null(gene_subset)) {
    feat_subset = gene_subset
  }

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))

  # average expression per cluster
  aggr_sc_clusters <- create_average_DT(gobject = gobject,
                                        spat_unit = spat_unit,
                                        feat_type = feat_type,
                                        meta_data_name = cluster_column,
                                        expression_values = values)
  aggr_sc_clusters_DT <- data.table::as.data.table(aggr_sc_clusters)
  aggr_sc_clusters_DT[, feats := rownames(aggr_sc_clusters)]
  aggr_sc_clusters_DT_melt <- data.table::melt.data.table(aggr_sc_clusters_DT,
                                                          variable.name = 'cluster',
                                                          id.vars = 'feats',
                                                          value.name = 'expression')

  # create matrix
  testmat = data.table::dcast.data.table(aggr_sc_clusters_DT_melt,
                                         formula = feats~cluster,
                                         value.var = 'expression')
  testmatrix = dt_to_matrix(testmat)

  # create subset if required
  if(!is.null(feat_subset)) {
    feat_subset_detected = feat_subset[feat_subset %in% rownames(testmatrix)]
    testmatrix = testmatrix[rownames(testmatrix) %in% feat_subset_detected, ]
  }

  return(testmatrix)

}





#' @title calculateMetaTable
#' @name calculateMetaTable
#' @description calculates the average gene expression for one or more (combined) annotation columns.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param metadata_cols annotation columns found in \code{pDataDT(gobject)}
#' @param selected_feats subset of features to use
#' @param selected_genes subset of genes to use (deprecated)
#' @return data.table with average expression values for each gene per (combined) annotation
#' @export
calculateMetaTable = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              expression_values =  c("normalized", "scaled", "custom"),
                              metadata_cols = NULL,
                              selected_feats = NULL,
                              selected_genes = NULL) {

  if(is.null(metadata_cols)) stop('\n You need to select one or more valid column names from pDataDT() \n')

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## deprecated arguments
  if(!is.null(selected_genes)) {
    selected_feats = selected_genes
    warning('selected_genes is deprecated, use selected_feats in the future \n')
  }

  # data.table variables
  uniq_ID = NULL

  ## get metadata and create unique groups
  metadata = data.table::copy(pDataDT(gobject, feat_type = feat_type, spat_unit = spat_unit))
  if(length(metadata_cols) > 1) {
    metadata[, uniq_ID := paste(.SD, collapse = '-'), by = 1:nrow(metadata), .SDcols = metadata_cols]
  } else {
    metadata[, uniq_ID := get(metadata_cols)]
  }

  ## possible groups
  possible_groups = unique(metadata[,metadata_cols, with = F])
  if(length(metadata_cols) > 1) {
    possible_groups[, uniq_ID := paste(.SD, collapse = '-'), by = 1:nrow(possible_groups), .SDcols = metadata_cols]
  } else {
    possible_groups[, uniq_ID := get(metadata_cols)]
  }

  ## get expression data
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      feat_type = feat_type,
                                      spat_unit = spat_unit,
                                      values = values)
  if(!is.null(selected_feats)) {
    expr_values = expr_values[rownames(expr_values) %in% selected_feats, ]
  }

  ## summarize unique groups (average)
  result_list = list()

  for(row in 1:nrow(possible_groups)) {

    uniq_identifiier = possible_groups[row][['uniq_ID']]
    selected_cell_IDs = metadata[uniq_ID == uniq_identifiier][['cell_ID']]
    sub_expr_values = expr_values[, colnames(expr_values) %in% selected_cell_IDs]

    if(is.vector(sub_expr_values) == FALSE) {
      subvec = rowMeans_flex(sub_expr_values)
    } else {
      subvec = sub_expr_values
    }
    result_list[[row]] = subvec
  }
  finaldt = data.table::as.data.table(do.call('rbind', result_list))
  possible_groups_res = cbind(possible_groups, finaldt)
  possible_groups_res_melt = data.table::melt.data.table(possible_groups_res, id.vars = c(metadata_cols, 'uniq_ID'))

  return(possible_groups_res_melt)

}


#' @title calculateMetaTableCells
#' @name calculateMetaTableCells
#' @description calculates the average metadata values for one or more (combined) annotation columns.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param value_cols metadata or enrichment value columns to use
#' @param metadata_cols annotation columns found in \code{pDataDT(gobject)}
#' @param spat_enr_names which spatial enrichment results to include
#' @return data.table with average metadata values per (combined) annotation
#' @export
calculateMetaTableCells = function(gobject,
                                   spat_unit = NULL,
                                   feat_type = NULL,
                                   value_cols = NULL,
                                   metadata_cols = NULL,
                                   spat_enr_names = NULL) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  if(is.null(metadata_cols)) stop('\n You need to select one or more valid column names from pDataDT() \n')
  if(is.null(value_cols)) stop('\n You need to select one or more valid value column names from pDataDT() \n')

  cell_metadata = combineMetadata(gobject = gobject,
                                  feat_type = feat_type,
                                  spat_unit = spat_unit,
                                  spat_enr_names = spat_enr_names)

  ## only keep columns that exist
  cell_metadata_cols = colnames(cell_metadata)

  if(!all(value_cols %in% cell_metadata_cols)) {
    missing_value_cols = value_cols[!value_cols %in% cell_metadata_cols]
    cat('These value columns were not found: ', missing_value_cols)
  }
  value_cols = value_cols[value_cols %in% cell_metadata_cols]

  if(!all(metadata_cols %in% cell_metadata_cols)) {
    missing_metadata_cols = metadata_cols[!metadata_cols %in% cell_metadata_cols]
    cat('These metadata columns were not found: ', missing_metadata_cols)
  }
  metadata_cols = metadata_cols[metadata_cols %in% cell_metadata_cols]

  if(!length(metadata_cols) > 0 | !length(value_cols) > 0) {
    stop('\n missing sufficient metadata or value columns \n')
  }

  workdt = cell_metadata[, lapply(.SD, mean), by = metadata_cols, .SDcols = value_cols]
  workdtmelt = data.table::melt.data.table(workdt, measure.vars = value_cols)

  return(workdtmelt)

}




#' @title combineMetadata
#' @name combineMetadata
#' @description This function combines the cell metadata with spatial locations and
#' enrichment results from \code{\link{runSpatialEnrich}}
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param spat_loc_name name of spatial locations to include
#' @param spat_enr_names names of spatial enrichment results to include
#' @return Extended cell metadata in data.table format.
#' @export
combineMetadata = function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL,
                           spat_loc_name = 'raw',
                           spat_enr_names = NULL) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # cell metadata
  metadata = pDataDT(gobject,
                     spat_unit = spat_unit,
                     feat_type = feat_type)

  # spatial locations
  if(!is.null(spat_loc_name)) {
    spatial_locs = get_spatial_locations(gobject = gobject,
                                         spat_unit = spat_unit,
                                         spat_loc_name = spat_loc_name)
  } else {
    spatial_locs = NULL
  }


  # data.table variables
  cell_ID = NULL

  if(!is.null(spatial_locs)) {
    metadata = cbind(metadata, spatial_locs[, cell_ID := NULL])
  }


  # cell/spot enrichment data
  available_enr = list_spatial_enrichments_names(gobject = gobject, spat_unit = spat_unit, feat_type = feat_type)

  # output warning if not found
  not_available = spat_enr_names[!spat_enr_names %in% available_enr]
  if(length(not_available) > 0) {
    cat('These spatial enrichment results have not been found: \n',
        not_available)
  }

  spat_enr_names = spat_enr_names[spat_enr_names %in% available_enr]

  if(!is.null(spat_enr_names) & length(spat_enr_names) > 0) {

    result_list = list()
    for(spatenr in 1:length(spat_enr_names)) {

      spatenr_name = spat_enr_names[spatenr]

      temp_spat = get_spatial_enrichment(gobject = gobject,
                                         spat_unit = spat_unit,
                                         feat_type = feat_type,
                                         enrichm_name = spatenr_name)

      temp_spat[, 'cell_ID' := NULL]

      result_list[[spatenr]] = temp_spat
    }
    final_meta = do.call('cbind', c(list(metadata), result_list))

    duplicates = sum(duplicated(colnames(final_meta)))
    if(duplicates > 0) cat('Some column names are not unique.
                           If you add results from multiple enrichments,
                           consider giving the signatures unique names')

  } else {

    final_meta = metadata

  }

  return(final_meta)

}





#' @title createMetafeats
#' @name createMetafeats
#' @description This function creates an average metafeat/metagene/module for clusters.
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param feat_clusters numerical vector with features as names
#' @param name name of the metagene results
#' @param return_gobject return giotto object
#' @return giotto object
#' @details An example for the 'gene_clusters' could be like this:
#' cluster_vector = c(1, 1, 2, 2); names(cluster_vector) = c('geneA', 'geneB', 'geneC', 'geneD')
#' @export
createMetafeats = function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           feat_clusters,
                           name = 'metafeat',
                           return_gobject = TRUE) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # expression values to be used
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values)


  ## calculate metafeat ##
  res_list = list()

  for(id in sort(unique(feat_clusters))) {

    clus_id = id

    selected_feats = names(feat_clusters[feat_clusters == clus_id])
    sub_mat = expr_values[rownames(expr_values) %in% selected_feats,]

    # calculate mean
    if(length(selected_feats) == 1) {
      mean_score = sub_mat
    } else{
      mean_score = colMeans_flex(sub_mat)
    }

    res_list[[id]] = mean_score
  }

  res_final = data.table::as.data.table(t(do.call('rbind', res_list)))
  colnames(res_final) = as.character(sort(unique(feat_clusters)))

  # data.table variables
  cell_ID = NULL

  res_final[, cell_ID := colnames(expr_values)]


  if(return_gobject == TRUE) {

    ## enrichment scores
    spenr_names = list_spatial_enrichments_names(gobject = gobject, spat_unit = spat_unit, feat_type = feat_type)

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }


    gobject = set_spatial_enrichment(gobject = gobject,
                                     spat_unit = spat_unit,
                                     feat_type = feat_type,
                                     enrichm_name = name,
                                     spatenrichment = res_final)

    ## update parameters used ##
    gobject = update_giotto_params(gobject, description = '_create_metafeat')
    return(gobject)

  } else {
    return(res_final)
  }

}



#' @title Create metagenes
#' @name createMetagenes
#' @description This function creates an average metagene for gene clusters.
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param gene_clusters numerical vector with genes as names
#' @param name name of the metagene results
#' @param return_gobject return giotto object
#' @return giotto object
#' @details An example for the 'gene_clusters' could be like this:
#' cluster_vector = c(1, 1, 2, 2); names(cluster_vector) = c('geneA', 'geneB', 'geneC', 'geneD')
#' @export
createMetagenes = function(gobject,
                           spat_unit = NULL,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           gene_clusters,
                           name = 'metagene',
                           return_gobject = TRUE) {

  warning("Deprecated and replaced by createMetafeats")

  createMetafeats(gobject = gobject,
                  spat_unit = spat_unit,
                  feat_type = NULL,
                  expression_values = expression_values,
                  feat_clusters = gene_clusters,
                  name = name,
                  return_gobject = return_gobject)

}


#' @title Find network neighbors
#' @name findNetworkNeighbors
#' @description Find the spatial neighbors for a selected group of cells within the selected spatial network.
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param spatial_network_name name of spatial network
#' @param source_cell_ids cell ids for which you want to know the spatial neighbors
#' @param name name of the results
#' @return data.table
#' @export
findNetworkNeighbors = function(gobject,
                                spat_unit = NULL,
                                spatial_network_name,
                                source_cell_ids = NULL,
                                name = 'nb_cells') {

  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  # get spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = get_spatialNetwork(gobject,
                                         spat_unit = spat_unit,
                                         name = spatial_network_name,
                                         return_network_Obj = FALSE)
  } else {
    stop('You need to select a spatial network')
  }

  # source cell ids that are found back
  all_cell_ids = gobject@cell_ID[[spat_unit]]
  source_cells = all_cell_ids[all_cell_ids %in% source_cell_ids]

  if(length(source_cells) == 0) {
    stop('No source cell ids were selected or found')
  }


  full_network_DT = convert_to_full_spatial_network(spatial_network)
  potential_target_cells = full_network_DT[source %in% source_cells][['target']]
  source_and_target_cells = potential_target_cells[potential_target_cells %in% source_cells]
  target_cells = potential_target_cells[!potential_target_cells %in% source_and_target_cells]

  cell_meta = pDataDT(gobject)

  # data.table variables
  nb_cells = cell_ID = NULL

  cell_meta[, nb_cells := ifelse(cell_ID %in% source_and_target_cells, 'both',
                                 ifelse(cell_ID %in% source_cells, 'source',
                                        ifelse(cell_ID %in% target_cells, 'neighbor', 'others')))]
  nb_annot = cell_meta[, c('cell_ID', 'nb_cells'), with = FALSE]
  data.table::setnames(nb_annot, 'nb_cells', name)

  return(nb_annot)

}



#' @title merge_spatial_locs_feat_info
#' @name merge_spatial_locs_feat_info
#' @description merge spatial cell and feature location information
#' @keywords internal
merge_spatial_locs_feat_info = function(spatial_info,
                                        feature_info) {

  # data.table variables
  cell_ID = used = NULL

  reslist = list()
  for(i in 1:length(unique(spatial_info$cell_ID))) {

    cell_i = unique(spatial_info$cell_ID)[i]

    temp = sp::point.in.polygon(point.x = feature_info$sdimx,
                                point.y = feature_info$sdimy,
                                pol.x = spatial_info[cell_ID == cell_i]$sdimx,
                                pol.y = spatial_info[cell_ID == cell_i]$sdimy)

    detected_feats = feature_info[temp == 1]
    detected_feats[, cell_ID := cell_i]

    reslist[[i]] = detected_feats

  }

  reslistfinal = do.call('rbind', reslist)

  # calculate how often a single transcript is used
  # > 1 means that a transcript was assigned to more than 1 cell
  reslistfinal[, used := .N, by = c('sdimx', 'sdimy', 'feat_ID')]

  return(reslistfinal)

}



#' @title combineSpatialCellFeatureInfo
#' @name combineSpatialCellFeatureInfo
#' @description Combine spatial cell information (e.g. polygon)
#' and spatial feature information (e.g. transcript locations)
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type(s)
#' @param selected_features select set of features
#' @return list of data.table(s)
#' @details
#' The returned data.table has the following columns: \cr
#' \itemize{
#'   \item{sdimx: spatial feature location on the x-axis}
#'   \item{sdimy: spatial feature location on the y-axis}
#'   \item{feat_ID: unique feature ID}
#'   \item{cell_ID: unique cell ID}
#'   \item{used: how often was the feature used/assigned to a cell}
#'   \item{feat: selected feature(s)}
#' }
#' @export
combineSpatialCellFeatureInfo = function(gobject,
                                         spat_unit = NULL,
                                         feat_type = NULL,
                                         selected_features = NULL) {

  # define for data.table
  feat_ID = NULL

  # combine
  # 1. spatial morphology information ( = polygon)
  # 2. spatial transcript location information

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  spatial_cell_info = gobject@spatial_info

  if(is.null(spatial_cell_info)) {
    stop('There is no available spatial segmentation/location information')
  }


  res_list = list()
  for(feat in unique(feat_type)) {

    spatial_feat_locs = gobject@feat_info[[feat]]

    if(!is.null(selected_features)) {
      spatial_feat_locs = spatial_feat_locs[feat_ID %in% selected_features]
    }

    if(is.null(spatial_feat_locs)) {
      stop('There is no available spatial feature location information for ', feat, '\n')
    }

    output = merge_spatial_locs_feat_info(spatial_info = spatial_cell_info,
                                          feature_info = spatial_feat_locs)
    output[, 'feat' := feat]

    res_list[[feat]] = output
  }

  return(res_list)

}



#' @title combineSpatialCellMetadataInfo
#' @name combineSpatialCellMetadataInfo
#' @description Combine cell metadata with spatial cell information (e.g. polygon)
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type(s)
#' @return list of data.table(s)
#' @details
#' The returned data.table has the following columns: \cr
#' \itemize{
#'   \item{sdimx: spatial feature location on the x-axis}
#'   \item{sdimy: spatial feature location on the y-axis}
#'   \item{cell_ID: unique cell ID}
#'   \item{feat: selected feature(s)}
#'   \item{other columns that are part of the cell metadata}
#' }
#' @export
combineSpatialCellMetadataInfo = function(gobject,
                                          spat_unit = NULL,
                                          feat_type = NULL) {


  # combine
  # 1. spatial morphology information ( = polygon)
  # 2. cell metadata

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # get spatial cell information
  spatial_cell_info = gobject@spatial_info

  res_list = list()
  for(feat in unique(feat_type)) {

    # get spatial cell metadata
    cell_meta = pDataDT(gobject,
                        spat_unit = spat_unit,
                        feat_type = feat)

    # merge
    spatial_cell_info_meta = merge.data.table(spatial_cell_info, cell_meta, by = 'cell_ID')

    spatial_cell_info_meta[, 'feat' := feat]

    res_list[[feat]] = spatial_cell_info_meta

  }

  return(res_list)

}




