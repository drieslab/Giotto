#' @title Remove hetero edges from igraph
#' @name .igraph_remove_hetero_edges
#' @description
#' Given an igraph `g` and set of node attributes `clus_att` that encode
#' different spatial clusters, remove edges that connect non-similar nodes.
#' This can be used when data is already clustered, but these clusters should
#' be further broken up based on whether they are spatially touching.
#' @param g igraph object
#' @param clus_attr character. A categorical node attribute
#' @md
#' @returns igraph
#' @keywords internal
.igraph_remove_hetero_edges <- function(g, clus_attr) {
    clus_attr_values <- igraph::vertex_attr(g, name = clus_attr)

    for (n in unique(clus_attr_values)) {
        # find all vertices of the attribute
        nv <- igraph::V(g)$name[clus_attr_values == n]

        # find edges that include these vertices
        n_all_edges <- igraph::E(g)[.inc(igraph::V(g)[nv])] %>%
            igraph::as_ids()

        # find edges associated with only these vertices
        n_internal_edges <- igraph::E(g)[nv %--% nv] %>%
            igraph::as_ids()

        het_edges <- n_all_edges[!n_all_edges %in% n_internal_edges]

        g <- igraph::delete_edges(g, edges = het_edges)
    }

    g
}




#' @title igraph vertex membership
#' @name .igraph_vertex_membership
#' @description
#' Get which weakly connected set of vertices each vertex is part of
#' @param g igraph
#' @param clus_name character. name to assign column of clustering info
#' @param all_ids (optional) character vector with all ids
#' @param missing_id_name character and name for vertices that were missing from g
#' @returns data.table
#' @keywords internal
.igraph_vertex_membership <- function(g, 
                                      clus_name, 
                                      all_ids = NULL, 
                                      missing_id_name) {
    
  # get membership
  membership <- igraph::components(g)$membership %>%
    data.table::as.data.table(keep.rownames = TRUE)
  data.table::setnames(membership, c("cell_ID", clus_name))
  
  # add vertices that were missing from g back 
  if(!is.null(all_ids)) {
    missing_ids = all_ids[!all_ids %in% igraph::V(g)$name]
    missing_membership = data.table::data.table('cell_ID' = missing_ids, 'cluster_name' = missing_id_name)
    data.table::setnames(missing_membership, c("cell_ID", clus_name))
    membership = data.table::rbindlist(list(membership, missing_membership))
  }
  
  return(membership)
  
}





#' @title Split cluster annotations based on a spatial network
#' @name spatialSplitCluster
#' @inheritParams data_access_params
#' @param spatial_network_name character. Name of spatial network to use
#' @param cluster_col character. Column in metadata containing original
#' clustering
#' @param split_clus_name character. Name to assign the split cluster results
#' @param include_all_ids Boolean. Include all ids, including vertex ids not found
#' in the spatial network
#' @param missing_id_name Character. Name for vertices that were missing from 
#' spatial network
#' @param return_gobject Boolean. Return giotto object
#' @returns giotto object with cluster annotations
#' @examples
#' library(Giotto)
#' g <- GiottoData::loadGiottoMini("vizgen")
#' activeSpatUnit(g) <- "aggregate"
#' spatPlot2D(g, cell_color = "leiden_clus")
#'
#' g <- spatialSplitCluster(g,
#'     cluster_col = "leiden_clus",
#'     split_clus_name = "new"
#' )
#' # don't show legend since there are too many categories generated
#' spatPlot2D(g, cell_color = "new", show_legend = FALSE)
#' @export
spatialSplitCluster <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spatial_network_name = "Delaunay_network",
        cluster_col,
        split_clus_name = paste0(cluster_col, "_split"), 
        include_all_ids = TRUE,
        missing_id_name = 'not_connected',
        return_gobject = TRUE) {
  
    # NSE vars
    cell_ID <- NULL

    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    cell_meta <- getCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = FALSE
    )

    sn <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "spatialNetworkObj",
        copy_obj = FALSE,
        verbose = FALSE,
    )

    clus_info <- cell_meta[, c("cell_ID", cluster_col), with = FALSE]
    # subset to needed cols
    g <- GiottoClass::spat_net_to_igraph(sn)
    # convert spatialNetworkObject to igraph

    # assign cluster info to igraph nodes
    clus_values <- clus_info[
        match(igraph::V(g)$name, cell_ID), get(cluster_col)
    ]
    igraph::V(g)$cluster <- clus_values

    # split cluster by spatial igraph
    g <- .igraph_remove_hetero_edges(
        g = g,
        clus_attr = "cluster"
    )

    # get new clusterings
    if(isTRUE(include_all_ids)) {
      # include all cell IDs
      all_ids = unique(cell_meta$cell_ID)
      new_clus_dt <- .igraph_vertex_membership(
        g = g,
        clus_name = split_clus_name,
        all_ids = all_ids, 
        missing_id_name = missing_id_name
      )
    } else {
      # only IDs present in graph
      new_clus_dt <- .igraph_vertex_membership(
        g = g,
        clus_name = split_clus_name,
        all_ids = NULL
      )
      
    }
   
    if(isTRUE(return_gobject)) {
      gobject <- addCellMetadata(
        gobject,
        spat_unit = spat_unit,
        new_metadata = new_clus_dt,
        by_column = TRUE,
        column_cell_ID = "cell_ID"
      )
      return(gobject)
    } else {
      new_clus_dt
    }
   
}





#' @title Split cluster annotations based on a spatial network
#' @name identifyTMAcores
#' @inheritParams data_access_params
#' @param spatial_network_name character. Name of spatial network to use
#' @param core_id_name metadata column name for the core information
#' @param include_all_ids Boolean. Include all ids, including vertex ids not found
#' in the spatial network
#' @param missing_id_name Character. Name for vertices that were missing from 
#' spatial network
#' @param return_gobject Boolean. Return giotto object
#' @returns cluster annotations
#' @export
identifyTMAcores <- function(gobject,
                             spat_unit = NULL,
                             feat_type = NULL,
                             spatial_network_name = "Delaunay_network",
                             core_id_name = 'core_id',
                             include_all_ids = TRUE,
                             missing_id_name = 'not_connected',
                             return_gobject = TRUE) {
  
  
  # NSE vars
  cell_ID <- NULL
  
  spat_unit <- set_default_spat_unit(
    gobject = gobject,
    spat_unit = spat_unit
  )
  feat_type <- set_default_feat_type(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type
  )
  
  cell_meta <- getCellMetadata(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type,
    output = "data.table",
    copy_obj = FALSE
  )
  
  sn <- getSpatialNetwork(
    gobject = gobject,
    spat_unit = spat_unit,
    name = spatial_network_name,
    output = "spatialNetworkObj",
    copy_obj = FALSE,
    verbose = FALSE,
  )
  
 
  g <- GiottoClass::spat_net_to_igraph(sn)
  # convert spatialNetworkObject to igraph
  
  
  # get new clusterings
  if(isTRUE(include_all_ids)) {
    # include all cell IDs
    all_ids = unique(cell_meta$cell_ID)
    new_clus_dt <- .igraph_vertex_membership(
      g = g,
      clus_name = core_id_name,
      all_ids = all_ids, 
      missing_id_name = missing_id_name
    )
  } else {
    # only IDs present in graph
    new_clus_dt <- .igraph_vertex_membership(
      g = g,
      clus_name = core_id_name,
      all_ids = NULL
    )
    
  }
  
  if(isTRUE(return_gobject)) {
    gobject <- addCellMetadata(
      gobject,
      spat_unit = spat_unit,
      new_metadata = new_clus_dt,
      by_column = TRUE,
      column_cell_ID = "cell_ID"
    )
    return(gobject)
  } else {
    new_clus_dt
  }
  
  
}





