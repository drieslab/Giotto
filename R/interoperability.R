
## sp classes ####



## gef object ####
#' @title Convert gef to Giotto
#' @name gefToGiotto
#' @description Converts .gef file (output stereo-seq pipeline) into
#' giotto subcellular object
#' @param gef_file path to .gef file
#' @param bin_size bin size to select from .gef file
#' @param verbose be verbose
#' @details Function in beta. Converts .gef object to Giotto object.
#'
#' There are six possible choices for bin_size: 1, 10, 20, 50, 100, 200.
#'
#' See SAW pipeline for additional information about the gef file.
#' @export
gefToGiotto = function(gef_file, bin_size = 'bin100', verbose = TRUE){

   # data.table vars
   genes = NULL

   # package check
   package_check(pkg_name = 'rhdf5', repository = 'Bioc')
   if(!file.exists(gef_file)) stop('File path to .gef file does not exist')

   # check if proper bin_size is selected
   bin_size_options = c('bin1', 'bin10', 'bin20', 'bin50', 'bin100', 'bin200')
   if(!(bin_size %in% bin_size_options)) stop('Please select valid bin size,
                                              see details for choices.')

   # step 1: read expression and gene data from gef file
   geneExpData = rhdf5::h5read(file = gef_file, name = 'geneExp')
   exprDT = data.table::as.data.table(geneExpData[[bin_size]][['expression']])
   geneDT = data.table::as.data.table(geneExpData[[bin_size]][['gene']])
   if(isTRUE(verbose)) wrap_msg('read in .gef file')

   # step 2: combine gene information from the geneDT to the exprDT
   exprDT[, genes := rep(x = geneDT$gene, geneDT$count)]

   # step 3: bin coordinates according to selected bin_size
   bin_size_int = as.integer(gsub("[^0-9.-]", "", bin_size))
   #TODO: update bin_shift for other shapes, not just rect_vertices
   bin_shift = ceiling(bin_size_int / 2) # ceiling catches bin_1
   bincoord = unique(exprDT[,.(x,y)])
   setorder(bincoord, x, y)
   setnames(bincoord, old = c('x', 'y'), new = c('sdimx', 'sdimy'))
   bincoord[, c('sdimx', 'sdimy') := list(sdimx+bin_shift, sdimy+bin_shift)]
   bincoord[, cell_ID := paste0('bin', 1:.N)]
   tx_data = exprDT[,.(genes, x, y, count)]
   tx_data[, c('x', 'y') := list(x+bin_shift, y+bin_shift)]
   if(isTRUE(verbose)) wrap_msg('shift and bin coordinates')

   # step 4: create rectangular polygons (grid) starting from the bin centroids
   x = polyStamp(stamp_dt = rectVertices(dims = c(x = (bin_size_int - 1),
                                                  y = (bin_size_int - 1))),
                 spatlocs = bincoord[,.(cell_ID, sdimx, sdimy)])
   pg = createGiottoPolygonsFromDfr(x)
   if(isTRUE(verbose)) wrap_msg('create polygon stamp')

   # step 5: create giotto subcellular object
   stereo = createGiottoObjectSubcellular(
      gpoints = list(rna = tx_data),
      gpolygons = list(cell = pg)
   )

   stereo = addSpatialCentroidLocations(gobject = stereo)
   if(isTRUE(verbose)) wrap_msg('giotto subcellular object created')

   return(stereo)
}



## anndata object ####

#' @title Convert anndata to Giotto
#' @name anndataToGiotto
#' @description Converts a spatial anndata (e.g. scanpy) .h5ad file into a Giotto object
#' @param anndata_path path to the .h5ad file
#' @param n_key_added equivalent of "key_added" argument from scanpy.pp.neighbors().
#'                    Cannot be "spatial". This becomes the name of the nearest network in the gobject.
#' @param spat_unit desired spatial unit for conversion, default NULL
#' @param feat_type desired feature type for conversion, default NULL
#' @param python_path path to python executable within a conda/miniconda environment
#' @return Giotto object
#' @details Function in beta. Converts a .h5ad file into a Giotto object.
#'    The returned Giotto Object will take default insructions with the
#'    exception of the python path, which may be customized.
#'    See \code{\link{changeGiottoInstructions}} to modify instructions after creation.
#' @export
anndataToGiotto = function(anndata_path = NULL,
                           n_key_added = NULL,
                           spatial_key_added = NULL,
                           spat_unit = NULL,
                           feat_type = NULL,
                           python_path = NULL) {

  # Preliminary file checks and guard clauses
  if (is.null(anndata_path)) {
    stop("Please provide a path to an AnnData .h5ad file for conversion.\n")
  }

  if(!file.exists(anndata_path)) {
    stop("The provided path to the AnnData .h5ad file does not exist.\n")
  }

  # Required step to properly initialize reticualte
  instrs = createGiottoInstructions(python_path = python_path)

  # test if scanpy is found
  module_test = reticulate::py_module_available('scanpy')
  py_path = reticulate::py_config()$python
  genv_in_use = grepl(pattern = "giotto_env", x = py_path)

  if(module_test == FALSE && !genv_in_use) {
    warning("scanpy python module is not installed:
            install in the environment or python path with:

            'pip install scanpy==1.9.0'

            Alternatively, install in the active python
            environment with reticulate:

            reticulate::py_install(packages = 'scanpy==1.9.0',
                                   pip = TRUE)
            \n
            ")
  } else if (module_test == FALSE && genv_in_use) {
    cat ("Python module scanpy is required for conversion.
          Installing scanpy now in the Giotto Miniconda Environment.\n")

    conda_path = reticulate::miniconda_path()
    py_ver = reticulate::py_config()$version_string
    py_ver = strsplit(py_ver,"|", fixed = T)[[1]][1]
    py_ver = gsub(" ","",py_ver, fixed = T)
    conda_full_path = paste0(conda_path,'/','bin/conda')
    full_envname = paste0(conda_path,'/envs/giotto_env')

    reticulate::py_install(packages = "scanpy==1.9.0",
                        envname = full_envname,
                        method = 'conda',
                        conda = conda_full_path,
                        pip = TRUE,
                        python_version = py_ver)
  } else cat ("Required Python module scanpy has been previously installed. Proceeding with conversion.\n")

  # Import ad2g, a python module for parsing anndata
  ad2g_path <- system.file("python","ad2g.py",package="Giotto")
  reticulate::source_python(ad2g_path)
  adata <- read_anndata_from_path(anndata_path)

  ### Set up expression matrix
  X <- extract_expression(adata)
  cID = extract_cell_IDs(adata)
  fID = extract_feat_IDs(adata)
  X@Dimnames[[1]] = fID
  X@Dimnames[[2]] = cID
  # Expression matrix X ready

  ### Set up spatial info
  sp = parse_obsm_for_spat_locs(adata)
  #Spatial locations sp ready

  ### Set up metadata
  cm = extract_cell_metadata(adata)
  cm = as.data.table(cm)
  if ('leiden' %in% names(cm)) {
    cm$leiden = as.numeric(cm$leiden)
  }

  fm = extract_feat_metadata(adata)
  fm = as.data.table(fm)
  # Metadata ready

  ### Create Minimal giottoObject
  gobject <- createGiottoObject(expression = X,
                                spatial_locs = sp,
                                instructions = instrs)

  ### Add metadata
  gobject = set_cell_metadata(gobject,
                              metadata = cm)
  gobject = set_feature_metadata(gobject,
                              metadata = fm)

  spat_unit = set_default_spat_unit(gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit)

  ### Set up PCA
  p = extract_pca(adata)
  if (!is.null(p)) {
    pca = p$pca
    evs = p$eigenvalues
    loads = p$loadings
    # Add PCA to giottoObject
    dobj = create_dimObject(name = 'pca.ad',
                            spat_unit = spat_unit,
                            feat_type = feat_type,
                            reduction_method = 'pca',
                            coordinates = pca,
                            misc = list(eigenvalues = evs,
                                        loadings = loads),
                            my_rownames = colnames(X))

    gobject = set_dimReduction(gobject = gobject,
                          spat_unit = spat_unit,
                          feat_type = feat_type,
                          reduction = 'cells',
                          reduction_method = 'pca',
                          name = 'pca.ad',
                          dimObject = dobj)

  }

  ### Set up UMAP
  u = extract_umap(adata)
  if (!is.null(u)) {
    # Add UMAP to giottoObject
    dobj = create_dimObject(name = 'umap.ad',
                            spat_unit = spat_unit,
                            feat_type = feat_type,
                            reduction_method = 'umap',
                            coordinates = u,
                            misc = NULL,
                            my_rownames = colnames(X))

    gobject = set_dimReduction(gobject = gobject,
                          spat_unit = spat_unit,
                          feat_type = feat_type,
                          reduction = 'cells',
                          reduction_method = 'umap',
                          name = 'umap.ad',
                          dimObject = dobj)

  }
  ### Set up TSNE
  t = extract_tsne(adata)
  if (!is.null(t)) {
    # Add UMAP to giottoObject
    dobj = create_dimObject(name = 'tsne.ad',
                            spat_unit = spat_unit,
                            feat_type = feat_type,
                            reduction_method = 'tsne',
                            coordinates = t,
                            misc = NULL,
                            my_rownames = colnames(X))

    gobject = set_dimReduction(gobject = gobject,
                          spat_unit = spat_unit,
                          feat_type = feat_type,
                          reduction = 'cells',
                          reduction_method = 'tsne',
                          name = 'tsne.ad',
                          dimObject = dobj)
  }

  ### NN Network

  # Need to create nnNetObj or igraph object to use with setter for NN

  weights_ad = NULL
  weights_ad = extract_NN_connectivities(adata, key_added = n_key_added)
  if (!is.null(weights_ad)) {
    distances_ad = extract_NN_distances(adata, key_added = n_key_added)
    ij_matrix = as(distances_ad, "TsparseMatrix")
    from_idx = ij_matrix@i + 1 #zero index!!!
    to_idx = ij_matrix@j + 1 #zero index!!!

    #pre-allocate DT variables
    from = to = weight = distance = from_cell_ID = to_cell_ID = uniq_ID = NULL
    nn_dt = data.table::data.table(from = from_idx,
                                  to = to_idx,
                                  weight = weights_ad@x,
                                  distance = distances_ad@x)

    nn_dt[, from_cell_ID := cID[from]]
    nn_dt[, to_cell_ID := cID[to]]
    nn_dt[, uniq_ID := paste0(from,to)]
    nn_dt[order(uniq_ID)]
    nn_dt[,uniq_ID := NULL]
    vert = unique(x = c(nn_dt$from_cell_ID, nn_dt$to_cell_ID))
    nn_network_igraph = igraph::graph_from_data_frame(nn_dt[,.(from_cell_ID, to_cell_ID, weight, distance)], directed = TRUE, vertices = vert)
    # og_nn = get_NearestNetwork(SS_seqfish, spat_unit = "cell", feat_type = "rna", nn_network_to_use = "kNN")
    # as written, no difference found between original and above from igraph::difference()

    nn_info = extract_NN_info(adata = adata, key_added = n_key_added)
    browser()
    net_type = "kNN" # anndata default
    if(("sNN" %in% n_key_added) & !is.null(n_key_added)){
      net_type = "sNN"
      net_name = paste0(n_key_added, ".", nn_info["method"])
    } else if (!("sNN" %in% n_key_added) & !is.null(n_key_added)) {
      net_name = paste0(n_key_added, ".", nn_info["method"])
    } else {
      net_name = paste0(net_type, ".", nn_info["method"])
    }

    gobject = set_NearestNetwork(gobject = gobject,
                                nn_network = nn_network_igraph,
                                spat_unit = spat_unit,
                                feat_type = feat_type,
                                nn_network_to_use = net_type,
                                network_name = net_name,
                                set_defaults = FALSE)
  }

  ### Layers
  lay_names = extract_layer_names(adata)
  if (!is.null(lay_names)) {
    for (l_n in lay_names){
      lay = extract_layered_data(adata, layer_name = l_n)
      if ("data.frame" %in% class(lay)){
        names(lay) = fID
        row.names(lay) = cID
      }
      else{
        lay@Dimnames[[1]] = fID
        lay@Dimnames[[2]] = cID
      }
      gobject = set_expression_values(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  name = l_n,
                                  values = lay)
    }
  }

  gobject <- update_giotto_params(gobject = gobject,
                                  description = "_AnnData_Conversion")

  wrap_msg("\nAnnData object successfully converted to Giotto.\n")
  return(gobject)

}


#' @title Convert Giotto to anndata
#' @name giottoToAnnData
#' @description Converts a Giotto object to a spatial anndata (e.g. scanpy) .h5ad file
#' @param gobject giotto object to be converted
#' @param spat_unit spatial unit which will be used in conversion.
#' @param feat_type feature type which will be used in conversion.
#' @param python_path path to python executable within a conda/miniconda environment
#' @param save_directory directory in which the file will be saved.
#' @return vector containing .h5ad file path(s)
#' @details Function in beta. Converts a Giotto object into .h5ad file(s).
#'
#' If there are multiple spatial units and/or feature types, but only
#' one spatial unit and/or feature type is specified, then only the
#' specified spatial unit and/or feature type will be used. If NULL,
#' by default, all spatial units will be used in conversion.
#'
#' If multiple spatial units or feature types are specified, multiple
#' AnnData object will be created and returned.
#'
#' The save_directory will be created if it does not already exist.
#' The default save_directory is the working directory.
#' @export
giottoToAnnData <- function(gobject = NULL,
                            spat_unit = NULL,
                            feat_type = NULL,
                            python_path = NULL,
                            save_directory = NULL){
  # Check gobject
  invalid_obj = !("giotto" %in% class(gobject))
  if (is.null(gobject) || invalid_obj) {
    stop(wrap_msg("Please provide a valid Giotto Object for conversion."))
  }

  # Python module import
  g2ad_path <- system.file("python","g2ad.py",package="Giotto")
  reticulate::source_python(g2ad_path)
  if (!is.null(save_directory)) dir_guard(save_directory)

  # Check directory, make it if it doesn't exist
  if (is.null(save_directory)) save_directory = paste0(getwd(),"/")
  else if (!dir.exists(save_directory)) {
    warning(wrap_msg("Provided save directory not found. Creating save directory at location:"))
    cat(save_directory)
    dir.create(save_directory, recursive = TRUE)
    if (dir.exists(save_directory)) cat("Created directory", save_directory)
    else stop(wrap_msg("Unable to create directory. Please change the provided path and try again."))
  }
  else {
    wrap_msg("Directory", save_directory,"found. The converted Giotto object will be saved here as a .h5ad file.")
  }

  # Expresion
  expr_dt <- list_expression(gobject)
  # ID spat_unit and feat_type if not already provided.
  if (is.null(spat_unit) && is.null(feat_type)) {
    spat_unit = unique(expr_dt$spat_unit)
    feat_type = unique(expr_dt$feat_type)
  } else if (is.null(spat_unit && !is.null(feat_type))) {
    spat_unit = unique(expr_dt$spat_unit)
  } else if (!is.null(spat_unit && is.null(feat_type))) {
    feat_type = unique(expr_dt$feat_type)
  }

  for (su in spat_unit) wrap_msg("Spatial unit(s)", su, "will be used in conversion.")
  for (ft in feat_type) wrap_msg("Feature type(s)", ft, "will be used in conversion.")

  # Iterate through spat_unit and feat_type to pull out expression data.
  # By default, the raw expression in the slot of the first spatial unit
  # and first feature type (if multiple are present) will be transferred to
  # the AnnData.anndata.X slot
  # Any other expression data will be inserted into AnnData.anndata.layers
  # By default, layer names are formed by "'spatial_unit'_'feature_type'_'value'"

  adata = NULL #scope
  su_ft_length = 0

  for (su in spat_unit) {
    for (ft in names(gobject@expression[[su]])) {
      su_ft_length = su_ft_length + 1
    }
  }


  adata_list = lapply(1:su_ft_length, function(i) adata)
  adata_pos = 1

  for (su in spat_unit) {
    for (ft in names(gobject@expression[[su]])) {
      expr_names = list_expression_names(gobject = gobject,
                                         spat_unit = su,
                                         feat_type = ft)

      for (en in expr_names) {
        if (en == "raw") {
          raw_x = get_expression_values(gobject = gobject,
                                values = en,
                                spat_unit = su,
                                feat_type = ft,
                                output = "matrix")

          adata = ad_obj(x = raw_x)
        } else {
          ad_layer_name = paste0(su,"_",ft,"_",en)

          x = get_expression_values(gobject = gobject,
                                    values = en,
                                    spat_unit = su,
                                    feat_type = ft,
                                    output = "matrix")

          if ("dgeMatrix" %in% class(x)) x = as(x,'array')

          adata = set_adg_layer_data(adata = adata,
                                     lay = x,
                                     lay_name = ad_layer_name)
        }

      }
      adata_list[[adata_pos]] = adata
      adata_pos = adata_pos + 1
      adata = NULL
    }
  }
  # Reset indexing variable
  adata_pos = 1

  # Spatial Locations
  for (su in spat_unit){
    for (ft_ in names(gobject@expression[[su]])) {
      sl = get_spatial_locations(gobject = gobject,
                                 output = "data.table",
                                 spat_unit = su)
      n_col_sl = dim(sl)[2]

      #preallocate data.table params
      sdimx = sdimy = sdimz = NULL

      if (n_col_sl == 3){
        sl = sl[, .(sdimx, sdimy)]
      } else {
        sl = sl[, .(sdimx, sdimy, sdimz)]
      }
      adata = adata_list[[adata_pos]]
      adata = set_adg_spat_locs(adata = adata,
                                spat_locs = sl)
      adata_pos = adata_pos + 1
    }
  }
  # Reset indexing variable
  adata_pos = 1

  # Spatial Info

  # Cell Metadata
  # Feat Metadata
  for (su in spat_unit){
    for (ft in names(gobject@expression[[su]])){

      cm = get_cell_metadata(gobject = gobject,
                             spat_unit = su,
                             feat_type = ft,
                             output = "data.table",
                             set_defaults = FALSE)

      fm = get_feature_metadata(gobject = gobject,
                                spat_unit = su,
                                feat_type = ft,
                                output = "data.table",
                                set_defaults = FALSE)

      adata_list[[adata_pos]] = set_adg_metadata(adata = adata_list[[adata_pos]],
                                                 cell_meta = cm,
                                                 feat_meta = fm)

      adata_pos = adata_pos + 1

    }
  }
  # Reset indexing variable
  adata_pos = 1

  # Dimension Reductions

  # error hanldling wrapper to get_dimReduction
  try_get_dimReduction = function(gobject,
                                  spat_unit,
                                  feat_type,
                                  reduction,
                                  reduction_method,
                                  name,
                                  output,
                                  set_defaults) {
  tryCatch(
    {
    dim_red = get_dimReduction(gobject = gobject,
                               spat_unit = spat_unit,
                               feat_type = feat_type,
                               reduction = reduction,
                               reduction_method = reduction_method,
                               name = name,
                               output = output,
                               set_defaults = set_defaults)
    return(dim_red)
    },
    error = function(e) {
      return(NULL)
    }
    )
  }

  ## PCA

  # pca on feats not supported by anndata because of dimensionality agreement reqs
  reduction_options = names(gobject@dimension_reduction)
  dim_red = NULL

  for (ro in reduction_options) {
    if (ro != "cells") {
      warning("AnnData does not support storing PCA by features. Skipping PCA data conversion.")
      break
    }
    for (su in spat_unit) {
      for (ft in names(gobject@expression[[su]])) {
        name = "pca"
        if (ft != "rna") name = paste0(ft,".pca")
        dim_red = try_get_dimReduction(gobject = gobject,
                                       spat_unit = su,
                                       feat_type = ft,
                                       reduction = ro,
                                       reduction_method = "pca",
                                       name = name,
                                       output = "dimObj",
                                       set_defaults = FALSE)
        if (is.null(dim_red)){
          adata_pos = adata_pos + 1
          next
        }
        pca_coord = dim_red[]
        pca_loadings = data.table(dim_red@misc$loadings)
        feats_used = dimnames(dim_red@misc$loadings)[[1]]
        evs = dim_red@misc$eigenvalues

        adata_list[[adata_pos]] = set_adg_pca(adata = adata_list[[adata_pos]],
                                              pca_coord = pca_coord,
                                              loadings = pca_loadings,
                                              eigenv = evs,
                                              feats_used = feats_used)
        adata_pos = adata_pos + 1
      }
    }
    adata_pos = 1

  }

  # Reset indexing variable
  adata_pos = 1

  ## UMAP
  for (ro in reduction_options) {
    for (su in spat_unit) {
      for (ft in names(gobject@expression[[su]])) {
        name = "umap"
        if (ft != "rna") name = paste0(ft,".umap")
        dim_red = try_get_dimReduction(gobject = gobject,
                                       spat_unit = su,
                                       feat_type = ft,
                                       reduction = ro,
                                       reduction_method = "umap",
                                       name = name,
                                       output = "dimObj",
                                       set_defaults = FALSE)

        if (is.null(dim_red)) {
          adata_pos = adata_pos + 1
          next
        }
        umap_data = dim_red[]
        adata_list[[adata_pos]] = set_adg_umap(adata = adata_list[[adata_pos]],
                                               umap_data = umap_data)
        adata_pos = adata_pos + 1
      }
    }
    # Reset indexing variable
    adata_pos = 1
  }

  # Reset indexing variable
  adata_pos = 1

  ## T-SNE
  for (ro in reduction_options) {
    for (su in spat_unit) {
      for (ft in names(gobject@expression[[su]])) {
        name = "tsne"
        if (ft != "rna") name = paste0(ft,".tsne")
        dim_red = try_get_dimReduction(gobject = gobject,
                                       spat_unit = su,
                                       feat_type = ft,
                                       reduction = ro,
                                       reduction_method = "tsne",
                                       name = name,
                                       output = "dimObj",
                                       set_defaults = FALSE)

        if (is.null(dim_red)) {
          adata_pos = adata_pos + 1
          next
        }
        tsne_data = dim_red[]
        adata_list[[adata_pos]] = set_adg_tsne(adata = adata_list[[adata_pos]],
                                               tsne_data = tsne_data)
        adata_pos = adata_pos + 1
      }
    }
    # Reset indexing variable
    adata_pos = 1
  }

  # Reset indexing variable
  adata_pos = 1

  # Nearest Neighbor Network

  # error hanldling wrapper to get_NearestNetwork
  try_get_NN = function(gobject,
                        spat_unit,
                        feat_type,
                        nn_network_to_use,
                        network_name,
                        output,
                        set_defaults) {
  tryCatch(
    {
    nearest_net = get_NearestNetwork(gobject = gobject,
                                     spat_unit = spat_unit,
                                     feat_type = feat_type,
                                     nn_network_to_use = nn_network_to_use,
                                     network_name = network_name,
                                     output = output,
                                     set_defaults = set_defaults)
    return(nearest_net)
    },
    error = function(e) {
      return(NULL)
    }
    )
  }

  for (su in spat_unit) {
    for (ft in names(gobject@expression[[su]])) {
      nn_network_to_use = c("sNN","kNN")
      for (nn_net_tu in nn_network_to_use) {
        network_name = list_nearest_networks_names(gobject = gobject,
                                                   spat_unit = su,
                                                   feat_type = ft,
                                                   nn_type = nn_net_tu)
        if (is.null(network_name)) {
          next
        }
        for (n_name in network_name) {
          gob_NN = try_get_NN(gobject = gobject,
                              spat_unit = su,
                              feat_type = ft,
                              nn_network_to_use = nn_net_tu,
                              network_name = n_name,
                              output = "nnNetObj",
                              set_defaults = FALSE)

          pidx = grep("nn_network", names(gobject@parameters))
          for (p in pidx) {
            if (gobject@parameters[[p]]["type"] == nn_net_tu) {
              kval = gobject@parameters[[p]]["k"]
              dim_red_used = gobject@parameters[[p]]["dim_red_to_use"]
            }
          }

          df_gob_NN = igraph::as_data_frame(gob_NN[])

          adata_list[[adata_pos]] = set_adg_nn(adata = adata_list[[adata_pos]],
                                               df_NN = df_gob_NN,
                                               net_name = n_name,
                                               n_neighbors = kval,
                                               dim_red_used = dim_red_used)

        }

      }
      adata_pos = adata_pos + 1
    }
  }

  # Reset indexing variable
  adata_pos = 1

  # Pipe non-expression data into AnnData object

  # Write AnnData object to .h5ad file
  # Verify it exists, and return upon success
  fname_list = lapply(1:su_ft_length, function(i) NULL)
  wrap_msg("\n")
  for (su in spat_unit) {
    for (ft in names(gobject@expression[[su]])) {
      adata = adata_list[[adata_pos]]
      path_adata = write_ad_h5ad(adata = adata,
                                 save_directory = save_directory,
                                 spat_unit = su,
                                 feat_type = ft)
      if (!is.null(path_adata)) {
        wrap_msg("Spatial unit", su, "and feature type", ft, "converted to:")
        wrap_msg(path_adata)
        fname_list[[adata_pos]] = path_adata
        adata_pos = adata_pos + 1
      } else {
        wrap_msg("Unable to convert spatial unit feature type pair", su, ft)
        stop(wrap_msg("File writing error. Please try again."))
      }
    }
  }

  wrap_msg("\nGiotto object successfully converted to .h5ad file(s)\n")

  return(fname_list)
}



## Seurat object ####


#' @title Convert Giotto to Seurat
#' @name giottoToSeurat
#' @description Converts Giotto object into a Seurat object. This functions extracts
#' specific sets of data belonging to specified spatial unit.
#' The default values are 'cell' and 'rna' respectively.
#' @param gobject Giotto object
#' @param obj_use Giotto object (deprecated, use gobject)
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param feat_type feature type (e.g. 'rna' or 'protein')
#' @param ... additional params to pass to \code{\link{get_spatial_locations}}
#' @return Seurat object
#' @export
giottoToSeurat <- function(gobject,
                           spat_unit = NULL,
                           obj_use = NULL,
                           ...){

  # data.table vars
  feat_type = name = dim_type = nn_type = NULL

  if(!is.null(obj_use)) {
    warning('obj_use param is deprecated. Please use "gobject')
    gobject = obj_use
  }

  # set default spat_unit and feat_type to be extracted as a Seurat assay
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  # verify if optional package is installed
  package_check(pkg_name = "Seurat", repository = "CRAN")
  requireNamespace('Seurat')

  # check whether any raw data exist -- required for Seurat
  avail_expr = list_expression(gobject = gobject, spat_unit = spat_unit)
  raw_exist = avail_expr[, 'raw' %in% name, by = feat_type]
  # raw_exist <- sapply(gobject@expression_feat,function(x)
  #   'raw' %in% names(gobject@expression[[x]]))
  if (nrow(raw_exist) > 0){
    assays_all = raw_exist[, feat_type]
    # assays_all <- names(raw_exist[1])
    # assays_all <- union(assays_all,gobject@expression_feat)
  } else {
    stop("Raw count data not found. Required for Seurat object.")
  }

  # create Seurat object when at least one raw data is available
  for (i in seq_along(assays_all)){
    assay_use <- assays_all[i]
    expr_use = lapply(avail_expr[feat_type == assay_use, name],
                      function(x) {
                        get_expression_values(gobject = gobject,
                                              spat_unit = spat_unit,
                                              feat_type = assay_use,
                                              values = x,
                                              output = 'exprObj')
                      })
    # expr_use <- gobject@expression[[assay_use]]
    names(expr_use) = unlist(lapply(expr_use, objName))
    slot_use <- names(expr_use)
    if (i == 1){
      data_raw <- expr_use[['raw']][]
      sobj <- Seurat::CreateSeuratObject(counts = data_raw,
                                         assay = assay_use)
      if ('normalized' %in% slot_use){
        sobj <- Seurat::SetAssayData(sobj,slot = 'data',
                                     new.data = expr_use[['normalized']][],
                                     assay = assay_use)
      }
      if ('scaled' %in% slot_use){
        sobj <- Seurat::SetAssayData(sobj,slot = 'scale.data',
                                     new.data = expr_use[['scaled']][],
                                     assay = assay_use)
      }
    } else {
      if ('raw' %in% slot_use){
        data_raw <- expr_use[['raw']][]
        flag_raw <- 1
      } else {
        flag_raw <- 0
      }
      if ('normalized' %in% slot_use){
        data_norm <- expr_use[['normalized']][]
        flag_norm <- 1
      } else {
        flag_norm <- 0
      }
      if (flag_raw == 1){
        assay_obj <- Seurat::CreateAssayObject(counts = data_raw)
      } else if (flag_raw == 0 & flag_norm == 1){
        assay_obj <- Seurat::CreateAssayObject(data = data_norm)
      } else {
        stop(paste0('Raw and normalized data not found for assay ',assay_use))
      }
      sobj[[assay_use]] <- assay_obj
      if ('scaled' %in% slot_use){
        data_scale <- expr_use[['scaled']][]
        sobj <- Seurat::SetAssayData(sobj,slot = "scale.data",
                                     new.data = data_scale,
                                     assay = assay_use)
      }
    }

    # add cell metadata
    meta_cells <- data.table::setDF(
      get_cell_metadata(
        gobject = gobject,
        feat_type = assay_use,
        output = 'data.table',
        copy_obj = TRUE
      )
    )
    rownames(meta_cells) <- meta_cells$cell_ID
    meta_cells <- meta_cells[,-which(colnames(meta_cells) == 'cell_ID')]
    if(ncol(meta_cells) > 0) {
      colnames(meta_cells) <- paste0(assay_use,"_",colnames(meta_cells))
    }
    sobj <- Seurat::AddMetaData(sobj,metadata = meta_cells[Seurat::Cells(sobj),])

    # add feature metadata
    meta_genes <- data.table::setDF(
      get_feature_metadata(
        gobject = gobject,
        feat_type = assay_use,
        output = 'data.table',
        copy_obj = TRUE
      )
    )
    rownames(meta_genes) <- meta_genes$feat_ID
    sobj[[assay_use]]@meta.features <- cbind(sobj[[assay_use]]@meta.features,meta_genes)


    # dim reduction
    # note: Seurat requires assay name specification for each dim reduc
    avail_dr = list_dim_reductions(gobject = gobject, spat_unit = spat_unit, feat_type = assay_use,)
    if (nrow(avail_dr) > 0){
      dr_use = avail_dr[, dim_type]
      # dr_use <- names(gobject@dimension_reduction[[1]])
      for (i in dr_use){
        dr_type = i
        dr_obj <- get_dimReduction(
          gobject = gobject,
          output = 'dimObj',
          spat_unit = spat_unit,
          feat_type = assay_use,
          reduction_method = dr_type  # allow default dim red name
        )
        emb_use <- dr_obj[][Seurat::Cells(sobj),]
        if (sum(c('loadings','eigenvalues') %in% names(slot(dr_obj, 'misc'))) == 2){
          loadings_use <- slot(dr_obj, 'misc')$loadings
          stdev_use <- slot(dr_obj, 'misc')$eigenvalues
          sobj[[dr_type]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(emb_use),
                                                          loadings = loadings_use,
                                                          key = dr_type,
                                                          stdev = stdev_use,
                                                          assay = assay_use)
        } else {
          sobj[[dr_type]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(emb_use),
                                                          key = dr_type,
                                                          assay = assay_use)
        }
      }
    }



    # network objects
    # expression network
    avail_nn = list_nearest_networks(gobject = gobject, spat_unit = spat_unit, feat_type = assay_use)
    nn_all = avail_nn[, nn_type]
    # nn_all <- names(gobject@nn_network)
    if (nrow(avail_nn) > 0){
      for (i in nn_all){
        nn_use <- get_NearestNetwork(
          gobject = gobject,
          spat_unit = spat_unit,
          feat_type = assay_use,
          nn_network_to_use = i,
          output = 'data.table' # allow default network name selection
        )
        idx1 <- match(nn_use$from,Seurat::Cells(sobj))
        idx2 <- match(nn_use$to,Seurat::Cells(sobj))
        edge_weight <- nn_use$weight
        nn_mtx <- Matrix::sparseMatrix(i = idx1,j = idx2,x = edge_weight,dims = c(ncol(sobj),ncol(sobj)))
        rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
        nn_name <- paste0('expr_',i)
        sobj[[nn_name]] <- Seurat::as.Graph(nn_mtx)
        sobj[[nn_name]]@assay.used <- assay_use
      }
    }

  }

  # spatial coordinates
  loc_use <- data.table::setDF(
    get_spatial_locations(
      gobject = gobject,
      spat_unit = spat_unit,
      output = 'data.table',
      copy_obj = TRUE,
      ... # allow setting of spat_loc_name through additional params
    )
  )
  rownames(loc_use) <- loc_use$cell_ID
  sobj <- Seurat::AddMetaData(sobj, metadata = loc_use)
  # add spatial coordinates as new dim reduct object
  loc_2 <- loc_use[,c('sdimx','sdimy')]
  colnames(loc_2) <- c('spatial_1','spatial_2')
  sobj[['spatial']] <- Seurat::CreateDimReducObject(embeddings = as.matrix(loc_2),
                                            assay = names(sobj@assays)[1],
                                            key = 'spatial_')



  # spatial network
  avail_sn = list_spatial_locations(gobject = gobject, spat_unit = spat_unit)
  sn_all = avail_sn[, name]
  # sn_all <- names(gobject@spatial_network)
  if (nrow(sn_all) > 0){
    for (i in sn_all){
      snt_use <- get_spatialNetwork(gobject = gobject, name = i)
      idx1 <- match(snt_use$from,Seurat::Cells(sobj))
      idx2 <- match(snt_use$to,Seurat::Cells(sobj))
      edge_weight <- snt_use$weight
      nn_mtx <- Matrix::sparseMatrix(i = idx1,j = idx2,x = edge_weight,dims = c(ncol(sobj),ncol(sobj)))
      rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
      nn_name <- paste0('spatial_',i)
      sobj[[nn_name]] <- Seurat::as.Graph(nn_mtx)
      sobj[[nn_name]]@assay.used <- names(sobj@assays)[1]
    }
  }

  return (sobj)
}


#' @title seuratToGiotto_OLD
#' @name seuratToGiotto_OLD
#' @description Converts Seurat object into a Giotto object. Deprecated, see \code{\link{giottoToSeurat}}
#' @param obj_use Seurat object
#' @param ... additional params to pass
#' @return Giotto object
#' @export
seuratToGiotto_OLD <- function(obj_use = NULL,
                               ...){
  requireNamespace('Seurat')
  requireNamespace('Giotto')

  # get general info in basic seurat structures
  obj_assays <- names(obj_use@assays)
  if ('Spatial' %in% obj_assays){
    obj_assays <- c('Spatial',obj_assays[-which(obj_assays == 'Spatial')])
  }

  obj_dimReduc <- names(obj_use@reductions)
  obj_dimReduc_assay <- sapply(obj_dimReduc,function(x)
    obj_use[[x]]@assay.used)

  obj_graph_expr <- names(obj_use@graphs)
  obj_graph_expr_assay <- sapply(obj_graph_expr,function(x)
    obj_use[[x]]@assay.used)

  obj_meta_cells <- obj_use@meta.data
  obj_meta_genes <- lapply(obj_assays,function(x)
    obj_use[[x]]@meta.features)
  names(obj_meta_genes) <- obj_assays

  obj_img <- obj_use@images
  obj_img_names <- names(obj_img)
  loc_use <- lapply(obj_img_names,function(x){
    temp <- obj_img[[x]]@coordinates
    temp <- as.data.frame(temp[,c('col','row')])
    # temp$region <- x
    return (temp)
  })
  loc_use <- Reduce(rbind,loc_use)

  # add assay data: raw, normalized & scaled
  for (i in 1:length(obj_assays)){
    data_raw <- Seurat::GetAssayData(obj_use,slot = 'counts',assay = obj_assays[i])
    data_norm <- Seurat::GetAssayData(obj_use,slot = 'data',assay = obj_assays[i])
    data_scale <- Seurat::GetAssayData(obj_use,slot = 'scale.data',assay = obj_assays[i])

    if (i == 1 & obj_assays[i] == 'Spatial'){
      feat_use <- 'rna'
      test <- createGiottoObject(expression = obj_use[[obj_assays[i]]]@counts,
                                 spatial_locs = loc_use,
                                 expression_feat = 'rna')
      test <- addCellMetadata(test,feat_type = feat_use,new_metadata = obj_meta_cells)
    } else {
      feat_use <- obj_assays[i]
      test@expression[[feat_use]][['raw']] <- data_raw
      test@feat_ID[[feat_use]] = rownames(data_raw)
      test@feat_metadata[[feat_use]] = data.table::data.table(feat_ID = test@feat_ID[[feat_use]])
    }
    if (nrow(data_norm) > 0){
      test@expression[[feat_use]][['normalized']] <- data_norm
    }
    if (nrow(data_scale) > 0){
      test@expression[[feat_use]][['scaled']] <- data_scale
    }

    # gene metadata
    if (length(obj_meta_genes[[i]]) > 0){
      test <- addFeatMetadata(test,feat_type = feat_use,
                              new_metadata = obj_meta_genes[[i]])
    }
  }

  # add dim reduction
  for (i in obj_dimReduc){
    if (!i %in% c('pca','umap','tsne')){
      next
    } else {
      dimReduc_name <- i
      dimReduc_method <- i
      dimReduc_coords <- obj_use[[i]]@cell.embeddings
      dimReduc_misc <- list(obj_use[[i]]@stdev,
                            obj_use[[i]]@feature.loadings,
                            obj_use[[i]]@feature.loadings.projected)
      names(dimReduc_misc) <- c('eigenvalues','loadings','loadings_projected')
      dimObject <- create_dimObject(name = dimReduc_name,reduction_method = dimReduc_method,
                                             coordinates = dimReduc_coords,misc = dimReduc_misc)
      test@dimension_reduction[['cells']][[dimReduc_method]][[dimReduc_name]] <- dimObject
    }
  }

  # add expr nearest neighbors
  for (i in obj_graph_expr){
    mtx_use <- obj_use[[i]]
    ig_use <- igraph::graph_from_adjacency_matrix(mtx_use,weighted = T)
    g_type <- unlist(strsplit(i,split = "_"))[2]
    g_val <- i
    test@nn_network[[g_type]][[g_val]][['igraph']] <- ig_use
  }
  return (test)

}




#' Convert a Seurat object to a Giotto object
#'
#' @param sobject Input Seurat object to convert to Giotto object
#' @param spatial_assay Specify name of the spatial assay slot in Seurat. Default is \code{"Spatial"}.
#' @param dim_reduction Specify which dimensional reduction computations to fetch from
#'  input Seurat object. Default is \code{"c('pca', 'umap')"}.
#' @param subcellular_assay Specify name of the subcellular assay in input object.
#'  Default is \code{"Vizgen"}.
#' @return A Giotto object converted from Seurat object with all computations stored in it.
#' @export
seuratToGiotto = function(sobject,
                          spatial_assay = 'Spatial',
                          dim_reduction = c('pca','umap'),
                          subcellular_assay = 'Vizgen'){
  package_check('Seurat')

  if(is.null(Seurat::GetAssayData(object = sobject, slot = "counts", assay = spatial_assay))) {
    cat('No raw expression values are provided in spatial_assay\n')
    return(sobject)

  } else {

    exp = Seurat::GetAssayData(object = sobject, slot = "counts", assay = spatial_assay)
    if(!is.null(sobject@assays$SCT)){
        normexp = Seurat::GetAssayData(object = sobject, slot = "counts", assay = 'SCT')
      }

      if(!is.null(sobject@assays$Spatial@data)){
        normexp = Seurat::GetAssayData(object = sobject, slot = "data", assay = 'Spatial')
      }

    # Cell Metadata
    cell_metadata = sobject@meta.data

    # Dimension Reduction
    if(sum(sapply(dim_reduction,function(x) length(sobject@reductions[[x]]))) == 0) {
      dim_reduc = NULL
    } else {
      dimReduc_list = lapply(dim_reduction,function(x){
        dim_coord = as.matrix(Seurat::Embeddings(object = sobject,reduction = x))
        dim_load = as.matrix(Seurat::Loadings(object = sobject, reduction = x))
        dim_eig = Seurat::Stdev(sobject, reduction = x)
        if (length(dim_eig) > 0){
          dim_eig = sapply(dim_eig, function(x) x ^ 2)
        }
        colnames(dim_coord) = paste0('Dim.',1:ncol(dim_coord))
        if (length(dim_load) > 0){
          colnames(dim_load) = paste0('Dim.',1:ncol(dim_load))
        }
        dimReduc = list(type = "cells",
                        spat_unit = "cell",
                        name = x,
                        reduction_method = x,
                        coordinates = dim_coord,
                        feat_type = 'rna',
                        misc =list(eigenvalues = dim_eig, loadings = dim_load))
        return (dimReduc)
      })
      # names(dimReduc_list) <- dim_reduction
    }

    # Spatial Locations
    if(length(sobject@assays[[spatial_assay]]) == 0) {
      spat_loc = NULL
    } else {

      spat_coord = Seurat::GetTissueCoordinates(sobject)
      spat_coord = cbind(rownames(spat_coord), data.frame(spat_coord, row.names=NULL))
      colnames(spat_coord) = c("cell_ID", "sdimy", "sdimx")
      spat_loc = spat_coord
    }


    # Subcellular
    name = names(sobject@images)
    if(length(sobject@assays[[subcellular_assay]]) == 1) {

      spat_coord = Seurat::GetTissueCoordinates(sobject)
      colnames(spat_coord) = c("sdimx", "sdimy", "cell_ID")
      exp = exp[  , c(intersect(spat_coord$cell_ID, colnames(exp)))]
      spat_loc = spat_coord
    }
    if (!length(sobject@images) == 0) {
      if ("molecules" %in% methods::slotNames(sobject@images[[name]]) == TRUE) {
        if(!length(sobject@images[[name]][["molecules"]]) == 0) {

          assay = names(sobject@assays)
          featnames = rownames(sobject@assays[[assay]]@meta.features)
          mol_spatlocs = data.table::data.table()

          for (x in featnames) {
            df = (Seurat::FetchData(sobject[[name]][["molecules"]], vars = x))
            mol_spatlocs = rbind(mol_spatlocs, df)
          }
          gpoints = createGiottoPoints(mol_spatlocs, feat_type = "rna")

        }
      }
    }
  }

  gobject = createGiottoObject(exp,
                               spatial_locs = spat_loc,
                               dimension_reduction = dimReduc_list)
  if (exists('normexp') == TRUE) {
    exprObj = create_expr_obj(name = 'normalized',
                              exprMat = normexp,
                              spat_unit = 'cell',
                              feat_type = 'rna',
                              provenance = 'cell')
    gobject = set_expression_values(gobject = gobject, values = exprObj, set_defaults = FALSE)
    # gobject@expression$cell$rna$normalized = normexp
  }
  gobject = addCellMetadata(gobject = gobject, new_metadata = cell_metadata)


  if (exists('gpoints') == TRUE) {
    gobject = addGiottoPoints(gobject = gobject,
                              gpoints = list(gpoints))
  }

  return (gobject)
}


## SpatialExperiment object ####

#' Utility function to convert a Giotto object to a SpatialExperiment object.
#'
#' @param giottoObj Input Giotto object to convert to a SpatialExperiment object.
#' @param verbose A boolean value specifying if progress messages should be displayed
#' or not. Default \code{TRUE}.
#'
#' @return A SpatialExperiment object that contains data from the input Giotto object.
#' @examples
#' mini_gobject <- GiottoData::loadGiottoMini('vizgen')
#' giottoToSpatialExperiment(mini_gobject)
#'
#' @export
giottoToSpatialExperiment <- function(giottoObj, verbose = TRUE){

  # Load required packages
  # package_check(pkg_name = "SummarizedExperiment", repository = 'Bioc') # SP should load this?
  # package_check(pkg_name = "SingleCellExperiment", repository = 'Bioc') # SP should load this?
  package_check(pkg_name = "SpatialExperiment", repository = 'Bioc')
  package_check(pkg_name = "S4Vectors", repository = 'Bioc')

  # SpatialExperiment objects list (one object for each spatial unit)
  speList <- list()

  # Expression Matrices
  giottoExpr <- list_expression(giottoObj)

  # Iterate over spatial units
  spatialUnits <- unique(giottoExpr$spat_unit) # a function to get spat units?
  for(su in seq(spatialUnits)){
    if(verbose) message("Processing spatial unit: '", spatialUnits[su], "'")
    # Check if expression matrices exist in input object
    if(!is.null(giottoExpr)){
      if(verbose) message("Copying expression matrix: '", giottoExpr[1]$name, "' for spatial unit: '", spatialUnits[su], "'")
      exprMat <- get_expression_values(
        gobject = giottoObj,
        spat_unit = spatialUnits[su],
        feat_type = giottoExpr[1]$feat_type,
        values = giottoExpr[1]$name,
        output = "matrix")
      names(rownames(exprMat)) <- NULL
      names(colnames(exprMat)) <- NULL
      exprMat <- list(exprMat)
      names(exprMat)[1] <- giottoExpr[1]$name
      # Creating SPE object with first expression matrix
      spe <- SpatialExperiment::SpatialExperiment(assays = exprMat)
      assayNames(spe) <- paste0(giottoExpr[1]$name, "_",
                                giottoExpr[1]$feat_type, "_",
                                spatialUnits[su])
      giottoExpr <- giottoExpr[-1, ]
    } else{
      stop("The input Giotto object must contain atleast one expression matrix.")
    }

    # Copying remaining expression matrices if they exist
    if(nrow(giottoExpr[spat_unit == spatialUnits[su]]) > 0){
      for(i in seq(nrow(giottoExpr))){
        if(verbose) message("Copying expression matrix: '", giottoExpr[i]$name, "' for spatial unit: '", spatialUnits[su], "'")
        # SPE does not have specific slots for different units, instead joining multiple unit names to identify them
        SummarizedExperiment::assay(
          spe,
          paste0(giottoExpr[i]$name, "_",
                 giottoExpr[i]$feat_type, "_",
                 spatialUnits[su]),
          withDimnames = FALSE) <- get_expression_values(
            gobject = giottoObj,
            spat_unit = spatialUnits[su],
            feat_type = giottoExpr[i]$feat_type,
            values = giottoExpr[i]$name,
            output = "matrix")
      }
    }

    # Cell Metadata to ColData
    pData <- pDataDT(gobject = giottoObj, spat_unit = spatialUnits[su])
    if(nrow(pData) > 0){
      if(verbose) message("Copying phenotype data for spatial unit: '", spatialUnits[su], "'")
      SummarizedExperiment::colData(spe) <- S4Vectors::DataFrame(pData, row.names = pData$cell_ID)
    } else{
      if(verbose)  message("No phenotype data found in input Giotto object")
    }

    # Feature Metadata to RowData
    fData <- fDataDT(gobject = giottoObj, spat_unit = spatialUnits[su])
    if(nrow(fData) > 0){
      if(verbose) message("Copying feature metadata for spatial unit: '", spatialUnits[su], "'")
      SummarizedExperiment::rowData(spe) <- fData
    } else{
      if(verbose) message("No feature metadata found in input Giotto object")
    }

    # Spatial Locations to Spatial Coordinates
    spatialLocs <- get_spatial_locations(gobject = giottoObj,
                                         spat_unit = spatialUnits[su],
                                         output = "data.table")
    if(!is.null(spatialLocs)){
      if(verbose) message("Copying spatial locations for spatial unit: '", spatialUnits[su], "'")
      SpatialExperiment::spatialCoords(spe) <- data.matrix(spatialLocs[, 1:2])
    } else{
      if(verbose) message("No spatial locations found in the input Giotto object")
    }

    # DimReductions
    giottoReductions <- list_dim_reductions(gobject = giottoObj, spat_unit = spatialUnits[su])
    if(!is.null(giottoReductions)){
      if(verbose) message("Copying reduced dimensions for spatial unit: '", spatialUnits[su], "'")
      for(i in seq(nrow(giottoReductions))){
        SingleCellExperiment::reducedDim(spe, giottoReductions[i]$name) <- get_dimReduction(
          gobject = giottoObj,
          reduction = "cells",
          spat_unit = spatialUnits[su],
          feat_type = giottoReductions[i]$feat_type,
          reduction_method = giottoReductions[i]$dim_type,
          name = giottoReductions[i]$name,
          output = "data.table")
      }
    } else{
      if(verbose) message("No reduced dimensions found in the input Giotto object")
    }


    # NN Graph
    giottoNearestNetworks <- list_nearest_networks(gobject = giottoObj, spat_unit = spatialUnits[su])
    if(!is.null(giottoNearestNetworks)){
      if(verbose) message("Copying nearest networks for spatial unit: '", spatialUnits[su], "'")
      for(i in seq(nrow(giottoNearestNetworks))){
        nn_network <- get_NearestNetwork(
          gobject = giottoObj,
          spat_unit = spatialUnits[su],
          nn_network_to_use = giottoNearestNetworks[i]$type,
          network_name = giottoNearestNetworks[i]$name,
          output = "data.table")

        # SPE stores in colpairs, with col indices instead of colnames
        cell1 <- match(nn_network$from, pData$cell_ID)
        cell2 <- match(nn_network$to, pData$cell_ID)

        SingleCellExperiment::colPair(spe, giottoNearestNetworks[i]$name) <- S4Vectors::SelfHits(
          cell1, cell2, nnode=ncol(spe),
          nn_network[, -1:-2] #removing from and to
          )
      }
    } else{
      if(verbose) message("No nearest networks found in the input Giotto object")
    }

    # Spatial Networks
    giottoSpatialNetworks <- list_spatial_networks(gobject = giottoObj, spat_unit = spatialUnits[su])
    if(!is.null(giottoSpatialNetworks)){
      if(verbose) message("Copying spatial networks for spatial unit: '", spatialUnits[su], "'")
      for(i in seq(nrow(giottoSpatialNetworks))){
        sp_network <- get_spatialNetwork(
          gobject = giottoObj,
          spat_unit = spatialUnits[su],
          name = giottoSpatialNetworks[i]$name,
          output = "networkDT")

        # spe stores in colpairs, with col indices instead of colnames
        cell1 <- match(sp_network$from, pData$cell_ID)
        cell2 <- match(sp_network$to, pData$cell_ID)

        SingleCellExperiment::colPair(spe, giottoSpatialNetworks[i]$name) <- S4Vectors::SelfHits(
          cell1, cell2, nnode=ncol(spe),
          sp_network[, -1:-2] #removing from and to
          )
      }
    } else{
      if(verbose) message("No spatial networks found in the input Giotto object")
    }

    # SpatialImages
    giottoImages <- list_images(gobject = giottoObj)
    if(!is.null(giottoImages)){
      for(i in seq(nrow(giottoImages))){

        img <- get_giottoImage(
          gobject = giottoObj,
          image_type = giottoImages[i]$img_type,
          name = giottoImages[i]$name)

        if(!is.null(img@file_path)){

          if(verbose) message("Copying spatial image: ", img@name)

          spe <- SpatialExperiment::addImg(spe,
                                           sample_id = spe$sample_id[i],
                                           image_id = img@name,
                                           imageSource = img@file_path,
                                           scaleFactor = mean(img@scale_factor),
                                           load = TRUE)
        }
        else{
          if(verbose) message("\t - Skipping image with NULL file path: ", img@name)
        }
        S4Vectors::metadata(spe)[[img@name]] <- img
      }
    } else{
      if(verbose) message("No spatial images found in the input Giotto object")
    }

    if(verbose) message("")

    # Add spe for current spatial unit to speList
    speList[[su]] <- spe
  }

  # return list of spe objects
  return(speList)
}


#' Utility function to convert a SpatialExperiment object to a Giotto object
#'
#' @param spe Input SpatialExperiment object to convert to a Giotto object.
#' @param nn_network Specify the name of the nearest neighbour network(s)
#' in the input SpatialExperiment object. Default \code{NULL} will use
#' all existing networks.
#' @param sp_network Specify the name of the spatial network(s) in the input
#' SpatialExperiment object. Default \code{NULL} will use all existing
#' networks.
#' @param verbose A boolean value specifying if progress messages should
#' be displayed or not. Default \code{TRUE}.
#' @return Giotto object
#' @examples
#' library(SpatialExperiment)
#' example(read10xVisium, echo = FALSE)
#' spatialExperimentToGiotto(spe)
#' @export
spatialExperimentToGiotto <- function(spe,
                                      nn_network = NULL,
                                      sp_network = NULL,
                                      verbose = TRUE){

  # Create Giotto object with first matrix
  exprMats <- assays(spe)
  exprMatsNames <- assayNames(spe)
  firstMatrix <- exprMats[[1]]

  #check unique colnames
  if(length(unique(colnames(firstMatrix))) != length(colnames(firstMatrix))){
    colnames(firstMatrix) <- make.names(colnames(firstMatrix), unique = TRUE)
  }

  if(verbose) message("Creating Giotto object with ", exprMatsNames[1], " matrix")
  suppressWarnings(suppressMessages(giottoObj <- createGiottoObject(expression = firstMatrix)))
  exprMats[[1]] <- NULL
  exprMatsNames <- exprMatsNames[-1]

  # Copying remaining matrices
  if(length(exprMats) > 0){
    for(i in seq(exprMats)){
      if(verbose) message("Copying expression matrix: ", exprMatsNames[i])
      giottoObj <- set_expression_values(gobject = giottoObj, name = exprMatsNames[i], values = exprMats[[i]])
    }
  }

  # Phenotype Data
  pData <- colData(spe)
  if(nrow(pData) > 0){
    if(verbose) message("Copying phenotype data")
    giottoObj <- addCellMetadata(gobject = giottoObj, new_metadata = as.data.table(pData))
  }

  # Feature Metadata
  fData <- rowData(spe)
  if(nrow(fData) > 0){
    if(verbose) message("Copying feature metadata")
    giottoObj <- addFeatMetadata(gobject = giottoObj, new_metadata = as.data.table(fData))
  }

  # Reduced Dimensions
  redDims <- reducedDims(spe)
  redDimsNames <- reducedDimNames(spe)
  if(length(redDims) > 0){
    for(i in seq(length(redDims))){
      if(verbose) message("Copying reduced dimensions")
      dimRedObj <- create_dimObject(name = redDimsNames[i],
                                             coordinates = redDims[[i]],
                                             reduction_method = redDimsNames[i])
      giottoObj <- set_dimReduction(gobject = giottoObj, dimObject = dimRedObj)
    }
  }

  # Spatial Locations
  spatialLocs <- spatialCoords(spe)
  if(ncol(spatialLocs) > 0){
    if(verbose) message("Copying spatial locations")
    spatialLocsDT <- data.table(sdimx = spatialLocs[, 1], sdimy = spatialLocs[, 2], cell_ID = rownames(spatialLocs))
    giottoObj <- set_spatial_locations(gobject = giottoObj, spatlocs = cbind(spatialLocsDT, cell_ID = colnames(spe)))
  }

  # Spatial Images
  spatialImages <- imgData(spe)
  if(nrow(spatialImages) > 0){
    for(i in seq(nrow(spatialImages))){
      if(verbose) message("Copying spatial images")
      spImg <- getImg(spe,
                      spatialImages[i, "sample_id"],
                      spatialImages[i, "image_id"])
      mObject <- magick::image_read(as.raster(spImg))
      giottoImage <- createGiottoImage(gobject = giottoObj,
                                                mg_object = mObject,
                                                scale_factor = spatialImages[i, "scaleFactor"])
      giottoObj <- addGiottoImage(gobject = giottoObj,
                                  images = list(giottoImage))

    }
  }

  # Networks
  networks <- colPairs(spe)
  # Spatial Networks
  if(!is.null(sp_network)){
    if(sp_network %in% names(networks)){
      for(i in seq(sp_network)){
        if(verbose) message("Copying spatial networks")
        spatNetObj <- create_spat_net_obj(
          networkDT = as.data.table(networks[[sp_network[i]]]))
        giottoObj <- set_spatialNetwork(gobject = giottoObj,
                                                 spatial_network = spatNetObj,
                                                 name = sp_network[i])
        networks[[sp_network[i]]] <- NULL
      }
    }
  }

  # Nearest Neighbour Networks
  if(!is.null(nn_network)){
    if(nn_network %in% names(networks)){
      for(i in seq(nn_network)){
        if(verbose) message("Copying nearest neighbour networks")
        giottoObj <- set_NearestNetwork(gobject = giottoObj,
                                                 nn_network = networks[[nn_network[i]]],
                                                 network_name = nn_network[i])
        networks[[nn_network[i]]] <- NULL
      }
    }
  }

  # if not specified, storing remaining as NN
  if(length(networks) > 0){
    for(i in seq(networks)){
      if(verbose) message("Copying additional networks")
      giottoObj <- set_NearestNetwork(gobject = giottoObj,
                                               nn_network = networks[[i]],
                                               network_name = names(networks)[i])
    }
  }

  return(giottoObj)
}



