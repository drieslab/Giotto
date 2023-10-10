
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
#' There are six possible choices for bin_size: bin1, bin10, bin20, bin50, bin100, bin200.
#'
#' See SAW pipeline for additional information about the gef file.
#' @export

gefToGiotto = function(gef_file, bin_size = 'bin100', verbose = FALSE){

  # data.table vars
  genes = y = sdimx = sdimy = cell_ID = count = NULL

  # package check
  package_check(pkg_name = 'rhdf5', repository = 'Bioc')
  if(!file.exists(gef_file)) stop('File path to .gef file does not exist')

  # check if proper bin_size is selected. These are determined in SAW pipeline.
  bin_size_options = c('bin1', 'bin10', 'bin20', 'bin50', 'bin100', 'bin200')
  if(!(bin_size %in% bin_size_options)) stop('Please select valid bin size,see details for choices.')

  # step 1: read expression and gene data from gef file
  if(isTRUE(verbose)) wrap_msg('reading in .gef file')
  geneExpData = rhdf5::h5read(file = gef_file, name = 'geneExp')
  exprDT = data.table::as.data.table(geneExpData[[bin_size]][['expression']])
  geneDT = data.table::as.data.table(geneExpData[[bin_size]][['gene']])

  # step 2: combine gene information from the geneDT to the exprDT
  exprDT[, genes := rep(x = geneDT$gene, geneDT$count)]

  # step 3: bin coordinates according to selected bin_size
  #TODO: update bin_shift for other shapes, not just rect_vertices
  bin_size_int = as.integer(gsub("[^0-9.-]", "", bin_size))
  bin_shift = ceiling(bin_size_int / 2) # ceiling catches bin_1
  bincoord = unique(exprDT[,.(x,y)])
  if(isTRUE(verbose)) wrap_msg('shifting and binning coordinates')
  data.table::setorder(bincoord, x, y)
  data.table::setnames(bincoord, old = c('x', 'y'), new = c('sdimx', 'sdimy'))
  bincoord[, c('sdimx', 'sdimy') := list(sdimx+bin_shift, sdimy+bin_shift)]
  bincoord[, cell_ID := paste0('bin', 1:.N)]
  tx_data = exprDT[,.(genes, x, y, count)]
  tx_data[, c('x', 'y') := list(x+bin_shift, y+bin_shift)]

  # step 4: create rectangular polygons (grid) starting from the bin centroids
  if(isTRUE(verbose)) wrap_msg('creating polygon stamp')
  x = polyStamp(stamp_dt = rectVertices(dims = c(x = (bin_size_int - 1),
                                                 y = (bin_size_int - 1))),
                spatlocs = bincoord[,.(cell_ID, sdimx, sdimy)])
  pg = createGiottoPolygonsFromDfr(x)

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


#' @title Check Scanpy Installation
#' @name check_py_for_scanpy
#' @description checks current python environment for scanpy 1.9.0
#' @keywords internal
check_py_for_scanpy = function(){
  # test if scanpy is found
  module_test = reticulate::py_module_available('scanpy')
  py_path = reticulate::py_config()$python
  genv_in_use = grepl(pattern = "giotto_env", x = py_path)

  if(module_test == FALSE && !genv_in_use) {
    stop(wrap_txt("scanpy python module is not installed:
            install in the environment or python path with:

            'pip install scanpy==1.9.0'

            Alternatively, install in the active python
            environment with reticulate:

            reticulate::py_install(packages = 'scanpy==1.9.0',
                                   pip = TRUE)
            \n
            ",errWidth = TRUE))
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

}


#' @title Convert anndata to Giotto
#' @name anndataToGiotto
#' @description Converts a spatial anndata (e.g. scanpy) .h5ad file into a Giotto object
#' @param anndata_path path to the .h5ad file
#' @param n_key_added equivalent of "key_added" argument from scanpy.pp.neighbors().
#'                    If multiple spatial networks are in the anndata object, a list of key_added
#'                            terms may be provided.
#'                            If converting an anndata object from giottoToAnnData, a .txt file may be
#'                            provided, which was generated in that function, 
#'                                i.e. {spat_unit}_{feat_type}_nn_network_keys_added.txt
#'                    Cannot be "spatial". This becomes the name of the nearest network in the gobject.
#' @param spatial_n_key_added equivalent of "key_added" argument from squidpy.gr.spatial_neighbors. 
#'                            If multiple spatial networks are in the anndata object, a list of key_added
#'                            terms may be provided.
#'                            If converting an anndata object from giottoToAnnData, a .txt file may be
#'                            provided, which was generated in that function, 
#'                                i.e. {spat_unit}_{feat_type}_spatial_network_keys_added.txt
#'                            Cannot be the same as n_key_added.
#' @param deluanay_spat_net binary parameter for spatial network. If TRUE, the spatial network is a deluanay network.
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
                           spatial_n_key_added = NULL,
                           deluanay_spat_net = TRUE,
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
  if (!is.null(n_key_added) && !is.null(spatial_n_key_added)){
    for (n in n_key_added){
      for (s in spatial_n_key_added) {
        if (n == s) stop("Arguments n_key_added and spatial_n_key_added may not take the same value.")
      }
    }
  }

  # Required step to properly initialize reticualte
  instrs = createGiottoInstructions(python_path = python_path)

  check_py_for_scanpy()

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
  cmeta = extract_cell_metadata(adata)
  cmeta = as.data.table(cmeta)
  if ('leiden' %in% names(cmeta)) {
    cmeta$leiden = as.numeric(cmeta$leiden)
  }

  fm = extract_feat_metadata(adata)
  fm = as.data.table(fm)
  # Metadata ready

  ### Create Minimal giottoObject
  gobject <- createGiottoObject(expression = X,
                                spatial_locs = sp,
                                instructions = instrs)

  ### Add metadata
  cmeta = readCellMetadata(cmeta)
  gobject = setCellMetadata(gobject,
                            x = cmeta)
  fm = readFeatMetadata(fm)
  gobject = setFeatureMetadata(gobject,
                               x = fm)

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
    dobj = create_dim_obj(name = 'pca.ad',
                          spat_unit = spat_unit,
                          feat_type = feat_type,
                          provenance = NULL,
                          reduction = 'cells',
                          reduction_method = 'pca',
                          coordinates = pca,
                          misc = list(eigenvalues = evs,
                                      loadings = loads),
                          my_rownames = colnames(X))

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject = set_dimReduction(gobject = gobject, dimObject = dobj)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  }

  ### Set up UMAP
  u = extract_umap(adata)
  if (!is.null(u)) {
    # Add UMAP to giottoObject
    dobj = create_dim_obj(name = 'umap.ad',
                          spat_unit = spat_unit,
                          feat_type = feat_type,
                          provenance = NULL,
                          reduction = 'cells',
                          reduction_method = 'umap',
                          coordinates = u,
                          misc = NULL,
                          my_rownames = colnames(X))

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject = set_dimReduction(gobject = gobject, dimObject = dobj)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  }
  ### Set up TSNE
  t = extract_tsne(adata)
  if (!is.null(t)) {
    # Add TSNE to giottoObject
    dobj = create_dim_obj(name = 'tsne.ad',
                          spat_unit = spat_unit,
                          feat_type = feat_type,
                          provenance = NULL,
                          reduction = 'cells',
                          reduction_method = 'tsne',
                          coordinates = t,
                          misc = NULL,
                          my_rownames = colnames(X))

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject = set_dimReduction(gobject = gobject, dimObject = dobj)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  }

  ### NN Network

  # Need to create nnNetObj or igraph object to use with setter for NN

  weights_ad = NULL
  num_NN_nets = length(n_key_added)

  if (is.null(n_key_added) && !is.null(extract_NN_connectivities(adata, key_added = n_key_added))) {
    num_NN_nets = 1
  }

  for (i in num_NN_nets){
    if (inherits(n_key_added, "list")){ 
      n_key_added_it = n_key_added[[i]]
    } else {
      n_key_added_it = n_key_added
    }
    
    weights_ad = extract_NN_connectivities(adata, key_added = n_key_added_it)
    #adw = methods::as(weights_ad, "TsparseMatrix")
    if (!is.null(weights_ad)) {
      distances_ad = extract_NN_distances(adata, key_added = n_key_added_it)

      nn_dt = align_network_data(distances = weights_ad, weights = distances_ad)

      #pre-allocate DT variables
      from = to = weight = distance = from_cell_ID = to_cell_ID = uniq_ID = NULL
      nn_dt = data.table::data.table(nn_dt)

      nn_dt[, from_cell_ID := cID[from]]
      nn_dt[, to_cell_ID := cID[to]]
      nn_dt[, uniq_ID := paste0(from,to)]
      nn_dt[order(uniq_ID)]
      nn_dt[,uniq_ID := NULL]
      vert = unique(x = c(nn_dt$from_cell_ID, nn_dt$to_cell_ID))
      nn_network_igraph = igraph::graph_from_data_frame(nn_dt[,.(from_cell_ID, to_cell_ID, weight, distance)], directed = TRUE, vertices = vert)

      nn_info = extract_NN_info(adata = adata, key_added = n_key_added_it)

      net_type = "kNN" # anndata default
      if(("sNN" %in% n_key_added_it) & !is.null(n_key_added_it)){
        net_type = "sNN"
        net_name = paste0(n_key_added_it, ".", nn_info["method"])
      } else if (!("sNN" %in% n_key_added_it) & !is.null(n_key_added_it)) {
        net_name = paste0(n_key_added_it, ".", nn_info["method"])
      } else {
        net_name = paste0(net_type, ".", nn_info["method"])
      }

      netObj = createNearestNetObj(name = net_name,
                                  network = nn_network_igraph,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type)

      gobject = set_NearestNetwork(gobject = gobject,
                                  nn_network = netObj,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  nn_network_to_use = net_type,
                                  network_name = net_name,
                                  set_defaults = FALSE)
    }
  }
  
  ## Spatial Network
  s_weights_ad = NULL
  num_SN_nets = length(spatial_n_key_added) 
  
  # Check for the case where NULL is provided, since the
  # anndata object takes the default value for SN

  if (is.null(spatial_n_key_added) && !is.null(extract_SN_connectivities(adata, key_added = spatial_n_key_added))) {
    num_SN_nets = 1
  }

  for (i in 1:num_SN_nets){

    if (inherits(spatial_n_key_added, "list")){ 
      spatial_n_key_added_it = spatial_n_key_added[[i]]
    } else {
      spatial_n_key_added_it = spatial_n_key_added
    }

    s_weights_ad = extract_SN_connectivities(adata, key_added = spatial_n_key_added_it)
    if (!is.null(s_weights_ad)){
      s_distances_ad = extract_SN_distances(adata, key_added = spatial_n_key_added_it)
      ij_matrix = methods::as(s_distances_ad, "TsparseMatrix")
      from_idx = ij_matrix@i + 1 #zero index!!!
      to_idx = ij_matrix@j + 1 #zero index!!!
      
      #pre-allocate DT variables
      from = to = weight = distance = from_cell_ID = to_cell_ID = uniq_ID = NULL
      sn_dt = data.table::data.table(from = from_idx,
                                    to = to_idx,
                                    weight = s_weights_ad@x,
                                    distance = s_distances_ad@x)
      
      sn_dt[, from_cell_ID := cID[from]]
      sn_dt[, to_cell_ID := cID[to]]

      sdimx = "sdimx"
      sdimy = "sdimy"
      xbegin_name = paste0(sdimx,'_begin')
      ybegin_name = paste0(sdimy,'_begin')
      xend_name = paste0(sdimx,'_end')
      yend_name = paste0(sdimy,'_end')

      network_DT = data.table::data.table(from = sn_dt$from_cell_ID,
                                                  to = sn_dt$to_cell_ID,
                                                  xbegin_name = sp[sn_dt$from, sdimx],
                                                  ybegin_name = sp[sn_dt$from, sdimy],
                                                  xend_name = sp[sn_dt$to, sdimx],
                                                  yend_name = sp[sn_dt$to, sdimy],
                                                  weight = s_weights_ad@x,
                                                  distance = s_distances_ad@x)
      data.table::setnames(network_DT,
                        old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                        new = c(xbegin_name, ybegin_name, xend_name, yend_name))
      data.table::setorder(network_DT, from, to)

      dist_mean = get_distance(network_DT, method = "mean")
      dist_median = get_distance(network_DT, method = "median")
      cellShapeObj = list("meanCellDistance" = dist_mean,
                        "medianCellDistance" = dist_median)

      #TODO filter network? 
      #TODO 3D handling?
      if (deluanay_spat_net){
        spatObj = create_spat_net_obj(name = "Spat_Net_from_AnnData", 
                                    method = "delaunay",
                                    networkDT=network_DT,
                                    cellShapeObj = cellShapeObj)
      } else {
        spatObj = create_spat_net_obj(name = "Spat_Net_from_AnnData", 
                                    method = "non-delaunay",
                                    networkDT=network_DT,
                                    cellShapeObj = cellShapeObj)
      }
      
      gobject = set_spatialNetwork(gobject = gobject,
                                  spatial_network = spatObj,
                                  name = "Spat_Net_from_AnnData")

    }
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
      layExprObj <- createExprObj(lay, name = l_n)
      gobject = set_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      name = l_n,
                                      values = layExprObj)
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
#' This function will create .txt files which will record any `key_added`
#' parameters for networks. They are named after the corresponding spatial unit
#' and feature type pair. 
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

          if ("dgeMatrix" %in% class(x)) x = methods::as(x,'array')

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

      cmeta = get_cell_metadata(gobject = gobject,
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
                                                 cell_meta = cmeta,
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

        fname_nn = paste0(su, "_", ft, "_nn_network_keys_added.txt")
        network_name = network_name[!grepl("kNN.", network_name)]
        append_n = FALSE
        if(length(network_name) != 0 ) {
          if (nn_net_tu == "kNN") append_n = TRUE
          write(network_name, fname_nn, append = append_n)
        }

      }
      adata_pos = adata_pos + 1
    }
  }

  # Reset indexing variable
  adata_pos = 1

  try_get_SN = function(gobject,
                        spat_unit,
                        name,
                        output,
                        set_defaults,
                        verbose) {
    tryCatch(
      {
        spatial_net = get_spatialNetwork(gobject = gobject,
                                         spat_unit = spat_unit,
                                         name = name,
                                         output = output,
                                         set_defaults = set_defaults,
                                         verbose = verbose)
        return(spatial_net)
      },
      error = function(e) {
        return(NULL)
      }
    )
  }


  for (su in spat_unit) {
    for (ft in names(gobject@expression[[su]])) { 
      # Spatial networks do not have a feature type slot. 
      # Iterate through anyways to properly assign to anndata objects
      network_name = list_spatial_networks_names(gobject = gobject,
                                                  spat_unit = su)
      if (is.null(network_name)) {
        next
      }
      for (sn_name in network_name) {
        gob_SN = try_get_SN(gobject = gobject,
                            spat_unit = su,
                            name = sn_name,
                            output = "networkDT",
                            set_defaults = FALSE)

        pidx = grep("spatial_network", names(gobject@parameters))
        for (p in pidx) {
          if (gobject@parameters[[p]]["name of spatial network"] == sn_name) {
            current_param = gobject@parameters[[p]]
            kval = current_param["k neighbours"]
            maxdist = current_param["maximum distance threshold"]
            dimused = current_param["dimensions used"]
          }
        }
        
        adata_list[[adata_pos]] = set_adg_sn(adata = adata_list[[adata_pos]],
                                              df_SN = gob_SN,
                                              net_name = sn_name,
                                              n_neighbors = kval,
                                              max_distance = maxdist,
                                              dim_used = dimused)
        
        }

      fname_sn = paste0(su,"_",ft, "_spatial_network_keys_added.txt")
      if(length(network_name) != 0) write(network_name, fname_sn)
      }

    adata_pos = adata_pos + 1
  }
  
  # Reset indexing variable
  adata_pos = 1

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
      dimObject <- create_dim_obj(name = dimReduc_name,
                                  reduction_method = dimReduc_method,
                                  coordinates = dimReduc_coords,
                                  misc = dimReduc_misc)
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


## SpatialExperiment object ####

#' Utility function to convert a Giotto object to a SpatialExperiment object.
#'
#' @param giottoObj Input Giotto object to convert to a SpatialExperiment object.
#' @param verbose A boolean value specifying if progress messages should be displayed
#' or not. Default \code{TRUE}.
#'
#' @return A SpatialExperiment object that contains data from the input Giotto object.
#' @examples
#' \dontrun{
#' mini_gobject <- GiottoData::loadGiottoMini('vizgen')
#' giottoToSpatialExperiment(mini_gobject)
#' }
#' @export
giottoToSpatialExperiment <- function(giottoObj, verbose = TRUE){

  spat_unit = NULL

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
      SummarizedExperiment::assayNames(spe) <- paste0(giottoExpr[1]$name, "_",
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

          tryCatch(
            expr = {
              spe <- SpatialExperiment::addImg(spe,
                                               sample_id = spe$sample_id[i],
                                               image_id = img@name,
                                               imageSource = img@file_path,
                                               scaleFactor = mean(img@scale_factor),
                                               load = TRUE)
            },
            error = function(e){
              message("Error copying spatial image: ", img@name, ". Please check if the image path is correct and the image exists at that path.")
            }
          )  
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
#' @param python_path Specify the path to python. 
#' @param nn_network Specify the name of the nearest neighbour network(s)
#' in the input SpatialExperiment object. Default \code{NULL} will use
#' all existing networks.
#' @param sp_network Specify the name of the spatial network(s) in the input
#' SpatialExperiment object. Default \code{NULL} will use all existing
#' networks. This can be a vector of multiple network names.
#' @param verbose A boolean value specifying if progress messages should
#' be displayed or not. Default \code{TRUE}.
#' @import data.table
#' @return Giotto object
#' @examples
#' \dontrun{
#' library(SpatialExperiment)
#' example(read10xVisium, echo = FALSE)
#' spatialExperimentToGiotto(spe)
#' }
#' @export
spatialExperimentToGiotto <- function(spe, 
                                      python_path,
                                      nn_network = NULL,
                                      sp_network = NULL,
                                      verbose = TRUE){

  # Create giotto instructions and set python path
  instrs <- createGiottoInstructions(python_path = python_path)
  
  # Create Giotto object with first matrix
  exprMats <- SummarizedExperiment::assays(spe)
  exprMatsNames <- SummarizedExperiment::assayNames(spe)
  firstMatrix <- exprMats[[1]]

  #check unique colnames
  if(length(unique(colnames(firstMatrix))) != length(colnames(firstMatrix))){
    colnames(firstMatrix) <- make.names(colnames(firstMatrix), unique = TRUE)
  }

  if(verbose) message("Creating Giotto object with ", exprMatsNames[1], " matrix")
  suppressWarnings(suppressMessages(giottoObj <- createGiottoObject(expression = firstMatrix, 
                                                                    instructions = instrs)))
  exprMats[[1]] <- NULL
  exprMatsNames <- exprMatsNames[-1]

  # Copying remaining matrices
  if(length(exprMats) > 0){
    for(i in seq(exprMats)){
      if(verbose) message("Copying expression matrix: ", exprMatsNames[i])
      exprObj <- create_expr_obj(name = exprMatsNames[i], exprMat = exprMats[[i]])
      giottoObj <- set_expression_values(gobject = giottoObj, values = exprObj)
    }
  }

  # Phenotype Data
  pData <- SummarizedExperiment::colData(spe)
  
  # To handle cell_ID duplication issues
  if ("cell_ID" %in% colnames(pData)) {pData$`cell_ID` <- NULL} 
  if(nrow(pData) > 0){
    if(verbose) message("Copying phenotype data")
    giottoObj <- addCellMetadata(gobject = giottoObj, new_metadata = as.data.table(pData))
  }

  # Feature Metadata
  fData <- SummarizedExperiment::rowData(spe)
  if(nrow(fData) > 0){
    if(verbose) message("Copying feature metadata")
    giottoObj <- addFeatMetadata(gobject = giottoObj, new_metadata = as.data.table(fData))
  }

  # Reduced Dimensions
  redDims <- SingleCellExperiment::reducedDims(spe)
  redDimsNames <- SingleCellExperiment::reducedDimNames(spe)
  if(length(redDims) > 0){
    for(i in seq(length(redDims))){
      if(verbose) message("Copying reduced dimensions")
      dimRedObj <- create_dim_obj(name = redDimsNames[i],
                                  coordinates = redDims[[i]],
                                  reduction_method = redDimsNames[i])
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      giottoObj <- set_dimReduction(gobject = giottoObj, dimObject = dimRedObj)
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    }
  }

  # Spatial Locations
  spatialLocs <- SpatialExperiment::spatialCoords(spe)
  if(ncol(spatialLocs) > 0){
    if(verbose) message("Copying spatial locations")
    spatialLocsDT <- data.table(sdimx = spatialLocs[, 1], sdimy = spatialLocs[, 2], cell_ID = colnames(spe))
    spatLocsObj <- Giotto:::create_spat_locs_obj(name = "spatLocs", coordinates = spatialLocsDT)
    giottoObj <- set_spatial_locations(gobject = giottoObj, spatlocs = spatLocsObj)
  }

  # Spatial Images
  spatialImages <- SpatialExperiment::imgData(spe)
  if(nrow(spatialImages) > 0){
    for(i in seq(nrow(spatialImages))){
      if(verbose) message("Copying spatial images")
      spImg <- SpatialExperiment::getImg(spe,
                                         spatialImages[i, "sample_id"],
                                         spatialImages[i, "image_id"])
      mObject <- magick::image_read(grDevices::as.raster(spImg))
      giottoImage <- createGiottoImage(gobject = giottoObj,
                                       mg_object = mObject,
                                       scale_factor = spatialImages[i, "scaleFactor"])
      giottoObj <- addGiottoImage(gobject = giottoObj,
                                  images = list(giottoImage))

    }
  }

  # Networks
  networks <- SingleCellExperiment::colPairs(spe)
  # Spatial Networks
  if(!is.null(sp_network)){
    if(all(sp_network %in% names(networks))){
      for(i in seq(sp_network)){
        if(verbose) message("Copying spatial networks")
        networkDT = as.data.table(networks[[sp_network[i]]])
        networkDT$to <-colnames(spe)[networkDT$to]
        networkDT$from <- colnames(spe)[networkDT$from]
        spatNetObj <- create_spat_net_obj(
          networkDT = networkDT)
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
        nnNetObj <- Giotto:::create_nn_net_obj(name = nn_network[i], 
                                               igraph = networks[[nn_network[i]]])
        giottoObj <- set_NearestNetwork(gobject = giottoObj,
                                        nn_network = nnNetObj)
        networks[[nn_network[i]]] <- NULL
      }
    }
  }

  # if not specified, storing remaining as NN
  if(length(networks) > 0){
    for(i in seq(networks)){
      if(verbose) message("Copying additional networks")
      nnNetObj <- Giotto:::create_nn_net_obj(name = names(networks)[i], 
                                             igraph = networks[[i]])
      giottoObj <- set_NearestNetwork(gobject = giottoObj,
                                      nn_network = nnNetObj)
    }
  }

  return(giottoObj)
}

#' Convert a master Giotto object to suite
#'
#' @param gobject A Giotto object created with master version
#' @param expression_feat available features (e.g. rna, protein, ...)
#'
#' @return A Giotto object compatible with suite version
#' @export
#'
giottoMasterToSuite <- function(gobject,
                                expression_feat = 'rna') {
  master_object = gobject
  
  spatial_locs = cell_metadata = feat_metadata = instructions = NULL
  
  if(!is.null( master_object@spatial_locs)) {
    spatial_locs = master_object@spatial_locs
  }
  
  if(!is.null(master_object@cell_metadata)) {
    cell_metadata = master_object@cell_metadata
  }
  
  if(!is.null(master_object@gene_metadata)) {
    feat_metadata = master_object@gene_metadata
    colnames(feat_metadata)[1] = 'feat_ID'
  }
  
  # instructions
  if(!is.null(master_object@instructions)) {
    instructions = master_object@instructions
  }
  
  # create Giotto object  
  gobject = createGiottoObject(expression = master_object@raw_exprs,
                               expression_feat = expression_feat,
                               spatial_locs = spatial_locs,
                               cell_metadata = cell_metadata,
                               feat_metadata = feat_metadata,
                               instructions = instructions)
  
  # add normalized expression matrix
  if(!is.null(master_object@norm_expr)){
    x = createExprObj(master_object@norm_expr,
                      name = 'normalized',
                      feat_type = expression_feat)
    gobject = setExpression(gobject,
                            x = x,
                            spat_unit = 'cell',
                            feat_type = expression_feat,
                            name = 'normalized')
  }
  
  # add scaled expression matrix
  if(!is.null(master_object@norm_scaled_expr)){
    x = createExprObj(master_object@norm_scaled_expr,
                      name = 'scaled',
                      feat_type = expression_feat)
    gobject = setExpression(gobject,
                            x = x,
                            spat_unit = 'cell',
                            feat_type = expression_feat,
                            name = 'scaled')
  }
  
  # dimension_reduction
  if(!is.null(master_object@dimension_reduction)) {
    dimension_reduction_master = master_object@dimension_reduction
    
    for (i in names(master_object@dimension_reduction$cells)) {
      j = names(dimension_reduction_master$cells[[i]])
      dimension_reduction = createDimObj(coordinates = dimension_reduction_master$cells[[i]][[j]]$coordinates,
                                         name = dimension_reduction_master$cells[[i]][[j]]$name,
                                         feat_type = expression_feat,
                                         method = dimension_reduction_master$cells[[i]][[j]]$reduction_method,
                                         misc = dimension_reduction_master$cells[[i]][[j]]$misc)
      gobject = setDimReduction(gobject,
                                dimension_reduction,
                                spat_unit = 'cell',
                                feat_type = expression_feat,
                                name = i,
                                reduction_method = dimension_reduction_master$cells[[i]][[j]]$reduction_method)
      
    }
  }
  
  # nn_network
  if(!is.null(master_object@nn_network)) {
    nn_network_master = master_object@nn_network
    
    for (i in names(nn_network_master)) {
      for (j in names(nn_network_master[[i]])) {
        k = names(nn_network_master[[i]][[j]])
        nn_network = createNearestNetObj(name = i,
                                         network = nn_network_master[[i]][[j]][[k]],
                                         feat_type = expression_feat)
        
        gobject = setNearestNetwork(gobject,
                                    nn_network,
                                    spat_unit = 'cell',
                                    feat_type = expression_feat,
                                    nn_type = i,
                                    name = j)
      }
      
    }
    
  }
  
  # spatial_network
  if(!is.null(master_object@spatial_network)) {
    spatial_network_master = master_object@spatial_network$spatial_network
    spatial_network = createSpatNetObj(network = spatial_network_master$networkDT,
                                       name = spatial_network_master$name,
                                       networkDT_before_filter = spatial_network_master$networkDT_before_filter,
                                       method = spatial_network_master$method,
                                       parameters = spatial_network_master$parameters,
                                       outputObj = spatial_network_master$outputObj,
                                       cellShapeObj = spatial_network_master$cellShapeObj,
                                       crossSectionObjects = spatial_network_master$crossSectionObjects,
                                       misc = spatial_network_master$misc
    )
    gobject = setSpatialNetwork(gobject,
                                spatial_network,
                                spat_unit = 'cell',
                                name = spatial_network_master$name)
  }
  
  # spatial_enrichment
  if(!is.null(master_object@spatial_enrichment)) {
    spatial_enrichment_master = master_object@spatial_enrichment
    
    for (i in names(spatial_enrichment_master)) {
      spatial_enrichment = createSpatEnrObj(spatial_enrichment_master[[i]],
                                            name = i,
                                            feat_type = expression_feat,
                                            method = i)
      
    }
    
    gobject = set_spatial_enrichment(gobject,
                                     spatial_enrichment,
                                     spat_unit = 'cell',
                                     feat_type = expression_feat,
                                     enrichm_name = i)
  }
  
  # spatial_grid
  if(!is.null(master_object@spatial_grid)) {
    spatial_grid_master = master_object@spatial_grid
    
    spatial_grid = new('spatialGridObj',
                       name = spatial_grid_master$spatial_grid$name,
                       method = spatial_grid_master$spatial_grid$method,
                       parameters = spatial_grid_master$spatial_grid$parameters,
                       gridDT = spatial_grid_master$spatial_grid$gridDT,
                       feat_type = expression_feat,
                       misc = spatial_grid_master$spatial_grid$misc)
    
    gobject = setSpatialGrid(gobject,
                             spatial_grid,
                             spat_unit = 'cell',
                             feat_type = expression_feat,
                             name = spatial_grid_master$spatial_grid$name)
  }
  
  return(gobject)
}


