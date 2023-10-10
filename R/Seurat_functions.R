## Seurat object ####


#' @title Convert Giotto to Seurat
#' @name giottoToSeurat
#' @description Converts Giotto object into a Seurat object. This functions extracts
#' specific sets of data belonging to specified spatial unit.
#' The default values are 'cell' and 'rna' respectively.
#' @param gobject Giotto object
#' @param obj_use Giotto object (deprecated, use gobject)
#' @param spat_unit spatial unit (e.g. 'cell')
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
                                     # does not accept 'dgeMatrix'
                                     new.data = as.matrix(expr_use[['scaled']][]),
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
        data_scale <- as.matrix(expr_use[['scaled']][])
        sobj <- Seurat::SetAssayData(sobj,slot = "scale.data",
                                     new.data = data_scale,
                                     assay = assay_use)
      }
    }
    
    # add cell metadata
    meta_cells <- data.table::setDF(
      get_cell_metadata(
        gobject = gobject,
        spat_unit = spat_unit,
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
        spat_unit = spat_unit,
        feat_type = assay_use,
        output = 'data.table',
        copy_obj = TRUE
      )
    )
    rownames(meta_genes) <- meta_genes$feat_ID
    sobj@assays$rna@meta.data<- meta_genes
    
    # dim reduction
    # note: Seurat requires assay name specification for each dim reduc
    avail_dr = list_dim_reductions(gobject = gobject, spat_unit = spat_unit, feat_type = assay_use,)
    if (nrow(avail_dr) > 0){
      dr_use = avail_dr[, name]
      for (i in seq(nrow(avail_dr))){
        dr_name = avail_dr[i, name]
        dr_type = avail_dr[i, dim_type]
        dr_obj <- get_dimReduction(
          gobject = gobject,
          output = 'dimObj',
          spat_unit = spat_unit,
          feat_type = assay_use,
          reduction_method = dr_type,
          name = dr_name
        )
        emb_use <- dr_obj[][Seurat::Cells(sobj),]
        if (sum(c('loadings','eigenvalues') %in% names(slot(dr_obj, 'misc'))) == 2){
          loadings_use <- slot(dr_obj, 'misc')$loadings
          stdev_use <- slot(dr_obj, 'misc')$eigenvalues
          sobj[[dr_name]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(emb_use),
                                                          loadings = loadings_use,
                                                          key = paste0(dr_name, '_'),
                                                          stdev = stdev_use,
                                                          assay = assay_use)
        } else {
          sobj[[dr_name]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(emb_use),
                                                          key = paste0(dr_name, '_'),
                                                          assay = assay_use)
        }
      }
    }
    
    
    
    # network objects
    # expression network
    avail_nn = list_nearest_networks(gobject = gobject, spat_unit = spat_unit, feat_type = assay_use)
    if (nrow(avail_nn) > 0){
      for (i in seq(nrow(avail_nn))){
        nn_name = avail_nn[i, name]
        nn_type = avail_nn[i, nn_type]
        nn_use <- get_NearestNetwork(
          gobject = gobject,
          spat_unit = spat_unit,
          feat_type = assay_use,
          nn_network_to_use = nn_type,
          network_name = nn_name,
          output = 'data.table'
        )
        idx1 <- match(nn_use$from,Seurat::Cells(sobj))
        idx2 <- match(nn_use$to,Seurat::Cells(sobj))
        edge_weight <- nn_use$weight
        nn_mtx <- Matrix::sparseMatrix(i = idx1,j = idx2,x = edge_weight,dims = c(ncol(sobj),ncol(sobj)))
        rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
        nn_name <- paste0('expr_',nn_name)
        sGraph <- Seurat::as.Graph(nn_mtx)
        sGraph@assay.used = assay_use
        sobj[[nn_name]] = sGraph
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
  avail_sn = list_spatial_networks(gobject = gobject, spat_unit = spat_unit)
  if (nrow(avail_sn) > 0){
    sn_all = avail_sn[, name]
    for (i in sn_all){
      snt_use <- get_spatialNetwork(gobject = gobject,
                                    spat_unit = spat_unit,
                                    name = i,
                                    output = 'networkDT')
      idx1 <- match(snt_use$from,Seurat::Cells(sobj))
      idx2 <- match(snt_use$to,Seurat::Cells(sobj))
      edge_weight <- snt_use$weight
      nn_mtx <- Matrix::sparseMatrix(i = idx1,j = idx2,x = edge_weight,dims = c(ncol(sobj),ncol(sobj)))
      rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
      nn_name <- paste0('spatial_',i)
      spatGraph = Seurat::as.Graph(nn_mtx)
      spatGraph@assay.used = names(sobj@assays)[1]
      sobj[[nn_name]] <- spatGraph
    }
  }
  
  return (sobj)
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
    wrap_msg('No raw expression values are provided in spatial_assay')
    return(sobject)
    
  } else {
    
    exp = Seurat::GetAssayData(object = sobject, slot = "counts", assay = spatial_assay)
    if(!is.null(sobject@assays$SCT)){
      normexp = Seurat::GetAssayData(object = sobject, slot = "counts", assay = 'SCT')
    }
    
    if(!is.null(slot(sobject, 'assays')[[spatial_assay]]@layers)){
      normexp = LayerData(object = sobject, assay = spatial_assay)
    }
    
    # Cell Metadata
    cell_metadata = sobject@meta.data
    # Feat Metadata 
    feat_metadata = sobject@assays$Vizgen@meta.data
    
    # Dimension Reduction
    if(sum(sapply(dim_reduction,function(x) length(sobject@reductions[[x]]))) == 0) {
      dimReduc_list = NULL
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
        dimReduc = create_dim_obj(name = x,
                                  reduction = 'cells',
                                  reduction_method = x,
                                  spat_unit = 'cell',
                                  feat_type = 'rna',
                                  provenance = NULL,
                                  coordinates = dim_coord,
                                  misc = list(eigenvalues = dim_eig, loadings = dim_load))
        
        return(dimReduc)
      })
      # names(dimReduc_list) <- dim_reduction
    }

    # Spatial Locations
    if(length(sobject@assays[[spatial_assay]]) == 0) {
      spat_loc = NULL
    } else {

      # !requires image objects!
      if (!is.null(Seurat::Images(object = sobject, assay = spatial_assay))) {
       spat_coord = Seurat::GetTissueCoordinates(sobject)
       # spat_coord = cbind(rownames(spat_coord), data.frame(spat_coord, row.names=NULL))
        colnames(spat_coord) = c("sdimx", "sdimy", "cell_ID")
        spat_loc = spat_coord
        length_assay <- length(sobject@assays$Vizgen@cells@.Data)
        spat_datatable <- data.table(cell_ID = character(length_assay),
                                     sdimx = rep(NA_real_, length_assay),
                                     sdimy = rep(NA_real_, length_assay))
        spat_datatable$cell_ID <- rownames(sobject@assays$Vizgen@cells@.Data)
        match_cell_ID <- match(spat_loc$cell_ID, spat_datatable$cell_ID)
        matching_indices <- match_cell_ID
        matching_indices <-  matching_indices[!is.na(matching_indices)]
        spat_datatable[matching_indices, c("sdimx", "sdimy") := list(spat_loc$sdimx, spat_loc$sdimy)]
        spat_loc <- spat_datatable
        
      } else {
        message("Images for RNA assay not found in the data. Skipping image processing.")
        spat_loc = NULL
      }

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
       if ("molecules" %in% names(sobject@images[["hippo"]]) == TRUE) {
         if(!length(sobject@images[["hippo"]][["molecules"]]) == 0) {

           assay = names(sobject@assays)
           featnames = rownames(sobject@assays$Vizgen@features@.Data)
           mol_spatlocs = data.table::data.table()

           for (x in featnames) {
             df = (Seurat::FetchData(sobject[["hippo"]][["molecules"]], vars = x))
             mol_spatlocs = rbind(mol_spatlocs, df)
           }
           gpoints = createGiottoPoints(mol_spatlocs, feat_type = "rna")

         }
       }
     }
  }
  
  #Find SueratImages, extract them, and pass to create seuratobj
  
  
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
  gobject = addFeatMetadata(gobject = gobject, new_metadata = feat_metadata)
  
  
  if (exists('gpoints') == TRUE) {
    gobject = addGiottoPoints(gobject = gobject,
                              gpoints = list(gpoints))
  }
  
  return (gobject)
}

