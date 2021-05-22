#' @title perform HMRF model to detect cell DT
#' @name doHMRF2
#' @description function to perform HMRF model to detect cell DT
#' @keywords external
doHMRF2 = function (
  gobject, 
  feat_type = NULL, 
  expression_values = c("normalized", "scaled", "custom"), 
  spatial_network_name = "Delaunay_network", 
  spat_loc_name = "raw", 
  spatial_genes = NULL, 
  spatial_dimensions = c("sdimx", "sdimy", "sdimz"), 
  dim_reduction_to_use = NULL, 
  dim_reduction_name = "pca", 
  dimensions_to_use = 1:10, 
  seed = 100, 
  name = "test", 
  k = 9, 
  beta = 50, 
  tolerance = 1e-10, 
  zscore = c("none"), 
  # numinit = 100, 
  # python_path = NULL, 
  # output_folder = NULL, 
  # overwrite_output = TRUE,
  
  ### additional para
  blocks = NULL,
  nstart = 100,
  factor.step = 1.05
) 
{
  ###### will remove this step
  if (!requireNamespace("smfishHmrf", quietly = TRUE)) {
    stop("\n package ", "smfishHmrf", " is not yet installed \n", 
         "To install: \n", "remotes::install_bitbucket(repo = 'qzhudfci/smfishhmrf-r', ref='master')", 
         "see http://spatial.rc.fas.harvard.edu/install.html for more information", 
         call. = FALSE)
  }
  
  to = from = NULL
  if (is.null(feat_type)) {
    feat_type = gobject@expression_feat[[1]]
  }
  
  # ### python_path might be removed
  # if (is.null(python_path)) {
  #   python_path = readGiottoInstructions(gobject, param = "python_path")
  # }
  # 
  # ### reader_path might be removed
  # reader_path = system.file("python", "reader2.py", package = "Giotto")
  # ###  
  
  ### output_folder might be removed
  # if (is.null(output_folder)) {
  #   output_folder = paste0(getwd(), "/", "HMRF_output")
  #   if (!file.exists(output_folder)) {
  #     dir.create(path = paste0(getwd(), "/", "HMRF_output"), 
  #                recursive = T)
  #   }
  # }
  # else if (!is.null(output_folder)) {
  #   if (!file.exists(output_folder)) {
  #     dir.create(path = output_folder, recursive = T)
  #   }
  # }
  
  
  if (!is.null(dim_reduction_to_use)) {
    expr_values = select_dimReduction(gobject = gobject, 
                                      # reduction = "genes", 
                                      reduction = "cells",
                                      reduction_method = dim_reduction_to_use, 
                                      name = dim_reduction_name, return_dimObj = FALSE)
    expr_values = expr_values[, dimensions_to_use]
    expr_values = t_flex(expr_values)
  }
  else {
    values = match.arg(expression_values, unique(c("normalized", 
                                                   "scaled", "custom", expression_values)))
    if(values=='custom')
    {
      zscore = match.arg(zscore, c("col", "row", "colrow", "rowcol"))
      expr_values = select_expression_values(gobject = gobject, 
                                             feat_type = feat_type, values = 'normalized')  
      if(zscore=='col'){expr_values = scale(expr_values)}
      if(zscore=='row'){expr_values = t(scale(t(expr_values)))}
      if(zscore=='colrow'){expr_values = t(scale(t(scale(expr_values))))}
      if(zscore=='rowcol'){expr_values = scale(t(scale(t(expr_values))))}
    }else{
      expr_values = select_expression_values(gobject = gobject, 
                                             feat_type = feat_type, values = values)  
      cat("\n customized zscore is not used \n")
    }
  }
  
  ### writing expression_file might be removed
  # expression_file = paste0(output_folder, "/", "expression_matrix.txt")
  # if (file.exists(expression_file) & overwrite_output == TRUE) {
  #   cat("\n expression_matrix.txt already exists at this location, will be overwritten \n")
  #   write.table(expr_values, file = expression_file, quote = F, 
  #               col.names = NA, row.names = T)
  # }
  # else if (file.exists(expression_file) & overwrite_output == 
  #          FALSE) {
  #   cat("\n expression_matrix.txt already exists at this location, will not be overwritten \n")
  # }
  # else {
  #   write.table(expr_values, file = expression_file, quote = F, 
  #               col.names = NA, row.names = T)
  # }
  
  ### expr_values is the expr matrix to use  
  ############################################  
  
  if (!is.null(dim_reduction_to_use)) {
    dimred_rownames = rownames(expr_values)
    spatial_genes_detected = dimred_rownames[dimensions_to_use]
    spatial_genes_detected = spatial_genes_detected[!is.na(spatial_genes_detected)]
  }
  else {
    if (is.null(spatial_genes)) {
      stop("\n you need to provide a vector of spatial genes (~500) \n")
    }
    spatial_genes_detected = spatial_genes[spatial_genes %in% 
                                             rownames(expr_values)]
  }
  
  ### writing spatial_genes might be removed
  # spatial_genes_file = paste0(output_folder, "/", "spatial_genes.txt")
  # if (file.exists(spatial_genes_file) & overwrite_output == 
  #     TRUE) {
  #   cat("\n spatial_genes.txt already exists at this location, will be overwritten \n")
  #   write.table(spatial_genes_detected, file = spatial_genes_file, 
  #               quote = F, col.names = F, row.names = F)
  # }
  # else if (file.exists(spatial_genes_file) & overwrite_output == 
  #          FALSE) {
  #   cat("\n spatial_genes.txt already exists at this location, will be used again \n")
  # }
  # else {
  #   write.table(spatial_genes_detected, file = spatial_genes_file, 
  #               quote = F, col.names = F, row.names = F)
  # }
  
  ##### spatial_genes_detected is the spatial gene set used in the model
  #############################################
  
  spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, 
                                          return_network_Obj = FALSE)
  # spatial_network = gobject@spatial_network$Delaunay_network$networkDT
  spatial_network = spatial_network[, .(to, from)]
  
  ### writing spatial_network might be removed
  # spatial_network_file = paste0(output_folder, "/", "spatial_network.txt")
  # if (file.exists(spatial_network_file) & overwrite_output == 
  #     TRUE) {
  #   cat("\n spatial_network.txt already exists at this location, will be overwritten \n")
  #   write.table(spatial_network, file = spatial_network_file, 
  #               row.names = F, col.names = F, quote = F, sep = "\t")
  # }
  # else if (file.exists(spatial_network_file) & overwrite_output == 
  #          FALSE) {
  #   cat("\n spatial_network.txt already exists at this location, will not be overwitten \n")
  # }
  # else {
  #   write.table(spatial_network, file = spatial_network_file, 
  #               row.names = F, col.names = F, quote = F, sep = "\t")
  # }
  
  ###  spatial_network is the spatial network, edges on each row to and from
  # ###########################################
  # spatial_location = select_spatial_locations(gobject = gobject, 
  #                                             spat_loc_name = spat_loc_name)
  # spatial_dimensions = spatial_dimensions[spatial_dimensions %in% 
  #                                           colnames(spatial_location)]
  # spatial_location = spatial_location[, c(spatial_dimensions, 
  #                                         "cell_ID"), with = F]
  
  ### writing spatial_location might be removed
  # spatial_location_file = paste0(output_folder, "/", "spatial_cell_locations.txt")
  # if (file.exists(spatial_location_file) & overwrite_output == 
  #     TRUE) {
  #   cat("\n spatial_cell_locations.txt already exists at this location, will be overwritten \n")
  #   write.table(spatial_location, file = spatial_location_file, 
  #               row.names = F, col.names = F, quote = F, sep = "\t")
  # }
  # else if (file.exists(spatial_location_file)) {
  #   cat("\n spatial_cell_locations.txt already exists at this location, will be used again \n")
  # }
  # else {
  #   write.table(spatial_location, file = spatial_location_file, 
  #               row.names = F, col.names = F, quote = F, sep = "\t")
  # }
  
  #### spatial_location   cell locations  with cell id in last column  
  #### cell_location      cell locations  with cell id in last column  
  # cell_location = spatial_location
  ##########################################################################
  
  
  ### expr_values is the expr matrix to use  
  ### spatial_genes
  ###  spatial_network is the spatial network, edges on each row to and from
  #### cell_location      cell locations  with cell id in last column  
  
  
  # dimensions_to_use = 1:10, 
  # seed = 100, 
  # name = "test", 
  # k = 10, 
  # betas = c(0, 2, 50), 
  # tolerance = 1e-10, 
  # zscore = c("none", "rowcol", "colrow"), 
  # numinit = 100, 
  # blocks = s$blocks, 
  # beta=28, 
  # mu=s$mu, 
  # sigma=s$sigma, 
  # err=1e-7, 
  # maxit=50, 
  # verbose=TRUE, dampFactor=s$damp, 
  # tolerance=1e-5
  # nstart = 100
  
  # spatial_genes = km_spatialgenes[1:150]$feats
  
  # expr_values = gobject@expression$rna$normalized[spatial_genes,]
  
  # spatial_network = gobject@spatial_network$Delaunay_network$networkDT
  # spatial_network = spatial_network[, .(to, from)]
  
  
  y0 = t(as.matrix(expr_values))
  
  cell.rm = setdiff(rownames(y0),unique(c(spatial_network$to,spatial_network$from)))
  
  y = y0[-match(cell.rm,rownames(y0)),]
  numcell<-dim(y)[1]
  m<-dim(y)[2]
  
  ### ncol of nei = max of neighbors
  ncol.nei = max(table(c(spatial_network$to,spatial_network$from)))
  nei = matrix(-1,ncol = ncol.nei,nrow = numcell)
  rownames(nei) = rownames(y)
  for(i in 1:numcell)
  {
    nei.i = c(spatial_network$from[spatial_network$to==rownames(nei)[i]],
              spatial_network$to[spatial_network$from==rownames(nei)[i]])
    if(length(nei.i)>0)nei[i,1:length(nei.i)] = sort(match(nei.i,rownames(y)))
  }
  numnei<-as.integer(rowSums(nei!=(-1)))
  
  
  kk = smfishHmrf.generate.centroid(y=y,par_k = k,par_seed=seed,nstart=nstart)
  mu<-t(kk$centers) #should be dimension (m,k)
  lclust<-lapply(1:k, function(x) which(kk$cluster == x))
  damp<-array(0, c(k)); 
  sigma<-array(0, c(m,m,k))
  for(i in 1:k){
    sigma[, , i] <- cov(matrix(y[lclust[[i]], ],ncol = ncol(y)))
    sigma[, , i][is.na(sigma[, , i])] = 0
    di<-findDampFactor(sigma[,,i], factor=factor.step, d_cutoff=1e-5, startValue=0.0001)
    ####
    di = max(di,0.01)
    ####
    damp[i]<-ifelse(is.null(di), 0, di)
  }
  
  # if(length(beta.h)>1){
  #   cat("\n only first element of beta.h is used \n")
  #   beta.h = beta.h[1]
  # }else if (is.null(beta.h)) {
  #   beta.h = 50
  #   cat("\n default beta.h = 50 is used \n")
  # }
  # beta.h = beta[1]
  
  if(!is.null(blocks))
  {
    if (mean(sort(unique(do.call(c,blocks)))==1:dim(y0)[1])<1) {
      stop("\n block elements and #cells do not match \n")
    }else{
      for(i in 1:length(blocks))
      {
        blocks[[i]] = match(setdiff(rownames(y0)[blocks[[i]]],cell.rm),rownames(y))
      }
    }
  }else{blocks = lclust}  
  
  result.hmrf = smfishHmrf.hmrfem.multi(y = y, neighbors = nei, numnei = numnei, blocks = blocks, beta=beta, 
                                        mu=mu, sigma=sigma, err=1e-7, maxit=50, verbose=TRUE, dampFactor=damp, 
                                        tolerance=tolerance)
  
  # res = smfishHmrf.hmrfem.multi(y = s$y, neighbors = s$nei, numnei = s$numnei, blocks = s$blocks, beta=beta, 
  #                               mu=s$mu, sigma=s$sigma, err=1e-7, maxit=50, verbose=TRUE, dampFactor=damp, 
  #                               tolerance=tolerance)
  
  # add class to hmrf result
  
  rownames(result.hmrf$prob) = rownames(y)
  colnames(result.hmrf$prob) = paste('Domain',1:10,sep = '_')
  
  rownames(result.hmrf$mu) = colnames(y)
  colnames(result.hmrf$mu) = paste('Domain',1:10,sep = '_')
  
  dimnames(result.hmrf$sigma)[[1]] = colnames(y)
  dimnames(result.hmrf$sigma)[[2]] = colnames(y)
  dimnames(result.hmrf$sigma)[[3]] = paste('Domain',1:10,sep = '_')
  
  rownames(result.hmrf$unnormprob) = rownames(y)
  colnames(result.hmrf$unnormprob) = paste('Domain',1:10,sep = '_')
  
  result.hmrf$feat_type = feat_type
  
  result.hmrf$k = k
  
  result.hmrf$beta = beta
  
  result.hmrf$cell.rm = cell.rm
  
  class(result.hmrf) <- append(class(result.hmrf),"HMRFoutput") 
  return(result.hmrf)
  
  # cell_location = paste0(output_folder, "/", "spatial_cell_locations.txt")
  # spatial_genes = paste0(output_folder, "/", "spatial_genes.txt")
  # spatial_network = paste0(output_folder, "/", "spatial_network.txt")
  # expression_data = paste0(output_folder, "/", "expression_matrix.txt")
  # output_data = paste0(output_folder, "/", "result.spatial.zscore")
  # if (!file.exists(output_data)) 
  #   dir.create(output_data)
  # cell_location = paste0("\"", cell_location, "\"")
  # spatial_genes = paste0("\"", spatial_genes, "\"")
  # spatial_network = paste0("\"", spatial_network, "\"")
  # expression_data = paste0("\"", expression_data, "\"")
  # output_data = paste0("\"", output_data, "\"")
  # zscore = match.arg(zscore, c("none", "rowcol", "colrow"))
  # betas_param = c("-b", betas)
  # betas_final = paste(betas_param, collapse = " ")
  # reader_command = paste0(python_path, " ", reader_path, " -l ", 
  #                         cell_location, " -g ", spatial_genes, " -n ", spatial_network, 
  #                         " -e ", expression_data, " -o ", output_data, " -a ", 
  #                         name, " -k ", k, " ", betas_final, " -t ", tolerance, 
  #                         " -z ", zscore, " -s ", seed, " -i ", numinit)
  # print(reader_command)
  # system(command = reader_command)
  # HMRFObj = list(name = name, output_data = output_data, k = k, 
  #                betas = betas, python_path = python_path)
  # class(HMRFObj) <- append(class(HMRFObj), "HMRFoutput")
  # return(HMRFObj)
}



#' @title get domain index of hmrf result
#' @name GetDT.HMRFoutput
#' @description function to select Domain type of cells in HMRFoutput
#' @keywords external
GetDT = function(x){UseMethod('GetDT',x)}
GetDT.HMRFoutput = function(x){
  cl = c(apply(x$prob,1,function(i){ifelse(test = sum(i>0.5)==1, yes = which(i>0.5), no = NA)}),
         rep(NA,length(x$cell.rm)));
  names(cl) = c(rownames(x$prob),x$cell.rm);
  return(cl)}

#which(is.na(GetDT(result.hmrf)))


#' @title add HMRF DT to cell meta data
#' @name addHMRF2
#' @description function to add HMRF Domain Type to cell meta data
#' @keywords external
addHMRF = 
  function (gobject, HMRFoutput, k = NULL, betas_to_add = NULL 
            #, hmrf_name = NULL
            ) 
  {
    if (!"HMRFoutput" %in% class(HMRFoutput)) {
      stop("\n HMRFoutput needs to be output from doHMRFextend \n")
    }
    feat_type = HMRFoutput$feat_type
    # get_result_path = system.file("python", "get_result2.py", 
    #                               package = "Giotto")
    # name = HMRFoutput$name
    # output_data = HMRFoutput$output_data
    # python_path = HMRFoutput$python_path
    # if (is.null(k)) {
    #   stop("\n you need to select a k that was used with doHMRFextend \n")
    # }
    k = HMRFoutput$k
    b = HMRFoutput$beta
    
    ordered_cell_IDs = gobject@cell_ID
    
    hmrf_DT = GetDT(HMRFoutput)
    hmrf_DT = hmrf_DT[match(ordered_cell_IDs,names(hmrf_DT))]
    
    gobject = addCellMetadata(gobject = gobject, feat_type = feat_type, 
                              column_cell_ID = "cell_ID", new_metadata = hmrf_DT,
                              vector_name = paste0('HMRF_k',k,'_b.',b), 
                              by_column = F)
    
    # possible_betas = seq(betas[1], to = betas[1] + (betas[2] * 
    #                                                   (betas[3] - 1)), by = betas[2])
    # betas_to_add_detected = betas_to_add[betas_to_add %in% possible_betas]
    # cell_metadata = pDataDT(gobject, feat_type = feat_type)
    # for (b in betas_to_add_detected) {
    #   result_command = paste0(python_path, " ", get_result_path, 
    #                           " -r ", output_data, " -a ", name, " -k ", k, " -b ", 
    #                           b)
    #   print(result_command)
    #   output = system(command = result_command, intern = T)
    #   annot_DT = data.table::data.table(temp_name = output)
    #   if (!is.null(hmrf_name)) {
    #     annot_name = paste0(hmrf_name, "_k", k, "_b.", b)
    #     setnames(annot_DT, old = "temp_name", new = annot_name)
    #   }
    #   else {
    #     annot_name = paste0("hmrf_k.", k, "_b.", b)
    #     data.table::setnames(annot_DT, old = "temp_name", 
    #                          new = annot_name)
    #   }
    # gobject = addCellMetadata(gobject = gobject, feat_type = feat_type, 
    #                           column_cell_ID = "cell_ID", new_metadata = annot_DT, 
    #                           by_column = F)
    # }
    return(gobject)
  }


