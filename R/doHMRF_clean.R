doHMRF = function (
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
  blocks,
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
  
  if (!is.null(dim_reduction_to_use)) {
    expr_values = select_dimReduction(gobject = gobject, 
                                      reduction = "genes", 
                                      # reduction = "cells", 
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
  

  ##### spatial_genes_detected is the spatial gene set used in the model
  #############################################
  
  spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, 
                                          return_network_Obj = FALSE)
  # spatial_network = gobject@spatial_network$Delaunay_network$networkDT
  spatial_network = spatial_network[, .(to, from)]
  
  
  spatial_genes = km_spatialgenes[1:150]$feats
  expr_values = my_giotto_object@expression$rna$normalized[spatial_genes,]
  spatial_network = my_giotto_object@spatial_network$Delaunay_network$networkDT
  spatial_network = spatial_network[, .(to, from)]
  
  y = t(as.matrix(expr_values))
  kk = smfishHmrf.generate.centroid(y=y,par_k = k,par_seed=seed,nstart=nstart)
  numcell<-dim(expr_values)[2]
  m<-dim(expr_values)[1]
  mu<-t(kk$centers) #should be dimension (m,k)
  lclust<-lapply(1:k, function(x) which(kk$cluster == x))
  damp<-array(0, c(k)); 
  sigma<-array(0, c(m,m,k))
  for(i in 1:k){
    sigma[, , i] <- cov(matrix(y[lclust[[i]], ],ncol = ncol(y)))
    sigma[, , i][is.na(sigma[, , i])] = 0
    di<-findDampFactor(sigma[,,i], factor=factor.step, d_cutoff=1e-5, startValue=0.0001)
    damp[i]<-ifelse(is.null(di), 0, di)
  }
  
  ### ncol of nei = max of neighbors
  ncol.nei = max(table(c(spatial_network$to,spatial_network$from)))
  nei = matrix(-1,ncol = ncol.nei,nrow = numcell)
  rownames(nei) = colnames(expr_values)
  for(i in 1:numcell)
  {
    nei.i = c(spatial_network$from[spatial_network$to==rownames(nei)[i]],
              spatial_network$to[spatial_network$from==rownames(nei)[i]])
    if(length(nei.i)>0)nei[i,1:length(nei.i)] = sort(match(nei.i,colnames(expr_values)))
  }
  numnei<-as.integer(rowSums(nei!=(-1)))
  
  result.hmrf = smfishHmrf.hmrfem.multi(y = y, neighbors = nei, numnei = numnei, blocks = blocks, beta=beta, 
                                        mu=mu, sigma=sigma, err=1e-7, maxit=50, verbose=TRUE, dampFactor=damp, 
                                        tolerance=tolerance)
  
  # add class to hmrf result
  result.hmrf$beta = beta
  class(result.hmrf) <- append(class(result.hmrf),"HMRFoutput") 
  return(result.hmrf)

}




GetClass = function(x){UseMethod('GetClass',x)}
GetClass.HMRFoutput = function(x){
  apply(x$prob,1,function(i){ifelse(test = sum((1-i)<1e-3)==1, yes = which((1-i)<1e-3), no = NA)})}

# GetClass(result.hmrf)


# additional function

select_dimReduction

t_flex

select_expression_values

select_spatialNetwork

select_spatial_locations


load('~/Downloads/gobject_working_for_hmrf.RData')

# y = t(as.matrix(gemSub1_test@dimension_reduction$cells$pca$pca$))
# expression$rna$normalized
y = (as.matrix(mini_smi@dimension_reduction$cells$pca$pca$coordinates))
# [-c(555,569),]

kk = smfishHmrf.generate.centroid(y=y,par_k = k,par_seed=seed,nstart=nstart)
numcell<-dim(y)[1]
m<-dim(y)[2]
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

# spatial_net = gemSub1_test@spatial_network$Delaunay_network$networkDT
spatial_net = mini_smi@spatial_network$Delaunay_network$networkDT
spatial_net = spatial_net[, .(to, from)]

### ncol of nei = max of neighbors
ncol.nei = max(table(c(spatial_net$to,spatial_net$from)))
nei = matrix(-1,ncol = ncol.nei,nrow = numcell)
rownames(nei) = rownames(y)
for(i in 1:numcell)
{
  nei.i = c(spatial_net$from[spatial_net$to==rownames(nei)[i]],
            spatial_net$to[spatial_net$from==rownames(nei)[i]])
  if(length(nei.i)>0)nei[i,1:length(nei.i)] = sort(match(nei.i,rownames(y)))
}
numnei<-as.integer(rowSums(nei!=(-1)))

blocks = lclust

which(numnei==0)
beta = 50


allgenes = gemSub1_test@feat_ID$rna

res = smfishHmrf.hmrfem.multi(y = y, neighbors = nei, numnei = numnei, blocks = blocks, beta=beta, 
                              mu=mu, sigma=sigma, err=1e-7, maxit=50, verbose=TRUE, dampFactor=damp, 
                              tolerance=tolerance)







##############################
### miscs


select_expression_values <- function(gobject,
                                     feat_type = 'rna',
                                     values) {
  
  potential_values = names(gobject@expression[[feat_type]])
  
  ## special cases for giotto standard pipeline
  if(values == 'scaled' & is.null(gobject@expression[[feat_type]][[values]])) {
    stop('run first scaling (& normalization) step')
  } else if(values == 'normalized' & is.null(gobject@expression[[feat_type]][[values]])) {
    stop('run first normalization step')
  } else if(values == 'custom' & is.null(gobject@expression[[feat_type]][[values]])) {
    stop('first add custom expression matrix')
  }
  
  
  if(values %in% potential_values) {
    expr_values = gobject@expression[[feat_type]][[values]]
    return(expr_values)
  } else {
    stop("The ", feat_type ," expression matrix with name ","'", values, "'"," can not be found \n")
  }
  
}


select_spatial_locations <- function(gobject,
                                     spat_loc_name = 'raw') {
  
  potential_names = names(gobject@spatial_locs)
  
  if(spat_loc_name %in% potential_names) {
    spatloc = data.table::copy(gobject@spatial_locs[[spat_loc_name]])
    return(spatloc)
  } else {
    stop("The spatial locations with name ","'", spat_loc_name, "'"," can not be found \n")
  }
  
}



