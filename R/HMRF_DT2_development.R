#' @title sampling_sp_genes
#' @name sampling_sp_genes
#' @description function to select a set of spatial genes
#' @keywords internal
sampling_sp_genes = function(clust, sample_rate=2, target=500){
  # clust = spat_cor_netw_DT$cor_clusters$spat_netw_clus
  # sample_rate=2
  # target=500
  tot=0
  num_cluster=length(unique(clust))
  gene_list = list()
  
  for(i in seq(1, num_cluster)){
    gene_list[[i]] = colnames(t(clust[which(clust==i)]))
  }
  for(i in seq(1, num_cluster)){
    num_g=length(gene_list[[i]])
    tot = tot+num_g/(num_g^(1/sample_rate))
  }
  factor=target/tot
  num_sample=c()
  genes=c()
  for(i in seq(1, num_cluster)){
    num_g=length(gene_list[[i]])
    genes[i] = num_g
    num_sample[i] = round(num_g/(num_g^(1/sample_rate)) * factor)
  }
  set.seed(10)
  samples=list()
  union_genes = c()
  for(i in seq(1, num_cluster)){
    if(length(gene_list[[i]])<num_sample[i]){
      samples[[i]] = gene_list[[i]]
    }else{
      samples[[i]] = sample(gene_list[[i]], num_sample[i])
    }
    union_genes = union(union_genes, samples[[i]])
  }
  union_genes = unique(union_genes)
  
  return(list(union_genes = union_genes, num_sample = num_sample, num_gene=genes, gene_list=gene_list))
  
}


#' @title addSpatialGenesTest
#' @name addSpatialGenesTest
#' @description function to add adjusted p value of spatial genes test to gene meta data
#' @keywords external
addSpatialGenesTest = function(
  gobject,
  method = 'kmeans',
  spatial_network_name = 'Delaunay_network',
  seed = 1234
){
  if(length(method)>1||!method%in%c('kmeans','rank'))
  {stop('\n method should be one of “kmeans”, “rank” \n')}
  
  sp_test = binSpect(gobject = gobject,bin_method = method,set.seed = seed,spatial_network_name = spatial_network_name)
  
  # setdiff(binSpect.obj$genes,gobject@gene_metadata$gene_ID)
  # if(!identical(sort(gobject@gene_metadata$gene_ID),sort(binSpect.obj$genes)))
  #   (stop('\n please provide a spatial test result with matched genes \n'))
  
  gobject@gene_metadata[,paste0(method,'_SpatialGenesTest')] = sp_test[match(gobject@gene_metadata$gene_ID,sp_test$genes),]$adj.p.value
}


#' @title initHMRF
#' @name initialize HMRF
#' @description function to initialize HMRF model
#' @keywords external
#' @param xxxxxxxxxxxxxxxxxxxxx
#' @param xxxxxxxxxxxxxxxxxxxxx
#' @param xxxxxxxxxxxxxxxxxxxxx
initHMRF = function(
  gobject, 
  # feat_type = NULL, 
  expression_values = c("scaled", "normalized", "custom"), 
  spatial_network_name = "Delaunay_network", 
  # input_features = 'spatial_genes',
  spatial_genes = NULL, 
  n_spatial_genes_select = 500, 
  # spat_loc_name = "raw", 
  # spatial_dimensions = c("sdimx", "sdimy", "sdimz"), 
  #disable using dimension reduction
  # dim_reduction_to_use = NULL, 
  # dim_reduction_name = "pca", 
  # dimensions_to_use = 1:10, 
  
  seed.init = 100, 
  # name = "test", 
  k = 9, 
  # # beta_init = 0, 
  # # beta_increment = 10,
  # # beta_num_iter = 5,
  tolerance = 1e-5,
  zscore = c("none"), 
  
  # numinit = 100, 
  ### additional para
  #blocks = NULL,
  nstart = 1000,
  factor.step = 1.05
){
  ###### will remove this step
  if (!requireNamespace("smfishHmrf", quietly = TRUE)) {
    stop("\n package ", "smfishHmrf", " is not yet installed \n", 
         "To install: \n", "devtools::install_bitbucket(\"qzhudfci/smfishHmrf-r\") \n", 
         "see http://spatial.rc.fas.harvard.edu for more information", 
         call. = FALSE)
  }

  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    stop("\n package ", "tidygraph", " is not yet installed \n", 
         call. = FALSE)
  }
  
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("\n package ", "ggraph", " is not yet installed \n", 
         call. = FALSE)
  }
  
  if (!requireNamespace("graphcoloring", quietly = TRUE)) {
    stop("\n package ", "graphcoloring", " is not yet installed \n", 
         "To install: \n", "devtools::install_bitbucket(\"qzhudfci/graphcoloring\") \n",  
         call. = FALSE)
    
  }
  
  library(smfishHmrf)
  library(tidygraph)
  library(ggraph)
  library(graphcoloring)
  
  # if (length(input_features)>1){stop("\n input_features needs to be either 'spatial_genes' or 'dim_reduction' \n")}
  # if (!input_features%in%c('spatial_genes','dim_reduction')){
  #   stop("\n input_features needs to be either 'spatial_genes' or 'dim_reduction' \n")
  # }
  
  spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, 
                                          return_network_Obj = FALSE)
  spatial_network = spatial_network[, .(to, from)]
  
  # 'spatial_genes' use spatial gene list or create spatial gene list
  values = match.arg(expression_values, unique(c("scaled","normalized", "custom", expression_values)))
  expr_values = select_expression_values(gobject = gobject, values = values)   
  
  if(zscore!="none"){
    zscore = match.arg(zscore, c("col", "row", "colrow", "rowcol"))
    expr_values = select_expression_values(gobject = gobject, 
                                           feat_type = feat_type, values = 'normalized')  
    if(zscore=='col'){expr_values = scale(expr_values)}
    if(zscore=='row'){expr_values = t(scale(t(expr_values)))}
    if(zscore=='colrow'){expr_values = t(scale(t(scale(expr_values))))}
    if(zscore=='rowcol'){expr_values = scale(t(scale(t(expr_values))))}
  }
  
  # if(!n_spatial_genes_select>0){stop("\n please provide a positive integer n_spatial_genes_select \n")}
  
  if(length(intersect(spatial_genes,rownames(expr_values)))>0){
    ### if input a sp gene list, use the list
    expr_values = expr_values[intersect(spatial_genes,rownames(expr_values)),]
    n0 = length(spatial_genes)
    n1 = length(intersect(spatial_genes,rownames(expr_values)))
    n = min(n_spatial_genes_select,500,n1)
    cat(paste0('\n total ',n0, 'spatial_genes imported \n'))
    cat(paste0('\n total ',n1, 'spatial_genes overlapped with expr data matrix \n'))
    cat(paste0('\n sample ',n, '(max of 500) spatial_genes selected for HMRF model \n'))
  }else{
    
    ### if no input any sp gene, pick significant genes from km test
    cat('\n none overlapped genes imported, use spatial genes from kmeans_SpatialGenesTest \n')
    
    if(!'kmeans_SpatialGenesTest'%in%names(gobject@gene_metadata))
    { 
      km_spatialgenes = binSpect(gobject = gobject,bin_method = 'kmeans',set.seed = 1234,spatial_network_name = spatial_network_name)
      spatial_genes = km_spatialgenes$genes[which(km_spatialgenes$adj.p.value<0.01)]
    }else{
      spatial_genes = gobject@gene_metadata$gene_ID[gobject@gene_metadata$kmeans_SpatialGenesTest<0.01]
    }
    
    n = min(n_spatial_genes_select,500,length(spatial_genes))
    cat(paste0('\n sample ',n,' (max of 500) genes out of total ',length(spatial_genes),' spatial genes for HMRF model \n'))
  }
  
  if(n<length(spatial_genes)){
    spat_cor_netw_DT = detectSpatialCorGenes(gobject = gobject ,method = 'network',
                                             spatial_network_name = spatial_network_name,
                                             subset_genes = spatial_genes,
                                             network_smoothing = 0)
    spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT,name = 'spat_netw_clus',k = 20)
    sample_genes = sampling_sp_genes(spat_cor_netw_DT$cor_clusters$spat_netw_clus, sample_rate=2, target=n)
    spatial_genes_selected = sample_genes$union_genes
  }else{spatial_genes_selected = spatial_genes}
  
  expr_values = expr_values[spatial_genes_selected,]
  
  y0 = t(as.matrix(expr_values))
  
  cell.rm = setdiff(rownames(y0),unique(c(spatial_network$to,spatial_network$from)))
  y = y0
  if(length(cell.rm)>0)
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
  
  nn<-nei
  #numcell <- dim(nn)[1]
  numedge <- 0
  for(i in 1:numcell){
    numedge<-numedge + length(nn[i,nn[i,]!=-1])
  }
  edgelist <- matrix(0, nrow=numedge, ncol=2)
  
  edge_ind <- 1
  for(i in 1:numcell){
    neighbors <- nn[i, nn[i,]!=-1]
    for(j in 1:length(neighbors)){
      edgelist[edge_ind,] <- c(i, neighbors[j]) 
      edge_ind<-edge_ind+1
    }
  }
  
  pp<-tbl_graph(edges=as.data.frame(edgelist), directed=F)
  yy<-pp%>%mutate(color=as.factor(color_dsatur()))
  colors<-as.list(yy)$nodes$color
  cl_color <- sort(unique(colors))
  blocks<-lapply(cl_color, function(cl){which(colors==cl)})
  
  kk = smfishHmrf.generate.centroid(y=y,par_k = k,par_seed=seed.init,nstart=nstart)
  mu<-t(kk$centers) #should be dimension (m,k)
  lclust<-lapply(1:k, function(x) which(kk$cluster == x))
  damp<-array(0, c(k)); 
  sigma<-array(0, c(m,m,k))
  for(i in 1:k){
    sigma[, , i] <- cov(y[lclust[[i]], ])
    di<-findDampFactor(sigma[,,i], factor=factor.step, d_cutoff=tolerance, startValue=0.0001)
    damp[i]<-ifelse(is.null(di), 0, di)
  }
  
  list(y=y, nei=nei, numnei=numnei, blocks=blocks, 
       # beta_init=beta_init, beta_increment=beta_increment, beta_num_iter=beta_num_iter, 
       damp=damp, mu=mu, sigma=sigma,k=k)
       # tolerance=tolerance, 
       # edgelist=edgelist)
}


#' @title doHMRF2
#' @name do HMRF
#' @description function to run HMRF model
#' @keywords external
#' @param xxxxxxxxxxxxxxxxxxxxx
#' @param xxxxxxxxxxxxxxxxxxxxx
#' @param xxxxxxxxxxxxxxxxxxxxx
doHMRF2 = function (HMRF_init_obj, betas = c(0,10,5))
  # y, nei, numnei, blocks, beta_init, beta_increment, beta_num_iter, damp, mu, sigma, k, tolerance) 
{
  if(!'y'%in%names(HMRF_init_obj))
  {stop('\n expression matrix \'y\' not in the intialization object \n')}
  if(!'nei'%in%names(HMRF_init_obj))
  {stop('\n neighbor matrix \'nei\' not in the intialization object \n')}
  if(!'numnei'%in%names(HMRF_init_obj))
  {stop('\n number of neighbors \'numnei\' not in the intialization object \n')}
  if(!'blocks'%in%names(HMRF_init_obj))
  {stop('\n iteration groups \'blocks\' not in the intialization object \n')}
  if(!'damp'%in%names(HMRF_init_obj))
  {stop('\n dampen factors \'damp\' not in the intialization object \n')}
  if(!'mu'%in%names(HMRF_init_obj))
  {stop('\n initial mean vector \'mu\' not in the intialization object \n')}
  if(!'sigma'%in%names(HMRF_init_obj))
  {stop('\n initial covariance matrix \'sigma\' not in the intialization object \n')}
  if(!'k'%in%names(HMRF_init_obj))
  {stop('\n cluster number \'k\' not in the intialization object \n')}
  
  y = HMRF_init_obj$y
  nei = HMRF_init_obj$nei
  numnei = HMRF_init_obj$numnei
  blocks = HMRF_init_obj$blocks
  damp = HMRF_init_obj$damp
  mu = HMRF_init_obj$mu
  sigma = HMRF_init_obj$sigma
  k = HMRF_init_obj$k
  
  if(length(betas)!=3||(sum(betas[1:3]<0)>0))
  {stop('\n please provide betas as a vector of 3 non-negative numbers (initial value, nicrement, total iteration number) \n')}
  beta_init = betas[1]
  beta_increment = betas[2]
  beta_num_iter = betas[3]
  # beta_current <- beta_init
  beta_seq = sequence(beta_num_iter,beta_init,beta_increment)
  beta_seq = sort(unique(c(0,beta_seq)))
  res <- c()
  for(beta_current in beta_seq){
    print(sprintf("Doing beta=%.3f", beta_current))
    tc.hmrfem<-smfishHmrf.hmrfem.multi(y=y, neighbors=nei, beta=beta_current, numnei=numnei, 
                                       blocks=blocks, mu=mu, sigma=sigma, verbose=T, err=1e-7, maxit=50, dampFactor=damp)
    #smfishHmrf.hmrfem.multi.save(name, outdir, beta_current, tc.hmrfem, k)
    #do_one(name, outdir, k, y, nei, beta_current, numnei, blocks, mu, sigma, damp)
    
    ### stop the loop if there is a samll maximum probablity (<0.51) of any gene
    if(sum(apply(tc.hmrfem$prob,1,max)<0.51)>0)
    {cat(paste0('\n HMRF is stopping at large beta >= ',beta_current,', numerical error occurs, results of smaller betas were stored\n'));
      break()}
    
    t_key <- sprintf("k=%d b=%.2f", k, beta_current)
    tc.hmrfem$sigma=NULL
    tc.hmrfem$mu=NULL
    rownames(tc.hmrfem$prob) = rownames(y)
    rownames(tc.hmrfem$unnormprob) = rownames(y)
    names(tc.hmrfem$class) = rownames(y)
    res[[t_key]] <- tc.hmrfem
    # beta_current <- beta_current + beta_increment
  }
  
  result.hmrf = res
  class(result.hmrf) <- append(class(result.hmrf),"HMRFoutput") 
  return(result.hmrf)
  
  #result.hmrf = smfishHmrf.hmrfem.multi(y = y, neighbors = nei, numnei = numnei, blocks = blocks, beta=beta, 
  #                                     mu=mu, sigma=sigma, err=1e-7, maxit=50, verbose=TRUE, dampFactor=damp, 
  #                                    tolerance=tolerance)
  
}



#' @title get domain index of hmrf result
#' @name GetDT.HMRFoutput
#' @description function to select Domain type of cells in HMRFoutput
#' @keywords external
GetDT = function(x){UseMethod('GetDT',x)}
GetDT.HMRFoutput = function(y){
  lapply(y,function(x){  cl = c(apply(x$prob,1,function(i){ifelse(test = sum(i>0.5)==1, yes = which(i>0.5), no = NA)}),
                                rep(NA,length(x$cell.rm)));
  names(cl) = c(rownames(x$prob),x$cell.rm);
  return(cl)})
}

#which(is.na(GetDT(result.hmrf)))


#' @title add HMRF DT to cell meta data
#' @name addHMRF2
#' @description function to add HMRF Domain Type to cell meta data
#' @keywords external
addHMRF2 = function (gobject, HMRFoutput
                     # , k = NULL, betas_to_add = NULL 
                     #, hmrf_name = NULL
) 
{
  if (!"HMRFoutput" %in% class(HMRFoutput)) {
    stop("\n HMRFoutput needs to be output from doHMRF extend \n")
  }
  # feat_type = HMRFoutput$feat_type
  # 
  # # get_result_path = system.file("python", "get_result2.py",
  # #                               package = "Giotto")
  # # name = HMRFoutput$name
  # # output_data = HMRFoutput$output_data
  # # python_path = HMRFoutput$python_path
  # # if (is.null(k)) {
  # #   stop("\n you need to select a k that was used with doHMRFextend \n")
  # # }
  # 
  # k = HMRFoutput$k
  # b = HMRFoutput$beta
  
  ordered_cell_IDs = gobject@cell_metadata$cell_ID
  
  # hmrf_DT = GetDT(HMRFoutput)
  # hmrf_DT = lapply(hmrf_DT,function(x){new.x = x[match(ordered_cell_IDs,names(x))];
  # names(new.x) = ordered_cell_IDs;return(new.x)})
  
  for(i in 1:length(HMRFoutput))
  {
    # if(names(hmrf_DT)[i]%in%names(gobject@cell_metadata))
    # {cat(paste0("\n cluster of ",names(hmrf_DT)[i]," exists, rewrite in meta data \n"))}
    gobject = addCellMetadata(gobject = gobject, 
                              # feat_type = feat_type, 
                              column_cell_ID = "cell_ID", new_metadata = HMRFoutputhmrf_DT[[i]]$class,
                              vector_name = names(HMRFoutput)[i], 
                              by_column = F)
  }
  
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



#' @title view HMRF result
#' @name viewHMRFresults2
#' @description function to plot spatial location with HMRF cluster of k and betas
#' @keywords external
viewHMRFresults2 = function (gobject, HMRFoutput, k, betas,
          third_dim = FALSE, ...) 
{
  if (!"HMRFoutput" %in% class(HMRFoutput)) {
    stop("\n HMRFoutput needs to be output from doHMRFextend \n")
  }
  beta_init = betas[1]
  beta_increment = betas[2]
  beta_num_iter = betas[3]
  # beta_current <- beta_init
  beta_seq = sequence(beta_num_iter,beta_init,beta_increment)
  # beta_seq = sort(unique(c(0,beta_seq)))
  
  t_key0 = sprintf("k=%d b=%.2f", k, 0)
  if(!t_key0%in%names(HMRFoutput))
  {stop(paste0('\n model of k = ',k,' was not calculated \n'))}
  
  dt_0 = HMRFoutput[[t_key0]]$class
  cat(paste0('\n plotting ',t_key0,' \n'))
  spatPlot2D(gobject = gobject, cell_color = dt_0, show_plot = T, 
             title = t_key0,...)
  if (third_dim == TRUE) {
    spatPlot3D(gobject = gobject, cell_color = dt_0, 
               show_plot = T, title = kk,...)
  }

  t_key = sprintf("k=%d b=%.2f", k, beta_seq) 
  if(sum(t_key%in%names(HMRFoutput))==0)
  {stop(paste0('\n model of k = ',k,' and betas = ',paste(beta_seq,collapse = ','),' was not calculated, only result of beta = 0 plotted \n'))}
  
  cat(paste0('\n plotting ',t_key[t_key%in%names(HMRFoutput)],' \n'))
  for(kk in intersect(t_key,names(HMRFoutput)))
  {
    dt_kk = HMRFoutput[[kk]]$class
    spatPlot2D(gobject = gobject, cell_color = dt_kk, show_plot = T, 
               title = kk,...)
    if (third_dim == TRUE) {
      spatPlot3D(gobject = gobject, cell_color = dt_kk, 
                 show_plot = T, title = kk,...)
    }
  }
  
}




