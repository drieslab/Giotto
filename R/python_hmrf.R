
sampling_sp_genes = function(clust, sample_rate=2, target=500, seed=10){
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
  set.seed(seed)
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

filterSpatialGenes <- function(gobject, spatial_genes, max=2500, name=c("binSpect", "silhouetteRank", "silhouetteRankTest"), method=c("none", "elbow")){
  name = match.arg(name, unique(c("binSpect", "silhouetteRank", "silhouetteRankTest", name)))
  method = match.arg(method, unique(c("none", "elbow", method)))

  gg = fDataDT(gobject)
  gx = gg[gene_ID %in% spatial_genes]

  if(name=="binSpect"){
    gx = gx[!is.na(binSpect.pval) & binSpect.pval<1]
    gx_sorted = gx[order(gx$binSpect.pval, decreasing=F),]
  }else if(name=="silhouetteRank"){
    gx = gx[!is.na(silhouetteRank.score)]
    gx_sorted = gx[order(gx$silhouetteRank.score, decreasing=T),]
  }else if(name=="silhouetteRankTest"){
    gx = gx[!is.na(silhouetteRankTest.pval) & silhouetteRankTest.pval<1]
    gx_sorted = gx[order(gx$silhouetteRankTest.pval, decreasing=F),]
  }

  if(method=="none"){
    gx = head(gx, n=max)
  }else if(method=="elbow"){
    
  }

  num_genes_removed = length(spatial_genes) - nrow(gx)

  return(list(genes=gx$gene_ID, num_genes_removed=num_genes_removed))
}


#' @title doHMRF
#' @name doHMRF
#' @description Run HMRF
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param spatial_network_name name of spatial network to use for HMRF
#' @param spatial_genes spatial genes to use for HMRF
#' @param spatial_dimensions select spatial dimensions to use, default is all possible dimensions
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name name of HMRF run
#' @param k  number of HMRF domains
#' @param seed seed to fix random number generator (for creating initialization of HMRF) (-1 if no fixing)
#' @param betas betas to test for. three numbers: start_beta, beta_increment, num_betas e.g. c(0, 2.0, 50)
#' @param tolerance tolerance
#' @param zscore zscore
#' @param numinit number of initializations
#' @param python_path python path to use
#' @param output_folder output folder to save results
#' @param overwrite_output overwrite output folder
#' @return Creates a directory with results that can be viewed with viewHMRFresults
#' @details Description of HMRF parameters ...
#' @export
initHMRF <- function(gobject,
                   expression_values = c('scaled', 'normalized', 'custom'),
                   spatial_network_name = 'Delaunay_network',
                   use_spatial_genes = c("binSpect", "silhouetteRank", "silhouetteRankTest"),
                   spatial_genes = NULL,
                   gene_samples = 500,
                   gene_sampling_rate = 2,
                   gene_sampling_seed = 10,
                   gene_sampling_from_top = 2500,
                   #gene_pval = "auto",
                   hmrf_seed = 100,
                   k = 10,
                   tolerance = 1e-5,
                   zscore = c('none','rowcol', 'colrow'),
                   nstart = 1000,
                   factor_step = 1.05) {


  if(!requireNamespace('smfishHmrf', quietly = TRUE)) {
    stop("\n package ", 'smfishHmrf' ," is not yet installed \n",
         "To install: \n",
         "devtools::install_bitbucket(repo = 'qzhudfci/smfishhmrf-r', ref='master')",
         "see http://spatial.rc.fas.harvard.edu/install.html for more information",
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

  if(spatial_genes !=NULL){
    cat(paste0("\n User supplied gene list detected.\n"))
    cat(paste0("\n Checking user gene list is spatial...\n"))
    if(!'binSpect.pval' %in% names(gobject@gene_metadata) ||
    !'silhouetteRank.score' %in% names(gobject@gene_metadata) ||
    !'silhouetteRankTest.pval' %in% names(gobject@gene_metadata)){
      stop(paste0("\n Cannot check user's gene list, because Giotto's spatial gene detection has not been run. Please run spatial gene detection first: binSpect, silhouetteRank.\n"), call.=FALSE)
    }
    filtered = filterSpatialGenes(gobject, spatial_genes, max=gene_sampling_from_top, name=use_spatial_genes, method="none")
    if(filtered$num_genes_removed>0){
      cat(paste0("\n Removed", filtered$num_genes_removed, "from user's input gene list due to being absent or non-spatial genes.\n"))
      cat(paste0("\n Kept", length(filtered$genes), "spatial genes for the sampling step next\n"))
    }
    spatial_genes = filtered$genes
    if(length(spatial_genes)==0){
      stop(paste0("\n No genes are remaining to do HMRF. Please give a larger gene list.\n"), call.=FALSE)
    }
  }

  if(!'binSpect.pval' %in% names(gobject@gene_metadata) ||
  !'silhouetteRank.score' %in% names(gobject@gene_metadata) ||
  !'silhouetteRankTest.pval' %in% names(gobject@gene_metadata)){
    stop(paste0('Giotto spatial gene detection has not been run. Please run spatial gene detection first: binSpect, silhouetteRank.\n'), call.=FALSE)
  }
  
  # if(!n_spatial_genes_select>0){stop("\n please provide a positive integer n_spatial_genes_select \n")}
  if(spatial_genes==NULL){
    all_genes = fDataDT(gobject)$gene_ID
    filtered = filterSpatialGenes(gobject, all_genes, max=gene_sampling_from_top, name=use_spatial_genes, method="none")
    spatial_genes = filtered$genes
  }

  n = min(spatial_gene_samples,500,length(spatial_genes))

  #if(length(intersect(spatial_genes,rownames(expr_values)))>0){
  #  ### if input a sp gene list, use the list
  #  expr_values = expr_values[intersect(spatial_genes,rownames(expr_values)),]
  #  n0 = length(spatial_genes)
  #  n1 = length(intersect(spatial_genes,rownames(expr_values)))
  #  n = min(spatial_gene_samples,500,n1)
  #  cat(paste0('\n total ',n0, 'spatial_genes imported \n'))
  #  cat(paste0('\n total ',n1, 'spatial_genes overlapped with expr data matrix \n'))
  #  cat(paste0('\n sample ',n, '(max of 500) spatial_genes selected for HMRF model \n'))
  #}else{
  #  
  #  ### if no input any sp gene, pick significant genes from km test
  #  cat('\n none overlapped genes imported, use spatial genes from kmeans_SpatialGenesTest \n')
  #  
  #  if(!'kmeans_SpatialGenesTest'%in%names(gobject@gene_metadata))
  #  { 
  #    km_spatialgenes = binSpect(gobject = gobject,bin_method = 'kmeans',set.seed = 1234,spatial_network_name = spatial_network_name)
  #    spatial_genes = km_spatialgenes$genes[which(km_spatialgenes$adj.p.value<0.01)]
  #  }else{
  #    spatial_genes = gobject@gene_metadata$gene_ID[gobject@gene_metadata$kmeans_SpatialGenesTest<0.01]
  #  }
  #  
  #  n = min(spatial_gene_samples,500,length(spatial_genes))
  #  cat(paste0('\n sample ',n,' (max of 500) genes out of total ',length(spatial_genes),' spatial genes for HMRF model \n'))
  #}
  
  if(n<length(spatial_genes)){
    spat_cor_netw_DT = detectSpatialCorGenes(gobject = gobject ,method = 'network',
                                             spatial_network_name = spatial_network_name,
                                             subset_genes = spatial_genes,
                                             network_smoothing = 0)
    spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT,name = 'spat_netw_clus',k = 20)
    sample_genes = sampling_sp_genes(spat_cor_netw_DT$cor_clusters$spat_netw_clus, sample_rate=spatial_gene_sampling_rate, target=n, seed=spatial_gene_sampling_seed)
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
  
  kk = smfishHmrf.generate.centroid(y=y,par_k = k,par_seed=hmrf_seed,nstart=nstart)
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
       damp=damp, mu=mu, sigma=sigma, k=k, genes=spatial_genes_selected, edgelist=edgelist)

}

#' @title doHMRF2
#' @name do HMRF
#' @description function to run HMRF model
#' @keywords external
#' @param xxxxxxxxxxxxxxxxxxxxx
#' @param xxxxxxxxxxxxxxxxxxxxx
#' @param xxxxxxxxxxxxxxxxxxxxx
doHMRF = function (HMRF_init_obj, betas = c(0,10,5))
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
  {stop('\n please provide betas as a vector of 3 non-negative numbers (initial beta, increment, number of betas) \n')}
  beta_init = betas[1]
  beta_increment = betas[2]
  beta_num_iter = betas[3]
  # beta_current <- beta_init
  #beta_seq = sequence(beta_num_iter,beta_init,beta_increment)
  #beta_seq = sort(unique(c(0,beta_seq)))
  res <- c()
  beta_current <- beta_init
  for(bx in 1:beta_num_iter){
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
    beta_current <- beta_current + beta_increment
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
addHMRF = function (gobject, HMRFoutput){
  if (!"HMRFoutput" %in% class(HMRFoutput)) {
    stop("\n HMRFoutput needs to be output from doHMRF extend \n")
  }
  
  ordered_cell_IDs = gobject@cell_metadata$cell_ID
  
  for(i in 1:length(HMRFoutput))
  {
    gobject = addCellMetadata(gobject = gobject, 
                              # feat_type = feat_type, 
                              column_cell_ID = "cell_ID", new_metadata = HMRFoutputhmrf_DT[[i]]$class,
                              vector_name = names(HMRFoutput)[i], 
                              by_column = F)
  }
  return(gobject)
}

#' @title view HMRF result
#' @name viewHMRFresults2
#' @description function to plot spatial location with HMRF cluster of k and betas
#' @keywords external
viewHMRFresults = function (gobject, HMRFoutput, k, betas,
          third_dim = FALSE, ...) {
  if (!"HMRFoutput" %in% class(HMRFoutput)) {
    stop("\n HMRFoutput needs to be output from doHMRFextend \n")
  }
  beta_init = betas[1]
  beta_increment = betas[2]
  beta_num_iter = betas[3]
  beta_current <- beta_init

  for(bx in 1:beta_num_iter){
    t_key <- sprintf("k=%d b=%.2f", k, beta_current)
    if(!t_key%in%names(HMRFoutput)){
      cat(paste0("\n", t_key, "was not calculated. Skipped.\n"))
    }else{
      cat(paste0('\n plotting ',t_key, '\n'))
      dt_kk = HMRFoutput[[t_key]]$class
      spatPlot2D(gobject = gobject, cell_color = dt_kk, show_plot = T, 
               title = t_key,...)
  
      if (third_dim == TRUE) {
        spatPlot3D(gobject = gobject, cell_color = dt_kk, 
                 show_plot = T, title = t_key,...)
      }
    }
    beta_current <- beta_current + beta_increment

  }

}

#' @title viewHMRFresults3D
#' @name viewHMRFresults3D
#' @description View results from doHMRF.
#' @param gobject giotto object
#' @param HMRFoutput HMRF output from doHMRF
#' @param k number of HMRF domains
#' @param betas_to_view results from different betas that you want to view
#' @param \dots additional parameters to spatPlot3D()
#' @return spatial plots with HMRF domains
#' @seealso \code{\link{spatPlot3D}}
#' @export
viewHMRFresults3D <- function(gobject,
                              HMRFoutput,
                              k = NULL,
                              betas_to_view = NULL,
                              ...) {


  if(!'HMRFoutput' %in% class(HMRFoutput)) {
    stop('\n HMRFoutput needs to be output from doHMRFextend \n')
  }

  ## reader.py and get_result.py paths
  # TODO: part of the package
  get_result_path = system.file("python", "get_result2.py", package = 'Giotto')

  # paths and name
  name = HMRFoutput$name
  output_data = HMRFoutput$output_data
  python_path = HMRFoutput$python_path

  # k-values
  if(is.null(k)) {
    stop('\n you need to select a k that was used with doHMRFextend \n')
  }
  k = HMRFoutput$k

  # betas
  betas = HMRFoutput$betas
  possible_betas = seq(betas[1], to = betas[1]+(betas[2]*(betas[3]-1)), by = betas[2])

  betas_to_view_detected = betas_to_view[betas_to_view %in% possible_betas]

  # plot betas
  for(b in betas_to_view_detected) {

    ## get results part ##
    result_command = paste0(python_path, ' ', get_result_path,
                            ' -r ', output_data,
                            ' -a ', name,
                            ' -k ', k,
                            ' -b ', b)

    print(result_command)

    output = system(command = result_command, intern = T)


    title_name = paste0('k = ', k, ' b = ',b)

    spatPlot3D(gobject = gobject, cell_color = output, show_plot = T, save_plot = F, title = title_name, ...)
  }
}
