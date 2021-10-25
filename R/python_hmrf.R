
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

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}

filterSpatialGenes <- function(gobject, spatial_genes, max=2500, name=c("binSpect", "silhouetteRank", "silhouetteRankTest"), method=c("none", "elbow")){
  name = match.arg(name, unique(c("binSpect", "silhouetteRank", "silhouetteRankTest", name)))
  method = match.arg(method, unique(c("none", "elbow", method)))

  #gg = fDataDT(gobject)
  #gx = gg[gene_ID %in% spatial_genes]
  #first determine how many spatial genes in this dataset
  gx = fDataDT(gobject)


  if(name=="binSpect"){
    gx = gx[!is.na(binSpect.pval) & binSpect.pval<1]
    gx_sorted = gx[order(gx$binSpect.pval, decreasing=F),]
  }else if(name=="silhouetteRank"){
    gx = gx[!is.na(silhouetteRank.score) & silhouetteRank.score>0]
    gx_sorted = gx[order(gx$silhouetteRank.score, decreasing=T),]
  }else if(name=="silhouetteRankTest"){
    gx = gx[!is.na(silhouetteRankTest.pval) & silhouetteRankTest.pval<1]
    gx_sorted = gx[order(gx$silhouetteRankTest.pval, decreasing=F),]
  }

  #print(gx_sorted)
  if(method=="none"){
    if(name=="binSpect"){
      gx_sorted = gx_sorted[binSpect.pval<0.01]
    }else if(name=="silhouetteRankTest"){
      gx_sorted = gx_sorted[silhouetteRankTest.pval<0.01]
    }
    gx_sorted = head(gx_sorted, n=max)

  }else if(method=="elbow"){
    y0 = c()
    if(name=="binSpect"){
      y0 = -log10(gx_sorted$binSpect.pval)
    }else if(name=="silhouetteRankTest"){
      y0 = -log10(gx_sorted$silhouetteRankTest.pval)
    }else if(name=="silhouetteRank"){
      y0 = gx_sorted$silhouetteRank.score
    }
    x0 = seq(1, nrow(gx_sorted))
    
    y0s<-sort(y0)
    y0s[y0s<0]<-0 #strictly positive
    #plot(x0, y0)
	  slope <- (max(y0s)-min(y0s))/length(y0s) #This is the slope of the line we want to slide. This is the diagonal.
	  #cat(paste0("slope is ", slope, ".\n"))
	  #tt<-optimize(numPts_below_line,lower=1,upper=length(y0s),myVector=y0s,slope=slope)
	  #print(tt)
	  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(y0s),myVector=y0s,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	  xPt <- length(y0s) - xPt
	  y_cutoff <- y0[xPt] #The y-value at this x point. This is our y_cutoff.
    gx_sorted = head(gx_sorted, n=xPt)
    cat(paste0("\nElbow method chosen to determine number of spatial genes.\n"))
    cat(paste0("\nElbow point determined to be at x=", xPt, " genes", " y=", y_cutoff, "\n"))
  }

  #filter user's gene list (spatial_genes)
  gx_sorted = gx_sorted[gene_ID %in% spatial_genes]

  num_genes_removed = length(spatial_genes) - nrow(gx_sorted)

  return(list(genes=gx_sorted$gene_ID, num_genes_removed=num_genes_removed))
}

chooseAvailableSpatialGenes <- function(gobject){
  eval1 = 'binSpect.pval' %in% names(gobject@gene_metadata)
  eval2 = 'silhouetteRankTest.pval' %in% names(gobject@gene_metadata)
  eval3 = 'silhouetteRank.score' %in% names(gobject@gene_metadata)
  if(eval1==TRUE){
    return("binSpect")
  }else if(eval2==TRUE){
    return("silhouetteRankTest")
  }else if(eval3==TRUE){
    return("silhouetteRank")
  }else{
    stop(paste0("No available spatial genes. Please run binSpect or silhouetteRank\n"), call.=FALSE)
  }
}

checkAndFixSpatialGenes <- function(gobject, use_spatial_genes, use_score=FALSE){
  if(use_spatial_genes=="silhouetteRank"){
    if(use_score==TRUE){
      use_spatial_genes = "silhouetteRank"
    }else{
      eval1 = 'silhouetteRank.score' %in% names(gobject@gene_metadata)
      eval2 = 'silhouetteRankTest.pval' %in% names(gobject@gene_metadata)
      if(eval1==TRUE && eval2==TRUE){
        #if both evaluate to true, then decide by use_score. 
        #silhouetteRank works only with score, silhouetteRankTest works only with pval
        if(use_score==TRUE){
          use_spatial_genes = "silhouetteRank"
        }else{
          use_spatial_genes = "silhouetteRankTest"
        }
      }else if(eval1==TRUE){
        use_spatial_genes = "silhouetteRank"
      }else if(eval2==TRUE){
        use_spatial_genes = "silhouetteRankTest"
      }else{
        stop(paste0("\n use_spatial_genes is set to silhouetteRank, but it has not been run yet. Run silhouetteRank first.\n"), call.=FALSE)
      }
    }
    return(use_spatial_genes)
  }
  else if(use_spatial_genes=="binSpect"){
    eval1 = 'binSpect.pval' %in% names(gobject@gene_metadata)
    if(eval1==FALSE){
      stop(paste0("\n use_spatial_genes is set to binSpect, but it has not been run yet. Run binSpect first.\n"), call.=FALSE)
    }
    return(use_spatial_genes)
  }else{
    stop(paste0("\n use_spatial_genes is set to one that is not supported.\n"), call.=FALSE)
  }
}


#' @title initHMRF
#' @name initHMRF
#' @description Initialize HMRF
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param spatial_network_name name of spatial network to use for HMRF
#' @param use_spatial_genes which of Giotto's spatial genes to use
#' @param gene_samples how many spatial genes to sample
#' @param gene_sampling_rate sampling rate parameter (1-50)
#' @param gene_sampling_seed sampling random generator seed
#' @param gene_sampling_from_top how many top spatial genes to perform sampling
#' @param user_gene_list user-specified genes (optional)
#' @param filter_method filter genes by top or by elbow method, prior to sampling
#' @param use_score use score as gene selection criterion (applies when use_spatial_genes=silhouetteRank)
#' @param hmrf_seed random generator seed for initialization of HMRF
#' @param k number of spatial clusters
#' @param tolerance error tolerance threshold
#' @param nstart number of Kmeans initializations from which to select the best initialization
#' @param factor_step dampened factor step
#' @return A list (see details)
#' @details There are two steps to running HMRF. This is the first step, the initialization. 
#' First, user specify which of Giotto's spatial genes to run, through use_spatial_genes. 
#' Spatial genes have been stored in the gene metadata table. A first pass of genes will filter genes 
#' that are not significantly spatial, as determined by filter_method. If filter_method is none, 
#' then top 2500 (gene_sampling_from_top) genes ranked by pvalue are considered spatial. 
#' If filter_method is elbow, then the exact cutoff is determined by the elbow in the 
#' -log10Pvalue vs. gene rank plot. 
#' Second, the filtered gene set is subject to sampling to select 500 
#' (controlled by gene_samples) genes for running HMRF. 
#' Third, once spatial genes are finalized, we are ready to initialize HMRF. 
#' This consists of running a K-means algorithm to determine initial centroids (nstart, hmrf_seed) of HMRF. 
#' The initialization is then finished. 
#' This function returns a list containing y (expression), nei (neighborhood structure), numnei (number of neighbors), blocks (graph colors), damp (dampened factor), mu (mean), sigma (covariance), k, genes, edgelist. This information is needed for the second step, doHMRF.
#'
#' @export
initHMRF <- function(gobject,
                   expression_values = c('scaled', 'normalized', 'custom'),
                   spatial_network_name = 'Delaunay_network',
                   use_spatial_genes = c("binSpect", "silhouetteRank"),
                   gene_samples = 500,
                   gene_sampling_rate = 2,
                   gene_sampling_seed = 10,
                   gene_sampling_from_top = 2500,
                   filter_method = c("none", "elbow"),
                   user_gene_list = NULL,
                   use_score = FALSE,
                   #gene_pval = "auto",
                   hmrf_seed = 100,
                   k = 10,
                   tolerance = 1e-5,
                   zscore = c('none','rowcol', 'colrow'),
                   nstart = 1000,
                   factor_step = 1.05) {

  message("\nIf used in published research, please cite:
  Q Zhu, S Shah, R Dries, L Cai, GC Yuan. 'Identification of spatially associated subpopulations by combining 
  scRNAseq and sequential fluorescence in situ hybridization data' Nature biotechnology 36 (12), 1183-1190. 2018\n")


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

  zscore = match.arg(zscore, unique(c('none','rowcol', 'colrow', zscore)))
  use_spatial_genes = match.arg(use_spatial_genes, 
    unique(c("binSpect", "silhouetteRank", use_spatial_genes)))
  filter_method = match.arg(filter_method, unique(c("none", "elbow", filter_method)))

  spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, 
                                          return_network_Obj = FALSE)
  spatial_network = spatial_network[, .(to, from)]
  
  # 'spatial_genes' use spatial gene list or create spatial gene list
  values = match.arg(expression_values, unique(c("scaled","normalized", "custom", expression_values)))
  expr_values = select_expression_values(gobject = gobject, values = values)   
  
  if(zscore!="none"){
    zscore = match.arg(zscore, c("none", "colrow", "rowcol"))
    #zscore = match.arg(zscore, c("col", "row", "colrow", "rowcol"))
    expr_values = select_expression_values(gobject = gobject, values = 'normalized')  
    #if(zscore=='col'){expr_values = scale(expr_values)}
    #if(zscore=='row'){expr_values = t(scale(t(expr_values)))}
    if(zscore=='colrow'){expr_values = t(scale(t(scale(expr_values))))}
    if(zscore=='rowcol'){expr_values = scale(t(scale(t(expr_values))))}
  }

  spatial_genes = c()

  if(! is.null(user_gene_list)){
    cat(paste0("\n User supplied gene list detected.\n"))
    cat(paste0("\n Checking user gene list is spatial...\n"))
    if(!'binSpect.pval' %in% names(gobject@gene_metadata) &&
    !'silhouetteRank.score' %in% names(gobject@gene_metadata) &&
    !'silhouetteRankTest.pval' %in% names(gobject@gene_metadata)){
      stop(paste0("\n Cannot check user's gene list, because Giotto's spatial gene detection has not been run. Please run spatial gene detection first: binSpect, silhouetteRank.\n"), call.=FALSE)
    }

    use_spatial_genes = chooseAvailableSpatialGenes(gobject)

    filtered = filterSpatialGenes(gobject, user_gene_list, max=gene_sampling_from_top, name=use_spatial_genes, method=filter_method)
    if(filtered$num_genes_removed>0){
      cat(paste0("\n Removed ", filtered$num_genes_removed, " from user's input gene list due to being absent or non-spatial genes.\n"))
      cat(paste0("\n Kept ", length(filtered$genes), " spatial genes for the sampling step next\n"))
    }
    spatial_genes = filtered$genes
    if(length(spatial_genes)==0){
      stop(paste0("\n No genes are remaining to do HMRF. Please give a larger gene list.\n"), call.=FALSE)
    }
  }

  if(!'binSpect.pval' %in% names(gobject@gene_metadata) &&
  !'silhouetteRank.score' %in% names(gobject@gene_metadata) &&
  !'silhouetteRankTest.pval' %in% names(gobject@gene_metadata)){
    stop(paste0('Giotto spatial gene detection has not been run. Please run spatial gene detection first: binSpect, silhouetteRank.\n'), call.=FALSE)
  }
  
  # if(!n_spatial_genes_select>0){stop("\n please provide a positive integer n_spatial_genes_select \n")}
  if(is.null(user_gene_list)){
    cat(paste0("\n Choosing spatial genes from the results of ", use_spatial_genes, "\n"))
    use_spatial_genes = checkAndFixSpatialGenes(gobject, use_spatial_genes, use_score=use_score)
    all_genes = fDataDT(gobject)$gene_ID
    filtered = filterSpatialGenes(gobject, all_genes, max=gene_sampling_from_top, name=use_spatial_genes, method=filter_method)
    cat(paste0("\n Kept ", length(filtered$genes), " top spatial genes for the sampling step next\n"))
    spatial_genes = filtered$genes
  }

  n = min(gene_samples,500,length(spatial_genes))
  
  if(n<length(spatial_genes)){
    cat(paste0("\n Computing spatial coexpression modules...\n"))
    spat_cor_netw_DT = detectSpatialCorGenes(gobject = gobject ,method = 'network',
                                             spatial_network_name = spatial_network_name,
                                             subset_genes = spatial_genes,
                                             network_smoothing = 0)
    spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT,name = 'spat_netw_clus',k = 20)
    cat(paste0("\n Sampling spatial genes from coexpression modules...\n"))
    sample_genes = sampling_sp_genes(spat_cor_netw_DT$cor_clusters$spat_netw_clus, sample_rate=gene_sampling_rate, target=n, seed=gene_sampling_seed)
    spatial_genes_selected = sample_genes$union_genes
    cat(paste0("\n Sampled ", length(spatial_genes_selected), " genes.\n"))
    
  }else{spatial_genes_selected = spatial_genes}

  cat(paste0("\n Will use ", length(spatial_genes_selected), " genes for init of HMRF.\n"))

  expr_values = expr_values[spatial_genes_selected,]
  
  y0 = t(Matrix::as.matrix(expr_values))
  
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
  
  cat(paste0("\n Parsing neighborhood graph...\n"))

  pp<-tbl_graph(edges=as.data.frame(edgelist), directed=F)
  yy<-pp%>%mutate(color=as.factor(color_dsatur()))
  colors<-as.list(yy)$nodes$color
  cl_color <- sort(unique(colors))
  blocks<-lapply(cl_color, function(cl){which(colors==cl)})

  cat(paste0("\n Kmeans initialization...\n"))

  kk = smfishHmrf.generate.centroid(y=y,par_k = k,par_seed=hmrf_seed,nstart=nstart)
  mu<-t(kk$centers) #should be dimension (m,k)
  lclust<-lapply(1:k, function(x) which(kk$cluster == x))
  damp<-array(0, c(k)); 
  sigma<-array(0, c(m,m,k))
  for(i in 1:k){
    sigma[, , i] <- cov(y[lclust[[i]], ])
    di<-findDampFactor(sigma[,,i], factor=factor_step, d_cutoff=tolerance, startValue=0.0001)
    damp[i]<-ifelse(is.null(di), 0, di)
  }

  cat(paste0("\n Done\n"))

  list(y=y, nei=nei, numnei=numnei, blocks=blocks, 
       damp=damp, mu=mu, sigma=sigma, k=k, genes=spatial_genes_selected, edgelist=edgelist)

}

#' @title doHMRF
#' @name do HMRF
#' @description function to run HMRF model
#' @keywords external
#' @param HMRF_init_obj return list of initHMRF() function
#' @param betas a vector of three values: initial beta, beta increment, and number of betas
#' @return A list (see details)
#' @details This function will run a HMRF model after initialization of HMRF.
#' Of note is the beta parameter, the smoothing parameter. We recommend running a range of betas,
#` hence betas specify what this range is.
#' For example, betas=c(0,10,5) will run for the following betas: 0, 10, 20, 30, 40.
#' betas=c(0,5,2) will run for betas: 0, 5, 10
#' Setting the beta can use the following guideline.
#' If number of genes N is 10<N<50, set betas=c(0, 1, 20)
#' If 50<N<100, set betas=c(0, 2, 25)
#' If 100<N<500, set betas=c(0, 5, 20)
#' Returns a list with class, probability, and model log-likelihood value.
doHMRF = function (HMRF_init_obj, betas = c(0,10,5))
  # y, nei, numnei, blocks, beta_init, beta_increment, beta_num_iter, damp, mu, sigma, k, tolerance) 
{
  message("\nIf used in published research, please cite:
  Q Zhu, S Shah, R Dries, L Cai, GC Yuan. 'Identification of spatially associated subpopulations by combining 
  scRNAseq and sequential fluorescence in situ hybridization data' Nature biotechnology 36 (12), 1183-1190. 2018\n")

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
  beta_seq = (1:beta_num_iter-1)*beta_increment+beta_init
  beta_seq = sort(unique(c(0,beta_seq)))
  
  res <- c()
  #beta_current <- beta_init
  for(beta_current in beta_seq){
    print(sprintf("Doing beta=%.3f", beta_current))
    tc.hmrfem<-smfishHmrf.hmrfem.multi(y=y, neighbors=nei, beta=beta_current, numnei=numnei, 
                                       blocks=blocks, mu=mu, sigma=sigma, verbose=T, err=1e-7, maxit=50, dampFactor=damp)
    #smfishHmrf.hmrfem.multi.save(name, outdir, beta_current, tc.hmrfem, k)
    #do_one(name, outdir, k, y, nei, beta_current, numnei, blocks, mu, sigma, damp)
    
    ### stop the loop if there is a samll maximum probablity (<0.51) of any gene
    if(sum(apply(tc.hmrfem$prob,1,max)<(1/k+0.05))>0)
    {cat(paste0('\n HMRF is stopping at large beta >= ',beta_current,', numerical error occurs, results of smaller betas were stored\n'));
      break()}
    
    t_key <- sprintf("k=%d b=%.2f", k, beta_current)
    tc.hmrfem$sigma=NULL
    tc.hmrfem$mu=NULL
    rownames(tc.hmrfem$prob) = rownames(y)
    rownames(tc.hmrfem$unnormprob) = rownames(y)
    names(tc.hmrfem$class) = rownames(y)
    res[[t_key]] <- tc.hmrfem
    #beta_current <- beta_current + beta_increment
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
                              column_cell_ID = "cell_ID", new_metadata = HMRFoutput[[i]]$class,
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
  #beta_current <- beta_init
  beta_seq = (1:beta_num_iter-1)*beta_increment+beta_init
  
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

  t_keys = sprintf("k=%d b=%.2f", k, beta_seq) 
  if(sum(t_keys%in%names(HMRFoutput))==0)
  {stop(paste0('\n model of k = ',k,' and betas = ',paste(beta_seq,collapse = ','),' was not calculated, only result of beta = 0 plotted \n'))}


  for(kk in intersect(t_keys, names(HMRFoutput)))
  {
    cat(paste0('\n plotting ', kk, '\n'))
    dt_kk = HMRFoutput[[kk]]$class
    spatPlot2D(gobject = gobject, cell_color = dt_kk, show_plot = T, 
               title = kk,...)
  
    if (third_dim == TRUE) {
      spatPlot3D(gobject = gobject, cell_color = dt_kk, 
                 show_plot = T, title = kk,...)
    }
   
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
