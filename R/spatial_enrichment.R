

#' @title makeSignMatrixPAGE
#' @description Function to convert a list of signature genes (e.g. for cell types or processes) into
#' a binary matrix format that can be used with the PAGE enrichment option. Each cell type or process should
#' have a vector of cell-type or process specific genes. These vectors need to be combined into a list (sign_list).
#' The names of the cell types or processes that are provided in the list need to be given (sign_names).
#' @param sign_names vector with names for each provided gene signature
#' @param sign_list list of genes (signature)
#' @return matrix
#' @seealso \code{\link{PAGEEnrich}}
#' @export
#' @examples
#'     makeSignMatrixPAGE()
makeSignMatrixPAGE = function(sign_names,
                              sign_list) {

  ## check input
  if(!is.list(sign_list)) {
    stop('\n sign_list needs to be a list of signatures for each cell type / process \n')
  }
  if(length(sign_names) != length(sign_list)) {
    stop('\n the number of names needs to match the number of signatures provided \n')
  }

  ## create genes and signatures
  genes = do.call('c', sign_list)
  types = lapply(1:length(sign_names), FUN = function(x) {

    subset = sign_list[[x]]
    name_subset = sign_names[[x]]

    res = rep(x = name_subset, length(subset))

  })
  mydt = data.table::data.table(genes = genes, types = unlist(types), value = 1)

  # convert data.table to signature matrix
  dtmatrix = dcast.data.table(mydt, formula = genes~types, value.var = 'value', fill = 0)
  final_sig_matrix = as.matrix(dtmatrix[,-1]); rownames(final_sig_matrix) = dtmatrix$genes

  return(final_sig_matrix)

}




#' @title makeSignMatrixRank
#' @description Function to convert a single-cell count matrix
#' and a corresponding single-cell cluster vector into
#' a rank matrix that can be used with the Rank enrichment option.
#' @param sc_matrix matrix of single-cell RNAseq expression data
#' @param sc_cluster_ids vector of cluster ids
#' @param gobject if giotto object is given then only genes present in both datasets will be considered
#' @return matrix
#' @seealso \code{\link{rankEnrich}}
#' @export
#' @examples
#'     makeSignMatrixRank()
makeSignMatrixRank <- function(sc_matrix,
                               sc_cluster_ids,
                               gobject = NULL) {

  # check input
  if(length(sc_cluster_ids) != ncol(sc_matrix)) {
    stop('Number of columns of the scRNAseq matrix needs to have the same length as the cluster ids')
  }

  mean_list = list()
  group_list = list()
  total_nr_genes = nrow(sc_matrix)

  # calculate means for each cluster group
  for(group in 1:length(unique(sc_cluster_ids))) {

    group_id = unique(sc_cluster_ids)[group]
    cell_ind = which(sc_cluster_ids == group_id)
    cluster_rowmeans = rowMeans(sc_matrix[,cell_ind])
    mean_list[[group_id]] = cluster_rowmeans
    group_list[[group]] = rep(group_id, total_nr_genes)
  }

  mean_list_res = as.data.table(do.call('c', mean_list))
  group_list_res = as.data.table(do.call('c', group_list))

  # average expression for all cells
  av_expression = rowMeans(sc_matrix)
  av_expression_res = rep(av_expression, length(unique(sc_cluster_ids)))

  gene_names = rownames(sc_matrix)
  gene_names_res = rep(gene_names, length(unique(sc_cluster_ids)))

  # create data.table with genes, mean expression per cluster, mean expression overall and cluster ids
  comb_dt = data.table(genes = gene_names_res, mean_expr = mean_list_res[[1]], av_expr = av_expression_res, clusters = group_list_res[[1]])

  # calculate fold change and rank of fold-change
  comb_dt[, fold := log2(mean_expr+1)-log2(av_expr+1)]
  comb_dt[, rankFold := frank(-fold, ties.method = 'first'), by = clusters]

  # create matrix
  comb_rank_mat = data.table::dcast.data.table(data = comb_dt, genes~clusters, value.var = 'rankFold')
  comb_rank_matrix = Giotto:::dt_to_matrix(comb_rank_mat)
  comb_rank_matrix = comb_rank_matrix[rownames(sc_matrix), unique(sc_cluster_ids)]


  #  # intersect rank matrix with giotto object if given
  #  if(!is.null(gobject) & class(gobject) %in% c('giotto')) {
  #    comb_rank_matrix = comb_rank_matrix[intersect(rownames(comb_rank_matrix), gobject@gene_ID), ]
  #  }

  return(comb_rank_matrix)

}


#' @title rankPermutation
#' @description creates permutation for the rankEnrich test
#' @examples
#'     rankPermutation()
rankPermutation <- function(sc_gene, n){
  random_df <- data.frame(matrix(ncol = n, nrow = length(sc_gene)))
  for (i in 1:n){
    set.seed(i)
    random_rank<-sample(1:length(sc_gene),length(sc_gene),replace=FALSE)
    random_df[,i]<-random_rank
  }
  rownames(random_df)<-sc_gene
  return(random_df)
}



#' @title pagePermutation
#' @description creates permutation for the PAGEEnrich test
#' @examples
#'     pagePermutation()
pagePermutation <- function(sc_gene, gene_number, n){
  all_sample_names<-NULL
  all_sample_list<-NULL
  for (i in 1:n){
    set.seed(i)
    random_gene=sample(sc_gene,gene_number,replace=FALSE)
    ct_name<-paste("ct",i,sep="")
    all_sample_names<-c(all_sample_names,ct_name)
    all_sample_list<-c(all_sample_list,list(random_gene))
  }
  random_sig<-makeSignMatrixPAGE(all_sample_names,all_sample_list)
  return(random_sig)
}




#' @title PAGEEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using PAGE.
#' @param gobject Giotto object
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param output_enrichment how to return enrichment output
#' @return data.table with enrichment results
#' @details
#' sign_matrix: a binary matrix with genes as row names and cell-types as column names.
#' Alternatively a list of signature genes can be provided to makeSignMatrixPAGE, which will create
#' the matrix for you. \cr
#'
#' The enrichment Z score is calculated by using method (PAGE) from
#' Kim SY et al., BMC bioinformatics, 2005 as \eqn{Z = ((Sm – mu)*m^(1/2)) / delta}.
#' For each gene in each spot, mu is the fold change values versus the mean expression
#' and delta is the standard deviation. Sm is the mean fold change value of a specific marker gene set
#' and  m is the size of a given marker gene set.
#' @seealso \code{\link{makeSignMatrixPAGE}}
#' @export
#' @examples
#'     PAGEEnrich(gobject)
PAGEEnrich <- function(gobject,
                       sign_matrix,
                       expression_values = c('normalized', 'scaled', 'custom'),
                       reverse_log_scale = TRUE,
                       logbase = 2,
                       output_enrichment = c('original', 'zscore')) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # output enrichment
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))

  # only continue with genes present in both datasets
  interGene = intersect(rownames(sign_matrix), rownames(expr_values))
  filterSig = sign_matrix[interGene,]
  signames = rownames(filterSig)[which(filterSig[,1]==1)]

  # calculate mean gene expression
  if(reverse_log_scale == TRUE) {
    mean_gene_expr = log(rowMeans(logbase^expr_values-1, dims = 1)+1)
  } else {
    mean_gene_expr = rowMeans(expr_values)
  }
  geneFold = expr_values - mean_gene_expr

  # calculate sample/spot mean and sd
  cellColMean = apply(geneFold,2,mean)
  cellColSd = apply(geneFold,2,sd)

  # get enrichment scores
  enrichment = matrix(data=NA,nrow = dim(filterSig)[2],ncol=length(cellColMean))
  for (i in (1:dim(filterSig)[2])){
    signames = rownames(filterSig)[which(filterSig[,i]==1)]
    sigColMean = apply(geneFold[signames,],2,mean)
    m = length(signames)
    vectorX = NULL
    for (j in(1:length(cellColMean))){
      Sm = sigColMean[j]
      u = cellColMean[j]
      sigma = cellColSd[j]
      zscore = (Sm - u)* m^(1/2) / sigma
      vectorX = append(vectorX,zscore)
    }
    enrichment[i,] = vectorX
  }

  rownames(enrichment) = colnames(filterSig)
  colnames(enrichment) = names(cellColMean)
  enrichment = t(enrichment)

  if(output_enrichment == 'zscore') {
    enrichment = scale(enrichment)
  }

  enrichmentDT = data.table::data.table(cell_ID = rownames(enrichment))
  enrichmentDT = cbind(enrichmentDT, as.data.table(enrichment))

  return(enrichmentDT)
}


#' @title rankEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using a rank based approach.
#' @param gobject Giotto object
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param output_enrichment how to return enrichment output
#' @return data.table with enrichment results
#' @details
#' sign_matrix: a rank-fold matrix with genes as row names and cell-types as column names.
#' Alternatively a scRNA-seq matrix and vector with clusters can be provided to makeSignMatrixRank, which will create
#' the matrix for you. \cr
#'
#' First a new rank is calculated as R = (R1*R2)^(1/2), where R1 is the rank of
#' fold-change for each gene in each spot and R2 is the rank of each marker in each cell type.
#' The Rank-Biased Precision is then calculated as: RBP = (1 - 0.99) * (0.99)^(R - 1)
#' and the final enrichment score is then calculated as the sum of top 100 RBPs.
#' @seealso \code{\link{makeSignMatrixRank}}
#' @export
#' @examples
#'     rankEnrich(gobject)
rankEnrich <- function(gobject,
                       sign_matrix,
                       expression_values = c('normalized', 'scaled', 'custom'),
                       reverse_log_scale = TRUE,
                       logbase = 2,
                       output_enrichment = c('original', 'zscore')) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # output enrichment
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))


  enrichment = matrix(data = NA,
                      nrow = dim(sign_matrix)[2],
                      ncol = dim(expr_values)[2])

  # calculate mean gene expression
  if(reverse_log_scale == TRUE) {
    mean_gene_expr = log(rowMeans(logbase^expr_values-1, dims = 1)+1)
  } else {
    mean_gene_expr = rowMeans(expr_values)
  }

  # fold change and ranking
  geneFold = expr_values - mean_gene_expr
  rankFold = t(matrixStats::colRanks(-geneFold, ties.method = "first"))


  rownames(rankFold) = rownames(expr_values)
  colnames(rankFold) = colnames(expr_values)

  for (i in (1:dim(sign_matrix)[2])){

    signames = rownames(sign_matrix)[which(sign_matrix[,i]>0)]
    interGene = intersect(signames, rownames(rankFold))
    filterSig = sign_matrix[interGene,]
    filterRankFold = rankFold[interGene,]

    multiplyRank = (filterRankFold*filterSig[,i])^(1/2)
    rpb = 0.01*(0.99^(multiplyRank-1))
    vectorX = NULL

    for (j in (1:dim(filterRankFold)[2])){
      toprpb = sort(rpb[,j],decreasing = T)
      zscore = sum(toprpb[1:100])
      vectorX = IRanges::append(vectorX, zscore)
    }
    enrichment[i,] = vectorX
  }

  rownames(enrichment) = colnames(sign_matrix)
  colnames(enrichment) = colnames(rankFold)

  enrichment = t(enrichment)

  if(output_enrichment == 'zscore') {
    enrichment = scale(enrichment)
  }

  enrichmentDT = data.table::data.table(cell_ID = rownames(enrichment))
  enrichmentDT = cbind(enrichmentDT, as.data.table(enrichment))

  return(enrichmentDT)

}


#' @title hyperGeometricEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using a hypergeometric test.
#' @param gobject Giotto object
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param top_percentage percentage of cells that will be considered to have gene expression with matrix binarization
#' @param output_enrichment how to return enrichment output
#' @return data.table with enrichment results
#' @details The enrichment score is calculated based on the p-value from the
#' hypergeometric test, -log10(p-value).
#' @export
#' @examples
#'     hyperGeometricEnrich(gobject)
hyperGeometricEnrich <- function(gobject,
                                 sign_matrix,
                                 expression_values = c('normalized', 'scaled', 'custom'),
                                 reverse_log_scale = TRUE,
                                 logbase = 2,
                                 top_percentage = 5,
                                 output_enrichment = c('original', 'zscore')) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # output enrichment
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))

  # calculate mean gene expression
  if(reverse_log_scale == TRUE) {
    expr_values = logbase^expr_values-1
  }

  interGene<-intersect(rownames(gobject@norm_expr),rownames(sign_matrix))
  inter_sign_matrix<-sign_matrix[interGene,]
  aveExp<-log2(2*(rowMeans(2^gobject@norm_expr-1, dims = 1))+1)
  foldChange<-gobject@norm_expr-aveExp
  top_q<-1-top_percentage/100
  quantilecut<-apply(foldChange, 2 , quantile , probs = top_q, na.rm = TRUE )
  expbinary<-t(ifelse(t(foldChange)>quantilecut,1,0))

  markerGenes = rownames(inter_sign_matrix)
  expbinaryOverlap = expbinary[markerGenes,]
  total = length(markerGenes)
  enrichment = matrix(data=NA,
                      nrow = dim(inter_sign_matrix)[2],
                      ncol=dim(expbinaryOverlap)[2])

  for (i in (1:dim(inter_sign_matrix)[2])){
    signames = rownames(inter_sign_matrix)[which(inter_sign_matrix[,i]==1)]
    vectorX  = NULL

    for (j in(1:dim(expbinaryOverlap)[2])){

      cellsiggene = names(expbinaryOverlap[which(expbinaryOverlap[,j]==1),j])
      x = length(intersect(cellsiggene,signames))
      m = length(rownames(inter_sign_matrix)[which(inter_sign_matrix[,i]==1)])
      n = total-m
      k = length(intersect(cellsiggene, markerGenes))
      enrich<-(0-log10(phyper(x, m, n, k, log = FALSE,lower.tail = FALSE)))
      vectorX = append(vectorX,enrich)
    }
    enrichment[i,] = vectorX
  }

  rownames(enrichment) = colnames(inter_sign_matrix)
  colnames(enrichment) = colnames(expbinaryOverlap)

  enrichment = t(enrichment)

  if(output_enrichment == 'zscore') {
    enrichment = scale(enrichment)
  }

  enrichmentDT = data.table::data.table(cell_ID = rownames(enrichment))
  enrichmentDT = cbind(enrichmentDT, as.data.table(enrichment))

  return(enrichmentDT)
}



#' @title createSpatialEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using a hypergeometric test.
#' @param gobject Giotto object
#' @param enrich_method method for gene signature enrichment calculation
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param p_value calculate p-value (default = TRUE)
#' @param n_times (page/rank) number of permutation iterations to calculate p-value
#' @param n_genes (page/rank) number of randomly selected genes for each permuation
#' @param top_percentage (hyper) percentage of cells that will be considered to have gene expression with matrix binarization
#' @param output_enrichment how to return enrichment output
#' @param name to give to spatial enrichment results, default = PAGE
#' @param return_gobject return giotto object
#' @return Giotto object or enrichment results if return_gobject = FALSE
#' @details For details see the individual functions:
#' \itemize{
#'   \item{PAGE: }{\code{\link{PAGEEnrich}}}
#'   \item{PAGE: }{\code{\link{rankEnrich}}}
#'   \item{PAGE: }{\code{\link{hyperGeometricEnrich}}}
#' }
#'
#'
#' @export
#' @examples
#'     createSpatialEnrich(gobject)
createSpatialEnrich = function(gobject,
                               enrich_method = c('PAGE', 'rank', 'hypergeometric'),
                               sign_matrix,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               reverse_log_scale = TRUE,
                               logbase = 2,
                               p_value = TRUE,
                               n_genes = 100,
                               n_times = 1000,
                               top_percentage = 5,
                               output_enrichment = c('original', 'zscore'),
                               name = 'PAGE',
                               return_gobject = TRUE) {


  enrich_method = match.arg(enrich_method, choices = c('PAGE', 'rank', 'hypergeometric'))
  expression_values = match.arg(expression_values, choices = c('normalized', 'scaled', 'custom'))
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))

  if(enrich_method == 'PAGE') {

    enrich_results_DT =  PAGEEnrich(gobject = gobject,
                                    sign_matrix = sign_matrix,
                                    expression_values = expression_values,
                                    reverse_log_scale = reverse_log_scale,
                                    logbase = logbase,
                                    output_enrichment = output_enrichment)
    # default name for page enrichment
    if(is.null(name)) name = 'PAGE'
    if (p_value==TRUE){
        random_sig = pagePermutation(rownames(gobject@norm_expr),n_genes,n_times)
        random_DT = PAGEEnrich(gobject,
                               sign_matrix = random_sig,
                               expression_values = expression_values,
                               reverse_log_scale = reverse_log_scale,
                               logbase = logbase,
                               output_enrichment = output_enrichment)
        background<-unlist(random_DT[,2:dim(random_DT)[2]])
        pvalue_DT<-enrich_results_DT
        pvalue_DT[,2:dim(pvalue_DT)[2]]<- lapply(pvalue_DT[,2:dim(pvalue_DT)[2]], function(x)
          {pnorm(x, mean = mean(background), sd = sd(background), lower.tail = FALSE,log.p = FALSE)})
    }

  } else if(enrich_method == 'rank') {

    enrich_results_DT =  rankEnrich(gobject = gobject,
                                    sign_matrix = sign_matrix,
                                    expression_values = expression_values,
                                    reverse_log_scale = reverse_log_scale,
                                    logbase = logbase,
                                    output_enrichment = output_enrichment)
    # default name for page enrichment
    if(is.null(name)) name = 'rank'
    if (p_value==TRUE){
        random_rank = rankPermutation(rownames(sign_matrix),n_times)
        random_DT = rankEnrich(gobject,
                               sign_matrix = random_rank,
                               expression_values = expression_values,
                               reverse_log_scale = reverse_log_scale,
                               logbase = logbase,
                               output_enrichment = output_enrichment)
        background<-unlist(random_DT[,2:dim(random_DT)[2]])
        fit.gamma <- fitdist(background, distr = "gamma", method = "mle")
        pvalue_DT<-enrich_results_DT
        pvalue_DT[,2:dim(pvalue_DT)[2]]<- lapply(pvalue_DT[,2:dim(pvalue_DT)[2]],function(x)
          {pgamma(x, fit.gamma$estimate[1], rate = fit.gamma$estimate[2], lower.tail = FALSE,log.p = FALSE)})
    }


  } else if(enrich_method == 'hypergeometric'){

    enrich_results_DT =  hyperGeometricEnrich(gobject = gobject,
                                              sign_matrix = sign_matrix,
                                              expression_values = expression_values,
                                              reverse_log_scale = reverse_log_scale,
                                              logbase = logbase,
                                              top_percentage = top_percentage,
                                              output_enrichment = output_enrichment)
    # default name for page enrichment
    if(is.null(name)) name = 'hypergeometric'
    if (p_value==TRUE){
        pvalue_DT<-enrich_results_DT
        pvalue_DT[,2:dim(pvalue_DT)[2]]<- lapply(pvalue_DT[,2:dim(pvalue_DT)[2]],function(x){10^(-x)})
    }
  }

  if(return_gobject == TRUE) {

    spenr_names = names(gobject@spatial_enrichment[[name]])

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_spatial_enrichment')

    # parameters to include
    parameters_list[[update_name]] = c('method used' = enrich_method,
                                       'expression values' = expression_values,
                                       'reverse log scale' = reverse_log_scale,
                                       'logbase' = logbase,
                                       'p-values calculated' = p_value,
                                       'output enrichment scores' = output_enrichment)
    gobject@parameters = parameters_list

    gobject@spatial_enrichment[[name]] = enrich_results_DT
    if (p_value==TRUE){
        pvalue_name=paste(name,"pvalue",sep="_")
        gobject@spatial_enrichment[[pvalue_name]] = pvalue_DT
    }

    return(gobject)

  } else {
    return(enrich_results_DT)
  }
}



