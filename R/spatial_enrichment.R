
## create spatial enrichment matrix ####

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
  dtmatrix = data.table::dcast.data.table(mydt, formula = genes~types, value.var = 'value', fill = 0)
  final_sig_matrix = Matrix::as.matrix(dtmatrix[,-1]); rownames(final_sig_matrix) = dtmatrix$genes

  return(final_sig_matrix)

}
## create spatialDWLS matrix ####

#' @title makeSignMatrixDWLS
#' @description Function to convert signature genes into an average expression matrix format that can be used for
#' spatialDWLS. gobject is the single cell Giotto object. sign_gene is in character format contains all DEGs.
#' cell_type are all the cell types used for single cell analysis.
#' @param gobject Giotto object of single cell
#' @param sign_gene all of DE genes (signature)
#' @param cell_type cell type for single cells
#' @return matrix
#' @seealso \code{\link{spatialDWLS}}
#' @export
makeSignMatrixDWLS = function(gobject,
                              sign_gene,
                              cell_type) {
  ## check input
  if(!is.character(sign_gene)) {
    stop('\n sign_gene needs to be a character of signatures for all cell types / process \n')
  }
  ## get un logged normalized expression value
  norm_exp<- 2^(gobject@norm_expr)-1
  id<- as.character(cell_type)
  intersect_sign_gene <- intersect(rownames(norm_exp), sign_gene)
  ExprSubset<-norm_exp[intersect_sign_gene,]
  sig_exp<-NULL
  for (i in unique(id)){
    sig_exp<-cbind(sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(sig_exp)<-unique(id)
  
  return(as.matrix(sig_exp))
}

#' @title makeSignMatrixRank
#' @description Function to convert a single-cell count matrix
#' and a corresponding single-cell cluster vector into
#' a rank matrix that can be used with the Rank enrichment option.
#' @param sc_matrix matrix of single-cell RNAseq expression data
#' @param sc_cluster_ids vector of cluster ids
#' @param ties_method how to handle rank ties
#' @param gobject if giotto object is given then only genes present in both datasets will be considered
#' @return matrix
#' @seealso \code{\link{rankEnrich}}
#' @export
makeSignMatrixRank <- function(sc_matrix,
                               sc_cluster_ids,
                               ties_method = c("random", "max"),
                               gobject = NULL) {

  if(methods::is(sc_matrix, "sparseMatrix")){
    sc_matrix = Matrix::as.matrix(sc_matrix)
  }

  # select ties_method
  ties_method =  match.arg(ties_method, choices = c("random", "max"))

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

  mean_list_res = data.table::as.data.table(do.call('c', mean_list))
  group_list_res = data.table::as.data.table(do.call('c', group_list))

  # average expression for all cells
  av_expression = rowMeans(sc_matrix)
  av_expression_res = rep(av_expression, length(unique(sc_cluster_ids)))

  gene_names = rownames(sc_matrix)
  gene_names_res = rep(gene_names, length(unique(sc_cluster_ids)))

  # create data.table with genes, mean expression per cluster, mean expression overall and cluster ids
  comb_dt = data.table::data.table(genes = gene_names_res,
                                   mean_expr = mean_list_res[[1]],
                                   av_expr = av_expression_res,
                                   clusters = group_list_res[[1]])

  # data.table variables
  fold = mean_expr = av_expr = rankFold = clusters = NULL

  # calculate fold change and rank of fold-change
  comb_dt[, fold := log2(mean_expr+1)-log2(av_expr+1)]
  comb_dt[, rankFold := data.table::frank(-fold, ties.method = ties_method), by = clusters]

  # create matrix
  comb_rank_mat = data.table::dcast.data.table(data = comb_dt, genes~clusters, value.var = 'rankFold')
  comb_rank_matrix = dt_to_matrix(comb_rank_mat)
  comb_rank_matrix = comb_rank_matrix[rownames(sc_matrix), unique(sc_cluster_ids)]


  #  # intersect rank matrix with giotto object if given
  #  if(!is.null(gobject) & class(gobject) %in% c('giotto')) {
  #    comb_rank_matrix = comb_rank_matrix[intersect(rownames(comb_rank_matrix), gobject@gene_ID), ]
  #  }

  return(comb_rank_matrix)

}



# * ####
## spatial enrichment functions ####


#' @title do_page_permutation
#' @description creates permutation for the PAGEEnrich test
#' @keywords internal
do_page_permutation<-function(gobject,
                          sig_gene,
                          ntimes){
  # check available gene
  available_ct<-c()
  for (i in colnames(sig_gene)){
    gene_i=rownames(sig_gene)[which(sig_gene[,i]==1)]
    overlap_i=intersect(gene_i,rownames(gobject@norm_expr))
    if (length(overlap_i)<=5){
      output<-paste0("Warning, ",i," only has ",length(overlap_i)," overlapped genes. Will remove it.")
      print(output)
    } else {
      available_ct<-c(available_ct,i)
    }
  }
  if (length(available_ct)==1){
    stop("Only one cell type available.")
  }
  # only continue with genes present in both datasets
  interGene = intersect(rownames(sig_gene), rownames(gobject@norm_expr))
  sign_matrix = sig_gene[interGene,available_ct]

  ct_gene_counts<-NULL
  for (i in 1:dim(sign_matrix)[2]){
    a<-length(which(sign_matrix[,i]==1))
    ct_gene_counts = c(ct_gene_counts,a)
  }
  uniq_ct_gene_counts = unique(ct_gene_counts)
  background_mean_sd = matrix(data=NA,nrow = length(uniq_ct_gene_counts)+1, ncol = 3)
  for (i in 1:length(uniq_ct_gene_counts)){
    gene_num<-uniq_ct_gene_counts[i]
    all_sample_names<-NULL
    all_sample_list<-NULL
    for (j in 1:ntimes){
      set.seed(j)
      random_gene=sample(rownames(gobject@norm_expr),gene_num,replace=FALSE)
      ct_name<-paste("ct",j,sep="")
      all_sample_names = c(all_sample_names,ct_name)
      all_sample_list = c(all_sample_list,list(random_gene))
    }
    random_sig = makeSignMatrixPAGE(all_sample_names,all_sample_list)
    random_DT = runPAGEEnrich(gobject, sign_matrix = random_sig, p_value = F)
    background = unlist(random_DT[,2:dim(random_DT)[2]])
    df_row_name = paste("gene_num_",uniq_ct_gene_counts[i],sep="")
    list_back_i = c(df_row_name,mean(background), stats::sd(background))
    background_mean_sd[i,] = list_back_i
  }
  background_mean_sd[length(uniq_ct_gene_counts)+1,] = c("temp","0","1")
  df_back = data.frame(background_mean_sd,row.names = 1)
  colnames(df_back) = c("mean","sd")
  return(df_back)
}



#' @title runPAGEEnrich_OLD
#' @description Function to calculate gene signature enrichment scores per spatial position using PAGE.
#' @param gobject Giotto object
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param output_enrichment how to return enrichment output
#' @param p_value calculate p-values (boolean, default = FALSE)
#' @param n_times number of permutations to calculate for p_value
#' @param name to give to spatial enrichment results, default = PAGE
#' @param return_gobject return giotto object
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
runPAGEEnrich_OLD <- function(gobject,
                          sign_matrix,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          reverse_log_scale = TRUE,
                          logbase = 2,
                          output_enrichment = c('original', 'zscore'),
                          p_value = FALSE,
                          n_times = 1000,
                          name = NULL,
                          return_gobject = TRUE) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # check parameters
  if(is.null(name)) name = 'PAGE'

  # check available gene
  available_ct<-c()
  for (i in colnames(sign_matrix)){
    gene_i=rownames(sign_matrix)[which(sign_matrix[,i]==1)]
    overlap_i=intersect(gene_i,rownames(expr_values))
    if (length(overlap_i)<=5){
      output<-paste0("Warning, ",i," only has ",length(overlap_i)," overlapped genes. Will remove it.")
      print(output)
    } else {
      available_ct<-c(available_ct,i)
    }
  }

  if (length(available_ct)==1){
    stop("Only one cell type available.")
  }

  # output enrichment
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))

  # only continue with genes present in both datasets
  interGene = intersect(rownames(sign_matrix), rownames(expr_values))
  filterSig = sign_matrix[interGene, available_ct]
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
  cellColSd = apply(geneFold,2,stats::sd)

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
  enrichmentDT = cbind(enrichmentDT, data.table::as.data.table(enrichment))



  ## calculate p-values if requested
  if (p_value==TRUE){

    # check available gene
    available_ct = c()
    for (i in colnames(sign_matrix)){
      gene_i = rownames(sign_matrix)[which(sign_matrix[,i]==1)]
      overlap_i = intersect(gene_i,rownames(gobject@norm_expr))

      if (length(overlap_i)<=5){
        output = paste0("Warning, ",i," only has ",length(overlap_i)," overlapped genes. It will be removed.")
        print(output)
      } else {
        available_ct = c(available_ct, i)
      }
    }

    if (length(available_ct) == 1){
      stop("Only one cell type available.")
    }

    # only continue with genes present in both datasets
    interGene = intersect(rownames(sign_matrix), rownames(gobject@norm_expr))
    filter_sign_matrix = sign_matrix[interGene,available_ct]

    background_mean_sd = do_page_permutation(gobject = gobject,
                                             sig_gene = filter_sign_matrix,
                                             ntimes = n_times)

    for (i in 1:dim(filter_sign_matrix)[2]){
      length_gene = length(which(filter_sign_matrix[,i] == 1))
      join_gene_with_length = paste("gene_num_", length_gene, sep = "")
      mean_i = as.numeric(as.character(background_mean_sd[join_gene_with_length,][[1]]))
      sd_i = as.numeric(as.character(background_mean_sd[join_gene_with_length,][[2]]))
      j = i+1
      enrichmentDT[[j]] = stats::pnorm(enrichmentDT[[j]], mean = mean_i, sd = sd_i, lower.tail = FALSE, log.p = FALSE)
    }
  }



  ## return object or results ##
  if(return_gobject == TRUE) {

    spenr_names = names(gobject@spatial_enrichment)

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_spatial_enrichment')

    # parameters to include
    parameters_list[[update_name]] = c('method used' = 'PAGE',
                                       'enrichment name' = name,
                                       'expression values' = expression_values,
                                       'reverse log scale' = reverse_log_scale,
                                       'logbase' = logbase,
                                       'p-values calculated' = p_value,
                                       'output enrichment scores' = output_enrichment,
                                       'p values calculated' = p_value,
                                       'nr permutations' = n_times)
    gobject@parameters = parameters_list

    gobject@spatial_enrichment[[name]] = enrichmentDT

    return(gobject)

  } else {
    return(enrichmentDT)
  }
}



#' @title PAGE_DT_method
#' @description PAGE data.table method
#' @keywords internal
PAGE_DT_method = function(sign_matrix,
                          expr_values,
                          min_overlap_genes = 5,
                          logbase = 2,
                          reverse_log_scale = TRUE,
                          output_enrichment = c('original', 'zscore'),
                          p_value = FALSE,
                          include_depletion = FALSE,
                          n_times = 1000,
                          max_block = 20e6,
                          verbose = TRUE) {


  # data.table variables
  Var1 = value = Var2 = V1 = marker = nr_markers = fc = cell_ID = zscore = colmean = colSd = pval = NULL
  mean_zscore = sd_zscore = pval_score = NULL

  # output enrichment
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))

  ## identify available cell types
  all_genes = rownames(expr_values)
  sign_matrix = as.matrix(sign_matrix)
  sign_matrix_DT = data.table::as.data.table(reshape2::melt(sign_matrix))
  sign_matrix_DT = sign_matrix_DT[Var1 %in% all_genes]
  detected_DT = sign_matrix_DT[, sum(value), by = Var2]

  lost_cell_types_DT = detected_DT[V1 <= min_overlap_genes]
  for(row in 1:nrow(lost_cell_types_DT)) {
    output<-paste0("Warning, ",lost_cell_types_DT[row][['Var2']]," only has ",lost_cell_types_DT[row][['V1']]," overlapping genes. Will be removed.")
    if(verbose) print(output)
  }
  available_ct = as.character(detected_DT[V1 > min_overlap_genes][['Var2']])

  if (length(available_ct) == 1){
    stop("Only one cell type available.")
  }

  # create subset of sinature matrix
  interGene = intersect(rownames(sign_matrix), rownames(expr_values))
  filterSig = sign_matrix[interGene, available_ct]

  # create fold expression for each gene in each spot
  # calculate mean gene expression
  if(reverse_log_scale == TRUE) {
    mean_gene_expr = log(rowMeans(logbase^expr_values-1, dims = 1)+1)
  } else {
    mean_gene_expr = rowMeans(expr_values)
  }
  geneFold = expr_values - mean_gene_expr

  # calculate sample/spot mean and sd
  cellColMean = colMeans(geneFold)
  cellColSd = apply(geneFold, 2, stats::sd)
  cellColMeanSd =  data.table::data.table(cell_ID = names(cellColMean),
                                          colmean = cellColMean,
                                          colSd = cellColSd)

  filterSig_DT = data.table::as.data.table(reshape2::melt(filterSig))
  colnames(filterSig_DT) = c('gene', 'cell_type', 'marker')
  sub_ct_DT = filterSig_DT[marker == 1]
  sub_ct_DT[, nr_markers := .N, by = cell_type]

  ## reshape gene fold-expression
  geneFold_DT = data.table::as.data.table(reshape2::melt(geneFold))
  colnames(geneFold_DT) = c('gene', 'cell_ID', 'fc')

  mergetest = data.table::merge.data.table(sub_ct_DT, geneFold_DT, by = 'gene')
  mergetest = mergetest[, mean(fc), by = .(cell_type, cell_ID, nr_markers)]
  mergetest = data.table::merge.data.table(mergetest, cellColMeanSd, by = 'cell_ID')
  mergetest[, zscore := ((V1 - colmean)* nr_markers^(1/2)) / colSd]

  if(output_enrichment == 'zscore') {
    mergetest[, zscore := scale(zscore), by = 'cell_type']
  }




  ## return p-values based on permutations ##
  if(p_value == TRUE) {

    ## 1. get number of markers instructions ##
    sample_intrs = unique(sub_ct_DT[,.(cell_type, nr_markers)])


    ## 2. first create the random samples all together ##
    cell_type_list = list()
    perm_type_list = list()
    for(row in 1:nrow(sample_intrs)) {

      cell_type = sample_intrs[row][['cell_type']]
      nr_genes = as.numeric(sample_intrs[row][['nr_markers']])

      gene_list = list()
      perm_list = list()
      for(i in 1:n_times) {
        sampled_genes = sample(rownames(expr_values), size = nr_genes)
        gene_list[[i]] = sampled_genes
        perm_list[[i]] = rep(paste0('p_',i), nr_genes)
      }

      gene_res = unlist(gene_list)
      names(gene_res) = rep(cell_type, length(gene_res))
      cell_type_list[[row]] = gene_res

      perm_res = unlist(perm_list)
      perm_type_list[[row]] = perm_res

    }

    cell_type_perm = unlist(cell_type_list)
    perm_round = unlist(perm_type_list)

    cell_type_perm_DT = data.table::data.table(cell_type = names(cell_type_perm),
                                               gene = cell_type_perm,
                                               round = perm_round)

    sample_intrs_vec = sample_intrs$nr_markers
    names(sample_intrs_vec) = sample_intrs$cell_type
    cell_type_perm_DT[, nr_markers := sample_intrs_vec[cell_type]]


    ## 3. decide on number of blocks to process ##
    nr_perm_lines = as.numeric(nrow(cell_type_perm_DT))
    nr_spots = as.numeric(ncol(expr_values))
    total_lines = nr_spots * nr_perm_lines
    nr_groups = round(total_lines / max_block)

    ## 4. create groups
    all_perms = unique(perm_round)
    all_perms_num = 1:length(all_perms)
    names(all_perms_num) = all_perms
    group_labels = paste0('group_',1:nr_groups)
    groups_vec = cut(all_perms_num, breaks = nr_groups, labels = group_labels)
    names(all_perms) = groups_vec


    ## 5. do random enrichment per block
    res_list = list()
    for(group_i in 1:length(group_labels)) {

      group = group_labels[group_i]
      sub_perms = all_perms[names(all_perms) == group]
      cell_type_perm_DT_sub = cell_type_perm_DT[round %in% sub_perms]

      mergetest_perm_sub = data.table::merge.data.table(cell_type_perm_DT_sub, geneFold_DT, allow.cartesian = TRUE)
      mergetest_perm_sub = mergetest_perm_sub[, mean(fc), by = .(cell_type, cell_ID, nr_markers, round)]
      mergetest_perm_sub = data.table::merge.data.table(mergetest_perm_sub, cellColMeanSd, by = 'cell_ID')
      mergetest_perm_sub[, zscore := ((V1 - colmean)* nr_markers^(1/2)) / colSd]

      res_list[[group_i]] = mergetest_perm_sub

    }

    res_list_comb = do.call('rbind', res_list)
    res_list_comb_average = res_list_comb[, .(mean_zscore = mean(zscore), sd_zscore = stats::sd(zscore)), by = c('cell_ID', 'cell_type')]
    mergetest_final = data.table::merge.data.table(mergetest, res_list_comb_average, by = c('cell_ID', 'cell_type'))

    ## calculate p.values based on normal distribution
    if(include_depletion == TRUE) {
      mergetest_final[, pval := stats::pnorm(abs(zscore), mean = mean_zscore, sd = sd_zscore, lower.tail = FALSE, log.p = FALSE)]
    } else {
      mergetest_final[, pval := stats::pnorm(zscore, mean = mean_zscore, sd = sd_zscore, lower.tail = FALSE, log.p = FALSE)]
    }

    data.table::setorder(mergetest_final, pval)

    ## calculate pval_score
    if(include_depletion == TRUE) {
      mergetest_final[, pval_score := sign(zscore)*-log10(pval)]
    } else {
      mergetest_final[, pval_score := -log10(pval)]
    }


    resultmatrix = data.table::dcast(mergetest_final, formula = cell_ID~cell_type, value.var = 'pval_score')
    return(list(DT = mergetest_final, matrix = resultmatrix))


  } else {

    resultmatrix = data.table::dcast(mergetest, formula = cell_ID~cell_type, value.var = 'zscore')
    return(list(DT = mergetest, matrix = resultmatrix))

  }

}



#' @title runPAGEEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using PAGE.
#' @param gobject Giotto object
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param min_overlap_genes minimum number of overlapping genes in sign_matrix required to calculate enrichment
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param output_enrichment how to return enrichment output
#' @param p_value calculate p-values (boolean, default = FALSE)
#' @param include_depletion calculate both enrichment and depletion
#' @param n_times number of permutations to calculate for p_value
#' @param max_block number of lines to process together (default = 20e6)
#' @param name to give to spatial enrichment results, default = PAGE
#' @param verbose be verbose
#' @param return_gobject return giotto object
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
runPAGEEnrich <- function(gobject,
                          sign_matrix,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          min_overlap_genes = 5,
                          reverse_log_scale = TRUE,
                          logbase = 2,
                          output_enrichment = c('original', 'zscore'),
                          p_value = FALSE,
                          include_depletion = FALSE,
                          n_times = 1000,
                          max_block = 20e6,
                          name = NULL,
                          verbose = TRUE,
                          return_gobject = TRUE) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # check parameters
  if(is.null(name)) name = 'PAGE'

  PAGE_results = PAGE_DT_method(sign_matrix = sign_matrix,
                                expr_values = expr_values,
                                min_overlap_genes = min_overlap_genes,
                                logbase = logbase,
                                reverse_log_scale = reverse_log_scale,
                                output_enrichment = c('original', 'zscore'),
                                p_value = p_value,
                                include_depletion = include_depletion,
                                n_times = n_times,
                                max_block = max_block,
                                verbose = verbose)


  ## return object or results ##
  if(return_gobject == TRUE) {

    spenr_names = names(gobject@spatial_enrichment)

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_spatial_enrichment')

    # parameters to include
    parameters_list[[update_name]] = c('method used' = 'PAGE',
                                       'enrichment name' = name,
                                       'expression values' = expression_values,
                                       'reverse log scale' = reverse_log_scale,
                                       'logbase' = logbase,
                                       'output enrichment scores' = output_enrichment,
                                       'p values calculated' = p_value,
                                       'include depletion' = include_depletion,
                                       'nr permutations' = n_times)

    gobject@parameters = parameters_list
    gobject@spatial_enrichment[[name]] = PAGE_results[['matrix']]

    return(gobject)

  } else {
    return(PAGE_results)
  }
}







#' @title PAGEEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using PAGE.
#' @inheritDotParams runPAGEEnrich
#' @seealso \code{\link{runPAGEEnrich}}
#' @export
PAGEEnrich <- function(...) {

  .Deprecated(new = "runPAGEEnrich")

  runPAGEEnrich(...)

}





#' @title do_rank_permutation
#' @description creates permutation for the rankEnrich test
#' @keywords internal
do_rank_permutation <- function(sc_gene, n){
  random_df <- data.frame(matrix(ncol = n, nrow = length(sc_gene)))
  for (i in 1:n){
    set.seed(i)
    random_rank = sample(1:length(sc_gene), length(sc_gene), replace=FALSE)
    random_df[,i]<-random_rank
  }
  rownames(random_df)<-sc_gene
  return(random_df)
}


#' @title runRankEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using a rank based approach.
#' @param gobject Giotto object
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param output_enrichment how to return enrichment output
#' @param ties_method how to handle rank ties
#' @param p_value calculate p-values (boolean, default = FALSE)
#' @param n_times number of permutations to calculate for p_value
#' @param rbp_p fractional binarization threshold (default = 0.99)
#' @param num_agg number of top genes to aggregate (default = 100)
#' @param name to give to spatial enrichment results, default = rank
#' @param return_gobject return giotto object
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
runRankEnrich <- function(gobject,
                       sign_matrix,
                       expression_values = c('normalized', "raw", 'scaled', 'custom'),
                       reverse_log_scale = TRUE,
                       logbase = 2,
                       output_enrichment = c('original', 'zscore'),
                       ties_method = c("random", "max"),
                       p_value = FALSE,
                       n_times = 1000,
                       rbp_p = 0.99,
                       num_agg = 100,
                       name = NULL,
                       return_gobject = TRUE) {

  # determine ties.method
  ties_method = match.arg(ties_method, choices = c("random", "max"))

  # expression values to be used
  values = match.arg(expression_values, c('normalized', "raw", 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  if(values=="raw"){
    expr_values = Matrix::as.matrix(expr_values)
  }

  # check parameters
  if(is.null(name)) name = 'rank'

  #check gene list
  interGene = intersect(rownames(sign_matrix), rownames(expr_values))
  if (length(interGene)<100){
    stop("Please check the gene numbers or names of scRNA-seq. The names of scRNA-seq should be consistent with spatial data.")
  }

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
  #geneFold = expr_values - mean_gene_expr
  #rankFold = t(matrixStats::colRanks(-geneFold, ties.method = "first"))

  ties_1 = ties_method
  ties_2 = ties_method
  if(ties_method == "max"){
    ties_1 = "min"
    ties_2 = "max"
  }
  #else ties_1=ties_2 is equal to random
  geneFold = expr_values
  geneFold = matrixStats::rowRanks(geneFold, ties.method = ties_1)
  rankFold = t(matrixStats::colRanks(-geneFold, ties.method = ties_2))

  rownames(rankFold) = rownames(expr_values)
  colnames(rankFold) = colnames(expr_values)

  for (i in (1:dim(sign_matrix)[2])){

    signames = rownames(sign_matrix)[which(sign_matrix[,i]>0)]
    interGene = intersect(signames, rownames(rankFold))
    filterSig = sign_matrix[interGene,]
    filterRankFold = rankFold[interGene,]

    multiplyRank = (filterRankFold*filterSig[,i])^(1/2)
    rpb = (1.0 - rbp_p)*(rbp_p^(multiplyRank-1))

    vectorX = rep(NA, dim(filterRankFold)[2])

    for (j in (1:dim(filterRankFold)[2])){
      toprpb = sort(rpb[,j],decreasing = T)
      zscore = sum(toprpb[1:num_agg])
      vectorX[j] = zscore
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
  enrichmentDT = cbind(enrichmentDT, data.table::as.data.table(enrichment))


  # default name for page enrichment

  if (p_value==TRUE){
    random_rank = do_rank_permutation(sc_gene = rownames(sign_matrix),
                                      n = n_times)
    random_DT = runRankEnrich(gobject = gobject,
                              sign_matrix = random_rank,
                              expression_values = expression_values,
                              reverse_log_scale = reverse_log_scale,
                              logbase = logbase,
                              output_enrichment = output_enrichment,
                              p_value = FALSE)

    background = unlist(random_DT[,2:dim(random_DT)[2]])
    fit.gamma = fitdistrplus::fitdist(background, distr = "gamma", method = "mle")
    pvalue_DT = enrichmentDT
    enrichmentDT[,2:dim(enrichmentDT)[2]] = lapply(enrichmentDT[,2:dim(enrichmentDT)[2]], function(x)
    {stats::pgamma(x, fit.gamma$estimate[1], rate = fit.gamma$estimate[2], lower.tail = FALSE,log.p = FALSE)})
  }


  ## return object or results ##
  if(return_gobject == TRUE) {

    spenr_names = names(gobject@spatial_enrichment)

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_spatial_enrichment')

    # parameters to include
    parameters_list[[update_name]] = c('method used' = 'rank',
                                       'enrichment name' = name,
                                       'expression values' = expression_values,
                                       'reverse log scale' = reverse_log_scale,
                                       'logbase' = logbase,
                                       'p-values calculated' = p_value,
                                       'output enrichment scores' = output_enrichment,
                                       'p values calculated' = p_value,
                                       'nr permutations' = n_times)
    gobject@parameters = parameters_list

    gobject@spatial_enrichment[[name]] = enrichmentDT

    return(gobject)

  } else {
    return(enrichmentDT)
  }

}



#' @title rankEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using a rank based approach.
#' @inheritDotParams runRankEnrich
#' @seealso \code{\link{runRankEnrich}}
#' @export
rankEnrich <- function(...) {

  .Deprecated(new = "runRankEnrich")

  runRankEnrich(...)

}



#' @title runHyperGeometricEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using a hypergeometric test.
#' @param gobject Giotto object
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param top_percentage percentage of cells that will be considered to have gene expression with matrix binarization
#' @param output_enrichment how to return enrichment output
#' @param p_value calculate p-values (boolean, default = FALSE)
#' @param name to give to spatial enrichment results, default = rank
#' @param return_gobject return giotto object
#' @return data.table with enrichment results
#' @details The enrichment score is calculated based on the p-value from the
#' hypergeometric test, -log10(p-value).
#' @export
runHyperGeometricEnrich <- function(gobject,
                                 sign_matrix,
                                 expression_values = c('normalized', 'scaled', 'custom'),
                                 reverse_log_scale = TRUE,
                                 logbase = 2,
                                 top_percentage = 5,
                                 output_enrichment = c('original', 'zscore'),
                                 p_value = FALSE,
                                 name = NULL,
                                 return_gobject = TRUE) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # check parameters
  if(is.null(name)) name = 'hypergeometric'

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
  quantilecut = apply(foldChange, 2 , stats::quantile , probs = top_q, na.rm = TRUE )
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
      enrich<-(0-log10(stats::phyper(x, m, n, k, log = FALSE,lower.tail = FALSE)))
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
  enrichmentDT = cbind(enrichmentDT, data.table::as.data.table(enrichment))


  ## calculate p-values ##
  if (p_value == TRUE){
    enrichmentDT[,2:dim(enrichmentDT)[2]] = lapply(enrichmentDT[,2:dim(enrichmentDT)[2]],function(x){10^(-x)})
  }



  ## return object or results ##
  if(return_gobject == TRUE) {

    spenr_names = names(gobject@spatial_enrichment)

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_spatial_enrichment')

    # parameters to include
    parameters_list[[update_name]] = c('method used' = 'rank',
                                       'enrichment name' = name,
                                       'expression values' = expression_values,
                                       'reverse log scale' = reverse_log_scale,
                                       'logbase' = logbase,
                                       'top percentage' = top_percentage,
                                       'p-values calculated' = p_value,
                                       'output enrichment scores' = output_enrichment,
                                       'p values calculated' = p_value)
    gobject@parameters = parameters_list

    gobject@spatial_enrichment[[name]] = enrichmentDT

    return(gobject)

  } else {
    return(enrichmentDT)
  }
}


#' @title hyperGeometricEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using a hypergeometric test.
#' @inheritDotParams runHyperGeometricEnrich
#' @seealso \code{\link{runHyperGeometricEnrich}}
#' @export
hyperGeometricEnrich <- function(...) {

  .Deprecated(new = "runHyperGeometricEnrich")

  runHyperGeometricEnrich(...)

}




#' @title runSpatialEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using an enrichment test.
#' @param gobject Giotto object
#' @param enrich_method method for gene signature enrichment calculation
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param min_overlap_genes minimum number of overlapping genes in sign_matrix required to calculate enrichment (PAGE)
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param p_value calculate p-value (default = FALSE)
#' @param n_times (page/rank) number of permutation iterations to calculate p-value
#' @param rbp_p (rank) fractional binarization threshold (default = 0.99)
#' @param num_agg (rank) number of top genes to aggregate (default = 100)
#' @param max_block number of lines to process together (default = 20e6)
#' @param top_percentage (hyper) percentage of cells that will be considered to have gene expression with matrix binarization
#' @param output_enrichment how to return enrichment output
#' @param name to give to spatial enrichment results, default = PAGE
#' @param verbose be verbose
#' @param return_gobject return giotto object
#' @return Giotto object or enrichment results if return_gobject = FALSE
#' @details For details see the individual functions:
#' \itemize{
#'   \item{PAGE: }{\code{\link{runPAGEEnrich}}}
#'   \item{Rank: }{\code{\link{runRankEnrich}}}
#'   \item{Hypergeometric: }{\code{\link{runHyperGeometricEnrich}}}
#' }
#'
#' @export
runSpatialEnrich = function(gobject,
                            enrich_method = c('PAGE', 'rank', 'hypergeometric'),
                            sign_matrix,
                            expression_values = c('normalized', 'scaled', 'custom'),
                            min_overlap_genes = 5,
                            reverse_log_scale = TRUE,
                            logbase = 2,
                            p_value = FALSE,
                            n_times = 1000,
                            rbp_p = 0.99,
                            num_agg = 100,
                            max_block = 20e6,
                            top_percentage = 5,
                            output_enrichment = c('original', 'zscore'),
                            name = NULL,
                            verbose = TRUE,
                            return_gobject = TRUE) {


  enrich_method = match.arg(enrich_method, choices = c('PAGE', 'rank', 'hypergeometric'))
  expression_values = match.arg(expression_values, choices = c('normalized', 'scaled', 'custom'))
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))


  if(enrich_method == 'PAGE') {

    results =  runPAGEEnrich(gobject = gobject,
                          sign_matrix = sign_matrix,
                          expression_values = expression_values,
                          reverse_log_scale = reverse_log_scale,
                          logbase = logbase,
                          output_enrichment = output_enrichment,
                          p_value = p_value,
                          n_times = n_times,
                          name = name,
                          return_gobject = return_gobject)

  } else if(enrich_method == 'rank') {

    results =  runRankEnrich(gobject = gobject,
                             sign_matrix = sign_matrix,
                             expression_values = expression_values,
                             reverse_log_scale = reverse_log_scale,
                             logbase = logbase,
                             output_enrichment = output_enrichment,
                             p_value = p_value,
                             n_times = n_times,
                             rbp_p = rbp_p,
                             num_agg = num_agg,
                             name = name,
                             return_gobject = return_gobject)


  } else if(enrich_method == 'hypergeometric'){

    results =  runHyperGeometricEnrich(gobject = gobject,
                                       sign_matrix = sign_matrix,
                                       expression_values = expression_values,
                                       reverse_log_scale = reverse_log_scale,
                                       logbase = logbase,
                                       top_percentage = top_percentage,
                                       output_enrichment = output_enrichment,
                                       p_value = p_value,
                                       name = name,
                                       return_gobject = return_gobject)
  }

  return(results)
}


#' @title createSpatialEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using an enrichment test.
#' @inheritDotParams runSpatialEnrich
#' @seealso \code{\link{runSpatialEnrich}}
#' @export
createSpatialEnrich <- function(...) {

  .Deprecated(new = "runSpatialEnrich")

  runSpatialEnrich(...)

}



# * ####
## spatial deconvolution functions ####


#' @title enrich_deconvolution
#' @description Rui to fill in
#' @keywords internal
enrich_deconvolution<-function(expr,
                               log_expr,
                               cluster_info,
                               ct_exp,
                               cutoff){
  #####generate enrich 0/1 matrix based on expression matrix
  ct_exp <- ct_exp[rowSums(ct_exp)>0,]
  enrich_matrix<-matrix(0,nrow=dim(ct_exp)[1],ncol=dim(ct_exp)[2])
  rowmax_col<-Rfast::rowMaxs(ct_exp)
  for (i in 1:length(rowmax_col)){
    enrich_matrix[i,rowmax_col[i]]=1
  }
  colsum_ct_binary <- colSums(enrich_matrix)
  for (i in 1:length(colsum_ct_binary)){
    if (colsum_ct_binary[i] <= 2){
      rank <- rank(-ct_exp[,i])
      enrich_matrix[rank <=2, i] =1
    }
  }
  rownames(enrich_matrix)<-rownames(ct_exp)
  colnames(enrich_matrix)<-colnames(ct_exp)
  # print(enrich_matrix)
  #####page enrich
  enrich_result<-enrich_analysis(log_expr,enrich_matrix)
  #####initialize dwls matrix
  dwls_results<-matrix(0,nrow =dim(enrich_matrix)[2],ncol = dim(expr)[2])
  rownames(dwls_results)<-colnames(enrich_matrix)
  colnames(dwls_results)<-colnames(expr)
  cluster_sort<-sort(unique(cluster_info))
  cluster_info<-cluster_info
  for (i in 1:length(cluster_sort)){
    cluster_i_enrich<-enrich_result[,which(cluster_info==cluster_sort[i])]
    row_i_max<-Rfast::rowMaxs(cluster_i_enrich,value = TRUE)
    ct<-rownames(enrich_result)[which(row_i_max>cutoff)]
    if (length(ct)<2){
      sort_rank<-sort(row_i_max,decreasing = T)
      ct<-rownames(enrich_result)[which(row_i_max>=sort_rank[2])]
    }
    ct_gene<-c()
    for (j in 1:length(ct)){
      sig_gene_j<-rownames(enrich_matrix)[which(enrich_matrix[,ct[j]]==1)]
      ct_gene<-c(ct_gene,sig_gene_j)
    }
    uniq_ct_gene<-intersect(rownames(expr),unique(ct_gene))
    select_sig_exp<-ct_exp[uniq_ct_gene,ct]
    cluster_i_cell<-which(cluster_info==cluster_sort[i])
    cluster_cell_exp<-expr[uniq_ct_gene,cluster_i_cell]
    # print(cluster_cell_exp)
    cluster_i_dwls<-optimize_deconvolute_dwls(cluster_cell_exp,select_sig_exp)
    dwls_results[ct,cluster_i_cell]<-cluster_i_dwls
  }
  #####remove negative values
  for (i in dim(dwls_results)[1]){
    negtive_index<-which(dwls_results[i,]<0)
    dwls_results[i,negtive_index]==0
  }
  return(dwls_results)
}

#' @title spot_deconvolution
#' @description Rui to fill in
#' @keywords internal
spot_deconvolution<-function(expr,
                             cluster_info,
                             ct_exp,
                             binary_matrix){
  #####generate enrich 0/1 matrix based on expression matrix
  enrich_matrix<-matrix(0,nrow=dim(ct_exp)[1],ncol=dim(ct_exp)[2])
  rowmax_col<-Rfast::rowMaxs(ct_exp)
  for (i in 1:length(rowmax_col)){
    enrich_matrix[i,rowmax_col[i]]=1
  }
  rownames(enrich_matrix)<-rownames(ct_exp)
  colnames(enrich_matrix)<-colnames(ct_exp)

    cluster_sort<-sort(unique(cluster_info))
  ####initialize dwls matrix
  dwls_results<-matrix(0,nrow =dim(ct_exp)[2],ncol = dim(expr)[2])
  rownames(dwls_results)<-colnames(ct_exp)
  colnames(dwls_results)<-colnames(expr)
  #print(binary_matrix)
  for (i in 1:length(cluster_sort)){
    cluster_i_matrix<-binary_matrix[,which(cluster_info==cluster_sort[i])]
    row_i_max<-Rfast::rowMaxs(cluster_i_matrix,value = TRUE)
    ct_i<-rownames(cluster_i_matrix)[which(row_i_max==1)]
    ########calculate proportion based on binarized deconvolution results at first step
    if (length(ct_i)==1){
      dwls_results[ct_i[1],which(cluster_info==cluster_sort[i])]==1
    } else {
      ct_gene<-c()
      for (j in 1:length(ct_i)){
        sig_gene_j<-rownames(enrich_matrix)[which(enrich_matrix[,ct_i[j]]==1)]
        ct_gene<-c(ct_gene,sig_gene_j)
      }
      uniq_ct_gene<-intersect(rownames(expr),unique(ct_gene))
      select_sig_exp<-ct_exp[uniq_ct_gene,ct_i]
      cluster_i_cell<-which(cluster_info==cluster_sort[i])
      cluster_cell_exp<-expr[uniq_ct_gene,cluster_i_cell]
      ######calculate
      ######overlap signature with spatial genes
      all_exp<-rowMeans(cluster_cell_exp)
      solution_all_exp<-solve_OLS_internal(select_sig_exp,all_exp)
      constant_J<-find_dampening_constant(select_sig_exp,all_exp,solution_all_exp)
      ######deconvolution for each spot
      for(k in 1:(dim(cluster_cell_exp)[2])){
        B<-Matrix::as.matrix(cluster_cell_exp[,k])
        ct_spot_k<-rownames(cluster_i_matrix)[which(cluster_i_matrix[,k]==1)]
        if (length(ct_spot_k)==1){
          dwls_results[ct_spot_k[1],colnames(cluster_cell_exp)[k]]<-1
        } else {
          ct_k_gene<-c()
          for (m in 1:length(ct_spot_k)){
            sig_gene_k<-rownames(enrich_matrix)[which(enrich_matrix[,ct_spot_k[m]]==1)]
            ct_k_gene<-c(ct_k_gene,sig_gene_k)
          }
          uniq_ct_k_gene<-intersect(rownames(ct_exp),unique(ct_k_gene))
          S_k<-Matrix::as.matrix(ct_exp[uniq_ct_k_gene,ct_spot_k])
          solDWLS<-optimize_solveDampenedWLS(S_k,B[uniq_ct_k_gene,],constant_J)
          dwls_results[names(solDWLS),colnames(cluster_cell_exp)[k]]<-solDWLS
        }
      }
    }
  }
  #####remove negative values
  for (i in dim(dwls_results)[1]){
    negtive_index<-which(dwls_results[i,]<0)
    dwls_results[i,negtive_index]==0
  }
  return(dwls_results)
}



#' @title cluster_enrich_analysis
#' @description Rui to fill in
#' @keywords internal
cluster_enrich_analysis <- function(exp_matrix,
                                    cluster_info,
                                    enrich_sig_matrix) {
  uniq_cluster<-sort(unique(cluster_info))
  if(length(uniq_cluster) == 1) {
    stop("Only one cluster identified, need at least two.")
  }
  cluster_exp<-NULL
  for (i in uniq_cluster){
    cluster_exp<-cbind(cluster_exp,(apply(exp_matrix,1,function(y) mean(y[which(cluster_info==i)]))))
  }
  log_cluster_exp<-log2(cluster_exp+1)
  colnames(log_cluster_exp)<-uniq_cluster
  cluster_enrich<-enrich_analysis(log_cluster_exp,enrich_sig_matrix)
  return(cluster_enrich)
}


#' @title enrich_analysis
#' @description Rui to fill in
#' @keywords internal
enrich_analysis <- function(expr_values,
                            sign_matrix) {

  # output enrichment
  # only continue with genes present in both datasets
  interGene = intersect(rownames(sign_matrix), rownames(expr_values))
  filterSig = sign_matrix[interGene,]
  signames = rownames(filterSig)[which(filterSig[,1]==1)]
  # calculate mean gene expression
  #mean_gene_expr = rowMeans(expr_values)
  mean_gene_expr = log2(rowMeans(2^expr_values-1, dims = 1)+1)
  geneFold = expr_values - mean_gene_expr
  # calculate sample/spot mean and sd
  cellColMean = apply(geneFold,2, mean)
  cellColSd = apply(geneFold,2, stats::sd)

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
  return(enrichment)
}



#' @title optimize_deconvolute_dwls
#' @description Rui to fill in
#' @keywords internal
optimize_deconvolute_dwls <- function(exp,
                                      Signature) {
  ######overlap signature with spatial genes
  Genes<-intersect(rownames(Signature),rownames(exp))
  S<-Signature[Genes,]
  S<-Matrix::as.matrix(S)
  Bulk<-Matrix::as.matrix(exp)
  subBulk = Bulk[Genes,]
  allCounts_DWLS<-NULL
  all_exp<-rowMeans(exp)

  solution_all_exp<-solve_OLS_internal(S,all_exp[Genes])

  constant_J<-find_dampening_constant(S,all_exp[Genes],solution_all_exp)
  for(j in 1:(dim(subBulk)[2])){
    B<-subBulk[,j]
    if (sum(B)>0){
      solDWLS<-optimize_solveDampenedWLS(S,B,constant_J)
    } else{
      solDWLS <- rep(0, length(B))
      names(solDWLS) <- names(B)
    }
    allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  }
  colnames(allCounts_DWLS)<-colnames(exp)
  return(allCounts_DWLS)
}


#' @title optimize_solveDampenedWLS
#' @description Rui to fill in
#' @keywords internal
optimize_solveDampenedWLS<-function(S,
                                    B,
                                    constant_J){
  #first solve OLS, use this solution to find a starting point for the weights
  solution<-solve_OLS_internal(S,B)
  #now use dampened WLS, iterate weights until convergence
  iterations<-0
  changes<-c()
  #find dampening constant for weights using cross-validation
  j<-constant_J
  change<-1
  while(change>.01 & iterations<1000){
    newsolution<-solve_dampened_WLSj(S,B,solution,j)
    #decrease step size for convergence
    solutionAverage<-rowMeans(cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
    change<-norm(Matrix::as.matrix(solutionAverage-solution))
    solution<-solutionAverage
    iterations<-iterations+1
    changes<-c(changes,change)
  }
  #print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}


#' @title find_dampening_constant
#' @description find a dampening constant for the weights using cross-validation
#' @keywords internal
find_dampening_constant<-function(S,
                                  B,
                                  goldStandard){

  solutionsSd<-NULL
  #goldStandard is used to define the weights
  sol = goldStandard
  ws = as.vector((1/(S%*%sol))^2)
  wsScaled = ws/min(ws)
  wsScaledMinusInf = wsScaled
  #ignore infinite weights
  if(max(wsScaled) == "Inf"){
    wsScaledMinusInf = wsScaled[-which(wsScaled == "Inf")]
  }
  #try multiple values of the dampening constant (multiplier)
  #for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
    multiplier = 1*2^(j-1)
    wsDampened = wsScaled
    wsDampened[which(wsScaled>multiplier)] = multiplier
    solutions<-NULL
    seeds<-c(1:100)
    for (i in 1:100){
      set.seed(seeds[i]) #make nondeterministic
      subset = sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
      #solve dampened weighted least squares for subset
      fit = stats::lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
      sol = fit$coef*sum(goldStandard)/sum(fit$coef)
      solutions = cbind(solutions,sol)
    }
    solutionsSd = cbind(solutionsSd,apply(solutions, 1, stats::sd))
  }
  #choose dampening constant that results in least cross-validation variance
  j = which.min(colMeans(solutionsSd^2))
  return(j)
}


#' @title solve_OLS_internal
#' @description basic functions for dwls
#' @keywords internal
solve_OLS_internal<-function(S,
                             B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  solution<-quadprog::solve.QP(D,d,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#
#' @title solve_dampened_WLSj
#' @description solve WLS given a dampening constant
#' @keywords internal
solve_dampened_WLSj<-function(S,
                              B,
                              goldStandard,
                              j){
  multiplier<-1*2^(j-1)
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsDampened<-wsScaled
  wsDampened[which(wsScaled>multiplier)]<-multiplier
  W<-diag(wsDampened)
  D<-t(S)%*%W%*%S
  d<- t(S)%*%W%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D,"2")
  solution<-quadprog::solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}


#' @title runDWLSDeconv
#' @description Function to perform DWLS deconvolution based on single cell expression data
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param logbase base used for log normalization
#' @param cluster_column name of cluster column
#' @param sign_matrix sig matrix for deconvolution
#' @param n_cell number of cells per spot
#' @param cutoff cut off (default = 2)
#' @param name name to give to spatial deconvolution results, default = DWLS
#' @param return_gobject return giotto object
#' @return giotto object or deconvolution results
#' @export
runDWLSDeconv <- function(gobject,
                          expression_values = c('normalized'),
                          logbase = 2,
                          cluster_column = 'leiden_clus',
                          sign_matrix,
                          n_cell = 50,
                          cutoff = 2,
                          name = NULL,
                          return_gobject = TRUE) {


  # verify if optional package is installed
  package_check(pkg_name = "quadprog", repository = "CRAN")


  ## check parameters ##

  # check parameters
  if(is.null(name)) name = 'DWLS'

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # #transform expression data to no log data
  nolog_expr = logbase^(expr_values)-1

  # cluster column
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n cluster column not found \n')
  }
  cluster = cell_metadata[[cluster_column]]


  #####getting overlapped gene lists
  sign_matrix <- as.matrix(sign_matrix)
  intersect_gene = intersect(rownames(sign_matrix), rownames(nolog_expr))
  filter_Sig = sign_matrix[intersect_gene,]
  filter_expr = nolog_expr[intersect_gene,]
  filter_log_expr = expr_values[intersect_gene,]

  #####first round spatial deconvolution ##spot or cluster
  enrich_spot_proportion <- enrich_deconvolution(expr = filter_expr,
                                                 log_expr = filter_log_expr,
                                                 cluster_info = cluster,
                                                 ct_exp = filter_Sig,
                                                 cutoff = cutoff)

  ####re-deconvolution based on spatial resolution
  resolution = (1/n_cell)
  binarize_proportion = ifelse(enrich_spot_proportion >= resolution, 1, 0)
  spot_proportion <- spot_deconvolution(expr = filter_expr,
                                        cluster_info = cluster,
                                        ct_exp = filter_Sig,
                                        binary_matrix = binarize_proportion)
  deconvolutionDT = data.table::data.table(cell_ID = colnames(spot_proportion))
  deconvolutionDT = cbind(deconvolutionDT, as.data.table(t(spot_proportion)))


  ## return object or results ##
  if(return_gobject == TRUE) {

    spenr_names = names(gobject@spatial_enrichment)

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_spatial_deconvolution')

    # parameters to include
    parameters_list[[update_name]] = c('method used' = 'DWLS',
                                       'deconvolution name' = name,
                                       'expression values' = expression_values,
                                       'logbase' = logbase,
                                       'cluster column used' = cluster_column,
                                       'number of cells per spot' = n_cell,
                                       'used cut off' = cutoff)

    gobject@parameters = parameters_list

    gobject@spatial_enrichment[[name]] = deconvolutionDT

    return(gobject)

  } else {
    return(deconvolutionDT)
  }

}



#' @title runSpatialDeconv
#' @description Function to perform deconvolution based on single cell expression data
#' @param gobject giotto object
#' @param deconv_method method to use for deconvolution
#' @param expression_values expression values to use
#' @param logbase base used for log normalization
#' @param cluster_column name of cluster column
#' @param sign_matrix signature matrix for deconvolution
#' @param n_cell number of cells per spot
#' @param cutoff cut off (default = 2)
#' @param name name to give to spatial deconvolution results
#' @param return_gobject return giotto object
#' @return giotto object or deconvolution results
#' @export
runSpatialDeconv <- function(gobject,
                             deconv_method = c('DWLS'),
                             expression_values = c('normalized'),
                             logbase = 2,
                             cluster_column = 'leiden_clus',
                             sign_matrix,
                             n_cell = 50,
                             cutoff = 2,
                             name = NULL,
                             return_gobject = TRUE) {


  deconv_method = match.arg(deconv_method, choices = c('DWLS'))


  if(deconv_method == 'DWLS') {

    results =  runDWLSDeconv(gobject = gobject,
                             expression_values = expression_values,
                             logbase = logbase,
                             cluster_column = cluster_column,
                             sign_matrix = sign_matrix,
                             n_cell = n_cell,
                             cutoff = cutoff,
                             name = name,
                             return_gobject = return_gobject)
  }

  return(results)

}