
#' @title convertSignListToMatrix
#' @description Function to convert list of signature genes (e.g. for cell types) into
#' a binary matrix format that can be used with the PAGE enrichment option.
#' @param sign_names vector with names for each provided gene signature
#' @param sign_list list of gene signatures
#' @return matrix
#' @export
#' @examples
#'     convertSignListToMatrix()
convertSignListToMatrix = function(sign_names,
                                   sign_list) {

  if(!is.list(sign_list)) {
    stop('\n sign_list needs to be a list of signatures for each cell type / process \n')
  }
  if(length(sign_names) != length(sign_list)) {
    stop('\n the number of names needs to match the number of signatures provided \n')
  }

  # create genes and signatures
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


#' @title PAGEEnrich
#' @description Function to calculate gene signature enrichment scores per spatial position using PAGE.
#' @param gobject Giotto object
#' @param sign_matrix Matrix of signature genes for each cell type / process
#' @param expression_values expression values to use
#' @param reverse_log_scale reverse expression values from log scale
#' @param logbase log base to use if reverse_log_scale = TRUE
#' @param output_enrichment how to return enrichment output
#' @return data.table with enrichment results
#' @details The enrichment Z score is calculated by using method (PAGE) from
#' Kim SY et al., BMC bioinformatics, 2005 as Z = ((Sm – μ)*m^(1/2)) / δ.
#' For each gene in each spot, μ is the fold change values versus the mean expression
#' and δ is the standard deviation. Sm is the mean fold change value of a specific marker gene set
#' and  m is the size of a given marker gene set.
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
#' @details First a new rank is calculated as R = (R1*R2)^(1/2), where R1 is the rank of
#' fold-change for each gene in each spot and R2 is the rank of each marker in each cell type.
#' The Rank-Biased Precision is then calculated as: RBP = (1 - 0.99) * (0.99)^(R - 1)
#' and the final enrichment score is then calculated as the sum of top 100 RBPs.
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
                                 output_enrichment = c('original', 'zscore')) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # output enrichment
  output_enrichment = match.arg(output_enrichment, choices = c('original', 'zscore'))

  # calculate mean gene expression
  if(reverse_log_scale == TRUE) {
    expr_values = log(logbase^expr_values-1)+1
  }

  expCutoff = (rowMeans(expr_values, dims = 1)) * 2
  expbinary = ifelse(expr_values > expCutoff, 1, 0)

  markerGenes = rownames(sign_matrix)
  expbinaryOverlap = expbinary[markerGenes,]
  total = length(markerGenes)
  enrichment = matrix(data=NA,
                      nrow = dim(sign_matrix)[2],
                      ncol=dim(expbinaryOverlap)[2])

  for (i in (1:dim(sign_matrix)[2])){
    signames = rownames(sign_matrix)[which(sign_matrix[,i]==1)]
    vectorX  = NULL

    for (j in(1:dim(expbinaryOverlap)[2])){

      cellsiggene = names(expbinaryOverlap[which(expbinaryOverlap[,j]==1),j])
      x = length(intersect(cellsiggene,signames))
      m = length(rownames(sign_matrix)[which(sign_matrix[,i]==1)])
      n = total-m
      k = length(intersect(cellsiggene, markerGenes))
      enrich = (0-log10(stats::dhyper(x, m, n, k, log = FALSE)))
      vectorX = append(vectorX,enrich)
    }
    enrichment[i,] = vectorX
  }

  rownames(enrichment) = colnames(sign_matrix)
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


  } else if(enrich_method == 'rank') {

    enrich_results_DT =  rankEnrich(gobject = gobject,
                                    sign_matrix = sign_matrix,
                                    expression_values = expression_values,
                                    reverse_log_scale = reverse_log_scale,
                                    logbase = logbase,
                                    output_enrichment = output_enrichment)
    # default name for page enrichment
    if(is.null(name)) name = 'rank'


  } else if(enrich_method == 'hypergeometric'){

    enrich_results_DT =  hyperGeometricEnrich(gobject = gobject,
                                              sign_matrix = sign_matrix,
                                              expression_values = expression_values,
                                              reverse_log_scale = reverse_log_scale,
                                              logbase = logbase,
                                              output_enrichment = output_enrichment)
    # default name for page enrichment
    if(is.null(name)) name = 'hypergeometric'
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
                                       'output enrichment scores' = output_enrichment)
    gobject@parameters = parameters_list

    gobject@spatial_enrichment[[name]] = enrich_results_DT


    return(gobject)

  } else {
    return(enrich_results_DT)
  }
}



