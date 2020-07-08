
#' @title createSpatialDeconvolution
#' @description Function to deconvolution based on single cell expression data
#' @param sign_matrix sig matrix for deconvolution
#' @param method list of genes (signature)
#' @return matrix
#' @export
#' @examples
#' createSpatialDeconvolution()

createSpatialDeconvolution<-function(gobject,
                                     sign_matrix,
                                     n_cell = 50){

  #########################################generate sig
  ######transform expression data to no log data
  nolog_expr<-2^(gobject@norm_expr)-1
  cluster<-gobject@cell_metadata$leiden_clus
  #####getting overlapped gene lists
  intersect_gene<-intersect(rownames(sign_matrix),rownames(nolog_expr))
  filter_Sig<-sign_matrix[intersect_gene,]
  filter_expr<-nolog_expr[intersect_gene,]
  filter_log_expr<-gobject@norm_expr[intersect_gene,]
  #####first round spatial deconvolution ##spot or cluster
  enrich_spot_proportion<-enrich_deconvolution(expr = filter_expr,
                                                    log_expr = filter_log_expr,
                                                    cluster_info = cluster,
                                                    ct_exp=filter_Sig,
                                                    cutoff=2)
  ####re-deconvolution based on spatial resolution
  resolution<-(1/n_cell)
  binarize_proportion<-ifelse(enrich_spot_proportion>=resolution, 1, 0)
  spot_proportion<-spot_deconvolution(expr = filter_expr,
                                      cluster_info = cluster,
                                      ct_exp=filter_Sig,
                                      binary_matrix = binarize_proportion)
  deconvolutionDT = data.table::data.table(cell_ID = colnames(spot_proportion))
  deconvolutionDT = cbind(deconvolutionDT, as.data.table(t(spot_proportion)))
  gobject@spatial_deconvolution = deconvolutionDT
  return(gobject)
}

#' @title enrich_deconvolution
enrich_deconvolution<-function(expr,log_expr,cluster_info,ct_exp,cutoff){
  #####generate enrich 0/1 matrix based on expression matrix
  enrich_matrix<-matrix(0,nrow=dim(ct_exp)[1],ncol=dim(ct_exp)[2])
  rowmax_col<-Rfast::rowMaxs(ct_exp)
  for (i in 1:length(rowmax_col)){
    enrich_matrix[i,rowmax_col[i]]=1
  }
  rownames(enrich_matrix)<-rownames(ct_exp)
  colnames(enrich_matrix)<-colnames(ct_exp)
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
spot_deconvolution<-function(expr,cluster_info,ct_exp,binary_matrix){
  #####generate enrich 0/1 matrix based on expression matrix
  enrich_matrix<-matrix(0,nrow=dim(ct_exp)[1],ncol=dim(ct_exp)[2])
  rowmax_col<-Rfast::rowMaxs(ct_exp)
  for (i in 1:length(rowmax_col)){
    enrich_matrix[i,rowmax_col[i]]=1
  }
  rownames(enrich_matrix)<-rownames(ct_exp)
  colnames(enrich_matrix)<-colnames(ct_exp)
  #########################################################
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
      solution_all_exp<-solveOLSInternal(select_sig_exp,all_exp)
      constant_J<-findDampeningConstant(select_sig_exp,all_exp,solution_all_exp)
      ######deconvolution for each spot
      for(k in 1:(dim(cluster_cell_exp)[2])){
        B<-as.matrix(cluster_cell_exp[,k])
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
          S_k<-as.matrix(ct_exp[uniq_ct_k_gene,ct_spot_k])
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
cluster_enrich_analysis <- function(exp_matrix,cluster_info,enrich_sig_matrix) {
  uniq_cluster<-sort(unique(cluster_info))
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
enrich_analysis <- function(expr_values,sign_matrix) {
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
  return(enrichment)
}

#########################dwls
#' @title optimize_deconvolute_dwls
optimize_deconvolute_dwls <- function(exp,Signature) {
  ######overlap signature with spatial genes
  Genes<-intersect(rownames(Signature),rownames(exp))
  S<-Signature[Genes,]
  S<-as.matrix(S)
  Bulk<-as.matrix(exp)
  subBulk = Bulk[Genes,]
  allCounts_DWLS<-NULL
  all_exp<-rowMeans(exp)
  solution_all_exp<-solveOLSInternal(S,all_exp[Genes])
  constant_J<-findDampeningConstant(S,all_exp[Genes],solution_all_exp)
  #print(constant_J)
  for(j in 1:(dim(subBulk)[2])){
    B<-subBulk[,j]
    solDWLS<-optimize_solveDampenedWLS(S,B,constant_J)
    allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)
  }
  colnames(allCounts_DWLS)<-colnames(exp)
  return(allCounts_DWLS)
}

#' @title optimize_solveDampenedWLS
optimize_solveDampenedWLS<-function(S,B,constant_J){
  #first solve OLS, use this solution to find a starting point for the weights
  solution<-solveOLSInternal(S,B)
  #now use dampened WLS, iterate weights until convergence
  iterations<-0
  changes<-c()
  #find dampening constant for weights using cross-validation
  j<-constant_J
  change<-1
  while(change>.01 & iterations<1000){
    newsolution<-solveDampenedWLSj(S,B,solution,j)
    #decrease step size for convergence
    solutionAverage<-rowMeans(cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
    change<-norm(as.matrix(solutionAverage-solution))
    solution<-solutionAverage
    iterations<-iterations+1
    changes<-c(changes,change)
  }
  #print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}

#find a dampening constant for the weights using cross-validation
#' @title findDampeningConstant
findDampeningConstant<-function(S,B,goldStandard){
  solutionsSd<-NULL
  #goldStandard is used to define the weights
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsScaledMinusInf<-wsScaled
  #ignore infinite weights
  if(max(wsScaled)=="Inf"){
    wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
  }
  #try multiple values of the dampening constant (multiplier)
  #for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
    multiplier<-1*2^(j-1)
    wsDampened<-wsScaled
    wsDampened[which(wsScaled>multiplier)]<-multiplier
    solutions<-NULL
    seeds<-c(1:100)
    for (i in 1:100){
      set.seed(seeds[i]) #make nondeterministic
      subset<-sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
      #solve dampened weighted least squares for subset
      fit = stats::lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
      sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
      solutions<-cbind(solutions,sol)
    }
    solutionsSd<-cbind(solutionsSd,apply(solutions,1,sd))
  }
  #choose dampening constant that results in least cross-validation variance
  j<-which.min(colMeans(solutionsSd^2))
  return(j)
}

########basic functions for dwls
#' @title solveOLSInternal
solveOLSInternal<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  solution<-quadprog::solve.QP(D,d,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#solve WLS given a dampening constant
#' @title solveDampenedWLSj
solveDampenedWLSj<-function(S,B,goldStandard,j){
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
