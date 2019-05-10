
#' @title calculateHVG
#' @name calculateHVG
#' @description compute highly variable genes
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param method method to calculate highly variable genes
#' @param reverse_log_scale reverse log-scale of expression values
#' @param logbase if reverse_log_scale is TRUE, which log base was used?
#' @param expression_threshold expression threshold to consider a gene detected
#' @param HVGname name for highly variable genes in cell metadata
#' @param difference_in_variance minimum difference in variance required
#' @param show_plot show plots
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object highly variable genes appended to gene metadata (fDataDT)
#' @details Description of how we compute highly variable genes.
#' @export
#' @examples
#'     calculateHVG(gobject)
calculateHVG <- function(gobject,
                         expression_values = c('normalized', 'scaled', 'custom'),
                         method = c('cov', 'gini'),
                         reverse_log_scale = T,
                         logbase = 2,
                         expression_threshold = 0,
                         HVGname = 'hvg',
                         difference_in_variance = 1,
                         show_plot = T,
                         return_gobject = T) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  if(reverse_log_scale == TRUE) {
    expr_values = (logbase^expr_values)-1
  }

  ## create data.table with relevant statistics ##
  gene_in_cells_detected <- data.table(genes = rownames(expr_values),
                                       nr_cells = rowSums(expr_values > expression_threshold),
                                       total_expr = rowSums(expr_values),
                                       mean_expr = rowMeans(expr_values))

  cov_level <- unlist(apply(expr_values, MARGIN = 1, FUN = function(x) {

    cov = sd(x)/mean(x)
    return(cov)

  }))
  gene_in_cells_detected[, cov := cov_level]

  gini_level <- unlist(apply(expr_values, MARGIN = 1, FUN = function(x) {

    gini = mygini_fun(x)
    return(gini)

  }))
  gene_in_cells_detected[, gini := gini_level]


  # predict variability with mean expression
  method = match.arg(method)
  if(method == 'cov') {
    loess_formula = paste0('cov~mean_expr')
  } else if(method == 'gini') {
    loess_formula = paste0('gini~mean_expr')
  }

  loess_model_sample <- loess(loess_formula, data = gene_in_cells_detected)
  gene_in_cells_detected$pred_var_genes <- predict(loess_model_sample, newdata = gene_in_cells_detected)
  gene_in_cells_detected[, var_diff := get(method)-pred_var_genes, by = 1:nrow(gene_in_cells_detected)]
  setorder(gene_in_cells_detected, -var_diff)

  # select HVG
  HVG <- gene_in_cells_detected[var_diff > difference_in_variance]$genes
  gene_in_cells_detected[, selected := ifelse(genes %in% HVG, 'yes', 'no')]
  print(table(gene_in_cells_detected$selected))



  if(show_plot == TRUE) {
    pl <- ggplot()
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + geom_point(data = gene_in_cells_detected, aes_string(x = 'mean_expr', y = method, color = 'selected'))
    pl <- pl + geom_line(data = gene_in_cells_detected, aes_string(x = 'mean_expr', y = 'pred_var_genes'), color = 'blue')
    hvg_line = paste0('pred_var_genes+',difference_in_variance)
    pl <- pl + geom_line(data = gene_in_cells_detected, aes_string(x = 'mean_expr', y = hvg_line), linetype = 2)
    pl <- pl + labs(x = 'mean expression', y = method)
    print(pl)
  }


  if(return_gobject == TRUE) {

    # add HVG metadata to gene_metadata
    gene_metadata = gobject@gene_metadata
    column_names_gene_metadata = colnames(gene_metadata)
    if(HVGname %in% column_names_gene_metadata) {
      cat('\n ', HVGname, ' has already been used, will be overwritten \n')
      gene_metadata[, eval(HVGname) := NULL]
      gobject@gene_metadata = gene_metadata
    }

    HVGgenes = gene_in_cells_detected[,.(genes, selected)]
    setnames(HVGgenes, 'selected', HVGname)
    gobject <- addGeneMetadata(gobject = gobject, new_metadata = HVGgenes, by_column = T, column_gene_ID = 'genes')



    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_hvg')
    # parameters to include
    parameters_list[[update_name]] = c('method used' = method,
                                       'expression values' = expression_values,
                                       'reversed log scale' = reverse_log_scale,
                                       'logbase' = logbase,
                                       'expression threshold' = expression_threshold,
                                       'difference in variance' = difference_in_variance,
                                       'name for HVGs' = HVGname)
    gobject@parameters = parameters_list

    return(gobject)

  } else {
    return(gene_in_cells_detected)
  }

}

