
#' @title calculateHVG
#' @name calculateHVG
#' @description compute highly variable genes
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param method method to calculate highly variable genes
#' @param reverse_log_scale reverse log-scale of expression values
#' @param logbase if reverse_log_scale is TRUE, which log base was used?
#' @param expression_threshold expression threshold to consider a gene detected
#' @param nr_expression_groups number of expression groups for cov_groups
#' @param zscore_threshold zscore to select hvg for cov_groups
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
                         method = c('cov_groups','cov_loess','gini_loess'),
                         reverse_log_scale = T,
                         logbase = 2,
                         expression_threshold = 0,
                         nr_expression_groups = 20,
                         zscore_threshold = 1.5,
                         HVGname = 'hvg',
                         difference_in_variance = 1,
                         show_plot = T,
                         return_gobject = T) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # method to use
  method = match.arg(method, choices = c('cov_groups', 'cov_loess', 'gini_loess'))

  if(reverse_log_scale == TRUE) {
    expr_values = (logbase^expr_values)-1
  }

  ## create data.table with relevant statistics ##
  gene_in_cells_detected <- data.table(genes = rownames(expr_values),
                                       nr_cells = rowSums(expr_values > expression_threshold),
                                       total_expr = rowSums(expr_values),
                                       mean_expr = rowMeans(expr_values),
                                       sd = unlist(apply(expr_values, 1, sd)))
  gene_in_cells_detected[, cov := (sd/mean_expr)]
  gini_level <- unlist(apply(expr_values, MARGIN = 1, FUN = function(x) {

    gini = Giotto:::mygini_fun(x)
    return(gini)

  }))
  gene_in_cells_detected[, gini := gini_level]


  if(method == 'cov_groups') {

    steps = 1/nr_expression_groups
    prob_sequence = seq(0, 1, steps)
    prob_sequence[length(prob_sequence)] = 1
    expr_group_breaks = stats::quantile(gene_in_cells_detected$mean_expr, probs = prob_sequence)

    ## remove zero's from cuts if there are too many and make first group zero
    if(any(duplicated(expr_group_breaks))) {
      m_expr_vector = gene_in_cells_detected$mean_expr
      expr_group_breaks = stats::quantile(m_expr_vector[m_expr_vector > 0], probs = prob_sequence)
      expr_group_breaks[[1]] = 0
    }

    expr_groups = cut(x = gene_in_cells_detected$mean_expr, breaks = expr_group_breaks,
                      labels = paste0('group_', 1:nr_expression_groups), include.lowest = T)
    gene_in_cells_detected[, expr_groups := expr_groups]
    gene_in_cells_detected[, cov_group_zscore := scale(cov), by =  expr_groups]
    gene_in_cells_detected[, selected := ifelse(cov_group_zscore > zscore_threshold, 'yes', 'no')]

    if(show_plot == TRUE) {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_point(data = gene_in_cells_detected, aes(x = mean_expr, y = cov, color = selected))
      pl <- pl + ggplot2::facet_wrap(~expr_groups, ncol = nr_expression_groups, scales = 'free_x')
      pl <- pl + ggplot2::theme(axis.text.x = element_blank(), strip.text = element_text(size = 4))
      pl <- pl + ggplot2::labs(x = 'expression groups', y = 'cov')
      print(pl)
    }
  } else {


    if(method == 'cov_loess') {
      loess_formula = paste0('cov~mean_expr')
      var_col = 'cov'
    } else if(method == 'gini_loess') {
      loess_formula = paste0('gini~mean_expr')
      var_col = 'gini'
    }

    loess_model_sample <- stats::loess(loess_formula, data = gene_in_cells_detected)
    gene_in_cells_detected$pred_var_genes <- predict(loess_model_sample, newdata = gene_in_cells_detected)
    gene_in_cells_detected[, var_diff := get(var_col)-pred_var_genes, by = 1:nrow(gene_in_cells_detected)]
    setorder(gene_in_cells_detected, -var_diff)
    gene_in_cells_detected[, selected := ifelse(var_diff > difference_in_variance, 'yes', 'no')]


    if(show_plot == TRUE) {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_point(data = gene_in_cells_detected, aes_string(x = 'mean_expr', y = var_col, color = 'selected'))
      pl <- pl + ggplot2::geom_line(data = gene_in_cells_detected, aes_string(x = 'mean_expr', y = 'pred_var_genes'), color = 'blue')
      hvg_line = paste0('pred_var_genes+',difference_in_variance)
      pl <- pl + ggplot2::geom_line(data = gene_in_cells_detected, aes_string(x = 'mean_expr', y = hvg_line), linetype = 2)
      pl <- pl + ggplot2::labs(x = 'mean expression', y = var_col)
      print(pl)
    }


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
                                       'number of expression groups ' = nr_expression_groups,
                                       'threshold for z-score ' = zscore_threshold,
                                       'difference in variance' = difference_in_variance,
                                       'name for HVGs' = HVGname)
    gobject@parameters = parameters_list

    return(gobject)

  } else {
    return(gene_in_cells_detected)
  }

}




