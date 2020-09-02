
#' @title calculateHVG
#' @name calculateHVG
#' @description compute highly variable genes
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param method method to calculate highly variable genes
#' @param reverse_log_scale reverse log-scale of expression values (default = FALSE)
#' @param logbase if reverse_log_scale is TRUE, which log base was used?
#' @param expression_threshold expression threshold to consider a gene detected
#' @param nr_expression_groups number of expression groups for cov_groups
#' @param zscore_threshold zscore to select hvg for cov_groups
#' @param HVGname name for highly variable genes in cell metadata
#' @param difference_in_cov minimum difference in coefficient of variance required
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object highly variable genes appended to gene metadata (fDataDT)
#' @details
#' Currently we provide 2 ways to calculate highly variable genes:
#'
#' \bold{1. high coeff of variance (COV) within groups: } \cr
#' First genes are binned (\emph{nr_expression_groups}) into average expression groups and
#' the COV for each gene is converted into a z-score within each bin. Genes with a z-score
#' higher than the threshold (\emph{zscore_threshold}) are considered highly variable.  \cr
#'
#' \bold{2. high COV based on loess regression prediction: } \cr
#' A predicted COV is calculated for each gene using loess regression (COV~log(mean expression))
#' Genes that show a higher than predicted COV (\emph{difference_in_cov}) are considered highly variable. \cr
#'
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # loads existing Giotto object
#'
#' # update a giotto object
#' mini_giotto_single_cell <- calculateHVG(gobject = mini_giotto_single_cell,
#'                                         zscore_threshold = 0.1,
#'                                         nr_expression_groups = 3)
#'
#' # return a data.table with the high variable genes annotated
#' hvg_dt <- calculateHVG(gobject = mini_giotto_single_cell,
#'                        zscore_threshold = 0.1, nr_expression_groups = 3,
#'                        return_plot = FALSE, return_gobject = FALSE)
#'
#' # return the ggplot object
#' hvg_plot <- calculateHVG(gobject = mini_giotto_single_cell,
#'                        zscore_threshold = 0.1, nr_expression_groups = 3,
#'                        return_plot = TRUE, return_gobject = FALSE)
#'
#'
calculateHVG <- function(gobject,
                         expression_values = c('normalized', 'scaled', 'custom'),
                         method = c('cov_groups','cov_loess'),
                         reverse_log_scale = FALSE,
                         logbase = 2,
                         expression_threshold = 0,
                         nr_expression_groups = 20,
                         zscore_threshold = 1.5,
                         HVGname = 'hvg',
                         difference_in_cov = 0.1,
                         show_plot = NA,
                         return_plot = NA,
                         save_plot = NA,
                         save_param = list(),
                         default_save_name = 'HVGplot',
                         return_gobject = TRUE) {

  sd = cov = mean_expr = gini = cov_group_zscore = selected = cov_diff = pred_cov_genes = genes = NULL

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # method to use
  method = match.arg(method, choices = c('cov_groups', 'cov_loess'))

  # not advised
  if(reverse_log_scale == TRUE) {
    expr_values = (logbase^expr_values)-1
  }

  ## create data.table with relevant statistics ##
  gene_in_cells_detected <- data.table::data.table(genes = rownames(expr_values),
                                       nr_cells = rowSums_giotto(expr_values > expression_threshold),
                                       total_expr = rowSums_giotto(expr_values),
                                       mean_expr = rowMeans_giotto(expr_values),
                                       sd = unlist(apply(expr_values, 1, sd)))
  gene_in_cells_detected[, cov := (sd/mean_expr)]
  gini_level <- unlist(apply(expr_values, MARGIN = 1, mygini_fun))
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

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                                                axis.text = ggplot2::element_text(size = 12))
    pl <- pl + ggplot2::geom_point(data = gene_in_cells_detected, ggplot2::aes(x = mean_expr, y = cov, color = selected))
    pl <- pl + ggplot2::scale_color_manual(values = c(no = 'lightgrey', yes = 'orange'), guide = ggplot2::guide_legend(title = 'HVG',
                                                                                                     override.aes = list(size=5)))
    pl <- pl + ggplot2::facet_wrap(~expr_groups, ncol = nr_expression_groups, scales = 'free_x')
    pl <- pl + ggplot2::theme(axis.text.x = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 4))
    pl <- pl + ggplot2::labs(x = 'expression groups', y = 'cov')

  } else {


    if(method == 'cov_loess') {
      loess_formula = paste0('cov~log(mean_expr)')
      var_col = 'cov'
    }

    loess_model_sample = stats::loess(loess_formula, data = gene_in_cells_detected)
    gene_in_cells_detected$pred_cov_genes = stats::predict(loess_model_sample, newdata = gene_in_cells_detected)
    gene_in_cells_detected[, cov_diff := get(var_col)-pred_cov_genes, by = 1:nrow(gene_in_cells_detected)]
    data.table::setorder(gene_in_cells_detected, -cov_diff)
    gene_in_cells_detected[, selected := ifelse(cov_diff > difference_in_cov, 'yes', 'no')]


    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                                                axis.text = ggplot2::element_text(size = 12))
    pl <- pl + ggplot2::geom_point(data = gene_in_cells_detected, ggplot2::aes_string(x = 'log(mean_expr)', y = var_col, color = 'selected'))
    pl <- pl + ggplot2::geom_line(data = gene_in_cells_detected, ggplot2::aes_string(x = 'log(mean_expr)', y = 'pred_cov_genes'), color = 'blue')
    hvg_line = paste0('pred_cov_genes+',difference_in_cov)
    pl <- pl + ggplot2::geom_line(data = gene_in_cells_detected, ggplot2::aes_string(x = 'log(mean_expr)', y = hvg_line), linetype = 2)
    pl <- pl + ggplot2::labs(x = 'log(mean expression)', y = var_col)
    pl <- pl + ggplot2::scale_color_manual(values = c(no = 'lightgrey', yes = 'orange'), guide = ggplot2::guide_legend(title = 'HVG',
                                                                                                     override.aes = list(size=5)))


  }


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    if(return_gobject == TRUE) {
      cat('return_plot = TRUE and return_gobject = TRUE \n
          plot will not be returned to object, but can still be saved with save_plot = TRUE or manually \n')
    } else {
      return(pl)
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
    data.table::setnames(HVGgenes, 'selected', HVGname)
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
                                       'difference in variance' = difference_in_cov,
                                       'name for HVGs' = HVGname)
    gobject@parameters = parameters_list

    return(gobject)

  } else {
    return(gene_in_cells_detected)
  }

  }




