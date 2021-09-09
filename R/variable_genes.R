
#' @name calc_cov_group_HVF
#' @keywords internal
calc_cov_group_HVF = function(feat_in_cells_detected,
                              nr_expression_groups = 20,
                              zscore_threshold = 1,
                              show_plot = NA,
                              return_plot = NA,
                              save_plot = NA) {


  steps = 1/nr_expression_groups
  prob_sequence = seq(0, 1, steps)
  prob_sequence[length(prob_sequence)] = 1
  expr_group_breaks = stats::quantile(feat_in_cells_detected$mean_expr, probs = prob_sequence)

  ## remove zero's from cuts if there are too many and make first group zero
  if(any(duplicated(expr_group_breaks))) {
    m_expr_vector = feat_in_cells_detected$mean_expr
    expr_group_breaks = stats::quantile(m_expr_vector[m_expr_vector > 0], probs = prob_sequence)
    expr_group_breaks[[1]] = 0
  }

  expr_groups = cut(x = feat_in_cells_detected$mean_expr, breaks = expr_group_breaks,
                    labels = paste0('group_', 1:nr_expression_groups), include.lowest = T)
  feat_in_cells_detected[, expr_groups := expr_groups]
  feat_in_cells_detected[, cov_group_zscore := scale(cov), by =  expr_groups]
  feat_in_cells_detected[, selected := ifelse(cov_group_zscore > zscore_threshold, 'yes', 'no')]

  if(show_plot == TRUE | return_plot == TRUE | save_plot == TRUE) {
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                                                         axis.text = ggplot2::element_text(size = 12))
    pl <- pl + ggplot2::geom_point(data = feat_in_cells_detected, ggplot2::aes(x = mean_expr, y = cov, color = selected))
    pl <- pl + ggplot2::scale_color_manual(values = c(no = 'lightgrey', yes = 'orange'), guide = ggplot2::guide_legend(title = 'HVF',
                                                                                                                       override.aes = list(size=5)))
    pl <- pl + ggplot2::facet_wrap(~expr_groups, ncol = nr_expression_groups, scales = 'free_x')
    pl <- pl + ggplot2::theme(axis.text.x = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 4))
    pl <- pl + ggplot2::labs(x = 'expression groups', y = 'cov')

    return(list(dt = feat_in_cells_detected, pl = pl))
  } else {
    return(list(dt = feat_in_cells_detected))
  }

}



#' @name calc_cov_loess_HVF
#' @keywords internal
calc_cov_loess_HVF = function(feat_in_cells_detected,
                              difference_in_cov = 0.1,
                              show_plot = NA,
                              return_plot = NA,
                              save_plot = NA) {

  # create loess regression
  loess_formula = paste0('cov~log(mean_expr)')
  var_col = 'cov'

  loess_model_sample = stats::loess(loess_formula, data = feat_in_cells_detected)
  feat_in_cells_detected$pred_cov_feats = stats::predict(loess_model_sample, newdata = feat_in_cells_detected)
  feat_in_cells_detected[, cov_diff := get(var_col)-pred_cov_feats, by = 1:nrow(feat_in_cells_detected)]
  data.table::setorder(feat_in_cells_detected, -cov_diff)
  feat_in_cells_detected[, selected := ifelse(cov_diff > difference_in_cov, 'yes', 'no')]

  if(show_plot == TRUE | return_plot == TRUE | save_plot == TRUE) {
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                                                         axis.text = ggplot2::element_text(size = 12))
    pl <- pl + ggplot2::geom_point(data = feat_in_cells_detected, ggplot2::aes_string(x = 'log(mean_expr)', y = var_col, color = 'selected'))
    pl <- pl + ggplot2::geom_line(data = feat_in_cells_detected, ggplot2::aes_string(x = 'log(mean_expr)', y = 'pred_cov_feats'), color = 'blue')
    hvg_line = paste0('pred_cov_feats+',difference_in_cov)
    pl <- pl + ggplot2::geom_line(data = feat_in_cells_detected, ggplot2::aes_string(x = 'log(mean_expr)', y = hvg_line), linetype = 2)
    pl <- pl + ggplot2::labs(x = 'log(mean expression)', y = var_col)
    pl <- pl + ggplot2::scale_color_manual(values = c(no = 'lightgrey', yes = 'orange'),
                                           guide = ggplot2::guide_legend(title = 'HVF',
                                                                         override.aes = list(size=5)))

    return(list(dt = feat_in_cells_detected, pl = pl))
  } else {
    return(list(dt = feat_in_cells_detected))
  }
}


#' @name calc_var_HVF
#' @keywords internal
calc_var_HVF = function(scaled_matrix,
                        var_threshold = 1.5,
                        var_number = NULL,
                        show_plot = NA,
                        return_plot = NA,
                        save_plot = NA) {

  test = apply(X = scaled_matrix, MARGIN = 1, FUN = function(x) var(x))
  test = sort(test, decreasing = T)

  dt_res = data.table::data.table(feats = names(test), var = test)

  if(!is.null(var_number) & is.numeric(var_number)) {
    dt_res[, selected := 1:.N]
    dt_res[, selected := ifelse(selected <= var_number, 'yes', 'no')]
  } else {
    dt_res[, selected := ifelse(var >= var_threshold, 'yes', 'no')]
  }


  if(show_plot == TRUE | return_plot == TRUE | save_plot == TRUE) {

    dt_res[, rank := 1:.N]

    pl = ggplot2::ggplot()
    pl = pl + ggplot2::geom_point(data = dt_res, aes(x = rank, y = var, color = selected))
    pl = pl + ggplot2::scale_x_reverse()
    pl = pl + ggplot2::theme_classic() + ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                                                         axis.text = ggplot2::element_text(size = 12))
    pl = pl + ggplot2::scale_color_manual(values = c(no = 'lightgrey', yes = 'orange'),
                                          guide = ggplot2::guide_legend(title = 'HVF',
                                                                        override.aes = list(size=5)))
    pl = pl + ggplot2::labs(x = 'feature rank', y = 'variance')

    dt_res_final = data.table::copy(dt_res)
    dt_res_final[, rank := NULL]

    return(list(dt = dt_res_final, pl = pl))

  } else {

    return(list(dt = dt_res))
  }

}

#' @title calculateHVF
#' @name calculateHVF
#' @description compute highly variable features
#' @param gobject giotto object
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param method method to calculate highly variable features
#' @param reverse_log_scale reverse log-scale of expression values (default = FALSE)
#' @param logbase if reverse_log_scale is TRUE, which log base was used?
#' @param expression_threshold expression threshold to consider a gene detected
#' @param nr_expression_groups [cov_groups] number of expression groups for cov_groups
#' @param zscore_threshold [cov_groups] zscore to select hvg for cov_groups
#' @param HVFname name for highly variable features in cell metadata
#' @param difference_in_cov [cov_loess] minimum difference in coefficient of variance required
#' @param var_threshold [var_p_resid] variance threshold for features for var_p_resid method
#' @param var_number [var_p_resid] number of top variance features for var_p_resid method
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object highly variable features appended to feature metadata (fDataDT)
#' @details
#' Currently we provide 2 ways to calculate highly variable genes:
#'
#' \bold{1. high coeff of variance (COV) within groups: } \cr
#' First genes are binned (\emph{nr_expression_groups}) into average expression groups and
#' the COV for each feature is converted into a z-score within each bin. Features with a z-score
#' higher than the threshold (\emph{zscore_threshold}) are considered highly variable.  \cr
#'
#' \bold{2. high COV based on loess regression prediction: } \cr
#' A predicted COV is calculated for each feature using loess regression (COV~log(mean expression))
#' Features that show a higher than predicted COV (\emph{difference_in_cov}) are considered highly variable. \cr
#'
#' @export
calculateHVF <- function(gobject,
                         feat_type = NULL,
                         expression_values = c('normalized', 'scaled', 'custom'),
                         method = c('cov_groups','cov_loess', 'var_p_resid'),
                         reverse_log_scale = FALSE,
                         logbase = 2,
                         expression_threshold = 0,
                         nr_expression_groups = 20,
                         zscore_threshold = 1.5,
                         HVFname = 'hvf',
                         difference_in_cov = 0.1,
                         var_threshold = 1.5,
                         var_number = NULL,
                         show_plot = NA,
                         return_plot = NA,
                         save_plot = NA,
                         save_param = list(),
                         default_save_name = 'HVFplot',
                         return_gobject = TRUE) {

  # set data.table variables to NULL
  sd = cov = mean_expr = gini = cov_group_zscore = selected = cov_diff = pred_cov_feats = feats = NULL

  # specify feat_type
  if(is.null(feat_type)) {
    feat_type = gobject@expression_feat[[1]]
  }

  # expression values to be used
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject, feat_type = feat_type, values = values)

  # not advised
  if(reverse_log_scale == TRUE) {
    expr_values = (logbase^expr_values)-1
  }


  # method to use
  method = match.arg(method, choices = c('cov_groups', 'cov_loess', 'var_p_resid'))


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)



  if(method == 'var_p_resid') {

    results = calc_var_HVF(scaled_matrix = expr_values,
                           var_threshold = var_threshold,
                           var_number = var_number,
                           show_plot = show_plot,
                           return_plot = return_plot,
                           save_plot = save_plot)

    feat_in_cells_detected = results[['dt']]
    pl = results[['pl']]

  } else {

    ## create data.table with relevant statistics ##
    feat_in_cells_detected <- data.table::data.table(feats = rownames(expr_values),
                                                     nr_cells = rowSums_flex(expr_values > expression_threshold),
                                                     total_expr = rowSums_flex(expr_values),
                                                     mean_expr = rowMeans_flex(expr_values),
                                                     sd = unlist(apply(expr_values, 1, sd)))
    feat_in_cells_detected[, cov := (sd/mean_expr)]
    gini_level <- unlist(apply(expr_values, MARGIN = 1, Giotto:::mygini_fun))
    feat_in_cells_detected[, gini := gini_level]


    if(method == 'cov_groups') {

      results = calc_cov_group_HVF(feat_in_cells_detected = feat_in_cells_detected,
                                   nr_expression_groups = nr_expression_groups,
                                   zscore_threshold = zscore_threshold,
                                   show_plot = show_plot,
                                   return_plot = return_plot,
                                   save_plot = save_plot)
      feat_in_cells_detected = results[['dt']]
      pl = results[['pl']]

    } else if(method == 'cov_loess') {

      results = calc_cov_loess_HVF(feat_in_cells_detected = feat_in_cells_detected,
                                   difference_in_cov = difference_in_cov,
                                   show_plot = show_plot,
                                   return_plot = return_plot,
                                   save_plot = save_plot)
      feat_in_cells_detected = results[['dt']]
      pl = results[['pl']]
    }


  }




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

    # add HVG metadata to feat_metadata
    feat_metadata = gobject@feat_metadata[[feat_type]]

    column_names_feat_metadata = colnames(feat_metadata)
    if(HVFname %in% column_names_feat_metadata) {
      cat('\n ', HVFname, ' has already been used, will be overwritten \n')
      feat_metadata[, eval(HVFname) := NULL]
      gobject@feat_metadata[[feat_type]] = feat_metadata
    }

    if(method == 'var_p_resid') {
      HVGfeats = feat_in_cells_detected[,.(feats, var, selected)]
      data.table::setnames(HVGfeats, 'selected', HVFname)
    } else {
      HVGfeats = feat_in_cells_detected[,.(feats, selected)]
      data.table::setnames(HVGfeats, 'selected', HVFname)
    }


    gobject = addFeatMetadata(gobject = gobject,
                              feat_type = feat_type,
                              new_metadata = HVGfeats,
                              by_column = TRUE,
                              column_feat_ID = 'feats')

    ## update parameters used ##
    gobject = update_giotto_params(gobject, description = '_hvf')

    return(gobject)

  } else {
    return(feat_in_cells_detected)
  }

}







#' @name calculateHVG
#' @description compute highly variable genes
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param method method to calculate highly variable genes
#' @param reverse_log_scale reverse log-scale of expression values (default = FALSE)
#' @param logbase if reverse_log_scale is TRUE, which log base was used?
#' @param expression_threshold expression threshold to consider a gene detected
#' @param nr_expression_groups [cov_groups] number of expression groups for cov_groups
#' @param zscore_threshold [cov_groups] zscore to select hvg for cov_groups
#' @param HVGname name for highly variable genes in cell metadata
#' @param difference_in_cov [cov_loess] minimum difference in coefficient of variance required
#' @param var_threshold [var_p_resid] variance threshold for features for var_p_resid method
#' @param var_number [var_p_resid] number of top variance features for var_p_resid method
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
calculateHVG <- function(gobject,
                         expression_values = c('normalized', 'scaled', 'custom'),
                         method = c('cov_groups','cov_loess', 'var_p_resid'),
                         reverse_log_scale = FALSE,
                         logbase = 2,
                         expression_threshold = 0,
                         nr_expression_groups = 20,
                         zscore_threshold = 1.5,
                         HVGname = 'hvf',
                         difference_in_cov = 0.1,
                         var_threshold = 1.5,
                         var_number = NULL,
                         show_plot = NA,
                         return_plot = NA,
                         save_plot = NA,
                         save_param = list(),
                         default_save_name = 'HVGplot',
                         return_gobject = TRUE) {


  warning("Deprecated and replaced by calculateHVF")

  result = calculateHVF(gobject = gobject,
                        feat_type = NULL,
                        expression_values = expression_values,
                        method = method,
                        reverse_log_scale = reverse_log_scale,
                        logbase = logbase,
                        expression_threshold = expression_threshold,
                        nr_expression_groups = nr_expression_groups,
                        zscore_threshold = zscore_threshold,
                        HVFname = HVGname,
                        difference_in_cov = difference_in_cov,
                        var_threshold = var_threshold,
                        var_number = var_number,
                        show_plot = show_plot,
                        return_plot = return_plot,
                        save_plot = save_plot,
                        save_param = save_param,
                        default_save_name = default_save_name,
                        return_gobject = return_gobject)

  return(result)
}





