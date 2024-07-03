## Giotto auxiliary functions ####




### matrix processing ####

#' @title Mean expression detected test
#' @param mymatrix matrix of expression info
#' @param detection_threshold detection threshold. Defaults to 1 count.
#' @returns numeric
#' @keywords internal
.mean_expr_det_test <- function(mymatrix, detection_threshold = 1) {
    unlist(apply(X = mymatrix, MARGIN = 1, FUN = function(x) {
        detected_x <- x[x > detection_threshold]
        mean(detected_x)
    }))
}

#' @title Normalize expression matrix for library size
#' @param mymatrix matrix object
#' @param scalefactor scalefactor
#' @returns matrix
#' @keywords internal
.lib_norm_giotto <- function(mymatrix, scalefactor) {
    libsizes <- colSums_flex(mymatrix)

    if (any(libsizes == 0)) {
        warning(wrap_txt("Total library size or counts for individual spat
                    units are 0.
                    This will likely result in normalization problems.
                    filter (filterGiotto) or impute (imputeGiotto) spatial
                    units."))
    }

    norm_expr <- t_flex(t_flex(mymatrix) / libsizes) * scalefactor
    return(norm_expr)
}

#' @title Log normalize expression matrix
#' @returns matrix
#' @keywords internal
.log_norm_giotto <- function(mymatrix, base, offset) {
    if (methods::is(mymatrix, "DelayedArray")) {
        mymatrix <- log(mymatrix + offset) / log(base)
        # } else if(methods::is(mymatrix, 'DelayedMatrix')) {
        #   mymatrix = log(mymatrix + offset)/log(base)
    } else if (methods::is(mymatrix, "dgCMatrix")) {
        mymatrix@x <- log(mymatrix@x + offset) / log(base)
        # replace with sparseMatrixStats
    } else if (methods::is(mymatrix, "Matrix")) {
        mymatrix@x <- log(mymatrix@x + offset) / log(base)
    } else {
        mymatrix <- log(as.matrix(mymatrix) + offset) / log(base)
    }

    return(mymatrix)
}












### Filter Values ####




#' @title filterDistributions
#' @name filterDistributions
#' @description show gene or cell distribution after filtering on expression
#' threshold
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param method method to create distribution (see details)
#' @param expression_threshold threshold to consider a gene expressed
#' @param detection consider features (e.g. genes) or cells
#' @param plot_type type of plot
#' @param scale_y scale y-axis (e.g. "log"), NULL = no scaling
#' @param nr_bins number of bins for histogram plot
#' @param fill_color fill color for plots
#' @param scale_axis ggplot transformation for axis (e.g. log2)
#' @param axis_offset offset to be used together with the scaling
#' transformation
#' @param show_plot logical. show plot
#' @param return_plot logical. return ggplot object
#' @param save_plot logical. directly save the plot
#' @param save_param list of saving parameters from
#' [GiottoVisuals::all_plots_save_function]
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @returns ggplot object
#' @details
#' There are 3 ways to create a distribution profile and summarize it for
#' either the features or the cells (spatial units) \cr
#' \itemize{
#'   \item{1. threshold: calculate features that cross a thresold (default)}
#'   \item{2. sum: summarize the features, i.e. total of a feature}
#'   \item{3. mean: calculate mean of the features, i.e. average expression}
#' }
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' filterDistributions(g)
#' @export
filterDistributions <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("raw", "normalized", "scaled", "custom"),
        method = c("threshold", "sum", "mean"),
        expression_threshold = 1,
        detection = c("feats", "cells"),
        plot_type = c("histogram", "violin"),
        scale_y = NULL,
        nr_bins = 30,
        fill_color = "lightblue",
        scale_axis = "identity",
        axis_offset = 0,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "filterDistributions") {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("raw", "normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )

    # plot distribution for feats or cells
    detection <- match.arg(detection, c("feats", "cells"))

    # method to calculate distribution
    method <- match.arg(method, c("threshold", "sum", "mean"))

    # plot type
    plot_type <- match.arg(plot_type, c("histogram", "violin"))

    # variables
    V1 <- NULL

    # for genes
    if (detection == "feats") {
        if (method == "threshold") {
            feat_detection_levels <- data.table::as.data.table(
                rowSums_flex(expr_values >= expression_threshold)
            )
            mytitle <- "feat detected in # of cells"
        } else if (method == "sum") {
            feat_detection_levels <- data.table::as.data.table(
                rowSums_flex(expr_values)
            )
            mytitle <- "total sum of feature detected in all cells"
        } else if (method == "mean") {
            feat_detection_levels <- data.table::as.data.table(
                rowMeans_flex(expr_values)
            )
            mytitle <- "average of feature detected in all cells"
        }

        y_title <- "count"
        if (!is.null(scale_y)) {
            feat_detection_levels[, V1 := do.call(what = scale_y, list(V1))]
            y_title <- paste(scale_y, y_title)
        }



        if (plot_type == "violin") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::geom_violin(
                data = feat_detection_levels,
                ggplot2::aes(x = "feats", y = V1 + axis_offset),
                fill = fill_color
            )
            pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
            pl <- pl + ggplot2::labs(y = mytitle, x = "")
        } else if (plot_type == "histogram") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::geom_histogram(
                data = feat_detection_levels,
                ggplot2::aes(x = V1 + axis_offset),
                color = "white", bins = nr_bins, fill = fill_color
            )
            pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
            pl <- pl + ggplot2::labs(x = mytitle, y = y_title)
        }

        # for cells
    } else if (detection == "cells") {
        if (method == "threshold") {
            cell_detection_levels <- data.table::as.data.table(
                colSums_flex(expr_values >= expression_threshold)
            )
            mytitle <- "feats detected per cell"
        } else if (method == "sum") {
            cell_detection_levels <- data.table::as.data.table(
                colSums_flex(expr_values)
            )
            mytitle <- "total features per cell"
        } else if (method == "mean") {
            cell_detection_levels <- data.table::as.data.table(
                colMeans_flex(expr_values)
            )
            mytitle <- "average number of features per cell"
        }

        y_title <- "count"
        if (!is.null(scale_y)) {
            cell_detection_levels[, V1 := do.call(what = scale_y, list(V1))]
            y_title <- paste(scale_y, y_title)
        }



        if (plot_type == "violin") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::geom_violin(
                data = cell_detection_levels,
                ggplot2::aes(x = "cells", y = V1 + axis_offset),
                fill = fill_color
            )
            pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
            pl <- pl + ggplot2::labs(y = mytitle, x = "")
        } else if (plot_type == "histogram") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::geom_histogram(
                data = cell_detection_levels,
                ggplot2::aes(x = V1 + axis_offset),
                color = "white", bins = nr_bins, fill = fill_color
            )
            pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
            pl <- pl + ggplot2::labs(x = mytitle, y = y_title)
        }
    }

    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}



#' @title filterCombinations
#' @name filterCombinations
#' @description Shows how many genes and cells are lost with combinations of
#' thresholds.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param expression_values expression values to use
#' @param expression_thresholds all thresholds to consider a gene expressed
#' @param feat_det_in_min_cells minimum # of cells that need to express a
#' feature
#' @param min_det_feats_per_cell minimum # of features that need to be
#' detected in a cell
#' @param scale_x_axis ggplot transformation for x-axis (e.g. log2)
#' @param x_axis_offset x-axis offset to be used together with the scaling
#' transformation
#' @param scale_y_axis ggplot transformation for y-axis (e.g. log2)
#' @param y_axis_offset y-axis offset to be used together with the scaling
#' transformation
#' @returns list of data.table and ggplot object
#' @details Creates a scatterplot that visualizes the number of genes and
#' cells that are lost with a specific combination of a gene and cell
#' threshold given an arbitrary cutoff to call a gene expressed. This function
#' can be used to make an informed decision at the filtering step with
#' filterGiotto.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' filterCombinations(g)
#' @export
filterCombinations <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("raw", "normalized", "scaled", "custom"),
        expression_thresholds = c(1, 2),
        feat_det_in_min_cells = c(5, 50),
        min_det_feats_per_cell = c(200, 400),
        scale_x_axis = "identity",
        x_axis_offset = 0,
        scale_y_axis = "identity",
        y_axis_offset = 0,
        show_plot = TRUE,
        return_plot = FALSE,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "filterCombinations") {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )


    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("raw", "normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values
    )[]

    # feat and cell minimums need to have the same length
    if (length(feat_det_in_min_cells) != length(min_det_feats_per_cell)) {
        stop("\n feat_det_in_min_cells and min_det_feats_per_cell need to be
            the same size \n")
    }

    # compute the number of removed feats and cells
    result_list <- list()
    for (thresh_i in seq_along(expression_thresholds)) {
        threshold <- expression_thresholds[thresh_i]

        det_feats_res <- list()
        det_cells_res <- list()
        for (combn_i in seq_along(feat_det_in_min_cells)) {
            min_cells_for_feat <- feat_det_in_min_cells[combn_i]
            min_feats_per_cell <- min_det_feats_per_cell[combn_i]


            # first remove feats
            filter_index_feats <- rowSums_flex(
                expr_values >= threshold
            ) >= min_cells_for_feat
            removed_feats <- length(filter_index_feats[
                filter_index_feats == FALSE
            ])
            det_cells_res[[combn_i]] <- removed_feats

            # then remove cells
            filter_index_cells <- colSums_flex(expr_values[
                filter_index_feats,
            ] >= threshold) >= min_feats_per_cell
            removed_cells <- length(filter_index_cells[
                filter_index_cells == FALSE
            ])
            det_feats_res[[combn_i]] <- removed_cells
        }

        temp_dt <- data.table::data.table(
            "threshold" = threshold,
            removed_feats = unlist(det_cells_res),
            removed_cells = unlist(det_feats_res)
        )

        result_list[[thresh_i]] <- temp_dt
    }

    result_DT <- do.call("rbind", result_list)

    # data.table variables
    feat_detected_in_min_cells <- min_detected_feats_per_cell <-
        combination <- NULL

    result_DT[["feat_detected_in_min_cells"]] <- feat_det_in_min_cells
    result_DT[["min_detected_feats_per_cell"]] <- min_det_feats_per_cell
    result_DT[["combination"]] <- paste0(
        result_DT$feat_detected_in_min_cells, "-",
        result_DT$min_detected_feats_per_cell
    )

    result_DT <- result_DT[, .(
        threshold,
        feat_detected_in_min_cells,
        min_detected_feats_per_cell,
        combination,
        removed_feats,
        removed_cells
    )]

    maximum_x_value <- max(result_DT[["removed_cells"]], na.rm = TRUE)
    maximum_y_value <- max(result_DT[["removed_feats"]], na.rm = TRUE)

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + ggplot2::geom_line(data = result_DT, aes(
        x = removed_cells + x_axis_offset,
        y = removed_feats + y_axis_offset,
        group = as.factor(threshold)
    ), linetype = 2)
    pl <- pl + ggplot2::geom_point(data = result_DT, aes(
        x = removed_cells + x_axis_offset,
        y = removed_feats + y_axis_offset,
        color = as.factor(threshold)
    ))
    pl <- pl + scale_color_discrete(
        guide = guide_legend(title = "threshold(s)")
    )
    pl <- pl + geom_text_repel(data = result_DT, aes(
        x = removed_cells + x_axis_offset,
        y = removed_feats + y_axis_offset,
        label = combination
    ))
    pl <- pl + ggplot2::scale_x_continuous(
        trans = scale_x_axis, limits = c(0, maximum_x_value)
    )
    pl <- pl + ggplot2::scale_y_continuous(
        trans = scale_y_axis, limits = c(0, maximum_y_value)
    )
    pl <- pl + ggplot2::labs(
        x = "number of removed cells", y = "number of removed feats"
    )


    return(plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = list(results = result_DT, ggplot = pl)
    ))
}


#' @title filterGiotto
#' @name filterGiotto
#' @description filter Giotto object based on expression threshold
#' @param gobject giotto object
#' @param spat_unit character. spatial unit. If more than one is provided then
#' the first will be filtered, the filtering results will be applied across the
#' other spat_units provided
#' @param feat_type character. feature type. If more than one is provided then
#' the first will be filtered, the filtering results will be applied across the
#' other feat_types provided.
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param feat_det_in_min_cells minimum # of cells that need to express a
#' feature
#' @param min_det_feats_per_cell minimum # of features that need to be detected
#' in a cell
#' @param all_spat_units deprecated. Use spat_unit_fsub = ":all:"
#' @param all_feat_types deprecated. Use feat_type_ssub = ":all:"
#' @param spat_unit_fsub character vector. (default = ':all:') limit features
#' to remove results to selected spat_units
#' @param feat_type_ssub character vector. (default = ':all:') limit cells to
#' remove results to selected feat_types
#' @param poly_info polygon information to use
#' @param tag_cells tag filtered cells in metadata vs. remove cells
#' @param tag_cell_name column name for tagged cells in metadata
#' @param tag_feats tag features in metadata vs. remove features
#' @param tag_feats_name column name for tagged features in metadata
#' @param verbose verbose
#'
#' @returns giotto object
#' @details The function \code{\link{filterCombinations}} can be used to
#' explore the effect of different parameter values.
#' Please note that this function filters data in a predefined order, features,
#' then cells.
#' After filtering in this order, certain features may be left over in the
#' metadata with a corresponding number of cells which is less than that of
#' the threshold value of cells,
#' feat_det_in_min_cells. This behavior is explained in detail here:
#' \url{https://github.com/drieslab/Giotto/issues/500#issuecomment-1396083446}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' filterGiotto(g)
#' @export
filterGiotto <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("raw", "normalized", "scaled", "custom"),
        expression_threshold = 1,
        feat_det_in_min_cells = 100,
        min_det_feats_per_cell = 100,
        spat_unit_fsub = ":all:",
        feat_type_ssub = ":all:",
        all_spat_units = NULL,
        all_feat_types = NULL,
        poly_info = NULL,
        tag_cells = FALSE,
        tag_cell_name = "tag",
        tag_feats = FALSE,
        tag_feats_name = "tag",
        verbose = TRUE) {
    # data.table vars
    cell_ID <- feat_ID <- NULL

    # handle deprecations
    if (!is.null(all_spat_units)) {
        if (all_spat_units) {
            spat_unit_fsub <- ":all:"
        } else {
            spat_unit_fsub <- spat_unit
        }

        warning(wrap_txt(
            'filterGiotto:
      all_spat_units param is deprecated.
      Please use spat_unit_fsub = \":all:\" instead. (this is the default)'
        ))
    }
    if (!is.null(all_feat_types)) {
        if (all_feat_types) {
            feat_type_ssub <- ":all:"
        } else {
            feat_type_ssub <- feat_type
        }

        warning(wrap_txt(
            'filterGiotto: all_feat_types param is deprecated.
            Please use feat_type_ssub = \":all:\" instead.
            (this is the default)'
        ))
    }


    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    # set poly_info
    if (is.null(poly_info)) {
        poly_info <- spat_unit
    }

    if (verbose && length(spat_unit) > 1L) {
        wrap_msg(
            "More than one spat_unit provided.\n",
            paste0("[", spat_unit[[1L]], "]"),
            "filtering results will be applied across spat_units:", spat_unit
        )
    }
    if (verbose && length(feat_type) > 1L) {
        wrap_msg(
            "More than one feat_type provided.\n",
            paste0("[", feat_type[[1L]], "]"),
            "filtering results will be applied across spat_units:", feat_type
        )
    }


    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("raw", "normalized", "scaled", "custom", expression_values))
    )

    # get expression values to perform filtering on
    # Only the first spat_unit and feat_type provided are filtered.
    # IF there are additional spat_units and feat_types provided, then the
    # filtering
    # results from this round will be applied to the other provided spat_units
    # and feat_types as well.
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit[[1L]],
        feat_type = feat_type[[1L]],
        values = values,
        output = "matrix"
    )

    # approach:
    # 1. first remove genes that are not frequently detected
    # 2. then remove cells that do not have sufficient detected genes

    ## filter features
    filter_index_feats <- rowSums_flex(
        expr_values >= expression_threshold
    ) >= feat_det_in_min_cells
    selected_feat_ids <- names(filter_index_feats[filter_index_feats == TRUE])



    ## filter cells
    filter_index_cells <- colSums_flex(expr_values[
        filter_index_feats,
    ] >= expression_threshold) >= min_det_feats_per_cell
    selected_cell_ids <- names(filter_index_cells[filter_index_cells == TRUE])



    # update cell metadata
    if (isTRUE(tag_cells)) {
        cell_meta <- getCellMetadata(gobject = gobject, copy_obj = TRUE)
        cell_meta[][, c(tag_cell_name) := ifelse(
            cell_ID %in% selected_cell_ids, 0, 1
        )]
        gobject <- setCellMetadata(
            gobject = gobject, x = cell_meta, initialize = FALSE
        )

        # set selected cells back to all cells
        selected_cell_ids <- names(filter_index_cells)
    }

    if (isTRUE(tag_feats)) {
        feat_meta <- getFeatureMetadata(gobject = gobject, copy_obj = TRUE)
        feat_meta[][, c(tag_feats_name) := ifelse(
            feat_ID %in% selected_feat_ids, 0, 1
        )]
        gobject <- setFeatureMetadata(
            gobject = gobject, x = feat_meta, initialize = FALSE
        )

        # set selected feats back to all feats
        selected_feat_ids <- names(filter_index_feats)
    }



    # update feature metadata
    newGiottoObject <- subsetGiotto(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        cell_ids = selected_cell_ids,
        feat_ids = selected_feat_ids,
        spat_unit_fsub = spat_unit_fsub,
        feat_type_ssub = feat_type_ssub,
        poly_info = poly_info,
        verbose = verbose
    )

    ## print output ##
    removed_feats <- length(filter_index_feats[filter_index_feats == FALSE])
    total_feats <- length(filter_index_feats)

    removed_cells <- length(filter_index_cells[filter_index_cells == FALSE])
    total_cells <- length(filter_index_cells)

    if (isTRUE(verbose)) {
        cat("\n")
        cat("Feature type: ", feat_type, "\n")

        if (isTRUE(tag_cells)) {
            cat(
                "Number of cells tagged: ", removed_cells, " out of ",
                total_cells, "\n"
            )
        } else {
            cat(
                "Number of cells removed: ", removed_cells, " out of ",
                total_cells, "\n"
            )
        }

        if (isTRUE(tag_feats)) {
            cat(
                "Number of feats tagged: ", removed_feats, " out of ",
                total_feats, "\n"
            )
        } else {
            cat(
                "Number of feats removed: ", removed_feats, " out of ",
                total_feats, "\n"
            )
        }
    }


    ## update parameters used ##

    # Do not update downstream of processGiotto
    # Parameters will be updated within processGiotto
    try(
        {
            upstream_func <- sys.call(-2)
            fname <- as.character(upstream_func[[1]])
            if (fname == "processGiotto") {
                return(newGiottoObject)
            }
        },
        silent = TRUE
    )


    # If this function call is not downstream of processGiotto, update normally
    newGiottoObject <- update_giotto_params(
        newGiottoObject,
        description = "_filter"
    )

    return(newGiottoObject)
}




### normalization ####


#' @title RNA standard normalization
#' @name .rna_standard_normalization
#' @description standard function for RNA normalization
#' @returns giotto object
#' @keywords internal
.rna_standard_normalization <- function(
        gobject,
        raw_expr,
        feat_type,
        spat_unit,
        library_size_norm = TRUE,
        scalefactor = 6e3,
        log_norm = TRUE,
        log_offset = 1,
        logbase = 2,
        scale_feats = TRUE,
        scale_cells = TRUE,
        scale_order = c("first_feats", "first_cells"),
        verbose = TRUE) {
    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: Standard normalization was developed for RNA data \n")
    }

    # evaluate provenance before modifying raw_expr in case h5_file exists
    if (isS4(raw_expr)) {
        provenance <- raw_expr@provenance
    } else {
        provenance <- NULL
    }


    feat_names <- rownames(raw_expr[])
    col_names <- colnames(raw_expr[])




    ## 1. library size normalize
    if (library_size_norm == TRUE) {
        norm_expr <- .lib_norm_giotto(
            mymatrix = raw_expr[],
            scalefactor = scalefactor
        )
    } else {
        norm_expr <- raw_expr[]
    }

    ## 2. lognormalize
    if (log_norm == TRUE) {
        norm_expr <- .log_norm_giotto(
            mymatrix = norm_expr,
            base = logbase,
            offset = log_offset
        )
    }

    ## 3. scale
    if (scale_feats == TRUE & scale_cells == TRUE) {
        scale_order <- match.arg(
            arg = scale_order, choices = c("first_feats", "first_cells")
        )

        if (scale_order == "first_feats") {
            if (isTRUE(verbose)) {
                wrap_msg("\n first scale feats and then cells \n")
            }

            norm_scaled_expr <- t_flex(standardise_flex(
                x = t_flex(norm_expr), center = TRUE, scale = TRUE
            ))
            norm_scaled_expr <- standardise_flex(
                x = norm_scaled_expr, center = TRUE, scale = TRUE
            )
        } else if (scale_order == "first_cells") {
            if (isTRUE(verbose)) {
                wrap_msg("\n first scale cells and then feats \n")
            }

            norm_scaled_expr <- standardise_flex(
                x = norm_expr, center = TRUE, scale = TRUE
            )
            norm_scaled_expr <- t_flex(standardise_flex(
                x = t_flex(norm_scaled_expr), center = TRUE, scale = TRUE
            ))
        } else {
            stop("\n scale order must be given \n")
        }
    } else if (scale_feats == TRUE) {
        norm_scaled_expr <- t_flex(standardise_flex(
            x = t_flex(norm_expr), center = TRUE, scale = TRUE
        ))
    } else if (scale_cells == TRUE) {
        norm_scaled_expr <- standardise_flex(
            x = norm_expr, center = TRUE, scale = TRUE
        )
    } else {
        norm_scaled_expr <- NULL
    }


    ## 4. add cell and gene names back
    if (!is.null(norm_expr)) {
        rownames(norm_expr) <- feat_names
        colnames(norm_expr) <- col_names
    }
    if (!is.null(norm_scaled_expr)) {
        rownames(norm_scaled_expr) <- feat_names
        colnames(norm_scaled_expr) <- col_names
    }

    ## 5. create and set exprObj
    norm_expr <- create_expr_obj(
        name = "normalized",
        exprMat = norm_expr,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = provenance,
        misc = NULL
    )

    norm_scaled_expr <- create_expr_obj(
        name = "scaled",
        exprMat = norm_scaled_expr,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = provenance,
        misc = NULL
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- set_expression_values(
        gobject = gobject,
        values = norm_expr
    )

    gobject <- set_expression_values(
        gobject = gobject,
        values = norm_scaled_expr
    )
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    ## 6. return Giotto object
    return(gobject)
}



#' @title RNA osmfish normalization
#' @name .rna_osmfish_normalization
#' @description function for RNA normalization according to osmFISH paper
#' @returns giotto object
#' @keywords internal
.rna_osmfish_normalization <- function(
        gobject,
        raw_expr,
        feat_type,
        spat_unit,
        name = "custom",
        verbose = TRUE) {
    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: osmFISH normalization was developed for RNA in situ
                data \n")
    }

    # 1. normalize per gene with scale-factor equal to number of genes
    norm_feats <- (raw_expr[] / rowSums_flex(raw_expr[])) * nrow(raw_expr[])
    # 2. normalize per cells with scale-factor equal to number of cells
    norm_feats_cells <- t_flex((t_flex(norm_feats) /
        colSums_flex(norm_feats)) * ncol(raw_expr[]))

    # return results to Giotto object
    if (verbose == TRUE) {
        message(
            "\n osmFISH-like normalized data will be returned to the",
            name, "Giotto slot \n"
        )
    }

    norm_feats_cells <- create_expr_obj(
        name = name,
        exprMat = norm_feats_cells,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = raw_expr@provenance
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- set_expression_values(
        gobject = gobject,
        values = norm_feats_cells
    )
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    return(gobject)
}


#' @title RNA pearson residuals normalization
#' @name .rna_pears_resid_normalization
#' @description function for RNA normalization according to Lause/Kobak et al
#' paper
#' Adapted from https://gist.github.com/hypercompetent/51a3c428745e1c06d826d76c3671797c#file-pearson_residuals-r
#' @returns giotto object
#' @keywords internal
.rna_pears_resid_normalization <- function(
        gobject,
        raw_expr,
        feat_type,
        spat_unit,
        theta = 100,
        name = "scaled",
        verbose = TRUE) {
    # print message with information #
    if (verbose) {
        message("using 'Lause/Kobak' method to normalize count matrix If used in
      published research, please cite:
      Jan Lause, Philipp Berens, Dmitry Kobak (2020).
      'Analytic Pearson residuals for normalization of single-cell RNA-seq UMI
      data' ")
    }


    # check feature type compatibility
    if (!feat_type %in% c("rna", "RNA")) {
        warning("Caution: pearson residual normalization was developed for RNA
                count normalization \n")
    }

    if (methods::is(raw_expr[], "HDF5Matrix")) {
        counts_sum0 <- methods::as(matrix(
            MatrixGenerics::colSums2(raw_expr[]),
            nrow = 1
        ), "HDF5Matrix")
        counts_sum1 <- methods::as(matrix(
            MatrixGenerics::rowSums2(raw_expr[]),
            ncol = 1
        ), "HDF5Matrix")
        counts_sum <- sum(raw_expr[])

        # get residuals
        mu <- (counts_sum1 %*% counts_sum0) / counts_sum
        z <- (raw_expr[] - mu) / sqrt(mu + mu^2 / theta)

        # clip to sqrt(n)
        n <- ncol(raw_expr[])
        z[z > sqrt(n)] <- sqrt(n)
        z[z < -sqrt(n)] <- -sqrt(n)
    } else {
        counts_sum0 <- methods::as(matrix(Matrix::colSums(
            raw_expr[]
        ), nrow = 1), "dgCMatrix")
        counts_sum1 <- methods::as(matrix(Matrix::rowSums(
            raw_expr[]
        ), ncol = 1), "dgCMatrix")
        counts_sum <- sum(raw_expr[])

        # get residuals
        mu <- (counts_sum1 %*% counts_sum0) / counts_sum
        z <- (raw_expr[] - mu) / sqrt(mu + mu^2 / theta)

        # clip to sqrt(n)
        n <- ncol(raw_expr[])
        z[z > sqrt(n)] <- sqrt(n)
        z[z < -sqrt(n)] <- -sqrt(n)
    }

    # return results to Giotto object
    if (verbose == TRUE) {
        message(
            "\n Pearson residual normalized data will be returned to the ",
            name, " Giotto slot \n"
        )
    }

    z <- create_expr_obj(
        name = name,
        exprMat = z,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = raw_expr@provenance
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- set_expression_values(
        gobject = gobject,
        values = z
    )
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
}




#' @title normalizeGiotto
#' @name normalizeGiotto
#' @description fast normalize and/or scale expresion values of Giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param norm_methods normalization method to use
#' @param library_size_norm normalize cells by library size
#' @param scalefactor scale factor to use after library size normalization
#' @param log_norm transform values to log-scale
#' @param log_offset offset value to add to expression matrix, default = 1
#' @param logbase log base to use to log normalize expression values
#' @param scale_feats z-score genes over all cells
#' @param scale_genes deprecated, use scale_feats
#' @param scale_cells z-score cells over all genes
#' @param scale_order order to scale feats and cells
#' @param theta theta parameter for the pearson residual normalization step
#' @param update_slot slot or name to use for the results from osmFISH and
#' pearson residual normalization
#' @param verbose be verbose
#' @returns giotto object
#' @details Currently there are two 'methods' to normalize your raw counts data.
#'
#' A. The standard method follows the standard protocol which can be adjusted
#' using the provided parameters and follows the following order: \cr
#' \itemize{
#'   \item{1. Data normalization for total library size and scaling by a custom scale-factor.}
#'   \item{2. Log transformation of data.}
#'   \item{3. Z-scoring of data by genes and/or cells.}
#' }
#' B. The normalization method as provided by the osmFISH paper is also implemented: \cr
#' \itemize{
#'   \item{1. First normalize genes, for each gene divide the counts by the total gene count and
#' multiply by the total number of genes.}
#'   \item{2. Next normalize cells, for each cell divide the normalized gene counts by the total
#' counts per cell and multiply by the total number of cells.}
#' }
#' C. The normalization method as provided by Lause/Kobak et al is also implemented: \cr
#' \itemize{
#'   \item{1. First calculate expected values based on Pearson correlations.}
#'   \item{2. Next calculate z-scores based on observed and expected values.}
#' }
#' By default the latter two results will be saved in the Giotto slot for
#' scaled expression, this can be changed by changing the update_slot parameters
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' normalizeGiotto(g)
#' @export
normalizeGiotto <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = "raw",
        norm_methods = c("standard", "pearson_resid", "osmFISH"),
        library_size_norm = TRUE,
        scalefactor = 6e3,
        log_norm = TRUE,
        log_offset = 1,
        logbase = 2,
        scale_feats = TRUE,
        scale_genes = NULL,
        scale_cells = TRUE,
        scale_order = c("first_feats", "first_cells"),
        theta = 100,
        update_slot = "scaled",
        verbose = TRUE) {
    ## deprecated arguments
    if (!is.null(scale_genes)) {
        scale_feats <- scale_genes
        warning("scale_genes is deprecated, use scale_feats in the future \n")
    }

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    ## default is to start from raw data
    values <- match.arg(expression_values, unique(c("raw", expression_values)))
    raw_expr <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "exprObj"
    )

    norm_methods <- match.arg(
        arg = norm_methods, choices = c("standard", "pearson_resid", "osmFISH")
    )

    # normalization according to standard methods
    if (norm_methods == "standard") {
        gobject <- .rna_standard_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            library_size_norm = library_size_norm,
            scalefactor = scalefactor,
            log_norm = log_norm,
            log_offset = log_offset,
            logbase = logbase,
            scale_feats = scale_feats,
            scale_cells = scale_cells,
            scale_order = scale_order,
            verbose = verbose
        )
    } else if (norm_methods == "osmFISH") {
        gobject <- .rna_osmfish_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            name = update_slot,
            verbose = verbose
        )
    } else if (norm_methods == "pearson_resid") {
        gobject <- .rna_pears_resid_normalization(
            gobject = gobject,
            raw_expr = raw_expr,
            feat_type = feat_type,
            spat_unit = spat_unit,
            theta = theta,
            name = update_slot,
            verbose = verbose
        )
    }

    ## update parameters used ##

    # Do not update downstream of processGiotto
    # Parameters will be updated within processGiotto
    try(
        {
            upstream_func <- sys.call(-2)
            fname <- as.character(upstream_func[[1]])
            if (fname == "processGiotto") {
                return(gobject)
            }
        },
        silent = TRUE
    )


    # If this function call is not downstream of processGiotto, update normally
    gobject <- update_giotto_params(gobject, description = "_normalize")

    return(gobject)
}



#' @title Adjust expression values
#' @name adjustGiottoMatrix
#' @description Adjust expression values to account for known batch effects or
#' technological covariates.
#' @inheritParams data_access_params
#' @param expression_values expression values to use
#' @param batch_columns metadata columns that represent different
#' batch (max = 2)
#' @param covariate_columns metadata columns that represent covariates to
#' regress out
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param update_slot expression slot that will be updated (default = custom)
#' @returns giotto object or exprObj
#' @details This function implements the \code{\link[limma]{removeBatchEffect}}
#' function to remove known batch effects and to adjust expression values
#' according to provided covariates.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' adjustGiottoMatrix(g, covariate_columns = "leiden_clus")
#' @export
adjustGiottoMatrix <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        batch_columns = NULL,
        covariate_columns = NULL,
        return_gobject = TRUE,
        update_slot = c("custom")) {
    # Catch for both batch and covariate being null
    if (is.null(batch_columns) & is.null(covariate_columns)) {
        stop("Metadata for either different batches or covariates must be
            provided.")
    }

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # metadata
    cell_metadata <- getCellMetadata(
        gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        output = "data.table",
        copy_obj = TRUE
    )

    if (!is.null(batch_columns)) {
        if (!all(batch_columns %in% colnames(cell_metadata))) {
            stop("batch column name(s) were not found in the cell metadata")
        }
    }

    if (!is.null(covariate_columns)) {
        if (!all(covariate_columns %in% colnames(cell_metadata))) {
            stop("covariate column name(s) were not found in the cell metadata")
        }
    }

    update_slot <- match.arg(
        update_slot, c("normalized", "scaled", "custom", update_slot)
    )

    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "exprObj"
    )


    # batch columns
    if (!is.null(batch_columns)) {
        batch_column_1 <- cell_metadata[[batch_columns[1]]]
        if (length(batch_columns) > 1) {
            batch_column_2 <- cell_metadata[[batch_columns[2]]]
        } else {
            batch_column_2 <- NULL
        }
    } else {
        batch_column_1 <- NULL
        batch_column_2 <- NULL
    }

    # covariate columns
    if (!is.null(covariate_columns)) {
        covariates <- as.matrix(
            cell_metadata[, covariate_columns, with = FALSE]
        )
    } else {
        covariates <- NULL
    }



    # TODO: implement ResidualMatrix to work with a delayed matrix
    adjusted_matrix <- limma::removeBatchEffect(
        x = expr_data[],
        batch = batch_column_1,
        batch2 = batch_column_2,
        covariates = covariates
    )

    if (return_gobject == TRUE) {
        adjusted_matrix <- create_expr_obj(
            name = update_slot,
            exprMat = adjusted_matrix,
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = expr_data@provenance
        )

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_expression_values(
            gobject = gobject,
            values = adjusted_matrix
        )
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        ## update parameters used ##

        # Do not update downstream of processGiotto
        # Parameters will be updated within processGiotto
        try(
            {
                test <- sys.call(-2)
                fname <- as.character(test[[1]])
                if (fname == "processGiotto") {
                    return(gobject)
                }
            },
            silent = TRUE
        )


        # If this function call is not downstream of processGiotto, update
        # normally
        gobject <- update_giotto_params(gobject, description = "_adj_matrix")

        return(gobject)
    } else {
        return(adjusted_matrix)
    }
}


#' @title processGiotto
#' @name processGiotto
#' @description Wrapper for the different Giotto object processing functions
#' @param gobject giotto object
#' @param filter_params additional parameters to filterGiotto
#' @param norm_params additional parameters to normalizeGiotto
#' @param stat_params additional parameters to addStatistics
#' @param adjust_params additional parameters to adjustGiottoMatrix; set to
#' NULL if not required
#' @param verbose be verbose (default is TRUE)
#' @returns giotto object
#' @details See \code{\link{filterGiotto}}, \code{\link{normalizeGiotto}},
#' \code{\link{addStatistics}}, and \code{\link{adjustGiottoMatrix}}. For more
#' information about the different parameters in each step. If you do not
#' provide them it will use the default values. If no adjustment is required,
#' adjust_params must be set to NULL
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' processGiotto(
#'     gobject = g,
#'     adjust_params = list(covariate_columns = "leiden_clus")
#' )
#' @export
processGiotto <- function(
        gobject,
        filter_params = list(),
        norm_params = list(),
        stat_params = list(),
        adjust_params = list(),
        verbose = TRUE) {
    # filter Giotto
    if (verbose == TRUE) message("1. start filter step")
    if (!inherits(filter_params, "list")) {
        stop("filter_params need to be a list of parameters for filterGiotto")
    }
    gobject <- do.call("filterGiotto", c(gobject = gobject, filter_params))

    # normalize Giotto
    if (verbose == TRUE) message("2. start normalization step")
    if (!inherits(norm_params, "list")) {
        stop("norm_params need to be a list of parameters for normalizeGiotto")
    }
    gobject <- do.call("normalizeGiotto", c(gobject = gobject, norm_params))

    # add Statistics
    if (verbose == TRUE) message("3. start cell and gene statistics step")
    if (!inherits(stat_params, "list")) {
        stop("stat_params need to be a list of parameters for addStatistics ")
    }
    stat_params[["return_gobject"]] <- TRUE # force this to be true
    gobject <- do.call("addStatistics", c(gobject = gobject, stat_params))

    # adjust Giotto, if applicable
    if (!is.null(adjust_params)) {
        if (verbose == TRUE) message("4. start adjusted matrix step")
        if (!inherits(adjust_params, "list")) {
            stop("adjust_params need to be a list of parameters for
                adjustGiottoMatrix")
        }
        adjust_params[["return_gobject"]] <- TRUE # force this to be true
        gobject <- do.call(
            "adjustGiottoMatrix", c(gobject = gobject, adjust_params)
        )
    }

    gobject <- update_giotto_params(gobject, description = "_process")

    return(gobject)
}










## * ####
## Feature & Cell statistics ####





#' @title Add feature statistics
#' @name addFeatStatistics
#' @description Adds feature statistics to the giotto object
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose be verbose
#' @returns giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to feature metadata:
#' \itemize{
#'   \item{\strong{nr_cells:} Denotes in how many cells the feature is
#'   detected}
#'   \item{\strong{per_cells:} Denotes in what percentage of cells the feature
#'   is detected}
#'   \item{\strong{total_expr:} Shows the total sum of feature expression in
#'   all cells}
#'   \item{\strong{mean_expr:} Average feature expression in all cells}
#'   \item{\strong{mean_expr_det:} Average feature expression in cells with
#'   detectable levels of the gene}
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addFeatStatistics(g)
#' @export
addFeatStatistics <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        detection_threshold = 0,
        return_gobject = TRUE,
        verbose = TRUE) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # expression values to be used
    expression_values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values,
        output = "exprObj",
        set_defaults = FALSE
    )

    # calculate stats
    feat_stats <- data.table::data.table(
        feats = rownames(expr_data[]),
        nr_cells = rowSums_flex(expr_data[] > detection_threshold),
        perc_cells = (rowSums_flex(expr_data[] > detection_threshold) /
            ncol(expr_data[])) * 100,
        total_expr = rowSums_flex(expr_data[]),
        mean_expr = rowMeans_flex(expr_data[])
    )

    # data.table variables
    mean_expr_det <- NULL

    mean_expr_detected <- .mean_expr_det_test(
        expr_data[],
        detection_threshold = detection_threshold
    )
    feat_stats[, mean_expr_det := mean_expr_detected]


    if (return_gobject == TRUE) {
        # remove previous statistics
        feat_metadata <- getFeatureMetadata(
            gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "featMetaObj",
            copy_obj = TRUE,
            set_defaults = FALSE
        )

        if (isS4(expr_data)) {
            if (!identical(expr_data@provenance, feat_metadata@provenance)) {
                warning("expression and feature metadata provenance mismatch")
            }
        }


        metadata_names <- colnames(feat_metadata[])

        if ("nr_cells" %in% metadata_names) {
            vmsg(
                .v = verbose, "feat statistics has already been applied",
                "once; overwriting"
            )
            feat_metadata[][, c(
                "nr_cells", "perc_cells", "total_expr", "mean_expr",
                "mean_expr_det"
            ) := NULL]
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- set_feature_metadata(gobject,
                metadata = feat_metadata,
                verbose = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }

        gobject <- addFeatMetadata(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            new_metadata = feat_stats,
            by_column = TRUE,
            column_feat_ID = "feats"
        )

        ## update parameters used ##

        # parent function name
        cl <- sys.call(-1)

        # Do not update downstream of processGiotto
        # Parameters will be updated within processGiotto
        try(
            {
                upstream_func <- sys.call(-3)
                fname <- as.character(upstream_func[[1]])
                if (fname == "processGiotto") {
                    return(gobject)
                }
            },
            silent = TRUE
        )


        # If this function call is not downstream of processGiotto, update
        # normally
        if (is.null(cl)) {
            gobject <- update_giotto_params(gobject,
                description = "_feat_stats"
            )
        } else {
            fname <- as.character(cl[[1]])
            if (fname == "addStatistics") {
                gobject <- update_giotto_params(gobject,
                    description = "_feat_stats",
                    toplevel = 3
                )
            } else {
                gobject <- update_giotto_params(gobject,
                    description = "_feat_stats"
                )
            }
        }


        return(gobject)
    } else {
        return(feat_stats)
    }
}





#' @title addCellStatistics
#' @name addCellStatistics
#' @description adds cells statistics to the giotto object
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose be verbose
#' @returns giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to cell metadata:
#' \itemize{
#'   \item{\strong{nr_feats:} Denotes in how many features are detected per
#'   cell}
#'   \item{\strong{perc_feats:} Denotes what percentage of features is
#'   detected per cell}
#'   \item{\strong{total_expr:} Shows the total sum of feature expression per
#'   cell}
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addCellStatistics(g)
#' @export
addCellStatistics <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        detection_threshold = 0,
        return_gobject = TRUE,
        verbose = TRUE) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # expression values to be used
    expression_values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values,
        output = "exprObj",
        set_defaults = FALSE
    )

    # calculate stats

    cell_stats <- data.table::data.table(
        cells = colnames(expr_data[]),
        nr_feats = colSums_flex(expr_data[] > detection_threshold),
        perc_feats = (colSums_flex(expr_data[] > detection_threshold) /
            nrow(expr_data[])) * 100,
        total_expr = colSums_flex(expr_data[])
    )

    if (return_gobject == TRUE) {
        # remove previous statistics
        cell_metadata <- getCellMetadata(
            gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "cellMetaObj",
            copy_obj = TRUE,
            set_defaults = FALSE
        )

        if (isS4(expr_data)) {
            if (!identical(expr_data@provenance, cell_metadata@provenance)) {
                warning("expression and feature metadata provenance mismatch")
            }
        }

        metadata_names <- colnames(cell_metadata[])
        if ("nr_feats" %in% metadata_names) {
            vmsg(
                .v = verbose, "cells statistics has already been applied",
                "once; overwriting"
            )
            cell_metadata[][, c("nr_feats", "perc_feats", "total_expr") := NULL]
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- set_cell_metadata(gobject,
                metadata = cell_metadata,
                verbose = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }




        gobject <- addCellMetadata(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            new_metadata = cell_stats,
            by_column = TRUE,
            column_cell_ID = "cells"
        )

        ## update parameters used ##

        # parent function name
        cl <- sys.call(-1)

        # Do not update downstream of processGiotto
        # Parameters will be updated within processGiotto
        try(
            {
                upstream_func <- sys.call(-3)
                fname <- as.character(upstream_func[[1]])
                if (fname == "processGiotto") {
                    return(gobject)
                }
            },
            silent = TRUE
        )

        # If this function call is not downstream of processGiotto, update
        # normally
        if (is.null(cl)) {
            gobject <- update_giotto_params(gobject,
                description = "_cell_stats"
            )
        } else {
            fname <- as.character(cl[[1]])
            if (fname == "addStatistics") {
                gobject <- update_giotto_params(gobject,
                    description = "_cell_stats",
                    toplevel = 3
                )
            } else {
                gobject <- update_giotto_params(gobject,
                    description = "_cell_stats"
                )
            }
        }


        return(gobject)
    } else {
        return(cell_stats)
    }
}


#' @title addStatistics
#' @name addStatistics
#' @description Adds feature and cell statistics to the giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a feature detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose be verbose
#' @returns giotto object if return_gobject = TRUE, else a list with results
#' @details See \code{\link{addFeatStatistics}} and
#' \code{\link{addCellStatistics}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addStatistics(g)
#' @export
addStatistics <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        detection_threshold = 0,
        return_gobject = TRUE,
        verbose = TRUE) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # get feats statistics
    feat_stats <- addFeatStatistics(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        expression_values = expression_values,
        detection_threshold = detection_threshold,
        return_gobject = return_gobject,
        verbose = verbose
    )

    if (return_gobject == TRUE) {
        gobject <- feat_stats
    }

    # get cell statistics
    cell_stats <- addCellStatistics(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        expression_values = expression_values,
        detection_threshold = detection_threshold,
        return_gobject = return_gobject,
        verbose = verbose
    )

    if (return_gobject == TRUE) {
        gobject <- cell_stats
        return(gobject)
    } else {
        return(feat_stats = feat_stats, cell_stats = cell_stats)
    }
}




#' @title addFeatsPerc
#' @name addFeatsPerc
#' @description Calculates the total percentage of (normalized) counts for a
#' subset of selected genes
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param feats vector of selected features
#' @param vector_name column name as seen in \code{\link{pDataDT}}
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @returns giotto object if \code{return_gobject = TRUE}, else a vector with %
#' results
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addFeatsPerc(g, feats = c("Gm19935", "9630013A20Rik", "2900040C04Rik"))
#' @export
addFeatsPerc <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        feats = NULL,
        vector_name = "feat_perc",
        return_gobject = TRUE) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # tests
    if (is.null(feats)) {
        stop("You need to provide a vector of feat names")
    }

    if (!methods::is(gobject, "giotto")) {
        stop("You need to provide a giotto object")
    }


    # expression values to be used
    expression_values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values,
        output = "matrix"
    )


    totalsum <- colSums_flex(expr_data)
    feat_sum <- colSums_flex(expr_data[rownames(expr_data) %in% feats, ])
    perc_feats <- round((feat_sum / totalsum) * 100, 2)

    if (return_gobject) {
        temp_gobj <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = perc_feats,
            vector_name = vector_name,
            by_column = TRUE
        )

        ## update parameters used ##
        temp_gobj <- update_giotto_params(temp_gobj,
            description = "_feats_perc"
        )

        return(temp_gobj)
    } else {
        return(perc_feats)
    }
}





## * ####
## Giotto auxiliary functions ####




#' @title Find network neighbors
#' @name findNetworkNeighbors
#' @description Find the spatial neighbors for a selected group of cells within
#' the selected spatial network.
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param spatial_network_name name of spatial network
#' @param source_cell_ids cell ids for which you want to know the spatial
#' neighbors
#' @param name name of the results
#' @returns data.table
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findNetworkNeighbors(
#'     gobject = g, spatial_network_name = "spatial_network",
#'     source_cell_ids = c("AACTCGATGGCGCAGT-1", "GGCTGGCTAGCTTAAA-1")
#' )
#' @export
findNetworkNeighbors <- function(
        gobject,
        spat_unit = NULL,
        spatial_network_name = NULL,
        source_cell_ids = NULL,
        name = "nb_cells") {
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )

    # get spatial network
    if (!is.null(spatial_network_name)) {
        spatial_network <- getSpatialNetwork(gobject,
            spat_unit = spat_unit,
            name = spatial_network_name,
            output = "networkDT"
        )
    } else {
        stop("You need to select a spatial network")
    }

    # source cell ids that are found back
    all_cell_ids <- gobject@cell_ID[[spat_unit]]
    source_cells <- all_cell_ids[all_cell_ids %in% source_cell_ids]

    if (length(source_cells) == 0) {
        stop("No source cell ids were selected or found")
    }


    full_network_DT <- convert_to_full_spatial_network(spatial_network)
    potential_target_cells <- full_network_DT[
        source %in% source_cells
    ][["target"]]
    source_and_target_cells <- potential_target_cells[
        potential_target_cells %in% source_cells
    ]
    target_cells <- potential_target_cells[
        !potential_target_cells %in% source_and_target_cells
    ]

    cell_meta <- pDataDT(gobject)

    # data.table variables
    nb_cells <- cell_ID <- NULL

    cell_meta[, nb_cells := ifelse(cell_ID %in% source_and_target_cells, "both",
        ifelse(cell_ID %in% source_cells, "source",
            ifelse(cell_ID %in% target_cells, "neighbor", "others")
        )
    )]
    nb_annot <- cell_meta[, c("cell_ID", "nb_cells"), with = FALSE]
    data.table::setnames(nb_annot, "nb_cells", name)

    return(nb_annot)
}
