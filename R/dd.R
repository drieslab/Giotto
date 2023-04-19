

# ------------------------------------------------------------------------- #
# This script contains reusable dummy documentation / templates for
# commonly used params.
#
# Use the @inheritParams tag to use these templates in documentation
# ------------------------------------------------------------------------- #

# No @title to prevent .Rd generation
# No @noRd tags can be used on these dummy documentations, otherwise @inheritParams
# cannot be used

# Note that these dummy documentations WILL be flagged with warnings when building
# the documentation, but this should be fine.








# Data Access ####


#' data_access_params
#' 
#' @name data_access_params
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param return_uniques return unique nesting names (ignores if final object exists/is correct class)
#' @param output what format in which to get information (e.g. "data.table")
#' @param set_defaults set default spat_unit and feat_type. Change to FALSE only when
#' expression and spat_info are not expected to exist.
#' @param copy_obj whether to deep copy/duplicate when getting the object (default = TRUE)
#' @param initialize (default = FALSE) whether to initialize the gobject before
#' returning
#' @keywords internal
NULL




# Read Functions ####

#' @name read_data_params
#' @param data_list (nested) list of input data to read
#' @param default_spat_unit (optional) default spat_unit to use
#' @param default_feat_type (optional) default feat_type to use
#' @param provenance (optional) provenance information
#' @param verbose be verbose
#' @keywords internal
NULL




# Plotting ####


#' plot_params
#' 
#' @name plot_params
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#'
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#'
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#'
#' @param show_other_cells display not selected cells
#' @param other_cell_color color for not selected cells
#' @param other_point_size point size for not selected cells
#'
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#'
#' @param center_point_size size of center points
#' @param center_point_border_col border color of center points
#' @param center_point_border_stroke border stroke size of center points
#'
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#'
#' @param spat_enr_names names of spatial enrichment results to include
#'
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#'
#' @param label_size  size of labels
#' @param label_fontface font of labels
#'
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_alpha transparency of points
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#'
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#'
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative heights of rows (e.g. c(1,2))
#' @param cow_rel_w cowplot param: relative widths of columns (e.g. c(1,2))
#' @param cow_align cowplot param: how to align
#'
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @keywords internal
NULL




















