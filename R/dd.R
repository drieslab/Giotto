

# ------------------------------------------------------------------------- #
# This script contains reusable dummy documentation for commonly used params
# Use @inheritParams data_access when documenting these parameters
# ------------------------------------------------------------------------- #

# No @title to prevent .Rd generation
# No @noRd so @inheritParams can be used

# Note that these dummy documentations will be flagged with warnings when building
# the documentation



#' @name data_access
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param return_uniques return unique nesting names (ignores if final object exists/is correct class)
#' @param output what format in which to get information (e.g. "data.table")
#' @param set_defaults set default spat_unit and feat_type. Change to FALSE only when
#' expression and spat_info are not expected to exist.
#' @param copy_obj whether to deep copy/duplicate when getting the object (default = TRUE)
#' @keywords internal
NULL





















