
# ** Spatial Method-Specific Convenience Functions for Giotto Object Creation ** #



# Common Utility Functions ####

#' @title Read a structured folder of exported data
#' @name read_data_folder
#' @description Framework function for reading the exported folder of a spatial
#' method and detecting the presence of needed files. NULL values denote missing
#' items.\cr
#' `.read_data_folder()` should not be called directly. Instead, specific
#' reader functions should be built using it as a base.
#' @param spat_method spatial method for which the data is being read
#' @param data_dir exported data directory to read from
#' @param dir_items named list of directory items to expect and keywords to match
#' @param data_to_use character. Which type(s) of expression data to build the
#' gobject with. Values should match with a *workflow* item in require_data_DT
#' (see details)
#' @param require_data_DT data.table detailing if expected data items are required
#' or optional for each \code{data_to_use} *workflow*
#' @param cores cores to use
#' @param verbose be verbose
#' @param toplevel stackframes back where the user-facing function was called.
#' default is one stackframe above `.read_data_folder`.
#' @details
#' **Steps performed:**
#' \itemize{
#'   \item{1. detection of items within \code{data_dir} by looking for keywords
#'   assigned through \code{dir_items}}
#'   \item{2. check of detected items to see if everything needed has been found.
#'   Dictionary of necessary vs optional items for each \code{data_to_use} *workflow*
#'   is provided through \code{require_data_DT}}
#'   \item{3. if multiple filepaths are found to be matching then select the first
#'   one. This function is only intended to find the first level subdirectories
#'   and files.}
#' }
#'
#' **Example reader implementation:**
#' \preformatted{
#'   foo <- function(x_dir,
#'                   data_to_use,
#'                   cores = NA,
#'                   verbose = NULL) {
#'     dir_items <- list(
#'       data1 = "regex_pattern1",
#'       data2 = "regex_pattern2",
#'       data3 = "regex_pattern3"
#'     )
#'
#'     # DT of info to check directory for. Has 3 cols
#'     require_data_DT <- data.table::data.table(
#'       workflow = "a", # data_to_use is matched against this
#'       item = c(
#'         "data1",
#'         "data2",
#'         "data3"
#'       ),
#'       needed = c(
#'         FALSE, # data1 optional for this workflow (if missing: warn)
#'         TRUE,  # data2 vital for this workflow (if missing: error)
#'         TRUE   # data3 vital for this workflow (if missing: error)
#'       )
#'     )
#'
#'     .read_data_folder(
#'       spat_method = "x_method",
#'       data_dir = x_dir,
#'       dir_items = dir_items,
#'       data_to_use = data_to_use,
#'       require_data_DT = require_data_DT,
#'       cores = cores,
#'       verbose = verbose
#'     )
#'   }
#' }
#'
#' @md
NULL

#' @describeIn read_data_folder Should not be used directly
#' @keywords internal
.read_data_folder <- function(spat_method = NULL,
                              data_dir = NULL,
                              dir_items,
                              data_to_use,
                              load_format = NULL,
                              require_data_DT,
                              cores = NA,
                              verbose = NULL,
                              toplevel = 2L) {

  ch = box_chars()

  # 0. check params
  if(is.null(data_dir) ||
     !dir.exists(data_dir)) {
    .gstop(.n = toplevel, 'The full path to a', spat_method, 'directory must be given.')
  }
  vmsg(.v = verbose, 'A structured', spat_method, 'directory will be used')
  if(!data_to_use %in% require_data_DT$workflow) {
    .gstop(.n = toplevel, 'Data requirements for data_to_use not found in require_data_DT')
  }

  # 1. detect items
  dir_items = lapply_flex(dir_items, function(x) {
    Sys.glob(paths = file.path(data_dir, x))
  }, cores = cores)
  # (length = 1 if present, length = 0 if missing)
  dir_items_lengths = lengths(dir_items)

  # 2. check directory contents
  vmsg(.v = verbose, 'Checking directory contents...')

  for(item in names(dir_items)) {

    # IF ITEM FOUND

    if(dir_items_lengths[[item]] > 0) {
      # print found items if verbose = "debug"
      if(isTRUE(verbose)) {
        vmsg(
          .v = verbose, .is_debug = TRUE,
          .initial = paste0(ch$s, '> '),
          item, ' found'
        )
        for(item_i in seq_along(dir_items[[item]])) { # print found item names
          subItem = gsub(pattern = '.*/', replacement = '', x = dir_items[[item]][[item_i]])
          vmsg(
            .v = verbose, .is_debug = TRUE,
            .initial = paste0(ch$s, ch$s, ch$l, ch$h, ch$h),
            subItem
          )
        }
      }
    } else {

      # IF ITEM MISSING
      # necessary (error)
      # optional (warning)

      # data.table variables
      workflow = needed = filetype = NULL


      require_data_DT = require_data_DT[workflow == data_to_use,]
      if(!is.null(load_format)) require_data_DT = require_data_DT[filetype == load_format,]

      if(item %in% require_data_DT[needed == TRUE, item]) stop(item, ' is missing\n')
      if(item %in% require_data_DT[needed == FALSE, item]) warning(item, 'is missing (optional)\n')

    }
  }

  # 3. select first path in list if multiple are detected
  if(any(dir_items_lengths > 1)) {
    warning(wrap_txt('Multiple matches for expected directory item(s).
                     First matching item selected'))

    multiples = which(dir_items_lengths > 1)
    for(mult_i in multiples) {
      message(names(dir_items)[[mult_i]], 'multiple matches found:')
      print(dir_items[[mult_i]])
      dir_items[[mult_i]] = dir_items[[mult_i]][[1]]
    }
  }
  vmsg(.v = verbose, 'Directory check done')

  return(dir_items)

}









# *---- object creation ----* ####






## Visium ####

#' @title Create a giotto object from 10x visium data
#' @name createGiottoVisiumObject
#' @description Create Giotto object directly from a 10X visium folder. Also accepts visium H5 outputs.
#'
#' @param visium_dir path to the 10X visium directory [required]
#' @param expr_data raw or filtered data (see details)
#' @param gene_column_index which column index to select (see details)
#' @param h5_visium_path path to visium 10X .h5 file
#' @param h5_gene_ids gene names as symbols (default) or ensemble gene ids
#' @param h5_tissue_positions_path path to tissue locations (.csv file)
#' @param h5_image_png_path path to tissue .png file (optional). Image autoscaling
#' looks for matches in the filename for either 'hires' or 'lowres'
#' @param h5_json_scalefactors_path path to .json scalefactors (optional)
#' @param png_name select name of png to use (see details)
#' @param do_manual_adj deprecated
#' @param xmax_adj deprecated
#' @param xmin_adj deprecated
#' @param ymax_adj deprecated
#' @param ymin_adj deprecated
#' @param instructions list of instructions or output result from \code{\link[GiottoClass]{createGiottoInstructions}}
#' @param cores how many cores or threads to use to read data if paths are provided
#' @param expression_matrix_class class of expression matrix to use (e.g. 'dgCMatrix', 'DelayedArray')
#' @param h5_file optional path to create an on-disk h5 file
#' @param verbose be verbose
#'
#' @return giotto object
#' @details
#' If starting from a Visium 10X directory:
#' \itemize{
#'   \item{expr_data: raw will take expression data from raw_feature_bc_matrix and filter from filtered_feature_bc_matrix}
#'   \item{gene_column_index: which gene identifiers (names) to use if there are multiple columns (e.g. ensemble and gene symbol)}
#'   \item{png_name: by default the first png will be selected, provide the png name to override this (e.g. myimage.png)}
#'   \item{the file scalefactors_json.json will be detected automatically and used to attempt to align the data}
#' }
#'
#' If starting from a Visium 10X .h5 file
#' \itemize{
#'   \item{h5_visium_path: full path to .h5 file: /your/path/to/visium_file.h5}
#'   \item{h5_tissue_positions_path: full path to spatial locations file: /you/path/to/tissue_positions_list.csv}
#'   \item{h5_image_png_path: full path to png: /your/path/to/images/tissue_lowres_image.png}
#'   \item{h5_json_scalefactors_path: full path to .json file: /your/path/to/scalefactors_json.json}
#' }
#'
#' @export
createGiottoVisiumObject = function(visium_dir = NULL,
                                    expr_data = c('raw', 'filter'),
                                    gene_column_index = 1,
                                    h5_visium_path = NULL,
                                    h5_gene_ids = c('symbols', 'ensembl'),
                                    h5_tissue_positions_path = NULL,
                                    h5_image_png_path = NULL,
                                    h5_json_scalefactors_path = NULL,
                                    png_name = NULL,
                                    do_manual_adj = FALSE, # deprecated
                                    xmax_adj = 0, # deprecated
                                    xmin_adj = 0, # deprecated
                                    ymax_adj = 0, # deprecated
                                    ymin_adj = 0, # deprecated
                                    instructions = NULL,
                                    expression_matrix_class = c("dgCMatrix", "DelayedArray"),
                                    h5_file = NULL,
                                    cores = NA,
                                    verbose = NULL) {

  # NSE vars
  barcode = row_pxl = col_pxl = in_tissue = array_row = array_col = NULL

  # handle deprecations
  img_dep_msg <- "The params 'do_manual_adj', 'xmax_adj', 'xmin_adj', 'ymax_adj,
  'ymin_adj' are no longer used. Please use the automated workflow."
  if(!isFALSE(do_manual_adj) ||
     xmax_adj != 0 ||
     xmin_adj != 0 ||
     ymax_adj != 0 ||
     ymin_adj != 0) {
    stop(wrap_txt(img_dep_msg))
  }

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)


  # get arguments list for object creation
  if(!is.null(h5_visium_path)) {

    argslist <- .visium_read_h5(
      h5_visium_path = h5_visium_path, # expression matrix file
      h5_gene_ids = h5_gene_ids, # symbol or ensembl
      h5_tissue_positions_path = h5_tissue_positions_path,
      h5_image_png_path = h5_image_png_path,
      h5_json_scalefactors_path = h5_json_scalefactors_path,
      verbose = verbose
    )

  } else {

    argslist <- .visium_read_folder(
      visium_dir = visium_dir,
      expr_data = expr_data, # type of expression matrix to load
      gene_column_index = gene_column_index, # symbol or ensembl
      png_name = png_name,
      verbose = verbose
    )

  }

  # additional args to pass to object creation
  argslist$verbose <- verbose
  argslist$expression_matrix_class <- expression_matrix_class
  argslist$h5_file <- h5_file
  argslist$instructions <- instructions

  giotto_object <- do.call(.visium_create, args = argslist)

  return(giotto_object)
}








.visium_create <- function(
    expr_counts_path,
    h5_gene_ids = NULL, # h5
    gene_column_index = NULL, # folder
    tissue_positions_path,
    image_path = NULL,
    scale_json_path = NULL,
    png_name = NULL,
    instructions = NULL,
    expression_matrix_class = c("dgCMatrix", "DelayedArray"),
    h5_file = NULL,
    verbose = NULL
) {

  # NSE vars
  barcode <- cell_ID <- row_pxl <- col_pxl <- in_tissue <- array_row <-
    array_col <- NULL

  # Assume path checking has been done

  # 1. expression
  if (!is.null(h5_gene_ids)) {
    expr_results <- get10Xmatrix_h5(path_to_data = expr_counts_path,
                                    gene_ids = h5_gene_ids)
  } else {
    expr_results <- get10Xmatrix(path_to_data = expr_counts_path,
                                 gene_column_index = gene_column_index)
  }

  # if expr_results is not a list, make it a list compatible with downstream
  if (!is.list(expr_results)) expr_results = list("Gene Expression" = expr_results)

  # format expected data into list to be used with readExprData()
  raw_matrix_list <- list("cell" = list("rna" = list("raw" = expr_results[["Gene Expression"]])))

  # add protein expression data to list if it exists
  if ('Antibody Capture' %in% names(expr_results)) {
    raw_matrix_list$cell$protein$raw <- expr_results[["Antibody Capture"]]
  }


  # 2. spatial locations
  spatial_results <- data.table::fread(tissue_positions_path)
  colnames(spatial_results) <- c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
  spatial_results <- spatial_results[match(colnames(raw_matrix_list$cell[[1]]$raw), barcode)]
  data.table::setnames(spatial_results, old = "barcode", new = "cell_ID")
  spatial_locs <- spatial_results[,.(cell_ID, row_pxl,-col_pxl)] # flip x and y
  colnames(spatial_locs) <- c("cell_ID", 'sdimx', 'sdimy')


  # 3. scalefactors (optional)
  json_info <- .visium_read_scalefactors(scale_json_path)


  # 4. image (optional)
  if (!is.null(image_path)) {

    visium_png_list <- .visium_image(
      image_path = image_path,
      json_info = json_info,
      verbose = verbose
    )
  }

  # 5. metadata
  meta_results <- spatial_results[,.(cell_ID, in_tissue, array_row, array_col)]
  expr_types <- names(raw_matrix_list$cell)
  meta_list <- list()
  for (etype in expr_types) {
    meta_list[[etype]] <- meta_results
  }


  # 6. giotto object
  giotto_object <- createGiottoObject(
    expression = raw_matrix_list,
    spatial_locs = spatial_locs,
    instructions = instructions,
    cell_metadata = meta_list,
    largeImages = visium_png_list
  )


  # 7. polygon information
  if(!is.null(json_info)){
    visium_polygons = .visium_spot_poly(
      spatlocs = spatial_locs,
      json_scalefactors = json_info
    )
    giotto_object = setPolygonInfo(
      gobject = giotto_object,
      x = visium_polygons,
      centroids_to_spatlocs = FALSE,
      verbose = FALSE,
      initialize = TRUE
    )
  }

  return(giotto_object)
}



# Find and check the filepaths within a structured visium directory
.visium_read_folder <- function(
    visium_dir = NULL,
    expr_data = c("raw", "filter"),
    gene_column_index = 1,
    png_name = NULL,
    verbose = NULL
)
{

  vmsg(.v = verbose, "A structured visium directory will be used")

  ## check arguments
  if(is.null(visium_dir)) .gstop('visium_dir needs to be a path to a visium directory \n')
  visium_dir = path.expand(visium_dir)
  if(!dir.exists(visium_dir)) .gstop(visium_dir, ' does not exist!')
  expr_data = match.arg(expr_data, choices = c('raw', 'filter'))


  ## 1. check expression
  expr_counts_path <- switch(
    expr_data,
    "raw" = paste0(visium_dir, '/', 'raw_feature_bc_matrix/'),
    "filter" = paste0(visium_dir, '/', 'filtered_feature_bc_matrix/')
  )
  if (!file.exists(expr_counts_path)) .gstop(expr_counts_path, "does not exist!")


  ## 2. check spatial locations
  spatial_dir = paste0(visium_dir, '/', 'spatial/')
  tissue_positions_path = Sys.glob(paths = file.path(spatial_dir, 'tissue_positions*'))


  ## 3. check spatial image
  if(is.null(png_name)) {
    png_list = list.files(spatial_dir, pattern = "*.png")
    png_name = png_list[1]
  }
  png_path = paste0(spatial_dir,'/',png_name)
  if(!file.exists(png_path)) .gstop(png_path, ' does not exist!')


  ## 4. check scalefactors
  scalefactors_path = paste0(spatial_dir,'/','scalefactors_json.json')
  if (!file.exists(scalefactors_path)) .gstop(scalefactors_path, "does not exist!")


  list(
    expr_counts_path = expr_counts_path,
    gene_column_index = gene_column_index,
    tissue_positions_path = tissue_positions_path,
    image_path = png_path,
    scale_json_path = scalefactors_path
  )
}



.visium_read_h5 <- function(
    h5_visium_path = h5_visium_path, # expression matrix
    h5_gene_ids = h5_gene_ids,
    h5_tissue_positions_path = h5_tissue_positions_path,
    h5_image_png_path = h5_image_png_path,
    h5_json_scalefactors_path = h5_json_scalefactors_path,
    verbose = NULL
) {

  # 1. filepaths
  vmsg(.v = verbose, "A path to an .h5 10X file was provided and will be used \n")
  if (!file.exists(h5_visium_path)) .gstop("The provided path ", h5_visium_path, " does not exist \n")
  if (is.null(h5_tissue_positions_path)) .gstop("A path to the tissue positions (.csv) needs to be provided to h5_tissue_positions_path \n")
  if (!file.exists(h5_tissue_positions_path)) .gstop("The provided path ", h5_tissue_positions_path, " does not exist \n")
  if (!is.null(h5_image_png_path)) {
    if (!file.exists(h5_image_png_path))
      .gstop("The provided h5 image path ", h5_image_png_path, "does not exist.
             Set to NULL to exclude or provide the correct path.\n")
  }
  if (!is.null(h5_json_scalefactors_path)) {
    if (!file.exists(h5_json_scalefactors_path)) {
      warning(wrap_txt(
        "No file found at h5_json_scalefactors_path.
      Scalefactors are needed for proper image alignment and polygon generation"
      ))
    }
  }

  list(
    expr_counts_path = h5_visium_path,
    h5_gene_ids = h5_gene_ids,
    tissue_positions_path = h5_tissue_positions_path,
    image_path = h5_image_png_path,
    scale_json_path = h5_json_scalefactors_path
  )
}









# Visium Polygon Creation

#' @title Add Visium Polygons to Giotto Object
#' @name addVisiumPolygons
#' @param gobject Giotto Object created with visium data, containing spatial
#' locations corresponding to spots
#' @param scalefactor_path path to scalefactors_json.json Visium output
#' @return Giotto Object with to-scale circular polygons added at each spatial
#' location
#' @details
#' Adds circular giottoPolygons to the spatial_info slot of a Giotto Object
#' for the "cell" spatial unit.
#' @export
addVisiumPolygons <- function(gobject,
                              scalefactor_path = NULL){
  assert_giotto(gobject)

  visium_spat_locs = getSpatialLocations(
    gobject = gobject,
    spat_unit = "cell"
  )

  scalefactors_list = .visium_read_scalefactors(
    json_path = scalefactor_path
  )

  visium_polygons = .visium_spot_poly(
    spatlocs = visium_spat_locs,
    json_scalefactors = scalefactors_list
  )

  gobject = addGiottoPolygons(
    gobject = gobject,
    gpolygons = list(visium_polygons)
  )

  return(gobject)
}





#' @title Read Visium ScaleFactors
#' @name .visium_read_scalefactors
#' @param json_path path to scalefactors_json.json for Visium experimental data
#' @return scalefactors within the provided json file as a named list,
#' or NULL if not discovered
#' @details asserts the existence of and reads in a .json file
#' containing scalefactors for Visium data in the expected format.
#' Returns NULL if no path is provided or if the file does not exist.
#' @keywords internal
.visium_read_scalefactors = function(json_path = NULL) {

  if (!checkmate::test_file_exists(json_path)) {
    if (!is.null(json_path)) {
      warning('scalefactors not discovered at: \n', json_path, call. = FALSE)
    }
    return(NULL)
  }

  json_scalefactors = jsonlite::read_json(json_path)

  # Intial assertion that json dimensions are appropriate
  checkmate::assert_list(
    x = json_scalefactors,
    types = 'numeric',
    min.len = 4L,
    max.len = 5L
  )

  expected_json_names = c(
    "regist_target_img_scalef", # NEW as of 2023
    "spot_diameter_fullres",
    "tissue_hires_scalef",
    "fiducial_diameter_fullres",
    "tissue_lowres_scalef"
  )

  # Visium assay with chemistry v2 contains an additional
  # keyword in the json file
  new_format_2023 = checkmate::test_list(
    x = json_scalefactors,
    types = 'numeric',
    len = 5L
  )

  # If the scalefactors are of size 4 (older assay), clip the new keyword
  if (!new_format_2023) expected_json_names = expected_json_names[2:5]

  if (!setequal(names(json_scalefactors), expected_json_names)) {
    warning(GiottoUtils::wrap_txt(
      'h5 scalefactors json names differ from expected.
       [Expected]:', expected_json_names, '\n',
      '[Actual]:', names(json_scalefactors)))
  }

  return (json_scalefactors)
}


#' @title Calculate Pixel to Micron Scalefactor
#' @name visium_micron_scalefactor
#' @param json_scalefactors list of scalefactors from .visium_read_scalefactors()
#' @return scale factor for converting pixel to micron
#' @details
#' Calculates pixel to micron scalefactor.
#' Visium xy coordinates are based on the fullres image
#' The values provided are directly usable for generating polygon information
#' or calculating the micron size relative to spatial coordinates for this set
#' of spatial information.
#' @keywords internal
.visium_micron_scale <- function(json_scalefactors) {
  # visium spots diameter                       : 55 micron
  # diameter of a spot at this spatial scaling  : scalefactor_list$spot_diameter_fullres
  px_to_micron <- 55 / json_scalefactors$spot_diameter_fullres
  return (px_to_micron)
}


#' @title Create Polygons for Visium Data
#' @name .visium_spot_poly
#' @param spatlocs spatial locations data.table or `spatLocsObj` containing
#' centroid locations of visium spots
#' @param json_scalefactors list of scalefactors from .visium_read_scalefactors()
#' @return giottoPolygon object
#' @details
#' Creates circular polygons for spatial representation of
#' Visium spots.
#' @keywords internal
#' @md
.visium_spot_poly <- function(spatlocs = NULL,
                             json_scalefactors) {

  if(inherits(spatlocs, "spatLocsObj")){
    spatlocs <- spatlocs[]
  }

  vis_spot_poly <- GiottoClass::circleVertices(
    radius = json_scalefactors$spot_diameter_fullres/2
  )

  GiottoClass::polyStamp(
    stamp_dt = vis_spot_poly,
    spatlocs = spatlocs,
    verbose = FALSE
  ) %>%
  createGiottoPolygonsFromDfr(calc_centroids = TRUE,
                              verbose = FALSE)

}






# json_info expects the list read output from .visium_read_scalefactors
# image_path should be expected to be full filepath
# should only be used when do_manual_adj (deprecated) is FALSE
.visium_image <- function(
    image_path,
    json_info = NULL,
    micron_scale = FALSE,
    verbose = NULL
) {

  # assume image already checked
  vmsg(.v = verbose, .initial = " - ", "found image")

  # 1. determine image scalefactor to use ----------------------------------- #
  if (!is.null(json_info)) checkmate::assert_list(json_info)
  png_name <- basename(image_path) # used for name pattern matching only

  if (is.null(json_info)) { # if none provided
    warning(wrap_txt(
      'No scalefactors json info provided.
      Visium image scale_factor defaulting to 1'
    ))
    scale_factor = 1

  } else { # if provided

    scale_factor <- NULL # initial value

    # determine type of visium image
    visium_img_type <- NULL
    possible_types <- c("lowres", "hires")
    for (img_type in possible_types) {
      if (grepl(img_type, png_name)) visium_img_type <- img_type
    }

    if (is.null(visium_img_type)) { # if not recognized visium image type
      .gstop(
        "\'image_path\' filename did not partial match either \'lowres\' or \'hires\'.
        Ensure specified image is either the Visium lowres or hires image and rename it accordingly"
      )
    }

    vmsg(.v = verbose, .initial = " - ",
         "found scalefactors.
         attempting automatic alignment for the",
         str_quote(visium_img_type), "image\n\n")

    scale_factor <- switch(
      visium_img_type,
      "lowres" = json_info[["tissue_lowres_scalef"]],
      "hires" = json_info[["tissue_hires_scalef"]]
    )
  }

  if (isTRUE(micron_scale)) {
    scale_factor <- scale_factor * .visium_micron_scale(json_info)
  }

  # 2. create image --------------------------------------------------------- #
  visium_img <- createGiottoLargeImage(
    raster_object = image_path,
    name = "image",
    negative_y = TRUE,
    scale_factor = (1/scale_factor)
  )

  visium_img_list = list(visium_img)
  names(visium_img_list) = c('image')

  return(visium_img_list)
}











## MERSCOPE ####


#' @title Create Vizgen MERSCOPE largeImage
#' @name createMerscopeLargeImage
#' @description
#' Read MERSCOPE stitched images as giottoLargeImage. Images will also be
#' transformed to match the spatial coordinate reference system of the paired
#' points and polygon data.
#' @param image_file character. Path to one or more MERSCOPE images to load
#' @param transforms_file character. Path to MERSCOPE transforms file. Usually
#' in the same folder as the images and named
#' 'micron_to_mosaic_pixel_transform.csv'
#' @param name character. name to assign the image. Multiple should be provided
#' if image_file is a list.
#' @export
createMerscopeLargeImage <- function(image_file,
                                     transforms_file,
                                     name = 'image') {

  checkmate::assert_character(transforms_file)
  tfsDT <- data.table::fread(transforms_file)
  if (inherits(image_file, "character")) {
    image_file <- as.list(image_file)
  }
  checkmate::assert_list(image_file)

  scalef <- c(1/tfsDT[[1,1]], 1/tfsDT[[2,2]])
  x_shift <- -tfsDT[[1,3]]/tfsDT[[1,1]]
  y_shift <- -tfsDT[[2,3]]/tfsDT[[2,2]]

  out <- lapply(seq_along(image_file), function(i) {
    gimg <- createGiottoLargeImage(
      raster_object = image_file[[i]],
      name = name[[i]],
      scale_factor = scalef,
      negative_y = FALSE
    )

    gimg <- spatShift(gimg, dx = x_shift, dy = y_shift)

    gimg@extent <- terra::ext(gimg@raster_object)
    return(gimg)
  })

  if (length(out) == 1L) {
    out <- unlist(out)
  }

  return(out)
}







#' @title Create Vizgen MERSCOPE Giotto Object
#' @name createGiottoMerscopeObject
#' @description Given the path to a MERSCOPE experiment directory, creates a Giotto
#' object.
#' @param merscope_dir full path to the exported merscope directory
#' @param data_to_use which of either the 'subcellular' or 'aggregate' information
#' to use for object creation
#' @param FOVs which FOVs to use when building the subcellular object. (default is NULL)
#' NULL loads all FOVs (very slow)
#' @param calculate_overlap whether to run \code{\link{calculateOverlapRaster}}
#' @param overlap_to_matrix whether to run \code{\link{overlapToMatrix}}
#' @param aggregate_stack whether to run \code{\link{aggregateStacks}}
#' @param aggregate_stack_param params to pass to \code{\link{aggregateStacks}}
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @return a giotto object
#' @export
#' @details
#' [\strong{Expected Directory}] This function generates a giotto object when given a
#' link to a MERSCOPE output directory. It expects the following items within the directory
#' where the \strong{bolded} portions are what this function matches against:
#' \itemize{
#'   \item{\strong{cell_boundaries} (folder .hdf5 files)}
#'   \item{\strong{images} (folder of .tif images and a scalefactor/transfrom table)}
#'   \item{\strong{cell_by_gene}.csv (file)}
#'   \item{cell_metadata\strong{fov_positions_file}.csv (file)}
#'   \item{detected_transcripts\strong{metadata_file}.csv (file)}
#' }
createGiottoMerscopeObject = function(merscope_dir,
                                      data_to_use = c('subcellular', 'aggregate'),
                                      FOVs = NULL,
                                      poly_z_indices = 1:7,
                                      calculate_overlap = TRUE,
                                      overlap_to_matrix = TRUE,
                                      aggregate_stack = TRUE,
                                      aggregate_stack_param = list(summarize_expression = 'sum',
                                                                   summarize_locations = 'mean',
                                                                   new_spat_unit = 'cell'),
                                      instructions = NULL,
                                      cores = NA,
                                      verbose = TRUE) {

  fovs = NULL

  # 0. setup
  merscope_dir = path.expand(merscope_dir)

  poly_z_indices = as.integer(poly_z_indices)
  if(any(poly_z_indices < 1)) stop(wrap_txt(
    'poly_z_indices is a vector of one or more integers starting from 1.',
    errWidth = TRUE
  ))

  # determine data to use
  data_to_use = match.arg(arg = data_to_use, choices = c('subcellular','aggregate'))

  # 1. test if folder structure exists and is as expected
  dir_items = .read_merscope_folder(merscope_dir = merscope_dir,
                                   data_to_use = data_to_use,
                                   cores = cores,
                                   verbose = verbose)

  # 2. load in directory items
  data_list = .load_merscope_folder(dir_items = dir_items,
                                   data_to_use = data_to_use,
                                   poly_z_indices = poly_z_indices,
                                   fovs = fovs,
                                   cores = cores,
                                   verbose = verbose)

  # 3. Create giotto object
  if(data_to_use == 'subcellular') {

    merscope_gobject = .createGiottoMerscopeObject_subcellular(data_list = data_list,
                                                              calculate_overlap = calculate_overlap,
                                                              overlap_to_matrix = overlap_to_matrix,
                                                              aggregate_stack = aggregate_stack,
                                                              aggregate_stack_param = aggregate_stack_param,
                                                              cores = cores,
                                                              verbose = verbose)

  } else if(data_to_use == 'aggregate') {

    merscope_gobject = .createGiottoMerscopeObject_aggregate(data_list = data_list,
                                                            cores = cores,
                                                            verbose = verbose)

  } else {
    stop(wrap_txt('data_to_use "', data_to_use, '" not implemented', sep = ''))
  }

  return(merscope_gobject)

}




#' @describeIn createGiottoMerscopeObject Create giotto object with 'subcellular' workflow
#' @param data_list list of loaded data from \code{\link{load_merscope_folder}}
#' @keywords internal
.createGiottoMerscopeObject_subcellular = function(data_list,
                                                  calculate_overlap = TRUE,
                                                  overlap_to_matrix = TRUE,
                                                  aggregate_stack = TRUE,
                                                  aggregate_stack_param = list(summarize_expression = 'sum',
                                                                               summarize_locations = 'mean',
                                                                               new_spat_unit = 'cell'),
                                                  cores = NA,
                                                  verbose = TRUE) {

  feat_coord = neg_coord = cellLabel_dir = instructions = NULL

  # unpack data_list
  poly_info = data_list$poly_info
  tx_dt = data_list$tx_dt
  micronToPixelScale = data_list$micronToPixelScale
  image_list = data_list$images

  # data.table vars
  gene = NULL

  # split tx_dt by expression and blank
  vmsg('Splitting detections by feature vs blank', .v = verbose)
  feat_id_all = tx_dt[, unique(gene)]
  blank_id = feat_id_all[grepl(pattern = 'Blank', feat_id_all)]
  feat_id = feat_id_all[!feat_id_all %in% blank_id]

  feat_dt = tx_dt[gene %in% feat_id,]
  blank_dt = tx_dt[gene %in% blank_id,]

  # extract transcript_id col and store as feature meta
  feat_meta = unique(feat_dt[, c('gene', 'transcript_id', 'barcode_id'), with = FALSE])
  blank_meta = unique(blank_dt[, c('gene', 'transcript_id', 'barcode_id'), with = FALSE])
  feat_dt[, c('transcript_id', 'barcode_id') := NULL]
  blank_dt[, c('transcript_id', 'barcode_id') := NULL]

  if(isTRUE(verbose)) {
    message('  > Features: ', feat_dt[, .N])
    message('  > Blanks: ', blank_dt[, .N])
  }

  # build giotto object
  vmsg('Building subcellular giotto object...', .v = verbose)
  z_sub = createGiottoObjectSubcellular(
    gpoints = list('rna' = feat_coord,
                   'neg_probe' = neg_coord),
    gpolygons = list('cell' = cellLabel_dir),
    polygon_mask_list_params = list(
      mask_method = 'guess',
      flip_vertical = TRUE,
      flip_horizontal = FALSE,
      shift_horizontal_step = FALSE
    ),
    instructions = instructions,
    cores = cores
  )

}




#' @describeIn createGiottoMerscopeObject Create giotto object with 'aggregate' workflow
#' @param data_list list of loaded data from \code{\link{load_merscope_folder}}
#' @keywords internal
.createGiottoMerscopeObject_aggregate = function(data_list,
                                                cores = NA,
                                                verbose = TRUE) {

  # unpack data_list
  micronToPixelScale = data_list$micronToPixelScale
  expr_dt = data_list$expr_dt
  cell_meta = data_list$expr_mat
  image_list = data_list$images

  # split expr_dt by expression and blank

  # feat_id_all =

}




## Spatial Genomics ####

#' @title Create Spatial Genomics Giotto Object
#' @name createSpatialGenomicsObject
#' @param sg_dir full path to the exported Spatial Genomics directory
#' @param instructions new instructions (e.g. result from createGiottoInstructions)
#' @description Given the path to a Spatial Genomics data directory, creates a
#' Giotto object.
#' @export
createSpatialGenomicsObject <- function(sg_dir = NULL,
                                        instructions = NULL) {
  # Find files in Spatial Genomics directory
  dapi = list.files(sg_dir, full.names = TRUE, pattern = 'DAPI')
  mask = list.files(sg_dir, full.names = TRUE, pattern = 'mask')
  tx = list.files(sg_dir, full.names = TRUE, pattern = 'transcript')
  # Create Polygons
  gpoly = createGiottoPolygonsFromMask(
    mask,
    shift_vertical_step = FALSE,
    shift_horizontal_step = FALSE,
    flip_horizontal = FALSE,
    flip_vertical = FALSE
  )
  # Create Points
  tx = data.table::fread(tx)
  gpoints = createGiottoPoints(tx)
  dim(tx)
  # Create object and add image
  gimg = createGiottoLargeImage(dapi, use_rast_ext = TRUE)
  sg = createGiottoObjectSubcellular(
    gpoints = list('rna' = gpoints),
    gpolygons = list('cell' = gpoly),
    instructions = instructions
  )
  sg = addGiottoLargeImage(sg, largeImages = list(image = gimg))
  # Return SG object
  return(sg)
}





## CosMx ####

#' @title Create Nanostring CosMx Giotto Object
#' @name createGiottoCosMxObject
#' @description Given the path to a CosMx experiment directory, creates a Giotto
#' object.
#' @param cosmx_dir full path to the exported cosmx directory
#' @param data_to_use which type(s) of expression data to build the gobject with
#' Default is \code{'all'} information available. \code{'subcellular'} loads the transcript
#' coordinates only. \code{'aggregate'} loads the provided aggregated expression matrix.
#' @param FOVs field of views to load (only affects subcellular data and images)
#' @param remove_background_polygon try to remove background polygon (default: FALSE)
#' @param background_algo algorithm to remove background polygon
#' @param remove_unvalid_polygons remove unvalid polygons (default: TRUE)
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @return a giotto object
#' @export
#' @details
#' [\strong{Expected Directory}] This function generates a giotto object when given a
#' link to a cosmx output directory. It expects the following items within the directory
#' where the \strong{bolded} portions are what this function matches against:
#' \itemize{
#'   \item{\strong{CellComposite} (folder of images)}
#'   \item{\strong{CellLabels} (folder of images)}
#'   \item{\strong{CellOverlay} (folder of images)}
#'   \item{\strong{CompartmentLabels} (folder of images)}
#'   \item{experimentname_\strong{exprMat_file}.csv (file)}
#'   \item{experimentname_\strong{fov_positions_file}.csv (file)}
#'   \item{experimentname_\strong{metadata_file}.csv (file)}
#'   \item{experimentname_\strong{tx_file}.csv (file)}
#' }
#'
#' [\strong{Workflows}] Workflow to use is accessed through the data_to_use param
#' \itemize{
#'   \item{'all' - loads and requires subcellular information from tx_file and fov_positions_file
#'   and also the existing aggregated information (expression, spatial locations, and metadata)
#'   from exprMat_file and metadata_file.}
#'   \item{'subcellular' - loads and requires subcellular information from tx_file and
#'   fov_positions_file only.}
#'   \item{'aggregate' - loads and requires the existing aggregate information (expression,
#'   spatial locations, and metadata) from exprMat_file and metadata_file.}
#' }
#'
#' [\strong{Images}] Images in the default CellComposite, CellLabels, CompartmentLabels, and CellOverlay
#' folders will be loaded as giotto largeImage objects in all workflows as long as they are available.
#' Additionally, CellComposite images will be converted to giotto image objects, making plotting with
#' these image objects more responsive when accessing them from a server.
#' \code{\link{showGiottoImageNames}} can be used to see the available images.
#'
#'
createGiottoCosMxObject = function(cosmx_dir = NULL,
                                   data_to_use = c('all','subcellular','aggregate'),
                                   remove_background_polygon = TRUE,
                                   background_algo = c('range'),
                                   remove_unvalid_polygons = TRUE,
                                   FOVs = NULL,
                                   instructions = NULL,
                                   cores = determine_cores(),
                                   verbose = TRUE) {

  # 0. setup
  cosmx_dir = path.expand(cosmx_dir)

  # determine data to use
  data_to_use = match.arg(arg = data_to_use, choices = c('all','subcellular','aggregate'))
  if(data_to_use %in% c('all', 'aggregate')) {
    stop(wrap_txt('Convenience workflows "all" and "aggregate" are not available yet'))
  }

  # Define for data.table
  fov = target = x_local_px = y_local_px = z = cell_ID = CenterX_global_px = CenterY_global_px =
    CenterX_local_px = CenterY_local_px = NULL


  # 1. test if folder structure exists and is as expected
  dir_items = .read_cosmx_folder(cosmx_dir = cosmx_dir,
                                verbose = verbose)


  # 2. load and create giotto object
  cosmx_gobject <- switch(data_to_use,
    "subcellular" = .createGiottoCosMxObject_subcellular(
      dir_items,
      FOVs = FOVs,
      remove_background_polygon = remove_background_polygon,
      background_algo = background_algo,
      remove_unvalid_polygons = remove_unvalid_polygons,
      cores = cores,
      verbose = verbose,
      instructions = instructions
    ),
    "aggregate" = .createGiottoCosMxObject_aggregate(
      dir_items,
      cores = cores,
      verbose = verbose,
      instructions = instructions
    ),
    "all" = .createGiottoCosMxObject_all(
      dir_items,
      FOVs = FOVs,
      remove_background_polygon = remove_background_polygon,
      background_algo = background_algo,
      remove_unvalid_polygons = remove_unvalid_polygons,
      cores = cores,
      verbose = verbose,
      instructions = instructions
    )
  )


  # load in subcellular information, subcellular FOV objects, then join


  # load in pre-generated aggregated expression matrix
  if(data_to_use == 'aggregate' | data_to_use == 'all') {

  }



  vmsg('done')
  return(cosmx_gobject)

}



#' @title Load and create a CosMx Giotto object from subcellular info
#' @name .createGiottoCosMxObject_subcellular
#' @inheritParams createGiottoCosMxObject
#' @keywords internal
.createGiottoCosMxObject_subcellular <- function(
    dir_items,
    FOVs = NULL,
    remove_background_polygon = TRUE,
    background_algo = c('range'),
    remove_unvalid_polygons = TRUE,
    cores,
    verbose = TRUE,
    instructions = NULL
) {

  target <- fov <- NULL

  # load tx detections and FOV offsets ------------------------------------- #
  data_list = .load_cosmx_folder_subcellular(
    dir_items = dir_items,
    FOVs = FOVs,
    cores = cores,
    verbose = verbose
  )

  # unpack data_list
  FOV_ID = data_list$FOV_ID
  fov_offset_file = data_list$fov_offset_file
  tx_coord_all = data_list$tx_coord_all

  # remove global xy values and cell_ID
  tx_coord_all[, c('x_global_px', 'y_global_px', 'cell_ID') := NULL]

  data.table::setcolorder(tx_coord_all, c('target', 'x_local_px', 'y_local_px', 'z', 'fov'))

  # feature detection type splitting --------------------------------------- #

  if(isTRUE(verbose)) wrap_msg('Splitting detections by feature vs neg probe')
  all_IDs = tx_coord_all[, unique(target)]
  neg_IDs = all_IDs[grepl(pattern = 'NegPrb', all_IDs)]
  feat_IDs = all_IDs[!all_IDs %in% neg_IDs]

  # split detections DT
  feat_coords_all = tx_coord_all[target %in% feat_IDs]
  neg_coords_all = tx_coord_all[target %in% neg_IDs]

  if(isTRUE(verbose)) {
    message('  > Features: ', feat_coords_all[, .N])
    message('  > NegProbes: ', neg_coords_all[, .N])
  }

  # FOV-based processing --------------------------------------------------- #

  fov_gobjects_list = lapply(FOV_ID, function(x) {

    # images --------------------------------------------------- #
    # build image paths
    if(isTRUE(verbose)) message('Loading image information...')

    composite_dir = Sys.glob(paths = file.path(dir_items$`CellComposite folder`, paste0('*',x, '*')))
    cellLabel_dir = Sys.glob(paths = file.path(dir_items$`CellLabels folder`, paste0('*',x, '*')))
    compartmentLabel_dir = Sys.glob(paths = file.path(dir_items$`CompartmentLabels folder`, paste0('*',x, '*')))
    cellOverlay_dir = Sys.glob(paths = file.path(dir_items$`CellOverlay folder`, paste0('*',x, '*')))

    # Missing warnings
    if(length(composite_dir) == 0) {
      warning('[ FOV ', x, ' ] No composite images found')
      composite_dir = NULL
    }
    if(length(cellLabel_dir) == 0) {
      stop('[ FOV ', x, ' ] No cell mask images found')
    } # cell masks are necessary
    if(length(compartmentLabel_dir) == 0) {
      warning('[ FOV ', x, ' ] No compartment label images found')
      compartmentLabel_dir = NULL
    }
    if(length(cellOverlay_dir) == 0) {
      warning('[ FOV ', x, ' ] No cell polygon overlay images found')
      cellOverlay_dir = NULL
    }

    if(isTRUE(verbose)) message('Image load done')

    if(isTRUE(verbose)) wrap_msg('[ FOV ', x, ']')


    # transcripts ---------------------------------------------- #
    # get FOV specific tx locations
    if(isTRUE(verbose)) wrap_msg('Assigning FOV feature detections...')


    # feature info
    coord_oldnames = c('target', 'x_local_px', 'y_local_px')
    coord_newnames = c('feat_ID', 'x', 'y')

    feat_coord = feat_coords_all[fov == as.numeric(x)]
    data.table::setnames(feat_coord, old = coord_oldnames, new = coord_newnames)
    # neg probe info
    neg_coord = neg_coords_all[fov == as.numeric(x)]
    data.table::setnames(neg_coord, old = coord_oldnames, new = coord_newnames)


    # build giotto object -------------------------------------- #
    if(isTRUE(verbose)) wrap_msg('Building subcellular giotto object...')
    fov_subset = createGiottoObjectSubcellular(
      gpoints = list('rna' = feat_coord,
                     'neg_probe' = neg_coord),
      gpolygons = list('cell' = cellLabel_dir),
      polygon_mask_list_params = list(
        mask_method = 'guess',
        flip_vertical = TRUE,
        flip_horizontal = FALSE,
        shift_horizontal_step = FALSE,
        remove_background_polygon = remove_background_polygon,
        background_algo = background_algo,
        remove_unvalid_polygons = remove_unvalid_polygons
      ),
      instructions = instructions,
      cores = cores
    )


    # find centroids as spatial locations ---------------------- #
    if(isTRUE(verbose)) wrap_msg('Finding polygon centroids as cell spatial locations...')
    fov_subset = addSpatialCentroidLocations(
      fov_subset,
      poly_info = 'cell',
      spat_loc_name = 'raw'
    )


    # create and add giotto image objects ---------------------- #
    if(isTRUE(verbose)) {
      message('Attaching image files...')
      print(composite_dir)
      print(cellOverlay_dir)
      print(compartmentLabel_dir)
    }

    gImage_list = list()

    # load image if files are found
    if(!is.null(composite_dir))
      gImage_list$composite = createGiottoLargeImage(
        raster_object = composite_dir,
        negative_y = FALSE,
        name = 'composite'
      )
    if(!is.null(cellOverlay_dir))
      gImage_list$overlay = createGiottoLargeImage(
        raster_object = cellOverlay_dir,
        negative_y = FALSE,
        name = 'overlay'
      )
    if(!is.null(compartmentLabel_dir))
      gImage_list$compartment = createGiottoLargeImage(
        raster_object = compartmentLabel_dir,
        negative_y = FALSE,
        name = 'compartment'
      ) #TODO



    if(length(gImage_list) > 0) {
      fov_subset = addGiottoImage(
        gobject = fov_subset,
        largeImages = gImage_list
      )

      # convert to MG for faster loading (particularly relevant for pulling from server)
      # TODO remove this
      fov_subset = convertGiottoLargeImageToMG(
        giottoLargeImage = gImage_list$composite,
        gobject = fov_subset,
        return_gobject = TRUE,
        verbose = FALSE
      )
      # fov_subset = convertGiottoLargeImageToMG(giottoLargeImage = gImage_list$overlay, gobject = fov_subset, return_gobject = TRUE)
      # fov_subset = convertGiottoLargeImageToMG(giottoLargeImage = gImage_list$compartment, gobject = fov_subset, return_gobject = TRUE)
    } else {
      message('No images found for fov')
    }


  }) #lapply end

  # returning -------------------------------------------------------------- #

  if(length(FOVs) == 1) {
    return(fov_gobjects_list[[1]])
  } else {
    # join giotto objects according to FOV positions file
    if(isTRUE(verbose)) message('Joining FOV gobjects...')
    new_gobj_names = paste0('fov', FOV_ID)
    id_match = match(as.numeric(FOV_ID), fov_offset_file$fov)
    x_shifts = fov_offset_file[id_match]$x_global_px
    y_shifts = fov_offset_file[id_match]$y_global_px

    # Join giotto objects
    cosmx_gobject = joinGiottoObjects(
      gobject_list = fov_gobjects_list,
      gobject_names = new_gobj_names,
      join_method = 'shift',
      x_shift = x_shifts,
      y_shift = y_shifts
    )
    return(cosmx_gobject)
  }

}



#' @title Load and create a CosMx Giotto object from aggregate info
#' @name .createGiottoCosMxObject_aggregate
#' @inheritParams createGiottoCosMxObject
#' @keywords internal
.createGiottoCosMxObject_aggregate = function(dir_items,
                                             cores,
                                             verbose = TRUE,
                                             instructions = NULL) {

  data_to_use = fov = NULL

  data_list = .load_cosmx_folder_aggregate(dir_items = dir_items,
                                          cores = cores,
                                          verbose = verbose)

  # unpack data_list
  spatlocs = data_list$spatlocs
  spatlocs_fov = data_list$spatlocs_fov
  metadata = data_list$metadata
  protM = data_list$protM
  spM = data_list$spM
  fov_shifts = data_list$fov_shifts


  # create standard gobject from aggregate matrix
  if(data_to_use == 'aggregate') {

    # Create aggregate gobject
    if(isTRUE(verbose)) message('Building giotto object...')
    cosmx_gobject = createGiottoObject(expression = list('raw' = spM, 'protein' = protM),
                                       cell_metadata = list('cell' = list('rna' = metadata,
                                                                          'protein' = metadata)),
                                       spatial_locs = spatlocs,
                                       instructions = instructions,
                                       cores = cores)


    # load in images
    img_ID = data.table::data.table(fov = fov_shifts[, fov],
                                    img_name = paste0('fov', sprintf('%03d', fov_shifts[, fov]), '-image'))

    if(isTRUE(verbose)) message('Attaching image files...')
    composite_dir = Sys.glob(paths = file.path(dir_items$`CellComposite folder`, paste0('/*')))
    cellLabel_dir = Sys.glob(paths = file.path(dir_items$`CellLabels folder`, paste0('/*')))
    compartmentLabel_dir = Sys.glob(paths = file.path(dir_items$`CompartmentLabels folder`, paste0('/*')))
    overlay_dir = Sys.glob(paths = file.path(dir_items$`CellOverlay folder`, paste0('/*')))

    if(length(cellLabel_imgList) > 0) cellLabel_imgList = lapply(cellLabel_dir, function(x) {createGiottoLargeImage(x,name = 'cellLabel',negative_y = TRUE)})
    if(length(composite_imgList) > 0) composite_imgList = lapply(composite_dir, function(x) {createGiottoLargeImage(x,name = 'composite',negative_y = TRUE)})
    if(length(compartmentLabel_dir) > 0) compartmentLabel_imgList = lapply(compartmentLabel_dir, function(x) {createGiottoLargeImage(x,name = 'composite',negative_y = TRUE)})
    if(length(overlay_dir) > 0) overlay_imgList = lapply(overlay_dir, function(x) {createGiottoLargeImage(x,name = 'composite',negative_y = TRUE)})



  }

}




#' @title Load and create a CosMx Giotto object from subcellular and aggregate info
#' @name .createGiottoCosMxObject_all
#' @param dir_items list of full directory paths from \code{.read_cosmx_folder}
#' @inheritParams createGiottoCosMxObject
#' @details Both \emph{subcellular} (subellular transcript detection information) and
#' \emph{aggregate} (aggregated detection count matrices by cell polygon from NanoString)
#' data will be loaded in. The two will be separated into 'cell' and 'cell_agg'
#' spatial units in order to denote the difference in origin of the two.
#' @seealso createGiottoCosMxObject .createGiottoCosMxObject_aggregate
#' .createGiottoCosMxObject_subcellular
#' @keywords internal
.createGiottoCosMxObject_all = function(dir_items,
                                       FOVs,
                                       remove_background_polygon = TRUE,
                                       background_algo = c('range'),
                                       remove_unvalid_polygons = TRUE,
                                       cores,
                                       verbose = TRUE,
                                       instructions = NULL,
                                       ...) {

  # 1. create subcellular giotto as spat_unit 'cell'
  cosmx_gobject = .createGiottoCosMxObject_subcellular(dir_items = dir_items,
                                                      FOVs = FOVs,
                                                      remove_background_polygon = remove_background_polygon,
                                                      background_algo = background_algo,
                                                      remove_unvalid_polygons = remove_unvalid_polygons,
                                                      cores = cores,
                                                      verbose = verbose,
                                                      instructions = instructions)

  # 2. load and append aggregated information in spat_unit 'cell_agg'
  agg_data = .load_cosmx_folder_aggregate(dir_items = dir_items,
                                         cores = cores,
                                         verbose = verbose)

  # unpack data_list
  spatlocs = agg_data$spatlocs
  spatlocs_fov = agg_data$spatlocs_fov
  metadata = agg_data$metadata
  protM = agg_data$protM
  spM = agg_data$spM

  # add in pre-generated aggregated expression matrix information for 'all' workflow

  # Add aggregate expression information
  if(isTRUE(verbose)) wrap_msg('Appending provided aggregate expression data as...
                               spat_unit: "cell_agg"
                               feat_type: "rna"
                               name: "raw"')
  # add expression data to expression slot
  s4_expr = createExprObj(
    name = 'raw',
    expression_data = spM,
    spat_unit = 'cell_agg',
    feat_type = 'rna',
    provenance = 'cell_agg'
  )

  cosmx_gobject = set_expression_values(cosmx_gobject, values = s4_expr)

  # Add spatial locations
  if(isTRUE(verbose)) wrap_msg('Appending metadata provided spatial locations data as...
                               --> spat_unit: "cell_agg" name: "raw"
                               --> spat_unit: "cell" name: "raw_fov"')
  if(isTRUE(verbose)) wrap_msg('Polygon centroid derived spatial locations assigned as...
                               --> spat_unit: "cell" name: "raw" (default)')

  locsObj = create_spat_locs_obj(name = 'raw',
                                 coordinates = spatlocs,
                                 spat_unit = 'cell_agg',
                                 provenance = 'cell_agg')
  locsObj_fov = create_spat_locs_obj(name = 'raw_fov',
                                     coordinates = spatlocs_fov,
                                     spat_unit = 'cell_agg',
                                     provenance = 'cell_agg')

  cosmx_gobject = set_spatial_locations(cosmx_gobject, spatlocs = locsObj)
  cosmx_gobject = set_spatial_locations(cosmx_gobject, spatlocs = locsObj_fov)

  # cosmx_gobject = set_spatial_locations(cosmx_gobject,
  #                                       spat_unit = 'cell_agg',
  #                                       spat_loc_name = 'raw',
  #                                       spatlocs = spatlocs)
  # cosmx_gobject = set_spatial_locations(cosmx_gobject,
  #                                       spat_unit = 'cell_agg',
  #                                       spat_loc_name = 'raw_fov',
  #                                       spatlocs = spatlocs_fov)

  # initialize cell and feat IDs and metadata slots for 'cell_agg' spat_unit
  agg_cell_ID = colnames(s4_expr[])
  agg_feat_ID = rownames(s4_expr[])

  sub_feat_ID <- featIDs(cosmx_gobject, feat_type = "rna")
  feat_ID_new = unique(c(agg_feat_ID, sub_feat_ID))

  # cosmx_gobject = set_cell_id(gobject = cosmx_gobject,
  #                             spat_unit = 'cell_agg',
  #                             cell_IDs = agg_cell_ID)
  # cosmx_gobject = set_feat_id(gobject = cosmx_gobject,
  #                             feat_type = 'rna',
  #                             feat_IDs = feat_ID_new)

  # cell metadata

  # Add metadata to both the given and the poly spat_units
  if(isTRUE(verbose)) message('Appending provided cell metadata...')
  cosmx_gobject = addCellMetadata(cosmx_gobject,
                                  spat_unit = 'cell',
                                  feat_type = 'rna',
                                  new_metadata = metadata,
                                  by_column = TRUE,
                                  column_cell_ID = 'cell_ID')
  cosmx_gobject = addCellMetadata(cosmx_gobject,
                                  spat_unit = 'cell_agg',
                                  feat_type = 'rna',
                                  new_metadata = metadata,
                                  by_column = TRUE,
                                  column_cell_ID = 'cell_ID')

  initialize(cosmx_gobject)
}










## Xenium ####

#' @title Create 10x Xenium Giotto Object
#' @name createGiottoXeniumObject
#' @description Given the path to a Xenium experiment output folder, creates a Giotto
#' object
#' @param xenium_dir full path to the exported xenium directory
#' @param data_to_use which type(s) of expression data to build the gobject with
#' (e.g. default: \strong{'subcellular'}, 'aggregate', or 'all')
#' @param load_format files formats from which to load the data. Either `csv` or
#' `parquet` currently supported.
#' @param h5_expression (boolean) whether to load cell_feature_matrix from .h5 file.
#' Default is \code{TRUE}
#' @param h5_gene_ids use gene symbols (default) or ensembl ids for the .h5 gene
#' expression matrix
#' @param bounds_to_load vector of boundary information to load (e.g. \code{'cell'}
#' or \code{'nucleus'} by themselves or \code{c('cell', 'nucleus')} to load both
#' at the same time.)
#' @param qv_threshold Minimum Phred-scaled quality score cutoff to be included as
#' a subcellular transcript detection (default = 20)
#' @param key_list (advanced) list of grep-based keywords to split the subcellular
#' feature detections by feature type. See details
#' @inheritParams get10Xmatrix
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @details
#'
#' [\strong{QC feature types}]
#' Xenium provides info on feature detections that include more than only the
#' Gene Expression specific probes. Additional probes for QC are included:
#' \emph{blank codeword}, \emph{negative control codeword}, and
#' \emph{negative control probe}. These additional QC probes each occupy and are treated
#' as their own feature types so that they can largely remain independent of the
#' gene expression information.
#'
#' [\strong{key_list}]
#' Related to \code{data_to_use = 'subcellular'} workflow only:
#' Additional QC probe information is in the subcellular feature detections information
#' and must be separated from the gene expression information during processing.
#' The QC probes have prefixes that allow them to be selected from the rest of the
#' feature IDs.
#' Giotto uses a named list of keywords (\code{key_list}) to select these QC probes,
#' with the list names being the names that will be assigned as the feature type
#' of these feature detections. The default list is used when \code{key_list} = NULL.
#'
#' Default list:
#' \preformatted{
#'  list(blank_code = 'BLANK_',
#'       neg_code = 'NegControlCodeword_',
#'       neg_probe = c('NegControlProbe_|antisense_'))
#' }
#'
#' The Gene expression subset is accepted as the subset of feat_IDs that do not
#' map to any of the keys.
#'
#' @export
createGiottoXeniumObject = function(xenium_dir,
                                    data_to_use = c('subcellular','aggregate'),
                                    load_format = 'csv',
                                    h5_expression = TRUE,
                                    h5_gene_ids = c('symbols', 'ensembl'),
                                    gene_column_index = 1,
                                    bounds_to_load = c('cell'),
                                    qv_threshold = 20,
                                    key_list = NULL,
                                    # include_analysis = FALSE,
                                    instructions = NULL,
                                    cores = NA,
                                    verbose = TRUE) {

  # 0. setup
  xenium_dir = path.expand(xenium_dir)

  # Determine data to load
  data_to_use = match.arg(arg = data_to_use, choices = c('subcellular','aggregate'))

  # Determine load formats
  load_format = 'csv' # TODO Remove this and add as param once other options are available
  load_format = match.arg(arg = load_format, choices = c('csv', 'parquet', 'zarr'))

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # 1. detect xenium folder and find filepaths to load

  # path_list contents:
  # tx_path
  # bound_paths
  # cell_meta_path
  # agg_expr_path
  # panel_meta_path
  path_list = .read_xenium_folder(xenium_dir = xenium_dir,
                                 data_to_use = data_to_use,
                                 bounds_to_load = bounds_to_load,
                                 load_format = load_format,
                                 h5_expression = h5_expression,
                                 verbose = verbose)


  # 2. load in data

  # data_list contents:
  # feat_meta
  # tx_dt
  # bound_dt_list
  # cell_meta
  # agg_expr
  data_list = .load_xenium_folder(path_list = path_list,
                                 load_format = load_format,
                                 data_to_use = data_to_use,
                                 h5_expression = h5_expression,
                                 h5_gene_ids = h5_gene_ids,
                                 gene_column_index = gene_column_index,
                                 cores = cores,
                                 verbose = verbose)


  # TODO load images


  # 3. Create giotto objects

  if(data_to_use == 'subcellular') {

    # ** feat type search keys **
    if(is.null(key_list)) {
      key_list = list(blank_code = 'BLANK_',
                      neg_code = 'NegControlCodeword_',
                      neg_probe = c('NegControlProbe_|antisense_'))
    }

    # needed:
    # feat_meta
    # tx_dt
    # bound_dt_list
    xenium_gobject = .createGiottoXeniumObject_subcellular(data_list = data_list,
                                                          qv_threshold = qv_threshold,
                                                          key_list = key_list,
                                                          instructions = instructions,
                                                          cores = cores,
                                                          verbose = verbose)

  }

  if(data_to_use == 'aggregate') {

    # needed:
    # feat_meta
    # cell_meta
    # agg_expr
    # optional?
    # tx_dt
    # bound_dt_list
    xenium_gobject = .createGiottoXeniumObject_aggregate(data_list = data_list,
                                                        instructions = instructions,
                                                        cores = cores,
                                                        verbose = verbose)

  }

  return(xenium_gobject)

}




#' @title Create a Xenium Giotto object from subcellular info
#' @name .createGiottoXeniumObject_subcellular
#' @description Subcellular workflow for createGiottoXeniumObject
#' @param data_list list of data loaded by \code{\link{.load_xenium_folder}}
#' @param key_list regex-based search keys for feature IDs to allow separation
#' into separate giottoPoints objects by feat_type
#' @param qv_threshold Minimum Phred-scaled quality score cutoff to be included as
#' a subcellular transcript detection (default = 20)
#' @inheritParams get10Xmatrix
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @seealso createGiottoXeniumObject .createGiottoXeniumObject_aggregate
#' @keywords internal
.createGiottoXeniumObject_subcellular = function(data_list,
                                                key_list = NULL,
                                                qv_threshold = 20,
                                                instructions = NULL,
                                                cores = NA,
                                                verbose = TRUE) {

  # data.table vars
  qv = NULL

  # Unpack data_list info
  feat_meta = data_list$feat_meta
  tx_dt = data_list$tx_dt
  bound_dt_list = data_list$bound_dt_list
  # cell_meta = data_list$cell_meta
  # agg_expr = data_list$agg_expr

  # define for data.table
  cell_id = feat_ID = feature_name = NULL

  vmsg('Building subcellular giotto object...', .v = verbose)
  # Giotto points object
  vmsg('> points data prep...', .v = verbose)

  # filter by qv_threshold
  vmsg('> filtering feature detections for Phred score >= ', qv_threshold, .v = verbose)
  n_before = tx_dt[,.N]
  tx_dt_filtered = tx_dt[qv >= qv_threshold]
  n_after = tx_dt_filtered[,.N]

  if(verbose) {
    cat('Number of feature points removed: ',
        n_before - n_after,
        ' out of ', n_before, '\n')
  }

  vmsg('> splitting detections by feat_type', .v = verbose)
  # discover feat_IDs for each feat_type
  all_IDs = tx_dt_filtered[, unique(feat_ID)]
  feat_types_IDs = lapply(key_list, function(x) all_IDs[grepl(pattern = x, all_IDs)])
  rna = list('rna' = all_IDs[!all_IDs %in% unlist(feat_types_IDs)])
  feat_types_IDs = append(rna, feat_types_IDs)

  # separate detections by feature type
  points_list = lapply(
    feat_types_IDs,
    function(types) {
      tx_dt_filtered[feat_ID %in% types]
    }
  )

  # Giotto polygons object
  vmsg('> polygons data prep...', .v = verbose)
  polys_list = lapply(
    bound_dt_list,
    function(bound_type) {
      bound_type[, cell_id := as.character(cell_id)]
    }
  )

  xenium_gobject = createGiottoObjectSubcellular(
    gpoints = points_list,
    gpolygons = polys_list,
    instructions = instructions,
    cores = cores,
    verbose = verbose
  )

  # generate centroids
  vmsg('Calculating polygon centroids...', .v = verbose)
  xenium_gobject = addSpatialCentroidLocations(
    xenium_gobject,
    poly_info = c(names(bound_dt_list)),
    provenance = as.list(names(bound_dt_list))
  )

  # add in feature metadata
  # xenium_gobject = addFeatMetadata(gobject = xenium_gobject,
  #                                  new_metadata = feat_meta,
  #                                  by_column = TRUE,
  #                                  column_feat_ID = 'feat_ID')

  return(xenium_gobject)

}





#' @title Create a Xenium Giotto object from aggregate info
#' @name .createGiottoXeniumObject_aggregate
#' @description Aggregate workflow for createGiottoXeniumObject
#' @param data_list list of data loaded by \code{.load_xenium_folder}
#' @inheritParams get10Xmatrix
#' @inheritParams GiottoClass::createGiottoObjectSubcellular
#' @seealso createGiottoXeniumObject .createGiottoXeniumObject_subcellular
#' @keywords internal
.createGiottoXeniumObject_aggregate = function(data_list,
                                              # include_analysis = FALSE,
                                              instructions = NULL,
                                              cores = NA,
                                              verbose = TRUE) {

  # Unpack data_list info
  feat_meta = data_list$feat_meta
  # tx_dt = data_list$tx_dt
  # bound_dt_list = data_list$bound_dt_list
  cell_meta = data_list$cell_meta
  agg_expr = data_list$agg_expr

  # define for data.table
  cell_ID = x_centroid = y_centroid = NULL

  # clean up names for aggregate matrices
  names(agg_expr) = gsub(pattern = ' ', replacement = '_' ,names(agg_expr))
  geneExpMat = which(names(agg_expr) == 'Gene_Expression')
  names(agg_expr)[[geneExpMat]] = 'raw'

  # set cell_id as character
  cell_meta = cell_meta[, data.table::setnames(.SD, 'cell_id', 'cell_ID')]
  cell_meta = cell_meta[, cell_ID := as.character(cell_ID)]

  # set up spatial locations
  agg_spatlocs = cell_meta[, .(x_centroid, y_centroid, cell_ID)]

  # set up metadata
  agg_meta = cell_meta[, !c('x_centroid','y_centroid')]

  vmsg('Building aggregate giotto object...', .v = verbose)
  xenium_gobject = createGiottoObject(expression = agg_expr,
                                      spatial_locs = agg_spatlocs,
                                      instructions = instructions,
                                      cores = cores,
                                      verbose = verbose)

  # append aggregate metadata
  xenium_gobject = addCellMetadata(gobject = xenium_gobject,
                                   new_metadata = agg_meta,
                                   by_column = TRUE,
                                   column_cell_ID = 'cell_ID')
  xenium_gobject = addFeatMetadata(gobject = xenium_gobject,
                                   new_metadata = feat_meta,
                                   by_column = TRUE,
                                   column_feat_ID = 'feat_ID')

  return(xenium_gobject)

}







# *---- folder reading and detection ----* ####


#' @describeIn read_data_folder Read a structured MERSCOPE folder
#' @keywords internal
  .read_merscope_folder <- function(merscope_dir,
                                    data_to_use,
                                    cores = NA,
                                    verbose = NULL) {

  # prepare dir_items list
  dir_items = list(`boundary info` = '*cell_boundaries*',
                   `image info` = '*images*',
                   `cell feature matrix` = '*cell_by_gene*',
                   `cell metadata` = '*cell_metadata*',
                   `raw transcript info` = '*transcripts*')

  # prepare require_data_DT
  sub_reqs = data.table::data.table(workflow = c('subcellular'),
                                    item = c('boundary info',
                                             'raw transcript info',
                                             'image info',
                                             'cell by gene matrix',
                                             'cell metadata'),
                                    needed = c(TRUE, TRUE, FALSE, FALSE, FALSE))

  agg_reqs = data.table::data.table(workflow = c('aggregate'),
                                    item = c('boundary info',
                                             'raw transcript info',
                                             'image info',
                                             'cell by gene matrix',
                                             'cell metadata'),
                                    needed = c(FALSE, FALSE, FALSE, TRUE, TRUE))

  require_data_DT = rbind(sub_reqs, agg_reqs)

  dir_items = .read_data_folder(spat_method = 'MERSCOPE',
                               data_dir = merscope_dir,
                               dir_items = dir_items,
                               data_to_use = data_to_use,
                               require_data_DT = require_data_DT,
                               cores = cores,
                               verbose = verbose)

  return(dir_items)

}



#' @title Read a structured CosMx folder
#' @name .read_cosmx_folder
#' @inheritParams createGiottoCosMxObject
#' @seealso createGiottoCosMxObject load_cosmx_folder
#' @return path_list a list of cosmx files discovered and their filepaths. NULL
#' values denote missing items
#' @keywords internal
.read_cosmx_folder = function(cosmx_dir,
                             verbose = TRUE) {

  ch = box_chars()

  if(is.null(cosmx_dir) | !dir.exists(cosmx_dir)) stop('The full path to a cosmx directory must be given.\n')
  vmsg('A structured CosMx directory will be used\n', .v = verbose)

  # find directories (length = 1 if present, length = 0 if missing)
  dir_items = list(`CellLabels folder` = '*CellLabels',
                   `CompartmentLabels folder` = '*CompartmentLabels',
                   `CellComposite folder` = '*CellComposite',
                   `CellOverlay folder` = '*CellOverlay',
                   `transcript locations file` = '*tx_file*',
                   `fov positions file` = '*fov_positions_file*',
                   `expression matrix file` = '*exprMat_file*',
                   `metadata file` = '*metadata_file*')
  dir_items = lapply(dir_items, function(x) Sys.glob(paths = file.path(cosmx_dir, x)))
  dir_items_lengths = lengths(dir_items)

  if(isTRUE(verbose)) {
    message('Checking directory contents...')
    for(item in names(dir_items)) {
      if(dir_items_lengths[[item]] > 0) {
        message(ch$s, '> ' ,item, ' found')
      } else {
        warning(item, ' is missing\n')
      }
    }
  }

  # select first directory in list if multiple are detected
  if(any(dir_items_lengths > 1)) {
    warning('Multiple matches for expected subdirectory item(s).\n First matching item selected')

    multiples = which(dir_items_lengths > 1)
    for(mult_i in multiples) {
      message(names(dir_items)[[mult_i]], 'multiple matches found:')
      print(dir_items[[mult_i]])
      dir_items[[mult_i]] = dir_items[[mult_i]][[1]]
    }
  }
  vmsg('Directory check done', .v = verbose)

  return(dir_items)

}




#' @title Read a structured xenium folder
#' @name .read_xenium_folder
#' @inheritParams createGiottoXeniumObject
#' @keywords internal
#' @return path_list a list of xenium files discovered and their filepaths. NULL
#' values denote missing items
.read_xenium_folder = function(xenium_dir,
                              data_to_use = 'subcellular',
                              bounds_to_load = c('cell'),
                              load_format = 'csv',
                              h5_expression = FALSE,
                              verbose = TRUE) {

  # Check needed packages
  if(load_format == 'parquet') {
    package_check(pkg_name = 'arrow', repository = 'CRAN')
    package_check(pkg_name = 'dplyr', repository = 'CRAN')
  }
  if(isTRUE(h5_expression)) {
    package_check(pkg_name = 'hdf5r', repository = 'CRAN')
  }

  ch = box_chars()


  # 0. test if folder structure exists and is as expected


  if(is.null(xenium_dir) | !dir.exists(xenium_dir)) stop('The full path to a xenium directory must be given.\n')
  vmsg('A structured Xenium directory will be used\n', .v = verbose)

  # find items (length = 1 if present, length = 0 if missing)
  dir_items = list(`analysis info` = '*analysis*',
                   `boundary info` = '*bound*',
                   `cell feature matrix` = '*cell_feature_matrix*',
                   `cell metadata` = '*cells*',
                   `image info` = '*tif',
                   `panel metadata` = '*panel*',
                   `raw transcript info` = '*transcripts*',
                   `experiment info (.xenium)` = '*.xenium')

  dir_items = lapply(dir_items, function(x) Sys.glob(paths = file.path(xenium_dir, x)))
  dir_items_lengths = lengths(dir_items)

  if(isTRUE(verbose)) {
    message('Checking directory contents...')
    for(item in names(dir_items)) {

      # IF ITEM FOUND

      if(dir_items_lengths[[item]] > 0) {
        message(ch$s, '> ' ,item, ' found')
        for(item_i in seq_along(dir_items[[item]])) { # print found item names
          subItem = gsub(pattern = '.*/', replacement = '', x = dir_items[[item]][[item_i]])
          message(ch$s, ch$s, ch$l,ch$h,ch$h, subItem)
        }
      } else {

        # IF ITEM MISSING
        # Based on workflow, determine if:
        # necessary (error)
        # optional (warning)

        if(data_to_use == 'subcellular') {
          # necessary items
          if(item %in% c('boundary info', 'raw transcript info')) stop(item, ' is missing\n')
          # optional items
          if(item %in% c('image info', 'experiment info (.xenium)', 'panel metadata')) warning(item, ' is missing (optional)\n')
          # items to ignore: analysis info, cell feature matrix, cell metadata
        } else if(data_to_use == 'aggregate') {
          # necessary items
          if(item %in% c('cell feature matrix', 'cell metadata')) stop(item, ' is missing\n')
          # optional items
          if(item %in% c('image info', 'experiment info (.xenium)', 'panel metadata', 'analysis info')) warning(item, ' is missing (optional)\n')
          # items to ignore: boundary info, raw transcript info
        }
      }
    }
  }


  # 1. Select data to load


  # **** transcript info ****
  tx_path = NULL
  tx_path = dir_items$`raw transcript info`[grepl(pattern = load_format, dir_items$`raw transcript info`)]
  # **** cell metadata ****
  cell_meta_path = NULL
  cell_meta_path = dir_items$`cell metadata`[grepl(pattern = load_format, dir_items$`cell metadata`)]

  # **** boundary info ****
  # Select bound load format
  if(load_format != 'zarr') { # No zarr available for boundary info
    dir_items$`boundary info` = dir_items$`boundary info`[grepl(pattern = load_format, dir_items$`boundary info`)]
  } else dir_items$`boundary info` = dir_items$`boundary info`[grepl(pattern = 'csv', dir_items$`boundary info`)]

  # Organize bound paths by type of bound (bounds_to_load param)
  bound_paths = NULL
  bound_names = bounds_to_load
  bounds_to_load = as.list(bounds_to_load)
  bound_paths = lapply(bounds_to_load, function(x) dir_items$`boundary info`[grepl(pattern = x, dir_items$`boundary info`)])
  names(bound_paths) = bound_names

  # **** aggregated expression info ****
  agg_expr_path = NULL
  if(isTRUE(h5_expression)) { # h5 expression matrix loading is default
    agg_expr_path = dir_items$`cell feature matrix`[grepl(pattern = 'h5', dir_items$`cell feature matrix`)]
  } else if(load_format == 'zarr') {
    agg_expr_path = dir_items$`cell feature matrix`[grepl(pattern = 'zarr', dir_items$`cell feature matrix`)]
  } else { # No parquet for aggregated expression - default to normal 10x loading
    agg_expr_path = dir_items$`cell feature matrix`[sapply(dir_items$`cell feature matrix`, function(x) file_test(op = '-d', x))]
    if(length(agg_expr_path) == 0) stop(wrap_txt(
      'Expression matrix cannot be loaded.\nHas cell_feature_matrix(.tar.gz) been unpacked into a directory?'
    ))
  }
  if(data_to_use == 'aggregate') {
    if(length(path_list$agg_expr_path) == 0) stop(wrap_txt(
      'Aggregated expression not found.\nPlease confirm h5_expression and load_format params are correct\n'
    ))
  }

  # **** panel info ****
  panel_meta_path = NULL
  panel_meta_path = dir_items$`panel metadata`


  vmsg('Directory check done', .v = verbose)

  path_list = list('tx_path' = tx_path,
                   'bound_paths' = bound_paths,
                   'cell_meta_path' = cell_meta_path,
                   'agg_expr_path' = agg_expr_path,
                   'panel_meta_path' = panel_meta_path)

  return(path_list)

}






# * ---- folder loading ---- * ####



## MERSCOPE ####

#' @title Load MERSCOPE data from folder
#' @name load_merscope_folder
#' @param dir_items list of full filepaths from \code{\link{.read_merscope_folder}}
#' @inheritParams createGiottoMerscopeObject
#' @return list of loaded-in MERSCOPE data
NULL

#' @rdname load_merscope_folder
#' @keywords internal
.load_merscope_folder = function(dir_items,
                                data_to_use,
                                fovs = NULL,
                                poly_z_indices = 1L:7L,
                                cores = NA,
                                verbose = TRUE) {

  # 1. load data_to_use-specific
  if(data_to_use == 'subcellular') {
    data_list = .load_merscope_folder_subcellular(dir_items = dir_items,
                                                 data_to_use = data_to_use,
                                                 fovs = fovs,
                                                 poly_z_indices = poly_z_indices,
                                                 cores = cores,
                                                 verbose = verbose)
  } else if(data_to_use == 'aggregate') {
    data_list = .load_merscope_folder_aggregate(dir_items = dir_items,
                                               data_to_use = data_to_use,
                                               cores = cores,
                                               verbose = verbose)
  } else {
    stop(wrap_txt('data_to_use "', data_to_use, '" not implemented', sep = ''))
  }

  # 2. Load images if available
  if(!is.null(dir_items$`image info`)) {
    ## micron to px scaling factor
    micronToPixelScale = Sys.glob(paths = file.path(dir_items$`image info`, '*micron_to_mosaic_pixel_transform*'))[[1]]
    micronToPixelScale = data.table::fread(micronToPixelScale, nThread = cores)
    # add to data_list
    data_list$micronToPixelScale = micronToPixelScale

    ## staining images
    ## determine types of stains
    images_filenames = list.files(dir_items$`image info`)
    bound_stains_filenames = images_filenames[grep(pattern = '.tif', images_filenames)]
    bound_stains_types = sapply(strsplit(bound_stains_filenames, '_'), `[`, 2)
    bound_stains_types = unique(bound_stains_types)

    img_list = lapply_flex(bound_stains_types, function(stype) {
      img_paths = Sys.glob(paths = file.path(dir_items$`image info`, paste0('*',stype,'*')))

      lapply_flex(img_paths, function(img) {
        createGiottoLargeImage(raster_object = img)
      }, cores = cores)
    }, cores = cores)
    # add to data_list
    data_list$images = img_list
  }



  return(data_list)

}



#' @describeIn load_merscope_folder Load items for 'subcellular' workflow
#' @keywords internal
.load_merscope_folder_subcellular = function(dir_items,
                                            data_to_use,
                                            cores = NA,
                                            poly_z_indices = 1L:7L,
                                            verbose = TRUE,
                                            fovs = NULL) {

  if(isTRUE(verbose)) wrap_msg('Loading transcript level info...')
  if(is.null(fovs)) {
    tx_dt = data.table::fread(dir_items$`raw transcript info`, nThread = cores)
  } else {
    vmsg('Selecting FOV subset transcripts')
    tx_dt = fread_colmatch(file = dir_items$`raw transcript info`,
                           col = 'fov',
                           values_to_match = fovs,
                           verbose = FALSE,
                           nThread = cores)
  }
  tx_dt[, c('x','y') := NULL] # remove unneeded cols
  data.table::setcolorder(tx_dt, c('gene', 'global_x', 'global_y', 'global_z'))

  if(isTRUE(verbose)) wrap_msg('Loading polygon info...')
  poly_info = readPolygonFilesVizgenHDF5(boundaries_path = dir_items$`boundary info`,
                                         z_indices = poly_z_indices,
                                         flip_y_axis = TRUE,
                                         fovs = fovs)

  data_list = list(
    'poly_info' = poly_info,
    'tx_dt' = tx_dt,
    'micronToPixelScale' = NULL,
    'expr_dt' = NULL,
    'cell_meta' = NULL,
    'images' = NULL
  )

}



#' @describeIn load_merscope_folder Load items for 'aggregate' workflow
#' @keywords internal
.load_merscope_folder_aggregate = function(dir_items,
                                          data_to_use,
                                          cores = NA,
                                          verbose = TRUE) {

  # metadata is polygon-related measurements
  vmsg('Loading cell metadata...', .v = verbose)
  cell_metadata_file = data.table::fread(dir_items$`cell metadata`, nThread = cores)

  vmsg('Loading expression matrix', .v = verbose)
  expr_dt = data.table::fread(dir_items$`cell feature matrix`, nThread = cores)


  data_list = list(
    'poly_info' = NULL,
    'tx_dt' = NULL,
    'micronToPixelScale' = NULL,
    'expr_dt' = expr_dt,
    'cell_meta' = cell_metadata_file,
    'images' = NULL
  )

}







## CosMx ####



.cosmx_transcript <- function(
        path,
        fovs = NULL,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        dropcols = c(
            "x_local_px",
            "y_local_px",
            "cell_ID",
            "cell"
        ),
        mm = FALSE,
        px2mm = 0.12028,
        cores = determine_cores(),
        verbose = NULL
    ) {

    if (missing(path)) {
        stop(wrap_txt(
            "No path to tx file provided or auto-detected"
        ), call. = FALSE)
    }

    checkmate::assert_file_exists(path)

    vmsg(.v = verbose, "loading feature detections...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    tx <- data.table::fread(input = path, nThread = cores, drop = dropcols)
    if (!is.null(fovs)) {
        # subset to only needed FOVs
        tx <- tx[fov %in% as.numeric(fovs),]
    }

    # mm scaling if desired
    if (mm) {
        tx[, x_global_px := x_global_px * px2mm]
        tx[, y_global_px := y_global_px * px2mm]
    }

    # giottoPoints ----------------------------------------------------- #

    # static gpoints params
    gpoints_params <- list()
    gpoints_params$feat_type <- feat_type
    gpoints_params$split_keyword <- split_keyword
    gpoints_params$x_colname <- "x_global_px"
    gpoints_params$y_colname <- "y_global_px"
    gpoints_params$feat_ID_colname <- "target"

    gpoints <- do.call(createGiottoPoints, c(list(x = tx), gpoints_params))
    # ensure output is always a list
    if (!is.list(gpoints)) {
        gpoints <- list(gpoints)
        names(gpoints) <- objName(gpoints[[1L]])
    }

    return(gpoints)
}

#' @returns data.table with three columns. 1. FOV (integer), xshift (numeric),
#' yshift (numeric). Values should always be in pixels
.cosmx_infer_fov_shifts <- function(tx_dt, meta_dt, flip_loc_y = NULL) {
    fov <- NULL # NSE vars

    if (!missing(tx_dt)) {
        flip_loc_y %null% TRUE # default = TRUE
        tx_head <- tx_dt[, head(.SD, 10L), by = fov]
        x <- tx_head[, mean(x_global_px - x_local_px), by = fov]
        if (flip_loc_y) {
            # use +y if local y values are flipped
            y <- tx_head[, mean(y_global_px + y_local_px), by = fov]
        } else {
            y <- tx_head[, mean(y_global_px - y_local_px), by = fov]
        }
    }

    if (!missing(meta_dt)) {
        flip_loc_y %null% FALSE # default = FALSE
        meta_head <- meta_dt[, head(.SD, 10L), by = fov]
        x <- meta_head[, mean(CenterX_global_px - CenterX_local_px), by = fov]
        if (flip_loc_y) {
            # use +y if local y values are flipped
            y <- meta_head[, mean(CenterY_global_px + CenterY_local_px),
                           by = fov]
        } else {
            y <- meta_head[, mean(CenterY_global_px - CenterY_local_px),
                           by = fov]
        }
    }

    res <- merge(x, y, by = "fov")
    data.table::setnames(res, new = c("fov", "x", "y"))

    return(res)
}

.cosmx_imgname_fovparser <- function(
        path
) {
    im_names <- list.files(path)
    fovs <- as.numeric(sub(".*F(\\d+)\\..*", "\\1", im_names))
    if (any(is.na(fovs))) {
        warning(wrap_txt(
            "Images to load should be sets of images/fov in subdirectories.
            No other files should be present."
        ))
    }
    return(fovs)
}

.cosmx_poly <- function(
        path,
        slide = 1,
        fovs = NULL,
        name = "cell",
        # VERTICAL FLIP + NO SHIFTS
        flip_vertical = TRUE,
        flip_horizontal = FALSE,
        shift_vertical_step = FALSE,
        shift_horizontal_step = FALSE,
        mm = FALSE,
        px2mm = 0.12028,
        offsets,
        verbose = NULL
) {
    # NSE params
    f <- x <- y <- NULL

    if (missing(path)) {
        stop(wrap_txt(
            "No path to polys subdirectory provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, "loading segmentation masks...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    mask_params <- list(
        # static params
        mask_method = "multiple",
        # if removal is TRUE, a real cell segmentation gets removed.
        # There is no background poly for nanostring masks
        remove_background_polygon = FALSE,
        fill_holes = TRUE,
        calc_centroids = TRUE,
        remove_unvalid_polygons = TRUE,
        # input params
        name = name,
        flip_vertical = flip_vertical,
        flip_horizontal = flip_horizontal,
        shift_vertical_step = shift_vertical_step,
        shift_horizontal_step = shift_horizontal_step,
        verbose = FALSE
    )

    fovs <- fovs %null% .cosmx_imgname_fovparser(path) # ALL if NULL
    progressr::with_progress({
        p <- progressr::progressor(along = fovs)

        gpolys <- lapply(fovs, function(f) {
            segfile <- Sys.glob(paths = sprintf("%s/*F%03d*", path, f))
            # naming format: c_SLIDENUMBER_FOVNUMBER_CELLID
            mask_params$ID_fmt = paste0(
                sprintf("c_%d_%d_", slide, f), "%d"
            )

            gpoly <- do.call(
                createGiottoPolygonsFromMask,
                args = c(list(maskfile = segfile), mask_params)
            )

            xshift <- offsets[fov == f, x]
            yshift <- offsets[fov == f, y]

            # if micron scale
            if (mm) {
                gpoly <- rescale(gpoly, fx = px2mm, fy = px2mm, x0 = 0, y0 = 0)
                xshift <- xshift * px2mm
                yshift <- yshift * px2mm
            }

            gpoly <- spatShift(x = gpoly, dx = xshift, dy = yshift)
            p(message = sprintf("F%03d", f))
            return(gpoly)
        })
    })

    if (length(gpolys) > 1L) {
        gpolys <- do.call(rbind, args = gpolys)
    }

    # never return lists. Only the single merged gpoly
    return(gpolys)
}

.cosmx_cellmeta <- function(
        path,
        slide = 1,
        fovs = NULL,
        dropcols = c(
            "CenterX_local_px",
            "CenterY_local_px",
            "CenterX_global_px",
            "CenterY_global_px"
        ),
        cores = determine_cores(),
        verbose = NULL
    ) {

    if (missing(path)) {
        stop(wrap_txt(
            "No path to metadata file provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, "loading cell metadata...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    verbose <- verbose %null% TRUE

    meta_dt <- data.table::fread(input = path, nThread = cores)

    # subset to needed fovs
    if (!is.null(fovs)) {
        fovs <- as.integer(fovs)
        meta_dt <- meta_dt[fov %in% fovs,]
    }

    dropcols <- dropcols[dropcols %in% meta_dt]
    meta_dt[, `:=`(dropcols, NULL)] # remove dropcols

    # create cell ID as `c_SLIDENUMBER_FOVNUMBER_CELLID`
    if ("cell" %in% colnames(meta_dt)) {
        # assume already formatted (current datasets Mar-27-2024)
        meta_dt[, c("fov", "cell_ID") := NULL]
        data.table::setnames(meta_dt, old = "cell", "cell_ID")
    } else {
        # older datasets
        meta_dt[, cell_ID := sprintf("c_%d_%d_%d", slide, fov, cell_ID)]
        # remove fov
        meta_dt[, fov := NULL]
    }


    # TODO figure out what to do about protein expression here.
    cx <- createCellMetaObj(
        metadata = meta_dt,
        spat_unit = "cell",
        feat_type = "rna",
        provenance = "cell",
        verbose = verbose
    )
    return(cx)
}

.cosmx_expression <- function(
        path,
        slide = 1,
        fovs = NULL,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        cores = determine_cores(),
        verbose = NULL
    ) {

    if (missing(path)) {
        stop(wrap_txt(
            "No path to exprMat file provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, "loading expression matrix...")
    vmsg(.v = verbose, .is_debug = TRUE, path)

    expr_dt <- data.table::fread(input = path, nThread = cores)

    # subset to needed fovs
    if (!is.null(fovs)) {
        fovs <- as.integer(fovs)
        expr_dt <- expr_dt[fov %in% fovs,]
    }

    # remove background values (cell 0)
    expr_dt <- expr_dt[cell_ID != 0L,]

    # create cell ID as `c_SLIDENUMBER_FOVNUMBER_CELLID`
    expr_dt[, cell_ID := sprintf("c_%d_%d_%d", slide, fov, cell_ID)]
    # remove fov
    expr_dt[, fov := NULL]

    # convert to Matrix
    expr_mat <- dt_to_matrix(expr_dt)
    expr_mat <- t_flex(expr_mat)
    feat_ids <- rownames(expr_mat)

    # split expression for rna / negprb if any split keywords provided.
    # Output of this chunk should always be a named list of 1 or more matrices
    if (length(split_keyword) > 0) {
        expr_list <- list()
        for (key_i in seq_along(split_keyword)) {
            bool <- grepl(pattern = split_keyword[[key_i]], x = feat_ids)
            # subset and store split matrix
            sub_mat <- expr_mat[bool,]
            expr_list[[feat_type[[key_i + 1L]]]] <- sub_mat
            # remaining matrix
            expr_mat <- expr_mat[!bool,]
        }
        expr_list[[feat_type[[1L]]]] <- expr_mat
    } else {
        expr_list <- list(expr_mat)
        names(expr_list) <- feat_type[[1L]]
    }

    expr_list <- lapply(seq_along(expr_list), function(expr_i) {
        createExprObj(expression_data = expr_list[[expr_i]],
                      spat_unit = "cell",
                      feat_type = names(expr_list)[[expr_i]],
                      name = "raw",
                      provenance = "cell")
    })

    return(expr_list)
}

.cosmx_image <- function(
        path,
        fovs = NULL,
        img_type = "composite",
        img_name_fmt = paste(img_type, "_fov%03d"),
        negative_y = FALSE,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        mm = FALSE,
        px2mm = 0.12028,
        offsets,
        verbose = NULL
    ) {

    if (missing(path)) {
        stop(wrap_txt(
            "No path to image subdirectory to load provided or auto-detected"
        ), call. = FALSE)
    }

    GiottoUtils::vmsg(.v = verbose, sprintf("loading %s images...", img_type))
    vmsg(.v = verbose, .is_debug = TRUE, path)

    fovs <- fovs %null% .cosmx_imgname_fovparser(path) # ALL if NULL
    verbose <- verbose %null% TRUE

    progressr::with_progress({
        p <- progressr::progressor(along = fovs)

        gimg_list <- lapply(fovs, function(f) {
            imgfile <- Sys.glob(paths = sprintf("%s/*F%03d*", path, f))
            img_name <- sprintf(img_name_fmt, f)

            gimg <- createGiottoLargeImage(
                raster_object = imgfile,
                name = img_name,
                negative_y = negative_y,
                flip_vertical = flip_vertical,
                flip_horizontal = flip_horizontal,
                verbose = verbose
            )

            xshift <- offsets[fov == f, x]
            yshift <- offsets[fov == f, y]

            if (mm) {
                gimg <- rescale(gimg, fx = px2mm, fy = px2mm, x0 = 0, y0 = 0)
                xshift <- xshift * px2mm
                yshift <- yshift * px2mm
            }

            gimg <- spatShift(x = gimg, dx = xshift, dy = yshift)
            p(message = sprintf("F%03d", f))
            return(gimg)
        })
    })


    return(gimg_list)
}



#' @title Load CosMx folder subcellular info
#' @name .load_cosmx_folder_subcellular
#' @description loads in the feature detections information. Note that the mask
#' images are still required for a working subcellular object, and those are loaded
#' in \code{\link{.createGiottoCosMxObject_subcellular}}
#' @inheritParams createGiottoCosMxObject
#' @keywords internal
.load_cosmx_folder_subcellular = function(dir_items,
                                         FOVs = NULL,
                                         cores,
                                         verbose = TRUE) {

  vmsg(.v = verbose, 'Loading subcellular information...')

  # subcellular checks
  if(!file.exists(dir_items$`transcript locations file`))
    stop(wrap_txt('No transcript locations file (.csv) detected'))
  if(!file.exists(dir_items$`fov positions file`))
    stop(wrap_txt('No fov positions file (.csv) detected'))

  # FOVs to load
  vmsg(.v = verbose, 'Loading FOV offsets...')
  fov_offset_file = fread(input = dir_items$`fov positions file`, nThread = cores)
  if(is.null(FOVs)) FOVs = fov_offset_file$fov # default to ALL FOVs
  FOV_ID = as.list(sprintf('%03d', FOVs))

  #TODO Load only relevant portions of file?

  vmsg(.v = verbose, 'Loading transcript level info...')
  tx_coord_all = fread(input = dir_items$`transcript locations file`, nThread = cores)
  vmsg(.v = verbose, 'Subcellular load done')

  data_list = list(
    'FOV_ID' = FOV_ID,
    'fov_offset_file' = fov_offset_file,
    'tx_coord_all' = tx_coord_all
  )

  return(data_list)

}



#' @title Load CosMx folder aggregate info
#' @name .load_cosmx_folder_aggregate
#' @inheritParams createGiottoCosMxObject
#' @keywords internal
.load_cosmx_folder_aggregate = function(dir_items,
                                       cores,
                                       verbose = TRUE) {

  # data.table vars
  fov = cell_ID = fov_cell_ID = CenterX_global_px = CenterY_global_px = CenterX_local_px =
    CenterY_local_px = x_shift = y_shift = NULL

  # load aggregate information
  vmsg(.v = verbose, 'Loading provided aggregated information...')

  # aggregate checks
  if(!file.exists(dir_items$`expression matrix file`)) stop(wrap_txt('No expression matrix file (.csv) detected'))
  if(!file.exists(dir_items$`metadata file`)) stop(wrap_txt('No metadata file (.csv) detected. Needed for cell spatial locations.'))

  # read in aggregate data
  expr_mat = fread(input = dir_items$`expression matrix file`, nThread = cores)
  metadata = fread(input = dir_items$`metadata file`, nThread = cores)

  # setorder expression and spatlocs
  data.table::setorder(metadata, fov, cell_ID)
  data.table::setorder(expr_mat, fov, cell_ID)


  # generate unique cell IDs
  expr_mat[, cell_ID := paste0('fov', sprintf('%03d', fov), '-', 'cell_', cell_ID)]
  # expr_mat$cell_ID = paste0('fov', sprintf('%03d', expr_mat$fov), '-', 'cell_', expr_mat$cell_ID)
  expr_mat = expr_mat[, fov := NULL]

  metadata[, fov_cell_ID := cell_ID]
  metadata[, cell_ID := paste0('fov', sprintf('%03d', fov), '-', 'cell_', cell_ID)]
  # metadata$cell_ID = paste0('fov', sprintf('%03d', metadata$fov), '-', 'cell_', metadata$cell_ID)
  # reorder
  data.table::setcolorder(x = metadata, c('cell_ID','fov','fov_cell_ID'))


  # extract spatial locations
  spatlocs = metadata[,.(CenterX_global_px, CenterY_global_px, cell_ID)]
  spatlocs_fov = metadata[,.(CenterX_local_px, CenterY_local_px, cell_ID)]
  # regenerate FOV shifts
  metadata[, x_shift := CenterX_global_px - CenterX_local_px]
  metadata[, y_shift := CenterY_global_px - CenterY_local_px]
  fov_shifts = metadata[, .(mean(x_shift), mean(y_shift)), fov]
  colnames(fov_shifts) = c('fov', 'x_shift', 'y_shift')


  # rename spatloc column names
  spatloc_oldnames = c('CenterX_global_px', 'CenterY_global_px', 'cell_ID')
  spatloc_oldnames_fov = c('CenterX_local_px', 'CenterY_local_px', 'cell_ID')
  spatloc_newnames = c('sdimx', 'sdimy', 'cell_ID')
  data.table::setnames(spatlocs, old = spatloc_oldnames, new = spatloc_newnames)
  data.table::setnames(spatlocs_fov, old = spatloc_oldnames_fov, new = spatloc_newnames)

  # cleanup metadata and spatlocs
  metadata = metadata[,c('CenterX_global_px', 'CenterY_global_px', 'CenterX_local_px', 'CenterY_local_px') := NULL]
  # find unique cell_IDs present in both expression and metadata
  giotto_cell_ID = unique(intersect(expr_mat$cell_ID, metadata$cell_ID))

  # subset to only unique cell_IDs
  expr_mat = expr_mat[cell_ID %in% giotto_cell_ID,]
  metadata = metadata[cell_ID %in% giotto_cell_ID,]


  # convert protein metadata to expr mat
  # take all mean intensity protein information except for MembraneStain and DAPI
  protein_meta_cols = colnames(metadata)
  protein_meta_cols = protein_meta_cols[grepl(pattern = 'Mean.*', x = protein_meta_cols)]
  protein_meta_cols = protein_meta_cols[!protein_meta_cols %in% c('Mean.MembraneStain', 'Mean.DAPI')]
  protein_meta_cols = c('cell_ID', protein_meta_cols)

  prot_expr = metadata[, protein_meta_cols, with = FALSE]
  prot_cell_ID = metadata[, cell_ID]
  protM = Matrix::Matrix(as.matrix(prot_expr[,-1]), dimnames = list(prot_expr[[1]], colnames(prot_expr[,-1])), sparse = FALSE)
  protM = t_flex(protM)

  # convert expression to sparse matrix
  spM = Matrix::Matrix(as.matrix(expr_mat[,-1]), dimnames = list(expr_mat[[1]], colnames(expr_mat[,-1])), sparse = TRUE)
  spM = t_flex(spM)

  ## Ready for downstream aggregate gobject creation or appending into existing subcellular Giotto object ##

  data_list = list(
    'spatlocs' = spatlocs,
    'spatlocs_fov' = spatlocs_fov,
    'metadata' = metadata,
    'protM' = protM,
    'spM' = spM,
    'fov_shifts' = fov_shifts
  )

  return(data_list)

}







## Xenium ####

#' @title Load xenium data from folder
#' @name load_xenium_folder
#' @param path_list list of full filepaths from .read_xenium_folder
#' @inheritParams createGiottoXeniumObject
#' @return list of loaded in xenium data
NULL

#' @rdname load_xenium_folder
#' @keywords internal
.load_xenium_folder = function(path_list,
                              load_format = 'csv',
                              data_to_use = 'subcellular',
                              h5_expression = 'FALSE',
                              h5_gene_ids = 'symbols',
                              gene_column_index = 1,
                              cores,
                              verbose = TRUE) {

  data_list = switch(
    load_format,
    "csv" = .load_xenium_folder_csv(
      path_list = path_list,
      data_to_use = data_to_use,
      h5_expression = h5_expression,
      h5_gene_ids = h5_gene_ids,
      gene_column_index = gene_column_index,
      cores = cores,
      verbose = verbose
    ),
    "parquet" = .load_xenium_folder_parquet(
      path_list = path_list,
      data_to_use = data_to_use,
      h5_expression = h5_expression,
      h5_gene_ids = h5_gene_ids,
      gene_column_index = gene_column_index,
      cores = cores,
      verbose = verbose
    ),
    "zarr" = stop("load_format zarr:\n Not yet implemented", call. = FALSE)
  )

  return(data_list)
}


#' @describeIn load_xenium_folder Load from csv files
#' @keywords internal
.load_xenium_folder_csv = function(path_list,
                                  cores,
                                  data_to_use = 'subcellular',
                                  h5_expression = FALSE,
                                  h5_gene_ids = 'symbols',
                                  gene_column_index = 1,
                                  verbose = TRUE) {

  # initialize return vars
  feat_meta = tx_dt = bound_dt_list = cell_meta = agg_expr = NULL

  vmsg("Loading feature metadata...", .v = verbose)
  # updated for pipeline v1.6 json format
  fdata_path <- path_list$panel_meta_path[[1]]
  fdata_ext <- GiottoUtils::file_extension(fdata_path)
  if ("json" %in% fdata_ext) {
    feat_meta <- .load_xenium_panel_json(path = fdata_path, gene_ids = h5_gene_ids)
  } else {
    feat_meta <- data.table::fread(fdata_path, nThread = cores)
    colnames(feat_meta)[[1]] <- 'feat_ID'
  }

  # **** subcellular info ****
  if(data_to_use == 'subcellular') {
    # append missing QC probe info to feat_meta
    if(isTRUE(h5_expression)) {
      h5 = hdf5r::H5File$new(path_list$agg_expr_path)
      tryCatch({
        root = names(h5)
        feature_id = h5[[paste0(root, "/features/id")]][]
        feature_info = h5[[paste0(root,"/features/feature_type")]][]
        feature_names = h5[[paste0(root, "/features/name")]][]
        features_dt = data.table::data.table(
          'id' = feature_id,
          'name' = feature_names,
          'feature_type' = feature_info
        )
      }, finally = {
        h5$close_all()
      })
    } else {
      features_dt <- data.table::fread(
        paste0(path_list$agg_expr_path, '/features.tsv.gz'),
        header = FALSE
      )
    }
    colnames(features_dt) = c('id', 'feat_ID', 'feat_class')
    feat_meta = merge(features_dt[,c(2,3)], feat_meta, all.x = TRUE, by = 'feat_ID')

    GiottoUtils::vmsg("Loading transcript level info...", .v = verbose)
    tx_dt = data.table::fread(path_list$tx_path[[1]], nThread = cores)
    data.table::setnames(x = tx_dt,
                         old = c('feature_name', 'x_location', 'y_location'),
                         new = c('feat_ID', 'x', 'y'))

    GiottoUtils::vmsg("Loading boundary info...", .v = verbose)
    bound_dt_list = lapply(
      path_list$bound_paths,
      function(x) data.table::fread(x[[1]], nThread = cores)
    )
  }

  # **** aggregate info ****
  GiottoUtils::vmsg("loading cell metadata...", .v = verbose)
  cell_meta = data.table::fread(path_list$cell_meta_path[[1]], nThread = cores)

  if(data_to_use == 'aggregate') {
    GiottoUtils::vmsg("Loading aggregated expression...", .v = verbose)
    if (isTRUE(h5_expression)) agg_expr = get10Xmatrix_h5(
      path_to_data = path_list$agg_expr_path,
      gene_ids = h5_gene_ids,
      remove_zero_rows = TRUE,
      split_by_type = TRUE
    )
    else agg_expr = get10Xmatrix(
      path_to_data = path_list$agg_expr_path,
      gene_column_index = gene_column_index,
      remove_zero_rows = TRUE,
      split_by_type = TRUE
    )
  }

  data_list = list(
    'feat_meta' = feat_meta,
    'tx_dt' = tx_dt,
    'bound_dt_list' = bound_dt_list,
    'cell_meta' = cell_meta,
    'agg_expr' = agg_expr
  )

  return(data_list)

}




#' @describeIn load_xenium_folder Load from parquet files
#' @keywords internal
.load_xenium_folder_parquet = function(path_list,
                                      cores,
                                      data_to_use = 'subcellular',
                                      h5_expression = FALSE,
                                      h5_gene_ids = 'symbols',
                                      gene_column_index = 1,
                                      verbose = TRUE) {

  # initialize return vars
  feat_meta = tx_dt = bound_dt_list = cell_meta = agg_expr = NULL
  # dplyr variable
  cell_id = NULL

  vmsg("Loading feature metadata...", .v = verbose)
  # updated for pipeline v1.6 json format
  fdata_path <- path_list$panel_meta_path[[1]]
  fdata_ext <- GiottoUtils::file_extension(fdata_path)
  if ("json" %in% fdata_ext) {
    feat_meta <- .load_xenium_panel_json(path = fdata_path, gene_ids = h5_gene_ids)
  } else {
    feat_meta <- data.table::fread(fdata_path, nThread = cores)
    colnames(feat_meta)[[1]] <- 'feat_ID'
  }

  # **** subcellular info ****
  if(data_to_use == 'subcellular') {

    # define for data.table
    transcript_id = feature_name = NULL

    # append missing QC probe info to feat_meta
    if(isTRUE(h5_expression)) {
      h5 = hdf5r::H5File$new(path_list$agg_expr_path)
      tryCatch({
        root = names(h5)
        feature_id = h5[[paste0(root, "/features/id")]][]
        feature_info = h5[[paste0(root,"/features/feature_type")]][]
        feature_names = h5[[paste0(root, "/features/name")]][]
        features_dt = data.table::data.table(
          'id' = feature_id,
          'name' = feature_names,
          'feature_type' = feature_info
        )
      }, finally = {
        h5$close_all()
      })
    } else {
      features_dt = arrow::read_tsv_arrow(paste0(path_list$agg_expr_path, '/features.tsv.gz'),
                                          col_names = FALSE) %>%
        data.table::setDT()
    }
    colnames(features_dt) = c('id', 'feat_ID', 'feat_class')
    feat_meta = merge(features_dt[,c(2,3)], feat_meta, all.x = TRUE, by = 'feat_ID')

    vmsg('Loading transcript level info...', .v = verbose)
    tx_dt = arrow::read_parquet(file = path_list$tx_path[[1]],
                                as_data_frame = FALSE) %>%
      dplyr::mutate(transcript_id = cast(transcript_id, arrow::string())) %>%
      dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
      dplyr::mutate(feature_name = cast(feature_name, arrow::string())) %>%
      as.data.frame() %>%
      data.table::setDT()
    data.table::setnames(x = tx_dt,
                         old = c('feature_name', 'x_location', 'y_location'),
                         new = c('feat_ID', 'x', 'y'))
    vmsg('Loading boundary info...', .v = verbose)
    bound_dt_list = lapply(path_list$bound_paths, function(x) {
      arrow::read_parquet(file = x[[1]], as_data_frame = FALSE) %>%
        dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
        as.data.frame() %>%
        data.table::setDT()})
  }
  # **** aggregate info ****
  if(data_to_use == 'aggregate') {
    vmsg('Loading cell metadata...', .v = verbose)
    cell_meta = arrow::read_parquet(file = path_list$cell_meta_path[[1]],
                                    as_data_frame = FALSE) %>%
      dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
      as.data.frame() %>%
      data.table::setDT()

    # NOTE: no parquet for agg_expr.
    vmsg('Loading aggregated expression...', .v = verbose)
    if(isTRUE(h5_expression)) agg_expr = get10Xmatrix_h5(
      path_to_data = path_list$agg_expr_path,
      gene_ids = h5_gene_ids,
      remove_zero_rows = TRUE,
      split_by_type = TRUE
    )
    else agg_expr = get10Xmatrix(
      path_to_data = path_list$agg_expr_path,
      gene_column_index = gene_column_index,
      remove_zero_rows = TRUE,
      split_by_type = TRUE
    )
  }

  data_list = list(
    'feat_meta' = feat_meta,
    'tx_dt' = tx_dt,
    'bound_dt_list' = bound_dt_list,
    'cell_meta' = cell_meta,
    'agg_expr' = agg_expr
  )

  return(data_list)

}



.load_xenium_panel_json <- function(path, gene_ids = "symbols") {
  gene_ids <- match.arg(gene_ids, c("symbols", "ensembl"))

  # tested on v1.6
  j <- jsonlite::fromJSON(path)
  # j$metadata # dataset meta
  # j$payload # main content
  # j$payload$chemistry # panel chemistry used
  # j$payload$customer # panel customer
  # j$payload$designer # panel designer
  # j$payload$spec_version # versioning
  # j$payload$panel # dataset panel stats

  panel_info <- j$payload$targets$type %>%
    data.table::as.data.table()

  switch(
    gene_ids,
    "symbols" = data.table::setnames(
      panel_info,
      old = c("data.id", "data.name", "descriptor"),
      new = c("ensembl", "feat_ID", "type")
    ),
    "ensembl" = data.table::setnames(
      panel_info,
      old = c("data.id", "data.name", "descriptor"),
      new = c("feat_ID", "symbol", "type")
    )
  )
  return(panel_info)
}



## ArchR ####

#' Create an ArchR project and run LSI dimension reduction
#'
#' @param fragmentsPath A character vector containing the paths to the input
#' files to use to generate the ArrowFiles.
#' These files can be in one of the following formats: (i) scATAC tabix files,
#' (ii) fragment files, or (iii) bam files.
#' @param genome A string indicating the default genome to be used for all ArchR
#' functions. Currently supported values include "hg19","hg38","mm9", and "mm10".
#' This value is stored as a global environment variable, not part of the ArchRProject.
#' This can be overwritten on a per-function basis using the given function's
#' geneAnnotationand genomeAnnotation parameter. For something other than one of
#' the currently supported, see createGeneAnnnotation() and createGenomeAnnnotation()
#' @param createArrowFiles_params list of parameters passed to `ArchR::createArrowFiles`
#' @param ArchRProject_params list of parameters passed to `ArchR::ArchRProject`
#' @param addIterativeLSI_params list of parameters passed to `ArchR::addIterativeLSI`
#' @param threads number of threads to use. Default = `ArchR::getArchRThreads()`
#' @param force Default = FALSE
#' @param verbose Default = TRUE
#'
#' @return An ArchR project with GeneScoreMatrix, TileMatrix, and TileMatrix-based LSI
#' @export
#'
createArchRProj <- function(fragmentsPath,
                            genome = c('hg19', 'hg38', 'mm9', 'mm10'),
                            createArrowFiles_params = list(sampleNames = 'sample1',
                                                           minTSS = 0,
                                                           minFrags = 0,
                                                           maxFrags = 1e+07,
                                                           minFragSize = 10,
                                                           maxFragSize = 2000,
                                                           offsetPlus = 0,
                                                           offsetMinus = 0,
                                                           TileMatParams = list(tileSize = 5000)),
                            ArchRProject_params = list(outputDirectory = getwd(),
                                                       copyArrows = FALSE),
                            addIterativeLSI_params = list(),
                            threads = ArchR::getArchRThreads(),
                            force = FALSE,
                            verbose = TRUE) {

  if(!requireNamespace('ArchR')) {
    wrap_msg('ArchR is needed. Install the package using remotes::install_github("GreenleafLab/ArchR")')
  }

  ## Add reference genome
  wrap_msg('Loading reference genome')
  ArchR::addArchRGenome(genome)

  # Creating Arrow Files
  wrap_msg('Creating Arrow files')
  ArrowFiles <- do.call(ArchR::createArrowFiles,
                        c(inputFiles = fragmentsPath,
                          verbose = verbose,
                          force = force,
                          createArrowFiles_params)
  )

  # Creating an ArchRProject
  wrap_msg('Creating ArchRProject')
  proj <- do.call(ArchR::ArchRProject,
                  c(list(ArrowFiles = ArrowFiles),
                    threads = threads,
                    ArchRProject_params)
  )

  # Data normalization and dimensionality reduction
  wrap_msg('Running dimension reduction')
  proj <- do.call(ArchR::addIterativeLSI,
                  c(ArchRProj = proj,
                    verbose = verbose,
                    name = "IterativeLSI",
                    threads = threads,
                    force = force,
                    addIterativeLSI_params)
  )
}

#' Create a Giotto object from an ArchR project
#'
#' @param archRproj ArchR project
#' @param expression expression information
#' @param expression_feat Giotto object available features (e.g. atac, rna, ...)
#' @param spatial_locs data.table or data.frame with coordinates for cell centroids
#' @param sampleNames A character vector containing the ArchR project sample name
#' @param ... additional arguments passed to `createGiottoObject`
#'
#' @return A Giotto object with at least an atac or epigenetic modality
#'
#' @export
#'
createGiottoObjectfromArchR <- function(archRproj,
                                        expression = NULL,
                                        expression_feat = 'atac',
                                        spatial_locs = NULL,
                                        sampleNames = 'sample1',
                                        ...) {
  # extract GeneScoreMatrix
  GeneScoreMatrix_summarizedExperiment = ArchR::getMatrixFromProject(archRproj)
  GeneScoreMatrix = slot(slot(GeneScoreMatrix_summarizedExperiment, 'assays'), 'data')[['GeneScoreMatrix']]

  ## get cell names
  cell_names = colnames(GeneScoreMatrix)
  cell_names = gsub(paste0(sampleNames,'#'),'',cell_names)
  cell_names = gsub('-1','',cell_names)

  ## get gene names
  gene_names = slot(GeneScoreMatrix_summarizedExperiment,'elementMetadata')[['name']]

  ## replace colnames with cell names
  colnames(GeneScoreMatrix) = cell_names

  ## replace rownames with gene names
  rownames(GeneScoreMatrix) = gene_names


  if(!is.null(expression)) {
    expression_matrix = data.table::fread(expression)

    expression_cell_names = colnames(expression_matrix)
    cell_names = intersect(cell_names, expression_cell_names)

    expression_matrix = Matrix::Matrix(as.matrix(expression_matrix[,-1]),
                                       dimnames = list(expression_matrix[[1]],
                                                       colnames(expression_matrix[,-1])),
                                       sparse = T)

    expression = expression_matrix[, cell_names]

    GeneScoreMatrix = GeneScoreMatrix[, cell_names]
  }


  ## filter spatial locations
  if(!is.null(spatial_locs)) {
    x = read.csv(spatial_locs)
    x = x[x$cell_ID %in% cell_names,]
    spatial_locs = x
  }

  # Creating GiottoObject
  wrap_msg('Creating GiottoObject')

  if(!is.null(expression)) {
    gobject <- createGiottoObject(expression = list(GeneScoreMatrix = GeneScoreMatrix,
                                                    raw = expression),
                                  expression_feat = expression_feat,
                                  spatial_locs = spatial_locs,
                                  ...)
  } else {
    gobject <- createGiottoObject(expression = list(GeneScoreMatrix = GeneScoreMatrix),
                                  expression_feat = expression_feat,
                                  spatial_locs = spatial_locs,
                                  ...)
  }

  # add LSI dimension reduction
  coordinates = slot(archRproj,'reducedDims')[['IterativeLSI']][['matSVD']]

  ## clean cell names
  lsi_cell_names = rownames(coordinates)
  lsi_cell_names = gsub(paste0(sampleNames,'#'),'',lsi_cell_names)
  lsi_cell_names = gsub('-1','',lsi_cell_names)

  rownames(coordinates) = lsi_cell_names

  coordinates = coordinates[cell_names,]

  dimension_reduction = Giotto::createDimObj(coordinates = coordinates,
                                             name = 'lsi',
                                             spat_unit = 'cell',
                                             feat_type = expression_feat[1],
                                             method = 'lsi')
  gobject <- setDimReduction(gobject,
                             dimension_reduction,
                             spat_unit = 'cell',
                             feat_type = expression_feat[1],
                             name = 'lsi',
                             reduction_method = 'lsi')

  return(gobject)
}



