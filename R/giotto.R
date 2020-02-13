#' @title S4 giotto Class
#' @description Framework of giotto object to store and work with spatial expression data
#' @keywords giotto, object
#' @slot raw_exprs raw expression counts
#' @slot norm_expr normalized expression counts
#' @slot norm_scaled_expr normalized and scaled expression counts
#' @slot custom_expr custom normalized counts
#' @slot spatial_locs spatial location coordinates for cells
#' @slot cell_metadata metadata for cells
#' @slot gene_metadata metadata for genes
#' @slot cell_ID unique cell IDs
#' @slot gene_ID unique gene IDs
#' @slot spatial_network spatial network in data.table/data.frame format
#' @slot spatial_grid spatial grid in data.table/data.frame format
#' @slot dimension_reduction slot to save dimension reduction coordinates
#' @slot nn_network nearest neighbor network in igraph format
#' @slot parameters slot to save parameters that have been used
#' @slot instructions slot for global function instructions
#' @slot offset_file offset file used to stitch together image fields
#' @slot OS_platform Operating System to run Giotto analysis on
#' @useDynLib Giotto
#' @importFrom Rcpp sourceCpp
#' @export
giotto <- setClass(
  "giotto",
  slots = c(
    raw_exprs = "ANY",
    norm_expr = "ANY",
    norm_scaled_expr = "ANY",
    custom_expr = "ANY",
    spatial_locs = "ANY",
    cell_metadata = "ANY",
    gene_metadata = "ANY",
    cell_ID = "ANY",
    gene_ID = "ANY",
    spatial_network = "ANY",
    spatial_grid = "ANY",
    spatial_enrichment = "ANY",
    dimension_reduction = 'ANY',
    nn_network = "ANY",
    parameters = "ANY",
    instructions = "ANY",
    offset_file = "ANY",
    OS_platform = "ANY"
  ),

  prototype = list(
    raw_exprs = NULL,
    norm_expr = NULL,
    norm_scaled_expr = NULL,
    custom_expr = NULL,
    spatial_locs = NULL,
    cell_metadata = NULL,
    gene_metadata = NULL,
    cell_ID = NULL,
    gene_ID = NULL,
    spatial_network = NULL,
    spatial_grid = NULL,
    spatial_enrichment = NULL,
    dimension_reduction = NULL,
    nn_network = NULL,
    parameters = NULL,
    instructions = NULL,
    offset_file = NULL,
    OS_platform = NULL
  ),

  validity = function(object) {

    if(any(lapply(list(object@raw_exprs), is.null) == TRUE)) {
      return('expression and spatial locations slots need to be filled in')
    }

    # check validity of instructions vector if provided by the user
    if(!is.null(object@instructions)) {

      instr_names = c("python_path", "save_dir", "plot_format", "dpi", "units", "height", "width")
      missing_names = instr_names[!instr_names %in% names(object@instructions)]
      if(length(missing_names) != 0) {
        return('\t Instruction parameters are missing for: ', missing_names, '\t')
      }
    }

    #if(any(lapply(list(object@raw_exprs, object@spatial_locs), is.null) == TRUE)) {
    #  return('expression and spatial locations slots need to be filled in')
    #}

    #if(ncol(object@raw_exprs) != nrow(object@spatial_locs)) {
    #  return('number of cells do not correspond between expression matrix and spatial locations')
    #}
    return(TRUE)
  }
)



#' @title show method for giotto class
#' @keywords giotto, object
#'
#' @export
setMethod(
  f = "show",
  signature = "giotto",
  definition = function(object) {
    cat(
      "An object of class",
      class(object),
      "\n",
      nrow(x = object@raw_exprs),
      "genes across",
      ncol(x = object@raw_exprs),
      "samples.\n \n"
    )
    cat('Steps and parameters used: \n \n')
    print(object@parameters)
    invisible(x = NULL)
  }
)


#' @title print method for giotto class
#' @description print method for giotto class.
#' Prints the chosen number of genes (rows) and cells (columns) from the raw count matrix.
#' Also print the spatial locations for the chosen number of cells.
#' @param nr_genes number of genes (rows) to print
#' @param nr_cells number of cells (columns) to print
#' @keywords giotto, object
#'
#' @export
setGeneric(name = "print.giotto",
           def = function(object, ...) {
             standardGeneric("print.giotto")
           })

setMethod(f = "print.giotto",
          signature = "giotto",
          definition = function(object, nr_genes = 5, nr_cells = 5) {
            print(object@raw_exprs[1:nr_genes, 1:nr_cells])
            cat('\n')
            print(object@spatial_locs[1:nr_cells,])
          })



#' @title createGiottoInstructions
#' @description Function to set global instructions for giotto functions
#' @param python_path path to python binary to use
#' @param show_plot print plot to console, default = TRUE
#' @param return_plot return plot as object, default = TRUE
#' @param save_plot automatically save plot, dafault = FALSE
#' @param save_dir path to directory where to save plots
#' @param dpi resolution for raster images
#' @param height height of plots
#' @param width width of  plots
#' @return named vector with giotto instructions
#' @export
#' @examples
#'     createGiottoInstructions()
createGiottoInstructions <- function(python_path =  NULL,
                                     show_plot = NULL,
                                     return_plot = NULL,
                                     save_plot = NULL,
                                     save_dir = NULL,
                                     plot_format = NULL,
                                     dpi = NULL,
                                     units = NULL,
                                     height = NULL,
                                     width = NULL) {

  # pyton path to use
  if(is.null(python_path)) {

    if(.Platform[['OS.type']] == 'unix') {
      python_path = try(system('which python', intern = T))
    } else if(.Platform[['OS.type']] == 'windows') {
      python_path = try(system('where python', intern = T))
      if(class(python_path) == "try-error") {
        cat('\n no python path found, set it manually when needed \n')
        python_path = '/set/your/python/path/manually/please/'
      }
    }
  }
  python_path = as.character(python_path)

  # print plot to console
  if(is.null(show_plot)) {
    show_plot = TRUE
  }

  # print plot to console
  if(is.null(return_plot)) {
    return_plot = TRUE
  }

  # print plot to console
  if(is.null(save_plot)) {
    save_plot = FALSE
  }

  # directory to save results to
  if(is.null(save_dir)) {
    save_dir = getwd()
  }
  save_dir = as.character(save_dir)

  # plot format
  if(is.null(plot_format)) {
    plot_format = "png"
  }
  plot_format = as.character(plot_format)

  # dpi of raster images
  if(is.null(dpi)) {
    dpi = 300
  }
  dpi = as.numeric(dpi)

  # units for height and width
  if(is.null(units)) {
    units = 'in'
  }
  units = as.character(units)

  # height of plot
  if(is.null(height)) {
    height = 9
  }
  height = as.numeric(height)

  # width of plot
  if(is.null(width)) {
    width = 9
  }
  width = as.numeric(width)

  instructions_list = list(python_path, show_plot, return_plot,
                           save_plot, save_dir, plot_format, dpi, units, height, width)
  names(instructions_list) = c('python_path', 'show_plot', 'return_plot',
                               'save_plot', 'save_dir', 'plot_format', 'dpi', 'units', 'height', 'width')
  return(instructions_list)

}

#' @title readGiottoInstrunctions
#' @description Retrieves the instruction associated with the provided parameter
#' @param giotto_instructions giotto object or result from createGiottoInstructions()
#' @param param parameter to retrieve
#' @return specific parameter
#' @export
#' @examples
#'     readGiottoInstrunctions()
readGiottoInstructions <- function(giotto_instructions,
                                   param = NULL) {

  # get instructions if provided the giotto object
  if(class(giotto_instructions) == 'giotto') {
    giotto_instructions = giotto_instructions@instructions
  }

  # stop if parameter is not found
  if(is.null(param)) {
    stop('\t readGiottoInstructions needs a parameter to work \t')
  } else if(!param %in% names(giotto_instructions)) {
    stop('\t parameter ', param, ' is not part of Giotto parameters \t')
  } else {
    specific_instruction = giotto_instructions[[param]]
  }
  return(specific_instruction)
}


#' @title showGiottoInstructions
#' @description Function to display all instructions from giotto object
#' @param gobject giotto object
#' @return named vector with giotto instructions
#' @export
#' @examples
#'     showGiottoInstructions()
showGiottoInstructions = function(gobject) {

  instrs = gobject@instructions
  return(instrs)
}


#' @title changeGiottoInstructions
#' @description Function to change one or more instructions from giotto object
#' @param gobject giotto object
#' @param params parameter(s) to change
#' @param new_values new value(s) for parameter(s)
#' @param return_gobject (boolean) return giotto object
#' @return named vector with giotto instructions
#' @export
#' @examples
#'     changeGiottoInstructions()
changeGiottoInstructions = function(gobject,
                                    params = NULL,
                                    new_values = NULL,
                                    return_gobject = TRUE) {

  instrs = gobject@instructions

  if(is.null(params) | is.null(new_values)) {
    stop('\t params and new_values can not be NULL \t')
  }

  if(length(params) != length(new_values)) {
    stop('\t length of params need to be the same as new values \t')
  }

  if(!all(params %in% names(instrs))) {
    stop('\t all params need to be part of Giotto instructions \t')
  }

  ## swap with new values
  instrs[params] = new_values

  ## make sure that classes remain consistent
  new_instrs = lapply(1:length(instrs), function(x) {

    if(names(instrs[x]) %in% c('dpi', 'height', 'width')) {
      instrs[[x]] = as.numeric(instrs[[x]])
    } else {
      instrs[[x]] = instrs[[x]]
    }

  })

  names(new_instrs) = names(instrs)



  if(return_gobject == TRUE) {
    gobject@instructions = new_instrs
    return(gobject)
  } else {
    return(new_instrs)
  }

}




#' @title replaceGiottoInstructions
#' @description Function to replace all instructions from giotto object
#' @param gobject giotto object
#' @param instructions new instructions (e.g. result from createGiottoInstructions)
#' @return named vector with giotto instructions
#' @export
#' @examples
#'     replaceGiottoInstructions()
replaceGiottoInstructions = function(gobject,
                                     instructions = NULL) {

  old_instrs_list = gobject@instructions

  # validate new instructions
  if(!all(names(instructions) %in% names(old_instrs_list)) | is.null(instructions)) {
    stop('\n You need to provide a named list for all instructions, like the outcome of createGiottoInstructions \n')
  } else {
    gobject@instructions = instructions
    return(gobject)
  }

}



#' @title create Giotto object
#' @description Function to create a giotto object
#' @param raw_exprs matrix with raw expression counts [required]
#' @param spatial_locs data.table or data.frame with coordinates for cell centroids
#' @param norm_expr normalized expression values
#' @param norm_scaled_expr scaled expression values
#' @param custom_expr custom expression values
#' @param cell_metadata cell annotation metadata
#' @param gene_metadata gene annotation metadata
#' @param spatial_network list of spatial network(s)
#' @param spatial_network_name list of spatial network name(s)
#' @param spatial_grid list of spatial grid(s)
#' @param spatial_grid_name list of spatial grid name(s)
#' @param spatial_enrichment list of spatial enrichment score(s) for each spatial region
#' @param spatial_enrichment_name list of spatial enrichment name(s)
#' @param dimension_reduction list of dimension reduction(s)
#' @param nn_network list of nearest neighbor network(s)
#' @param offset_file file used to stitch fields together (optional)
#' @param instructions list of instructions or output result from createGiottoInstructions
#' @return giotto object
#' @details
#' [\strong{Requirements}] To create a giotto object you need to provide at least a matrix with genes as
#' row names and cells as column names. To include spatial information about cells (or regions)
#' you need to provide a data.table or data.frame with coordinates for all spatial dimensions.
#' This can be 2D (x and y) or 3D (x, y, x).The row order for the cell coordinates should be
#' the same as the column order for the provided expression data.
#'
#' [\strong{Instructions}] Additionally an instruction file, generated manually or with \code{\link{createGiottoInstructions}}
#' can be provided to instructions, if not a default instruction file will be created
#' for the Giotto object.
#'
#' [\strong{Multiple fields}] In case a dataset consists of multiple fields, like seqFISH+ for example,
#' an offset file can be provided to stitch the different fields together. \code{\link{stitchFieldCoordinates}}
#' can be used to generate such an offset file.
#'
#' [\strong{Processed data}] Processed count data, such as normalized data, can be provided using
#' one of the different expression slots (norm_expr, norm_scaled_expr, custom_expr).
#'
#' [\strong{Metadata}] Cell and gene metadata can be provided using the cell and gene metadata slots.
#' This data can also be added afterwards using the \code{\link{addGeneMetadata}} or \code{\link{addCellMetadata}} functions.
#'
#' [\strong{Other information}] Additional information can be provided through the appropriate slots:
#' \itemize{
#'   \item{spatial networks}
#'   \item{spatial girds}
#'   \item{spatial enrichments}
#'   \item{dimensions reductions}
#'   \item{nearest neighbours networks}
#' }
#'
#' @keywords giotto
#' @export
#' @examples
#'     createGiottoObject(raw_exprs, spatial_locs)
createGiottoObject <- function(raw_exprs,
                               spatial_locs = NULL,
                               norm_expr = NULL,
                               norm_scaled_expr = NULL,
                               custom_expr = NULL,
                               cell_metadata = NULL,
                               gene_metadata = NULL,
                               spatial_network = NULL,
                               spatial_network_name = NULL,
                               spatial_grid = NULL,
                               spatial_grid_name = NULL,
                               spatial_enrichment = NULL,
                               spatial_enrichment_name = NULL,
                               dimension_reduction = NULL,
                               nn_network = NULL,
                               offset_file = NULL,
                               instructions = NULL) {

  # create minimum giotto
  gobject = giotto(raw_exprs = raw_exprs,
                   spatial_locs = spatial_locs,
                   norm_expr = NULL,
                   norm_scaled_expr = NULL,
                   custom_expr = NULL,
                   cell_metadata = cell_metadata,
                   gene_metadata = gene_metadata,
                   cell_ID = NULL,
                   gene_ID = NULL,
                   spatial_network = NULL,
                   spatial_grid = NULL,
                   spatial_enrichment = NULL,
                   dimension_reduction = NULL,
                   nn_network = NULL,
                   parameters = NULL,
                   offset_file = offset_file,
                   instructions = instructions,
                   OS_platform = .Platform[['OS.type']])


  # check if all optional packages are installed
  extra_packages = c("scran", "MAST", "png", "tiff", "biomaRt", "trendsceek", "multinet")
  pack_index = extra_packages %in% rownames(installed.packages())
  extra_installed_packages = extra_packages[pack_index]
  extra_not_installed_packages = extra_packages[!pack_index]

  if(any(pack_index == FALSE) == TRUE) {
    cat("Consider to install these (optional) packages to run all possible Giotto commands: ",
        extra_not_installed_packages)
    cat("\n Giotto does not automatically install all these packages as they are not absolutely required and this reduces the number of dependencies")
  }


  # check input of raw_exprs & force it as matrix
  if(!any(c('matrix','data.frame') %in% class(raw_exprs)) | 'data.table' %in% class(raw_exprs)) {
    stop("raw_exprs needs to be of class 'matrix', check class(raw_exprs)")
  }

  raw_exprs = as.matrix(raw_exprs)
  gobject@raw_exprs = raw_exprs

  if(any(duplicated(rownames(raw_exprs)))) {
    stop("row names contain duplicates, please remove or rename")
  }

  if(any(duplicated(colnames(raw_exprs)))) {
    stop("column names contain duplicates, please remove or rename")
  }

  # prepare other slots
  gobject@cell_ID = colnames(raw_exprs)
  gobject@gene_ID = rownames(raw_exprs)
  gobject@parameters = list()

  if(is.null(instructions)) {
    # create all default instructions
    gobject@instructions = createGiottoInstructions()
  }

  ## if no spatial information is given; create dummy spatial data
  if(is.null(spatial_locs)) {
    cat('\n spatial locations are not given, dummy 3D data will be created \n')
    spatial_locs = data.table::data.table(x = 1:ncol(raw_exprs),
                                          y = 1:ncol(raw_exprs),
                                          z = 1:ncol(raw_exprs))
    gobject@spatial_locs = spatial_locs
  }


  ## spatial
  if(nrow(spatial_locs) != ncol(raw_exprs)) {
    stop('\n Number of rows of spatial location must equal number of columns of expression matrix \n')
  } else {
    spatial_dimensions = c('x', 'y', 'z')
    colnames(gobject@spatial_locs) <- paste0('sdim', spatial_dimensions[1:ncol(gobject@spatial_locs)])
    gobject@spatial_locs = data.table::as.data.table(gobject@spatial_locs)
    gobject@spatial_locs[, cell_ID := colnames(raw_exprs)]
  }


  ## OPTIONAL:
  # add other normalized expression data
  if(!is.null(norm_expr)) {

    if(all(dim(norm_expr) == dim(raw_exprs)) &
       all(colnames(norm_expr) == colnames(raw_exprs)) &
       all(rownames(norm_expr) == rownames(raw_exprs))) {

      gobject@norm_expr = as.matrix(norm_expr)
    } else {
      stop('\n dimensions, row or column names are not the same between normalized and raw expression \n')
    }
  }

  # add other normalized and scaled expression data
  if(!is.null(norm_scaled_expr)) {

    if(all(dim(norm_scaled_expr) == dim(raw_exprs)) &
       all(colnames(norm_scaled_expr) == colnames(raw_exprs)) &
       all(rownames(norm_scaled_expr) == rownames(raw_exprs))) {

      gobject@norm_scaled_expr = as.matrix(norm_scaled_expr)
    } else {
      stop('\n dimensions, row or column names are not the same between normalized + scaled and raw expression \n')
    }
  }

  # add other custom normalized expression data
  if(!is.null(custom_expr)) {

    if(all(dim(custom_expr) == dim(raw_exprs)) &
       all(colnames(custom_expr) == colnames(raw_exprs)) &
       all(rownames(custom_expr) == rownames(raw_exprs))) {

      gobject@custom_expr = as.matrix(custom_expr)
    } else {
      stop('\n dimensions, row or column names are not the same between custom normalized and raw expression \n')
    }
  }

  # cell metadata
  if(is.null(cell_metadata)) {
    gobject@cell_metadata = data.table::data.table(cell_ID = colnames(raw_exprs))
  } else {
    gobject@cell_metadata = data.table::as.data.table(gobject@cell_metadata)
    gobject@cell_metadata[, cell_ID := colnames(raw_exprs)]
    # put cell_ID first
    all_colnames = colnames(gobject@cell_metadata)
    other_colnames = grep('cell_ID', all_colnames, invert = T, value = T)
    gobject@cell_metadata = gobject@cell_metadata[, c('cell_ID', other_colnames), with = FALSE]
  }

  # gene metadata
  if(is.null(gene_metadata)) {
    gobject@gene_metadata = data.table::data.table(gene_ID = rownames(raw_exprs))
  } else {
    gobject@gene_metadata = data.table::as.data.table(gobject@gene_metadata)
    gobject@gene_metadata[, gene_ID := rownames(raw_exprs)]
  }


  ### OPTIONAL:
  ## spatial network
  if(!is.null(spatial_network)) {
    if(is.null(spatial_network_name) | length(spatial_network) != length(spatial_network_name)) {
      stop('\n each spatial network must be given a unique name \n')
    } else {

      for(network_i in 1:length(spatial_network)) {

        networkname = spatial_network_name[[network_i]]
        network     = spatial_network[[network_i]]

        if(any(c('data.frame', 'data.table') %in% class(network))) {
          if(all(c('to', 'from', 'weight', 'sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') %in% colnames(network))) {
            gobject@spatial_network[[networkname]] = network
          } else {
            stop('\n network ', networkname, ' does not have all necessary column names, see details \n')
          }
        } else {
          stop('\n network ', networkname, ' is not a data.frame or data.table \n')
        }
      }
    }
  }


  ## spatial grid
  if(!is.null(spatial_grid)) {
    if(is.null(spatial_grid_name) | length(spatial_grid) != length(spatial_grid_name)) {
      stop('\n each spatial grid must be given a unique name \n')
    } else {

      for(grid_i in 1:length(spatial_grid)) {

        gridname = spatial_grid_name[[grid_i]]
        grid     = spatial_grid[[grid_i]]

        if(any(c('data.frame', 'data.table') %in% class(grid))) {
          if(all(c('x_start', 'y_start', 'x_end', 'y_end', 'gr_name') %in% colnames(grid))) {
            gobject@spatial_grid[[gridname]] = grid
          } else {
            stop('\n grid ', gridname, ' does not have all necessary column names, see details \n')
          }
        } else {
          stop('\n grid ', gridname, ' is not a data.frame or data.table \n')
        }
      }
    }
  }

  ## spatial enrichment
  if(!is.null(spatial_enrichment)) {
    if(is.null(spatial_enrichment_name) | length(spatial_enrichment) != length(spatial_enrichment_name)) {
      stop('\n each spatial enrichment data.table or data.frame must be given a unique name \n')
    } else {

      for(spat_enrich_i in 1:length(spatial_enrichment)) {

        spatenrichname = spatial_enrichment_name[[spat_enrich_i]]
        spatenrich     = spatial_enrichment[[spat_enrich_i]]

        if(nrow(spatenrich) != nrow(goject@cell_metadata)) {
          stop('\n spatial enrichment ', spatenrichname, ' does not have the same number of rows as spots/cells, see details \n')
        } else {

          gobject@spatial_enrichment[[spatenrichname]] = spatenrich

        }
      }
    }
  }


  ## dimension reduction
  if(!is.null(dimension_reduction)) {

    for(dim_i in 1:length(dimension_reduction)) {

      dim_red = dimension_reduction[[dim_i]]

      if(all(c('type', 'name', 'reduction_method', 'coordinates', 'misc') %in% names(dim_red))) {

        coord_data = dim_red[['coordinates']]

        if(all(rownames(coord_data) %in% gobject@cell_ID)) {

          type_value = dim_red[['type']] # cells or genes
          reduction_meth_value = dim_red[['reduction_method']] # e.g. umap, tsne, ...
          name_value = dim_red[['name']]  # uniq name
          misc_value = dim_red[['misc']]  # additional data

          gobject@dimension_reduction[[type_value]][[reduction_meth_value]][[name_value]] = dim_red[c('name', 'reduction_method', 'coordinates', 'misc')]
        } else {
          stop('\n rownames for coordinates are not found in gobject IDs \n')
        }

      } else {
        stop('\n each dimension reduction list must contain all required slots, see details. \n')
      }

    }

  }

  # NN network
  if(!is.null(nn_network)) {

    for(nn_i in 1:length(nn_network)) {

      nn_netw = nn_network[[nn_i]]

      if(all(c('type', 'name', 'igraph') %in% names(nn_netw))) {

        igraph_data = nn_netw[['igraph']]

        if(all(names(V(igraph_data)) %in% gobject@cell_ID)) {

          type_value = nn_netw[['type']] # sNN or kNN
          name_value = nn_netw[['name']]  # uniq name

          gobject@nn_network[[type_value]][[name_value]][['igraph']] = igraph_data
        } else {
          stop('\n igraph vertex names are not found in gobject IDs \n')
        }

      } else {
        stop('\n each nn network list must contain all required slots, see details. \n')
      }

    }

  }

  # other information
  # TODO

  return(gobject)

}






