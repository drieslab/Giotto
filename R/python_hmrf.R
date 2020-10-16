
#' @title doHMRF
#' @name doHMRF
#' @description Run HMRF
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param spatial_network_name name of spatial network to use for HMRF
#' @param spatial_genes spatial genes to use for HMRF
#' @param spatial_dimensions select spatial dimensions to use, default is all possible dimensions
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name name of HMRF run
#' @param k  number of HMRF domains
#' @param seed seed to fix random number generator (for creating initialization of HMRF) (-1 if no fixing)
#' @param betas betas to test for. three numbers: start_beta, beta_increment, num_betas e.g. c(0, 2.0, 50)
#' @param tolerance tolerance
#' @param zscore zscore
#' @param numinit number of initializations
#' @param python_path python path to use
#' @param output_folder output folder to save results
#' @param overwrite_output overwrite output folder
#' @return Creates a directory with results that can be viewed with viewHMRFresults
#' @details Description of HMRF parameters ...
#' @export
doHMRF <- function(gobject,
                   expression_values = c('normalized', 'scaled', 'custom'),
                   spatial_network_name = 'Delaunay_network',
                   spatial_genes = NULL,
                   spatial_dimensions = c('sdimx', 'sdimy', 'sdimz'),
                   dim_reduction_to_use = NULL,
                   dim_reduction_name = 'pca',
                   dimensions_to_use = 1:10,
                   seed = 100,
                   name = 'test',
                   k = 10,
                   betas = c(0, 2, 50),
                   tolerance = 1e-10,
                   zscore = c('none','rowcol', 'colrow'),
                   numinit = 100,
                   python_path = NULL,
                   output_folder = NULL,
                   overwrite_output = TRUE) {


  if(!requireNamespace('smfishHmrf', quietly = TRUE)) {
    stop("\n package ", 'smfishHmrf' ," is not yet installed \n",
         "To install: \n",
         "remotes::install_bitbucket(repo = 'qzhudfci/smfishhmrf-r', ref='master')",
         "see http://spatial.rc.fas.harvard.edu/install.html for more information",
         call. = FALSE)
  }


  # data.table set global variable
  to = from = NULL

  ## check or make paths
  # python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
    #python_path = system('which python')
  }

  ## reader.py and get_result.py paths
  reader_path = system.file("python", "reader2.py", package = 'Giotto')

  ## output folder
  # no folder path specified
  if(is.null(output_folder)) {
    output_folder = paste0(getwd(),'/','HMRF_output')
    if(!file.exists(output_folder)) {
      dir.create(path = paste0(getwd(),'/','HMRF_output'), recursive = T)
    }
  }
  # folder path specified
  else if(!is.null(output_folder)) {
    if(!file.exists(output_folder)) {
      dir.create(path = output_folder, recursive = T)
    }
  }


  ## first write necessary txt files to output folder ##
  # cell location / spatial network / expression data and selected spatial genes

  ## 1. expression values
  if(!is.null(dim_reduction_to_use)) {
    expr_values = gobject@dimension_reduction[['cells']][[dim_reduction_to_use]][[dim_reduction_name]][['coordinates']][, dimensions_to_use]
    expr_values = t_giotto(expr_values)
  } else {
    values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
    expr_values = select_expression_values(gobject = gobject, values = values)
  }
  expression_file = paste0(output_folder,'/', 'expression_matrix.txt')

  # overwrite if exists
  if(file.exists(expression_file) & overwrite_output == TRUE) {
    cat('\n expression_matrix.txt already exists at this location, will be overwritten \n')
    write.table(expr_values,
                file = expression_file,
                quote = F, col.names = NA, row.names = T)
  } else if(file.exists(expression_file) & overwrite_output == FALSE) {
    cat('\n expression_matrix.txt already exists at this location, will be used again \n')
  } else {
    write.table(expr_values,
                file = expression_file,
                quote = F, col.names = NA, row.names = T)
  }






  ## 2. spatial genes
  if(!is.null(dim_reduction_to_use)) {
    dimred_rownames = rownames(expr_values)
    spatial_genes_detected = dimred_rownames[dimensions_to_use]
    spatial_genes_detected = spatial_genes_detected[!is.na(spatial_genes_detected)]
  } else {
    if(is.null(spatial_genes)) {
      stop('\n you need to provide a vector of spatial genes (~500) \n')
    }
    spatial_genes_detected = spatial_genes[spatial_genes %in% rownames(expr_values)]
  }
  spatial_genes_file = paste0(output_folder,'/', 'spatial_genes.txt')

  # overwrite if exists
  if(file.exists(spatial_genes_file) & overwrite_output == TRUE) {
    cat('\n spatial_genes.txt already exists at this location, will be overwritten \n')
    write.table(spatial_genes_detected,
                file = spatial_genes_file,
                quote = F, col.names = F, row.names = F)
  } else if(file.exists(spatial_genes_file) & overwrite_output == FALSE) {
    cat('\n spatial_genes.txt already exists at this location, will be used again \n')
  } else {
    write.table(spatial_genes_detected,
                file = spatial_genes_file,
                quote = F, col.names = F, row.names = F)
  }




  ## 3. spatial network
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  spatial_network = spatial_network[,.(to,from)]
  spatial_network_file = paste0(output_folder,'/', 'spatial_network.txt')

  if(file.exists(spatial_network_file) & overwrite_output == TRUE) {
    cat('\n spatial_network.txt already exists at this location, will be overwritten \n')
    write.table(spatial_network,
                file = spatial_network_file,
                row.names = F, col.names = F, quote = F, sep = '\t')
  } else if(file.exists(spatial_network_file) & overwrite_output == FALSE) {
    cat('\n spatial_network.txt already exists at this location, will be used again \n')
  } else {
    write.table(spatial_network,
                file = spatial_network_file,
                row.names = F, col.names = F, quote = F, sep = '\t')
  }




  ## 4. cell location
  spatial_location = gobject@spatial_locs

  # select spatial dimensions that are available #
  spatial_dimensions = spatial_dimensions[spatial_dimensions %in% colnames(spatial_location)]
  spatial_location = spatial_location[, c(spatial_dimensions,'cell_ID'), with = F]
  spatial_location_file = paste0(output_folder,'/', 'spatial_cell_locations.txt')

  if(file.exists(spatial_location_file) & overwrite_output == TRUE) {
    cat('\n spatial_cell_locations.txt already exists at this location, will be overwritten \n')
    write.table(spatial_location,
                file = spatial_location_file,
                row.names = F, col.names = F, quote = F, sep = '\t')
  } else if(file.exists(spatial_location_file)) {
    cat('\n spatial_cell_locations.txt already exists at this location, will be used again \n')
  } else {
    write.table(spatial_location,
                file = spatial_location_file,
                row.names = F, col.names = F, quote = F, sep = '\t')
  }




  # prepare input paths
  cell_location = paste0(output_folder,'/','spatial_cell_locations.txt')
  spatial_genes = paste0(output_folder,'/','spatial_genes.txt')
  spatial_network = paste0(output_folder,'/','spatial_network.txt')
  expression_data = paste0(output_folder,'/', 'expression_matrix.txt')

  # create output subfolder for HMRF
  output_data = paste0(output_folder,'/', 'result.spatial.zscore')
  if(!file.exists(output_data)) dir.create(output_data)

  # encapsulate to avoid path problems
  # python code also needs to be updated internally
  cell_location =  paste0('"', cell_location, '"')
  spatial_genes =  paste0('"', spatial_genes, '"')
  spatial_network =  paste0('"', spatial_network, '"')
  expression_data =  paste0('"', expression_data, '"')
  output_data =  paste0('"', output_data, '"')

  # process other params
  zscore = match.arg(zscore, c('none','rowcol', 'colrow'))
  betas_param = c('-b', betas)
  betas_final = paste(betas_param, collapse = ' ')

  ## reader part ##
  reader_command = paste0(python_path, ' ', reader_path,
                          ' -l ', cell_location,
                          ' -g ', spatial_genes,
                          ' -n ', spatial_network,
                          ' -e ', expression_data,
                          ' -o ', output_data,
                          ' -a ', name,
                          ' -k ', k,
                          ' ', betas_final,
                          ' -t ', tolerance,
                          ' -z ', zscore,
                          ' -s ', seed,
                          ' -i ', numinit)

  print(reader_command)
  system(command = reader_command)


  # store parameter results in HMRF S3 object
  HMRFObj = list(name = name,
                 output_data = output_data,
                 k = k,
                 betas = betas,
                 python_path = python_path)

  class(HMRFObj) <- append(class(HMRFObj), 'HMRFoutput')


  return(HMRFObj)

}



#' @title loadHMRF
#' @name loadHMRF
#' @description load previous HMRF
#' @param name_used name of HMRF that was run
#' @param k_used  number of HMRF domains that was tested
#' @param betas_used betas that were tested
#' @param python_path_used python path that was used
#' @param output_folder_used output folder that was used
#' @return reloads a previous ran HMRF from doHRMF
#' @details Description of HMRF parameters ...
#' @export
loadHMRF = function(name_used = 'test',
                    output_folder_used,
                    k_used = 10,
                    betas_used,
                    python_path_used) {

  output_data = paste0(output_folder_used,'/', 'result.spatial.zscore')
  if(!file.exists(output_data)) {
    stop('\n doHMRF was not run in this output directory \n')
  }

  # check if it indeed exists

  HMRFObj = list(name = name_used,
                 output_data = output_data,
                 k = k_used,
                 betas = betas_used,
                 python_path = python_path_used)

  class(HMRFObj) <- append(class(HMRFObj), 'HMRFoutput')


  return(HMRFObj)

}



#' @title viewHMRFresults
#' @name viewHMRFresults
#' @description View results from doHMRF.
#' @param gobject giotto object
#' @param HMRFoutput HMRF output from doHMRF
#' @param k number of HMRF domains
#' @param betas_to_view results from different betas that you want to view
#' @param third_dim 3D data (boolean)
#' @param \dots additional paramters (see details)
#' @return spatial plots with HMRF domains
#' @seealso \code{\link{spatPlot2D}} and \code{\link{spatPlot3D}}
#' @export
viewHMRFresults <- function(gobject,
                            HMRFoutput,
                            k = NULL,
                            betas_to_view = NULL,
                            third_dim = FALSE,
                            ...) {


  if(!'HMRFoutput' %in% class(HMRFoutput)) {
    stop('\n HMRFoutput needs to be output from doHMRFextend \n')
  }

  ## reader.py and get_result.py paths
  # TODO: part of the package
  get_result_path = system.file("python", "get_result2.py", package = 'Giotto')

  # paths and name
  name = HMRFoutput$name
  output_data = HMRFoutput$output_data
  python_path = HMRFoutput$python_path

  # k-values
  if(is.null(k)) {
    stop('\n you need to select a k that was used with doHMRFextend \n')
  }
  k = HMRFoutput$k

  # betas
  betas = HMRFoutput$betas
  possible_betas = seq(betas[1], to = betas[1]+(betas[2]*(betas[3]-1)), by = betas[2])

  betas_to_view_detected = betas_to_view[betas_to_view %in% possible_betas]

  # plot betas
  for(b in betas_to_view_detected) {

    ## get results part ##
    result_command = paste0(python_path, ' ', get_result_path,
                            ' -r ', output_data,
                            ' -a ', name,
                            ' -k ', k,
                            ' -b ', b)

    print(result_command)

    output = system(command = result_command, intern = T)


    title_name = paste0('k = ', k, ' b = ',b)

    spatPlot2D(gobject = gobject, cell_color = output, show_plot = T, title = title_name, ...)

    if(third_dim == TRUE) {
      spatPlot3D(gobject = gobject, cell_color = output, show_plot = T, ...)
    }
    #visPlot(gobject = gobject, sdimz = third_dim, cell_color = output, show_plot = T, title = title_name,...)
  }
}



#' @title writeHMRFresults
#' @name writeHMRFresults
#' @description write results from doHMRF to a data.table.
#' @param gobject giotto object
#' @param HMRFoutput HMRF output from doHMRF
#' @param k k to write results for
#' @param betas_to_view results from different betas that you want to view
#' @param print_command see the python command
#' @return data.table with HMRF results for each b and the selected k
#' @export
writeHMRFresults <- function(gobject,
                             HMRFoutput,
                             k = NULL,
                             betas_to_view = NULL,
                             print_command = F) {


  if(!'HMRFoutput' %in% class(HMRFoutput)) {
    stop('\n HMRFoutput needs to be output from doHMRFextend \n')
  }

  ## reader.py and get_result.py paths
  # TODO: part of the package
  get_result_path = system.file("python", "get_result2.py", package = 'Giotto')

  # paths and name
  name = HMRFoutput$name
  output_data = HMRFoutput$output_data
  python_path = HMRFoutput$python_path

  # k-values
  if(is.null(k)) {
    stop('\n you need to select a k that was used with doHMRFextend \n')
  }
  k = HMRFoutput$k

  # betas
  betas = HMRFoutput$betas
  possible_betas = seq(betas[1], to = betas[1]+(betas[2]*(betas[3]-1)), by = betas[2])

  betas_to_view_detected = betas_to_view[betas_to_view %in% possible_betas]

  result_list = list()

  # plot betas
  for(i in 1:length(betas_to_view_detected)) {

    b = betas_to_view_detected[i]

    ## get results part ##
    result_command = paste0(python_path, ' ', get_result_path,
                            ' -r ', output_data,
                            ' -a ', name,
                            ' -k ', k,
                            ' -b ', b)

    if(print_command == TRUE) {
      print(result_command)
    }

    output = system(command = result_command, intern = T)
    title_name = paste0('k.', k, '.b.',b)
    result_list[[title_name]] = output

  }

  result_DT = data.table::as.data.table(do.call('cbind', result_list))
  result_DT = cbind(data.table::data.table('cell_ID' = gobject@cell_ID), result_DT)
  return(result_DT)

}




#' @title addHMRF
#' @name addHMRF
#' @description Add selected results from doHMRF to the giotto object
#' @param gobject giotto object
#' @param HMRFoutput HMRF output from doHMRF()
#' @param k number of domains
#' @param betas_to_add results from different betas that you want to add
#' @param hmrf_name specify a custom name
#' @return giotto object
#' @export
addHMRF <- function(gobject,
                    HMRFoutput,
                    k = NULL,
                    betas_to_add = NULL,
                    hmrf_name = NULL) {


  if(!'HMRFoutput' %in% class(HMRFoutput)) {
    stop('\n HMRFoutput needs to be output from doHMRFextend \n')
  }

  ## reader.py and get_result.py paths
  # TODO: part of the package
  get_result_path = system.file("python", "get_result2.py", package = 'Giotto')

  # paths and name
  name = HMRFoutput$name
  output_data = HMRFoutput$output_data
  python_path = HMRFoutput$python_path

  # k-values
  if(is.null(k)) {
    stop('\n you need to select a k that was used with doHMRFextend \n')
  }
  k = HMRFoutput$k

  # betas
  betas = HMRFoutput$betas
  possible_betas = seq(betas[1], to = betas[1]+(betas[2]*(betas[3]-1)), by = betas[2])

  betas_to_add_detected = betas_to_add[betas_to_add %in% possible_betas]


  # get cell metadata for object
  cell_metadata = pDataDT(gobject)


  # plot betas
  for(b in betas_to_add_detected) {

    ## get results part ##
    result_command = paste0(python_path, ' ', get_result_path,
                            ' -r ', output_data,
                            ' -a ', name,
                            ' -k ', k,
                            ' -b ', b)
    print(result_command)
    output = system(command = result_command, intern = T)

    # create unique name
    annot_DT = data.table::data.table(temp_name = output)

    if(!is.null(hmrf_name)) {
      annot_name = paste0(hmrf_name,'_k', k, '_b.',b)
      setnames(annot_DT, old = 'temp_name', new = annot_name)
    } else {
      annot_name = paste0('hmrf_k.', k, '_b.',b)
      data.table::setnames(annot_DT, old = 'temp_name', new = annot_name)
    }


    gobject = addCellMetadata(gobject = gobject, column_cell_ID = 'cell_ID',
                              new_metadata = annot_DT,
                              by_column = F)


  }

  return(gobject)

}





#' @title viewHMRFresults2D
#' @name viewHMRFresults2D
#' @description View results from doHMRF.
#' @param gobject giotto object
#' @param HMRFoutput HMRF output from doHMRF
#' @param k number of HMRF domains
#' @param betas_to_view results from different betas that you want to view
#' @param \dots additional parameters to spatPlot2D()
#' @return spatial plots with HMRF domains
#' @seealso \code{\link{spatPlot2D}}
#' @export
viewHMRFresults2D <- function(gobject,
                            HMRFoutput,
                            k = NULL,
                            betas_to_view = NULL,
                            ...) {


  if(!'HMRFoutput' %in% class(HMRFoutput)) {
    stop('\n HMRFoutput needs to be output from doHMRFextend \n')
  }

  ## reader.py and get_result.py paths
  # TODO: part of the package
  get_result_path = system.file("python", "get_result2.py", package = 'Giotto')

  # paths and name
  name = HMRFoutput$name
  output_data = HMRFoutput$output_data
  python_path = HMRFoutput$python_path

  # k-values
  if(is.null(k)) {
    stop('\n you need to select a k that was used with doHMRFextend \n')
  }
  k = HMRFoutput$k

  # betas
  betas = HMRFoutput$betas
  possible_betas = seq(betas[1], to = betas[1]+(betas[2]*(betas[3]-1)), by = betas[2])

  betas_to_view_detected = betas_to_view[betas_to_view %in% possible_betas]

  # plot betas
  for(b in betas_to_view_detected) {

    ## get results part ##
    result_command = paste0(python_path, ' ', get_result_path,
                            ' -r ', output_data,
                            ' -a ', name,
                            ' -k ', k,
                            ' -b ', b)

    print(result_command)

    output = system(command = result_command, intern = T)


    title_name = paste0('k = ', k, ' b = ',b)

    spatPlot2D(gobject = gobject, cell_color = as.factor(output), show_plot = T, save_plot = F, title = title_name, ...)
  }
}


#' @title viewHMRFresults3D
#' @name viewHMRFresults3D
#' @description View results from doHMRF.
#' @param gobject giotto object
#' @param HMRFoutput HMRF output from doHMRF
#' @param k number of HMRF domains
#' @param betas_to_view results from different betas that you want to view
#' @param \dots additional parameters to spatPlot3D()
#' @return spatial plots with HMRF domains
#' @seealso \code{\link{spatPlot3D}}
#' @export
viewHMRFresults3D <- function(gobject,
                              HMRFoutput,
                              k = NULL,
                              betas_to_view = NULL,
                              ...) {


  if(!'HMRFoutput' %in% class(HMRFoutput)) {
    stop('\n HMRFoutput needs to be output from doHMRFextend \n')
  }

  ## reader.py and get_result.py paths
  # TODO: part of the package
  get_result_path = system.file("python", "get_result2.py", package = 'Giotto')

  # paths and name
  name = HMRFoutput$name
  output_data = HMRFoutput$output_data
  python_path = HMRFoutput$python_path

  # k-values
  if(is.null(k)) {
    stop('\n you need to select a k that was used with doHMRFextend \n')
  }
  k = HMRFoutput$k

  # betas
  betas = HMRFoutput$betas
  possible_betas = seq(betas[1], to = betas[1]+(betas[2]*(betas[3]-1)), by = betas[2])

  betas_to_view_detected = betas_to_view[betas_to_view %in% possible_betas]

  # plot betas
  for(b in betas_to_view_detected) {

    ## get results part ##
    result_command = paste0(python_path, ' ', get_result_path,
                            ' -r ', output_data,
                            ' -a ', name,
                            ' -k ', k,
                            ' -b ', b)

    print(result_command)

    output = system(command = result_command, intern = T)


    title_name = paste0('k = ', k, ' b = ',b)

    spatPlot3D(gobject = gobject, cell_color = output, show_plot = T, save_plot = F, title = title_name, ...)
  }
}
