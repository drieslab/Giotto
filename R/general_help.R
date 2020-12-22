


#' @title determine_cores
#' @description guesses how many cores to use
#' @return numeric
#' @keywords internal
determine_cores = function(cores, min_cores = 1, max_cores = 10) {

  if(is.na(cores) | !is.numeric(cores) | (is.numeric(cores) & cores <= 0)) {
    cores = parallel::detectCores()

    if(cores <= 2) {
      cores = ifelse(cores < min_cores, cores, min_cores)
    } else {
      cores = cores - 2
      cores = ifelse(cores > max_cores, max_cores, cores)
    }
    return(cores)
  } else {
    cores = cores
    return(cores)
  }
}

#' @title getDistinctColors
#' @description Returns a number of distint colors based on the RGB scale
#' @param n number of colors wanted
#' @return number of distinct colors
#' @export
getDistinctColors <- function(n) {
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))));

  if(n > length(col_vector)) {

    # get all possible colors
    all_colors = grDevices::colors()
    all_colors_no_grey = grep(x = all_colors, pattern = 'grey|gray', value = T, invert = T)
    grey_colors = grep(x = all_colors, pattern = 'grey', value = T, invert = F)
    admitted_grey_colors = grey_colors[seq(1, 110, 10)]
    broad_colors = c(all_colors_no_grey, admitted_grey_colors)

    # if too many colors stop
    if(n > length(broad_colors)) {
      warning('\n not enough unique colors in R, maximum = 444 \n')
      col_vector = sample(x = broad_colors, size = n, replace = T)
    } else {
      col_vector = sample(x = broad_colors, size = n, replace = F)
    }

  } else {

    xxx <- grDevices::col2rgb(col_vector);
    dist_mat <- as.matrix(stats::dist(t(xxx)));
    diag(dist_mat) <- 1e10;
    while (length(col_vector) > n) {
      minv <- apply(dist_mat,1,function(x)min(x));
      idx <- which(minv==min(minv))[1];
      dist_mat <- dist_mat[-idx, -idx];
      col_vector <- col_vector[-idx]
    }

  }
  return(col_vector)
}


#' @title get_os
#' @description return the type of operating system, see https://conjugateprior.org/2015/06/identifying-the-os-from-r/
#' @return character osx, linux or windows
#' @keywords internal
get_os <- function(){

  if(.Platform[['OS.type']] == 'windows') {
    os = 'windows'
  } else {

    sysinf <- Sys.info()
    if (!is.null(sysinf)){
      os = sysinf['sysname']
      if (os == 'Darwin')
        os = "osx"
    } else { ## mystery machine
      os = .Platform$OS.type
      if (grepl("^darwin", R.version$os))
        os = "osx"
      if (grepl("linux-gnu", R.version$os))
        os = "linux"
    }

  }
  return(tolower(os))
}



#' @title dt_to_matrix
#' @description converts data.table to matrix
#' @param x data.table object
#' @keywords internal
dt_to_matrix <- function(x) {
  rownames = as.character(x[[1]])
  mat = methods::as(as.matrix(x[,-1]), 'Matrix')
  rownames(mat) = rownames
  return(mat)
}


#' @title mygini_fun
#' @description calculate gini coefficient
#' @keywords internal
#' @return gini coefficient
mygini_fun <- function(x,
                       weights = rep(1,length(x))) {

  # adapted from R package GiniWegNeg
  dataset = cbind(x, weights)
  ord_x = order(x)
  dataset_ord = dataset[ord_x,]
  x       = dataset_ord[,1]
  weights = dataset_ord[,2]
  N  = sum(weights)
  xw = x*weights
  C_i = cumsum(weights)
  num_1 = sum(xw*C_i)
  num_2 = sum(xw)
  num_3 = sum(xw*weights)
  G_num = (2/N^2)*num_1-(1/N)*num_2-(1/N^2)*num_3
  t_neg = subset(xw, xw<=0)
  T_neg = sum(t_neg)
  T_pos = sum(xw)+abs(T_neg)
  n_RSV = (2*(T_pos+(abs(T_neg)))/N)
  mean_RSV = (n_RSV/2)
  G_RSV = (1/mean_RSV)*G_num
  return(G_RSV)
}


#' @title extended_gini_fun
#' @description calculate gini coefficient on a minimum length vector
#' @keywords internal
#' @return gini coefficient
extended_gini_fun <- function(x,
                              weights = rep(1, length = length(x)),
                              minimum_length = 16) {

  if(length(x) < minimum_length) {
    difference = minimum_length - length(x)
    min_value = min(x)
    x = c(x,rep(min_value, difference))
  }

  result <- mygini_fun(x = x, weights = weights)
  return(result)
}


#' @title stitchFieldCoordinates
#' @description Helper function to stitch field coordinates together to form one complete picture
#' @param location_file location dataframe with X and Y coordinates
#' @param offset_file dataframe that describes the offset for each field (see details)
#' @param cumulate_offset_x (boolean) Do the x-axis offset values need to be cumulated?
#' @param cumulate_offset_y (boolean) Do the y-axis offset values need to be cumulated?
#' @param field_col column that indicates the field within the location_file
#' @param X_coord_col column that indicates the x coordinates
#' @param Y_coord_col column that indicates the x coordinates
#' @param reverse_final_x (boolean) Do the final x coordinates need to be reversed?
#' @param reverse_final_y (boolean) Do the final y coordinates need to be reversed?
#' @return Updated location dataframe with new X ['X_final'] and Y ['Y_final'] coordinates
#' @details Stitching of fields:
#' \itemize{
#'   \item{1. have cell locations: }{at least 3 columns: field, X, Y}
#'   \item{2. create offset file: }{offset file has 3 columns: field, x_offset, y_offset}
#'   \item{3. create new cell location file by stitching original cell locations with stitchFieldCoordinates}
#'   \item{4. provide new cell location file to \code{\link{createGiottoObject}}}
#' }
#'
#' @export
stitchFieldCoordinates <- function(location_file,
                                   offset_file,
                                   cumulate_offset_x = F,
                                   cumulate_offset_y = F,
                                   field_col = 'Field of View',
                                   X_coord_col = 'X',
                                   Y_coord_col = 'Y',
                                   reverse_final_x = F,
                                   reverse_final_y = T) {


  # data.table variables
  x_offset_final = x_offset = y_offset_final = y_offset = field = NULL


  # cumulate offset values or not for offset file
  if(cumulate_offset_x == TRUE) {
    offset_file[, x_offset_final := cumsum(x_offset)]
  } else {
    offset_file[, x_offset_final := x_offset]
  }

  if(cumulate_offset_y == TRUE) {
    offset_file[, y_offset_final := cumsum(y_offset)]
  } else {
    offset_file[, y_offset_final := y_offset]
  }

  copy_loc_file = data.table::copy(location_file)

  new_x_coord = rep(0, nrow(copy_loc_file))
  new_y_coord = rep(0, nrow(copy_loc_file))

  for(row in 1:nrow(copy_loc_file)) {

    myrow = copy_loc_file[row,]

    field_select = myrow[[field_col]]
    X_select = myrow[[X_coord_col]]
    Y_select = myrow[[Y_coord_col]]

    X_offset = offset_file[field == field_select][['x_offset_final']]
    Y_offset = offset_file[field == field_select][['y_offset_final']]

    final_x = X_select+X_offset
    final_y = Y_select+Y_offset

    new_x_coord[row] = final_x
    new_y_coord[row] = final_y

  }

  if(reverse_final_x == TRUE) new_x_coord = new_x_coord*-1
  if(reverse_final_y == TRUE) new_y_coord = new_y_coord*-1

  copy_loc_file = data.table(copy_loc_file)

  copy_loc_file[, c('X_final', 'Y_final') := list(new_x_coord, new_y_coord)]

  return(copy_loc_file)
}


#' @title stitchTileCoordinates
#' @description Helper function to stitch tile coordinates together to form one complete picture
#' @param location_file location dataframe with X and Y coordinates
#' @param Xtilespan numerical value specifying the width of each tile
#' @param Ytilespan numerical value specifying the height of each tile
#' @export
stitchTileCoordinates <- function (location_file,
                                   Xtilespan,
                                   Ytilespan) {

  # data.table variables
  Xcoord = X.X = XtileIndex = Ycoord = Y.Y = YtileIndex = NULL

  if (is.null(location_file$X.X)){
    print("X coordinates missing in input file.")
  }else if (is.null(location_file$Y.Y)){
    print("Y coordinates missing in input file.")
  } else if (is.null(location_file$XtileIndex)){
    print("X tile index missing in input file.")
  }else if (is.null(location_file$YtileIndex)){
    print("Y tile index missing in input file.")
  }else{
    copy_loc_file = data.table::copy(location_file)
    copy_loc_file[,Xcoord := X.X + Xtilespan*(XtileIndex-1)]
    copy_loc_file[,Ycoord := Y.Y + Ytilespan*(YtileIndex-1)]
    return(copy_loc_file)
  }
}





#' @title my_arowMeans
#' @description arithmic rowMeans that works for a single column
#' @keywords internal
my_arowMeans = function(x) {
  if(is.null(nrow(x))) {
    x # if only one column is selected
    #mean(x)
  } else {
    rowMeans_giotto(x)
  }
}

#' @title my_growMeans
#' @description geometric rowMeans that works for a single column
#' @keywords internal
my_growMeans = function(x, offset = 0.1) {
  if(is.null(nrow(x))) {
    x # if only one column is selected
    #exp(mean(log(x+offset)))-offset
  } else {
    exp(rowMeans_giotto(log(x+offset)))-offset
  }
}

#' @title my_rowMeans
#' @description arithmic or geometric rowMeans that works for a single column
#' @keywords internal
my_rowMeans = function(x, method = c('arithmic', 'geometric'), offset = 0.1) {
  method = match.arg(method, c('arithmic', 'geometric'))
  if(method == 'arithmic') return(my_arowMeans(x))
  if(method == 'geometric') return(my_growMeans(x))
}



## matrix binarization methods ####

#' @title kmeans_binarize
#' @name kmeans_binarize
#' @description create binarized scores from a vector using kmeans
#' @keywords internal
kmeans_binarize = function(x, nstart = 3, iter.max = 10, set.seed = NULL) {

  if(!is.null(set.seed)) set.seed(1234)
  sel_gene_km = stats::kmeans(x, centers = 2, nstart = nstart, iter.max = iter.max)$cluster
  mean_1 = mean(x[sel_gene_km == 1])
  mean_2 = mean(x[sel_gene_km == 2])

  if(mean_1 > mean_2) {
    mean_1_value = 1
    mean_2_value = 0
  } else {
    mean_1_value = 0
    mean_2_value = 1
  }

  sel_gene_bin = x
  sel_gene_bin[sel_gene_km == 1] = mean_1_value
  sel_gene_bin[sel_gene_km == 2] = mean_2_value

  return(sel_gene_bin)

}

#' @title kmeans_arma_binarize
#' @name kmeans_arma_binarize
#' @description create binarized scores from a vector using kmeans_arma
#' @keywords internal
kmeans_arma_binarize = function(x, n_iter = 5, set.seed = NULL) {


  if(!is.null(set.seed)) set.seed(1234)
  sel_gene_km_res = ClusterR::KMeans_arma(data = as.matrix(x),
                                          clusters = 2,
                                          n_iter = n_iter)
  sel_gene_km = ClusterR::predict_KMeans(data = as.matrix(x),
                                         CENTROIDS = sel_gene_km_res)

  mean_1 = mean(x[sel_gene_km == 1])
  mean_2 = mean(x[sel_gene_km == 2])

  if(mean_1 > mean_2) {
    mean_1_value = 1
    mean_2_value = 0
  } else {
    mean_1_value = 0
    mean_2_value = 1
  }

  sel_gene_bin = x
  sel_gene_bin[sel_gene_km == 1] = mean_1_value
  sel_gene_bin[sel_gene_km == 2] = mean_2_value

  return(sel_gene_bin)

}

#' @title kmeans_arma_subset_binarize
#' @name kmeans_arma_subset_binarize
#' @description create binarized scores from a subsetted vector using kmeans_arma
#' @keywords internal
kmeans_arma_subset_binarize = function(x, n_iter = 5, extreme_nr = 20, sample_nr = 200, set.seed = NULL) {

  length_x = length(x)

  vector_x = sort(x)
  first_set = vector_x[1:extreme_nr]
  last_set = vector_x[(length_x-(extreme_nr-1)):length_x]
  random_set = sample(vector_x[(extreme_nr+1):(length_x-extreme_nr)], size = sample_nr)
  testset = c(first_set, last_set, random_set)

  if(!is.null(set.seed)) set.seed(1234)
  sel_gene_km_res = ClusterR::KMeans_arma(data = as.matrix(testset),
                                          clusters = 2,
                                          n_iter = n_iter)
  sel_gene_km = ClusterR::predict_KMeans(data = as.matrix(x),
                                         CENTROIDS = sel_gene_km_res)

  mean_1 = mean(x[sel_gene_km == 1])
  mean_2 = mean(x[sel_gene_km == 2])

  if(mean_1 > mean_2) {
    mean_1_value = 1
    mean_2_value = 0
  } else {
    mean_1_value = 0
    mean_2_value = 1
  }

  sel_gene_bin = x
  sel_gene_bin[sel_gene_km == 1] = mean_1_value
  sel_gene_bin[sel_gene_km == 2] = mean_2_value

  return(sel_gene_bin)

}


#' @title kmeans_binarize_wrapper
#' @name kmeans_binarize_wrapper
#' @description wrapper for different binarization functions
#' @keywords internal
kmeans_binarize_wrapper = function(gobject,
                                   expression_values = c('normalized', 'scaled', 'custom'),
                                   subset_genes = NULL,
                                   kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                                   nstart = 3,
                                   iter_max = 10,
                                   extreme_nr = 50,
                                   sample_nr = 50,
                                   set.seed = NULL) {


  # expression
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
  }

  # check parameter
  kmeans_algo = match.arg(arg = kmeans_algo, choices = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'))

  if(kmeans_algo == 'kmeans') {
    bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = kmeans_binarize,
                                nstart = nstart, iter.max = iter_max, set.seed = set.seed))
  } else if(kmeans_algo == 'kmeans_arma') {
    bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = kmeans_arma_binarize,
                                n_iter = iter_max, set.seed = set.seed))
  } else if(kmeans_algo == 'kmeans_arma_subset') {
    bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = kmeans_arma_subset_binarize,
                                n_iter = iter_max,
                                extreme_nr = extreme_nr,
                                sample_nr = sample_nr,
                                set.seed = set.seed))
  }

  return(bin_matrix)

}




#' @title rank_binarize
#' @name rank_binarize
#' @description create binarized scores from a vector using arbitrary rank
#' @keywords internal
rank_binarize = function(x, max_rank = 200) {

  sel_gene_rank = rank(-x, ties.method = 'average')

  sel_gene_bin = x
  sel_gene_bin[sel_gene_rank <= max_rank] = 1
  sel_gene_bin[sel_gene_rank > max_rank] = 0

  return(sel_gene_bin)

}


## data.table helper functions ####

#' @title DT_removeNA
#' @name DT_removeNA
#' @description set NA values to 0 in a data.table object
#' @keywords internal
DT_removeNA = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0]
  return(DT)
}


#' @title sort_combine_two_DT_columns
#' @name sort_combine_two_DT_columns
#' @description fast sorting and pasting of 2 character columns in a data.table
#' @keywords internal
sort_combine_two_DT_columns = function(DT,
                                       column1,
                                       column2,
                                       myname = 'unif_gene_gene') {

  # data.table variables
  values_1_num = values_2_num = scolumn_1 = scolumn_2 = unif_sort_column = NULL

  # maybe faster with converting to factors??

  # make sure columns are character
  selected_columns = c(column1, column2)
  DT[,(selected_columns):= lapply(.SD, as.character), .SDcols = selected_columns]

  # convert characters into numeric values
  uniq_values = sort(unique(c(DT[[column1]], DT[[column2]])))
  uniq_values_num = 1:length(uniq_values)
  names(uniq_values_num) = uniq_values


  DT[,values_1_num := uniq_values_num[get(column1)]]
  DT[,values_2_num := uniq_values_num[get(column2)]]


  DT[, scolumn_1 := ifelse(values_1_num < values_2_num, get(column1), get(column2))]
  DT[, scolumn_2 := ifelse(values_1_num < values_2_num, get(column2), get(column1))]

  DT[, unif_sort_column := paste0(scolumn_1,'--',scolumn_2)]
  DT[, c('values_1_num', 'values_2_num', 'scolumn_1', 'scolumn_2') := NULL]
  data.table::setnames(DT, 'unif_sort_column', myname)

  return(DT)
}





## package checks ####

#' @title package_check
#' @name package_check
#' @param pkg_name name of package
#' @param repository where is the package
#' @param github_repo name of github repository if needed
#' @description check if package is available and provide installation instruction if not available
#' @keywords internal
package_check = function(pkg_name,
                         repository = c('CRAN', 'Bioc', 'github'),
                         github_repo = NULL) {

  repository = match.arg(repository, choices = c('CRAN', 'Bioc', 'github'))

  if(repository == 'CRAN') {

    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "install.packages('",pkg_name,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }


  } else if(repository == 'Bioc') {

    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager');
         BiocManager::install('",pkg_name,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }

  } else if(repository == 'github') {

    if(is.null(github_repo)) stop("provide the github repo of package, e.g. 'johndoe/cooltool' ")

    if(!requireNamespace(pkg_name, quietly = TRUE)) {
      stop("\n package ", pkg_name ," is not yet installed \n",
           "To install: \n",
           "devtools::install_github('",github_repo,"')",
           call. = FALSE)
    } else {
      return(TRUE)
    }

  }

}



## dataset helpers ####

#' @title getSpatialDataset
#' @name getSpatialDataset
#' @param dataset dataset to download
#' @param directory directory to save the data to
#' @param \dots additional parameters to \code{\link[utils]{download.file}}
#' @description This package will automatically download the spatial locations and
#' expression matrix for the chosen dataset. These files are already in the right format
#' to create a Giotto object. If wget is installed on your machine, you can add
#' 'method = wget' to the parameters to download files faster.
#' @export
getSpatialDataset = function(dataset = c('ST_OB1',
                                         'ST_OB2',
                                         'codex_spleen',
                                         'cycif_PDAC',
                                         'starmap_3D_cortex',
                                         'osmfish_SS_cortex',
                                         'merfish_preoptic',
                                         'seqfish_SS_cortex',
                                         'seqfish_OB',
                                         'slideseq_cerebellum'),
                             directory = getwd(),
                             ...) {

  sel_dataset = match.arg(dataset, choices = c('ST_OB1',
                                               'ST_OB2',
                                               'codex_spleen',
                                               'cycif_PDAC',
                                               'starmap_3D_cortex',
                                               'osmfish_SS_cortex',
                                               'merfish_preoptic',
                                               'seqfish_SS_cortex',
                                               'seqfish_OB',
                                               'slideseq_cerebellum'))

  # check operating system first
  os_specific_system = get_os()

  #if(os_specific_system == 'windows') {
  #  stop('This function is currently not supported on windows systems,
  #       please visit https://github.com/RubD/spatial-datasets and manually download your files')
  #}


  # check directory
  if(!file.exists(directory)) {
    warning('The output directory does not exist and will be created \n')
    dir.create(directory, recursive = T)
  }

  datasets_file = system.file("extdata", "datasets.txt", package = 'Giotto')
  datasets_file = data.table::fread(datasets_file)



  ## check if wget is installed
  #message = system("if ! command -v wget &> /dev/null
  #                  then
  #                  echo 'wget could not be found, please install wget first'
  #                  exit
  #                  fi", intern = TRUE)

  #if(identical(message, character(0))) {
  #  print('wget was found, start downloading datasets: ')
  #} else {
  #  stop(message)
  #}

  ## alternative
  #wget_works = try(system('command -v wget', intern = T))

  #if(class(wget_works) == 'try-error' | is.na(wget_works[1])) {
  #  stop('wget was not found, please install wget first \n')
  #} else {
  #  print('wget was found, start downloading datasets: \n')
  #}



  # get url to spatial locations and download
  spatial_locs_url = datasets_file[dataset == sel_dataset][['spatial_locs']]
  myfilename = basename(spatial_locs_url)
  mydestfile = paste0(directory,'/', myfilename)

  print(spatial_locs_url)
  print(mydestfile)

  utils::download.file(url = spatial_locs_url, destfile = mydestfile, ...)

  #system(paste0("wget -P ", "'",directory,"'"," ", spatial_locs_url))


  # get url to expression matrix and download
  expr_matrix_url = datasets_file[dataset == sel_dataset][['expr_matrix']]
  myfilename = basename(expr_matrix_url)
  mydestfile = paste0(directory,'/', myfilename)
  utils::download.file(url = expr_matrix_url, destfile = mydestfile, ...)

  #system(paste0("wget -P ", "'",directory,"'"," ", expr_matrix_url))

  # get url(s) to additional metadata files and download
  metadata_url = datasets_file[dataset == sel_dataset][['metadata']][[1]]
  metadata_url = unlist(strsplit(metadata_url, split = '\\|'))

  if(identical(metadata_url, character(0))) {
    NULL
  } else {
    for(url in metadata_url) {
      myfilename = basename(url)
      mydestfile = paste0(directory,'/', myfilename)
      utils::download.file(url = url, destfile = mydestfile, ...)
      #system(paste0("wget -P ", "'",directory,"'"," ", url))
    }
  }

}


#' @title get10Xmatrix
#' @description This function creates an expression matrix from a 10X structured folder
#' @param path_to_data path to the 10X folder
#' @param gene_column_index which column from the features or genes .tsv file to use for row ids
#' @return sparse expression matrix from 10X
#' @details A typical 10X folder is named raw_feature_bc_matrix or raw_feature_bc_matrix and it has 3 files:
#' \itemize{
#'   \item{barcodes.tsv(.gz)}
#'   \item{features.tsv(.gz) or genes.tsv(.gz)}
#'   \item{matrix.mtx(.gz)}
#' }
#' By default the first column of the features or genes .tsv file will be used, however if multiple
#' annotations are provided (e.g. ensembl gene ids and gene symbols) the user can select another column.
#' @export
get10Xmatrix = function(path_to_data, gene_column_index = 1) {

  # data.table variables
  total = gene_symbol = gene_id = gene_id_num = cell_id = cell_id_num = sort_gene_id_num = NULL

  # data directory
  files_10X = list.files(path_to_data)

  # get barcodes and create vector
  barcodes_file = grep(files_10X, pattern = 'barcodes', value = T)
  barcodesDT = fread(input = paste0(path_to_data,'/',barcodes_file), header = F)
  barcodes_vec = barcodesDT$V1
  names(barcodes_vec) = 1:nrow(barcodesDT)

  # get features and create vector
  features_file = grep(files_10X, pattern = 'features|genes', value = T)
  featuresDT = fread(input = paste0(path_to_data,'/',features_file), header = F)

  g_name = colnames(featuresDT)[gene_column_index]
  ## convert ensembl gene id to gene symbol ##
  ## TODO

  featuresDT[, total := .N, by = get(g_name)]
  featuresDT[, gene_symbol := ifelse(total > 1, paste0(get(g_name),'--',1:.N), get(g_name)), by = get(g_name)]
  features_vec = featuresDT$gene_symbol
  names(features_vec) = 1:nrow(featuresDT)

  # get matrix
  matrix_file = grep(files_10X, pattern = 'matrix', value = T)
  matrixDT = fread(input = paste0(path_to_data,'/',matrix_file), header = F, skip = 3)
  colnames(matrixDT) = c('gene_id_num', 'cell_id_num', 'umi')

  # convert barcodes and features
  matrixDT[, gene_id := features_vec[gene_id_num]]
  matrixDT[, cell_id := barcodes_vec[cell_id_num]]

  # make sure that gene id are consecutive
  sort_gene_id_vec = 1:length(unique(matrixDT$gene_id))
  names(sort_gene_id_vec) = unique(matrixDT$gene_id)
  matrixDT[, sort_gene_id_num := sort_gene_id_vec[gene_id]]

  sparsemat = Matrix::sparseMatrix(i = matrixDT$sort_gene_id_num, j = matrixDT$cell_id_num, x = matrixDT$umi,
                                   dimnames = list(unique(matrixDT$gene_id), unique(matrixDT$cell_id)))

  return(sparsemat)

  # create a final matrix
  #matrix_ab = data.table::dcast.data.table(data = matrixDT, gene_id~cell_id, value.var = 'umi')
  #matrix_ab_mat = dt_to_matrix(matrix_ab)
  #matrix_ab_mat[is.na(matrix_ab_mat)] = 0

}



#' @title get10Xmatrix_h5
#' @description This function creates an expression matrix from a 10X h5 file path
#' @param path_to_data path to the 10X .h5 file
#' @param gene_ids use gene symbols (default) or ensembl ids for the gene expression matrix
#' @return (list of) sparse expression matrix from 10X
#' @details If the .h5 10x file has multiple modalities (e.g. RNA and protein),
#'  multiple matrices will be returned
#' @export
get10Xmatrix_h5 = function(path_to_data, gene_ids = c('symbols', 'ensembl')) {

  ## function inspired by and modified from the VISION package
  ## see read_10x_h5_v3 in https://github.com/YosefLab/VISION/blob/master/R/Utilities.R

  # verify if optional package is installed
  package_check(pkg_name = "hdf5r", repository = "CRAN")

  # select parameter
  gene_ids = match.arg(gene_ids, choices = c('symbols', 'ensembl'))

  h5 = hdf5r::H5File$new(path_to_data)

  tryCatch({

    # list objects part of the h5 file
    # hdf5r::list.objects(h5)

    # get root folder name e.g. 'matrix'
    root <- names(h5)
    root <- root[1]

    # extraction information
    data <- h5[[paste0(root, "/data")]][]
    data <- as.numeric(data)

    genome = unique(h5[[paste0(root, "/features/genome")]][])
    barcodes = h5[[paste0(root, "/barcodes")]][]
    feature_id = h5[[paste0(root, "/features/id")]][]
    data_shape = h5[[paste0(root, "/shape")]][]
    feature_names = h5[[paste0(root, "/features/name")]][]
    feature_tag_keys = h5[[paste0(root, "/features/_all_tag_keys")]][]
    feature_types = unique(h5[[paste0(root, "/features/feature_type")]][])
    indices = h5[[paste0(root, "/indices")]][]
    indptr = h5[[paste0(root, "/indptr")]][]

    # create a feature data.table
    features_dt = data.table::data.table(
      'id' = feature_id,
      'name' = feature_names,
      'feature_type' = feature_types,
      'genome' = genome
    )

    # create uniq name symbols
    # duplicate gene symbols will be given a suffix '_1', '_2', ...

    # data.table variables
    nr_name = name = uniq_name = NULL

    features_dt[, nr_name := 1:.N, by = name]
    features_dt[, uniq_name := ifelse(nr_name == 1, name, paste0(name, '_', (nr_name-1)))]


    # dimension names
    dimnames = list(feature_id, barcodes)

    sparsemat = Matrix::sparseMatrix(i = indices + 1,
                                     p = indptr,
                                     x = data,
                                     dims = data_shape,
                                     dimnames = dimnames)

    # multiple modalities? add for future improvement
    result_list = list()

    for(modality in unique(feature_types)) {

      result_list[[modality]] = sparsemat[features_dt$feature_type == modality, ]

      # change names to gene symbols if it's expression
      if(modality == 'Gene Expression' & gene_ids == 'symbols') {

        conv_vector = features_dt$uniq_name
        names(conv_vector) = features_dt$id

        current_names = rownames(result_list[[modality]])
        new_names = conv_vector[current_names]
        rownames(result_list[[modality]]) = new_names
      }
    }

  },
  finally = {
    h5$close_all()
  })

  return(result_list)

}



#' @title convertEnsemblToGeneSymbol
#' @description This function convert ensembl gene IDs from a matrix to official gene symbols
#' @param matrix an expression matrix with ensembl gene IDs as rownames
#' @param species species to use for gene symbol conversion
#' @return expression matrix with gene symbols as rownames
#' @details This function requires that the biomaRt library is installed
#' @export
convertEnsemblToGeneSymbol = function(matrix,
                                      species = c('mouse', 'human')) {

  # data.table: set global variable
  dupes = mgi_symbol = gene_symbol = ensembl_gene_id = hgnc_symbol = NULL

  if("biomaRt" %in% rownames(installed.packages()) == FALSE) {
    cat("\n package 'biomaRt' is not yet installed and is required for this function \n")
  }

  species = match.arg(species, choices = c('mouse', 'human'))

  if(species == 'mouse') {

    # ensembl IDs to change
    ensemblsIDS = rownames(matrix)

    # prepare ensembl database
    ensembl = biomaRt::useMart("ensembl",
                               dataset = "mmusculus_gene_ensembl")
    gene_names = biomaRt::getBM(attributes= c('mgi_symbol', 'ensembl_gene_id'),
                                filters = 'ensembl_gene_id',
                                values = ensemblsIDS,
                                mart = ensembl)
    gene_names_DT = data.table::as.data.table(gene_names)
    gene_names_DT[, dupes := duplicated(mgi_symbol)]
    gene_names_DT[, gene_symbol := ifelse(any(dupes) == FALSE, mgi_symbol,
                                          ifelse(mgi_symbol == "", ensembl_gene_id, 'temporary')), by = mgi_symbol]
    gene_names_DT[, gene_symbol := ifelse(mgi_symbol == '', ensembl_gene_id, gene_symbol)]
    gene_names_DT[, gene_symbol := ifelse(gene_symbol == 'temporary', paste0(mgi_symbol,'--', 1:.N), gene_symbol), by = mgi_symbol]

    # filter
    matrix = matrix[rownames(matrix) %in% gene_names_DT$ensembl_gene_id, ]

    # create swapping vector
    new_symbols = gene_names_DT[['gene_symbol']]
    names(new_symbols) = gene_names_DT[['ensembl_gene_id']]

    # replace
    new_rownames = new_symbols[rownames(matrix)]
    rownames(matrix) = new_rownames

    return(matrix)

  }

  if(species == 'human') {

    # ensembl IDs to change
    ensemblsIDS = rownames(matrix)

    # prepare ensembl database
    ensembl = biomaRt::useMart("ensembl",
                               dataset = "hsapiens_gene_ensembl")
    gene_names = biomaRt::getBM(attributes= c('hgnc_symbol', 'ensembl_gene_id'),
                                filters = 'ensembl_gene_id',
                                values = ensemblsIDS,
                                mart = ensembl)
    gene_names_DT = data.table::as.data.table(gene_names)
    gene_names_DT[, dupes := duplicated(hgnc_symbol)]
    gene_names_DT[, gene_symbol := ifelse(any(dupes) == FALSE, hgnc_symbol,
                                          ifelse(hgnc_symbol == "", ensembl_gene_id, 'temporary')), by = hgnc_symbol]
    gene_names_DT[, gene_symbol := ifelse(hgnc_symbol == '', ensembl_gene_id, gene_symbol)]
    gene_names_DT[, gene_symbol := ifelse(gene_symbol == 'temporary', paste0(hgnc_symbol,'--', 1:.N), gene_symbol), by = hgnc_symbol]

    # filter
    matrix = matrix[rownames(matrix) %in% gene_names_DT$ensembl_gene_id, ]

    # create swapping vector
    new_symbols = gene_names_DT[['gene_symbol']]
    names(new_symbols) = gene_names_DT[['ensembl_gene_id']]

    # replace
    new_rownames = new_symbols[rownames(matrix)]
    rownames(matrix) = new_rownames

    return(matrix)

  }


}




