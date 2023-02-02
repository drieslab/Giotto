


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

  if(n < 1) stop('Error in getDistinctColors: number of colors wanted must be at least 1')

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


#' @title Over-allocation for giotto DT-based info
#' @description Finds DT based objects, overallocates the data.tables, then sets
#' the objects back in the giotto object
#' @param gobject giotto object
#' @keywords internal
giotto_alloc_dt_slots = function(gobject) {

  # data.table vars
  spat_unit = feat_type = name = NULL

  # metadata
  avail_cm = list_cell_metadata(gobject)
  if(!is.null(avail_cm)) {
    for(cm_i in seq(nrow(avail_cm))) {
      cm = get_cell_metadata(gobject = gobject,
                             spat_unit = avail_cm[cm_i, spat_unit],
                             feat_type = avail_cm[cm_i, feat_type],
                             output = 'cellMetaObj',
                             copy_obj = FALSE)
      if(!is.null(cm[])) {
        cm[] = data.table::setalloccol(cm[])
        gobject = set_cell_metadata(gobject = gobject,
                                    metadata = cm,
                                    set_defaults = FALSE,
                                    verbose = FALSE)
      }
    }
  }

  avail_fm = list_feat_metadata(gobject)
  if(!is.null(avail_fm)) {
    for(fm_i in seq(nrow(avail_fm))) {
      fm = get_feature_metadata(gobject = gobject,
                                spat_unit = avail_fm[fm_i, spat_unit],
                                feat_type = avail_fm[fm_i, feat_type],
                                output = 'featMetaObj',
                                copy_obj = FALSE)
      if(!is.null(fm[])) {
        fm[] = data.table::setalloccol(fm[])
        gobject = set_feature_metadata(gobject = gobject,
                                       metadata = fm,
                                       set_defaults = FALSE,
                                       verbose = FALSE)
      }
    }
  }

  # spatlocs
  avail_sl = list_spatial_locations(gobject)
  if(!is.null(avail_sl)) {
    for(sl_i in seq(nrow(avail_sl))) {
      sl = get_spatial_locations(gobject = gobject,
                                 spat_unit = avail_sl[sl_i, spat_unit],
                                 spat_loc_name = avail_sl[sl_i, name],
                                 output = 'spatLocsObj',
                                 copy_obj = FALSE)
      if(!is.null(sl[])) {
        sl[] = data.table::setalloccol(sl[])
        gobject = set_spatial_locations(gobject = gobject,
                                        spatlocs = sl,
                                        set_defaults = FALSE,
                                        verbose = FALSE)
      }
    }
  }

  # spatial enrichment
  avail_se = list_spatial_enrichments(gobject)
  if(!is.null(avail_se)) {
    for(se_i in seq(nrow(avail_se))) {
      se = get_spatial_enrichment(gobject = gobject,
                                  spat_unit = avail_se[se_i, spat_unit],
                                  feat_type = avail_se[se_i, feat_type],
                                  enrichm_name = avail_se[se_i, name],
                                  output = 'spatEnrObj',
                                  copy_obj = FALSE)
      if(!is.null(se[])) {
        se[] = data.table::setalloccol(se[])
        gobject = set_spatial_enrichment(gobject = gobject,
                                         spatenrichment = se,
                                         set_defaults = FALSE,
                                         verbose = FALSE)
      }
    }
  }

  # spatial network
  avail_sn = list_spatial_networks(gobject)
  if(!is.null(avail_sn)) {
    for(sn_i in seq(nrow(avail_sn))) {
      sn = get_spatialNetwork(gobject = gobject,
                              spat_unit = avail_sn[sn_i, spat_unit],
                              name = avail_sn[sn_i, name],
                              output = 'spatialNetworkObj')
      if(!is.null(slot(sn, 'networkDT_before_filter'))) {
        slot(sn, 'networkDT_before_filter') = data.table::setalloccol(slot(sn, 'networkDT_before_filter'))
      }
      if(!is.null(sn[])) {
        sn[] = data.table::setalloccol(sn[])
        gobject = set_spatialNetwork(gobject = gobject,
                                     spatial_network = sn,
                                     verbose = FALSE,
                                     set_defaults = FALSE)
      }
    }
  }

  # spatial grid
  avail_sg = list_spatial_grids(gobject)
  if(!is.null(avail_sg)) {
    for(sg_i in seq(nrow(avail_sg))) {
      sg = get_spatialGrid(gobject = gobject,
                           spat_unit = avail_sg[sg_i, spat_unit],
                           feat_type = avail_sg[sg_i, feat_type],
                           name = avail_sg[sg_i, feat_type],
                           return_grid_Obj = TRUE)
      if(!is.null(sg[])) {
        sg[] = data.table::setalloccol(sg[])
        gobject = set_spatialGrid(gobject = gobject,
                                  spatial_grid = sg,
                                  verbose = FALSE,
                                  set_defaults = FALSE)
      }
    }
  }
  return(gobject)
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
#' @name stitchFieldCoordinates
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
#' @name stitchTileCoordinates
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
    rowMeans_flex(x)
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
    exp(rowMeans_flex(log(x+offset)))-offset
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
kmeans_binarize_wrapper = function(expr_values,
                                   subset_feats = NULL,
                                   kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                                   nstart = 3,
                                   iter_max = 10,
                                   extreme_nr = 50,
                                   sample_nr = 50,
                                   set.seed = NULL) {



  # expression values
  if(!is.null(subset_feats)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_feats, ]
  }

  # check parameter
  kmeans_algo = match.arg(arg = kmeans_algo, choices = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'))

  if(kmeans_algo == 'kmeans') {
    bin_matrix = t_flex(apply(X = expr_values, MARGIN = 1, FUN = kmeans_binarize,
                              nstart = nstart, iter.max = iter_max, set.seed = set.seed))
  } else if(kmeans_algo == 'kmeans_arma') {
    bin_matrix = t_flex(apply(X = expr_values, MARGIN = 1, FUN = kmeans_arma_binarize,
                              n_iter = iter_max, set.seed = set.seed))
  } else if(kmeans_algo == 'kmeans_arma_subset') {
    bin_matrix = t_flex(apply(X = expr_values, MARGIN = 1, FUN = kmeans_arma_subset_binarize,
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



#' @title rank_binarize_wrapper
#' @name rank_binarize_wrapper
#' @description wrapper for rank binarization function
#' @keywords internal
rank_binarize_wrapper = function(expr_values,
                                 subset_feats = NULL,
                                 percentage_rank = 30) {

  # expression values
  if(!is.null(subset_feats)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_feats, ]
  }

  max_rank = (ncol(expr_values)/100)*percentage_rank
  bin_matrix = t_flex(apply(X = expr_values, MARGIN = 1, FUN = rank_binarize, max_rank = max_rank))

  return(bin_matrix)
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


#' @title Check for updates to Giotto Suite
#' @name check_github_suite_ver
#' @description Checks the Giotto Suite github repository and compares the version
#' number to the currently installed.
#' @keywords internal
check_github_suite_ver = function() {
  current_ver = utils::packageVersion('Giotto')
  url = paste0('https://raw.githubusercontent.com/drieslab/Giotto/suite/DESCRIPTION')
  # suppress warnings and errors if inaccessible
  x = suppressWarnings(try(readLines(url), silent = TRUE))
  if(!inherits(x, 'try-error')) {
    gh_ver = x[grep(pattern = 'Version:', x)]
    gh_ver = gsub(pattern = 'Version: ', replacement = '', gh_ver)
    ver_compare = utils::compareVersion(gh_ver, as.character(current_ver))

    if(ver_compare == 1) wrap_msg('Newer devel version of Giotto on GitHub:', gh_ver)
  }
}






#' @title package_check
#' @name package_check
#' @param pkg_name name of package
#' @param repository where is the package
#' @param github_repo name of github repository if needed
#' @param optional whether the package is optional. \code{stop()} is used if TRUE
#' and only \code{message()} will be sent if FALSE.
#' @param custom_msg custom message to be sent instead of default error or message
#' @description check if package is available and provide installation instruction if not available
#' @keywords internal
package_check = function(pkg_name,
                         repository = c('CRAN', 'Bioc', 'github', 'pip'),
                         github_repo = NULL,
                         optional = FALSE,
                         custom_msg = NULL) {

  repository = match.arg(repository, choices = c('CRAN', 'Bioc', 'github', 'pip'))

  check_message = function(default_msg, custom_msg, optional) {
    if(!isTRUE(optional)) {
      if(is.null(custom_msg)) stop(default_msg, call. = FALSE)
      else stop(custom_msg, call. = FALSE)
    } else {
      if(is.null(custom_msg)) message(default_msg)
      else message(custom_msg)
    }
  }

  if(repository == 'CRAN') {

    default_msg = c("\n package ", pkg_name ," is not yet installed \n",
    "To install: \n",
    "install.packages('",pkg_name,"')")

    if(!requireNamespace(pkg_name, quietly = TRUE)) {

      check_message(default_msg = default_msg,
                    custom_msg = custom_msg,
                    optional = optional)

    } else {
      return(TRUE)
    }


  } else if(repository == 'Bioc') {

    default_msg = c("\n package ", pkg_name ," is not yet installed \n",
                    "To install: \n",
                    "if(!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager');\nBiocManager::install('",pkg_name,"')")

    if(!requireNamespace(pkg_name, quietly = TRUE)) {

      check_message(default_msg = default_msg,
                    custom_msg = custom_msg,
                    optional = optional)

    } else {
      return(TRUE)
    }

  } else if(repository == 'github') {

    if(is.null(github_repo)) stop(wrap_txt("provide the github repo of package, e.g. 'johndoe/cooltool' ", sep = ''))

    default_msg = c("\n package ", pkg_name ," is not yet installed \n",
                    "To install: \n",
                    "devtools::install_github('",github_repo,"')")

    if(!requireNamespace(pkg_name, quietly = TRUE)) {

      check_message(default_msg = default_msg,
                    custom_msg = custom_msg,
                    optional = optional)

    } else {
      return(TRUE)
    }

  } else if(repository == 'pip') {

    default_msg = c("\n package ", pkg_name ," is not yet installed \n",
                    "To install for default Giotto miniconda environment: \n",
                    "reticulate::conda_install(envname = 'giotto_env',packages = '",pkg_name,"',pip = TRUE)")

    if(!reticulate::py_module_available(pkg_name)) {

      check_message(default_msg = default_msg,
                    custom_msg = custom_msg,
                    optional = optional)

    }
  }

}



## I/O helpers ####

#' @title saveGiotto
#' @name saveGiotto
#' @description Saves a Giotto object to a specific folder structure
#' @param gobject Giotto object
#' @param foldername Folder name
#' @param dir Directory where to create the folder
#' @param method method to save main object
#' @param method_params additional method parameters for RDS or qs
#' @param overwrite Overwrite existing folders
#' @param image_filetype the image filetype to use, see \code{\link[terra]{writeRaster}}
#' @param verbose be verbose
#' @param ... additional parameters for \code{\link[terra]{writeRaster}}
#' @return Creates a directory with Giotto object information
#' @details Works together with \code{\link{loadGiotto}} to save and re-load
#' Giotto objects. Additional method_params need to be provided as a list and will
#' go to \code{\link[base]{saveRDS}} or \code{\link[qs]{qsave}}
#' @export
saveGiotto = function(gobject,
                      foldername = 'saveGiottoDir',
                      dir = getwd(),
                      method = c('RDS', 'qs'),
                      method_params = list(),
                      overwrite = FALSE,
                      image_filetype = 'PNG',
                      verbose = TRUE,
                      ...) {


  final_dir = paste0(dir,'/', foldername)

  if(dir.exists(final_dir)) {
    if(overwrite == FALSE) {
      stop('Folder already exist and overwrite = FALSE, abort saving \n')
    } else {
      wrap_msg('Folder already exist and overwrite = TRUE, overwrite folder \n')
      unlink(x = final_dir, recursive = TRUE)
      dir.create(final_dir)
    }
  } else {
    dir.create(final_dir, recursive = TRUE)
  }

  ## save spatVector objects related to feature information
  if(verbose) wrap_msg('1. Start writing feature information \n')
  feat_info_names = list_feature_info_names(gobject)

  if(!is.null(feat_info_names)) {
    feat_dir = paste0(final_dir,'/','Features')
    dir.create(feat_dir)
    for(feat in feat_info_names) {
      if(verbose) wrap_msg('For feature: ', feat, '\n')

      # original spatvector
      if(!is.null(gobject@feat_info[[feat]]@spatVector)) {
        filename = paste0(feat_dir, '/', feat, '_feature_spatVector.shp')
        terra::writeVector(gobject@feat_info[[feat]]@spatVector, filename = filename)
      }

      # network
      # ? data.table object

    }
  }


  ## save spatVector objects related to spatial information
  if(verbose) wrap_msg('2. Start writing spatial information \n')
  spat_info_names = list_spatial_info_names(gobject)

  if(!is.null(spat_info_names)) {
    spatinfo_dir = paste0(final_dir,'/','SpatialInfo')
    dir.create(spatinfo_dir)
    for(spatinfo in spat_info_names) {

      if(verbose) wrap_msg('For spatial information: ', spatinfo, '\n')

      # original spatVectors
      if(!is.null(gobject@spatial_info[[spatinfo]]@spatVector)) {
        filename = paste0(spatinfo_dir, '/', spatinfo, '_spatInfo_spatVector.shp')
        terra::writeVector(gobject@spatial_info[[spatinfo]]@spatVector, filename = filename)
      }

      # spatVectorCentroids
      if(!is.null(gobject@spatial_info[[spatinfo]]@spatVectorCentroids)) {
        filename = paste0(spatinfo_dir, '/', spatinfo, '_spatInfo_spatVectorCentroids.shp')
        terra::writeVector(gobject@spatial_info[[spatinfo]]@spatVectorCentroids, filename = filename)
      }

      # overlap information
      if(!is.null(gobject@spatial_info[[spatinfo]]@overlaps)) {

        for(feature in names(gobject@spatial_info[[spatinfo]]@overlaps)) {
          filename = paste0(spatinfo_dir, '/', feature, '_', spatinfo, '_spatInfo_spatVectorOverlaps.shp')
          terra::writeVector(gobject@spatial_info[[spatinfo]]@overlaps[[feature]], filename = filename)
        }
      }

    }
  }



  ## save spatVector objects related to spatial information
  if(verbose) wrap_msg('3. Start writing image information \n')
  image_names = list_images_names(gobject, img_type = 'largeImage')

  if(!is.null(image_names)) {
    image_dir = paste0(final_dir,'/','Images')
    dir.create(image_dir)
    for(image in image_names) {
      if(verbose) wrap_msg('For image information: ', image, '\n')

      if(!is.null(gobject@largeImages[[image]]@raster_object)) {
        # save extent info just in case
        gobject@largeImages[[image]]@extent = terra::ext(gobject@largeImages[[image]]@raster_object)[]

        # save raster
        filename = paste0(image_dir, '/', image, '_spatRaster')
        terra::writeRaster(x = gobject@largeImages[[image]]@raster_object, filename = filename, filetype = image_filetype)
      }
    }
  }


  ## save whole Giotto object
  method = match.arg(arg = method, choices = c('RDS', 'qs'))

  if(method == 'RDS') {
    do.call('saveRDS', c(object = gobject, file = paste0(final_dir, '/', 'gobject.RDS'), method_params))
  } else if(method == 'qs') {
    package_check(pkg_name = 'qs', repository = 'CRAN')
    qsave_fun = get("qsave", asNamespace("qs"))
    do.call(qsave_fun, c(x = gobject, file = paste0(final_dir, '/', 'gobject.qs'), method_params))
  }

}


#' @title loadGiotto
#' @name loadGiotto
#' @description Saves a Giotto object to a specific folder structure
#' @param path_to_folder path to folder where Giotto object was stored with \code{\link{saveGiotto}}
#' @param load_params additional parameters for loading or reading giotto object
#' @param reconnect_giottoImage (default = TRUE) whether to attempt reconnection of magick based image objects
#' @param python_path (optional) manually set your python path
#' @param verbose be verbose
#' @return Giotto object
#' @details Works together with \code{\link{saveGiotto}} to save and re-load
#' Giotto objects.
#' Additional load_params need to be provided as a list and will
#' go to \code{\link[base]{readRDS}} or \code{\link[qs]{qread}}
#' You can set the python path, alternatively it will look for an existing
#' Giotto python environment.
#' @export
loadGiotto = function(path_to_folder,
                      load_params = list(),
                      reconnect_giottoImage = TRUE,
                      python_path = NULL,
                      verbose = TRUE) {

  # data.table vars
  img_type = NULL

  if(!file.exists(path_to_folder)) {
    stop('path_to_folder does not exist \n')
  }

  ## 1. load giotto object
  if(verbose) wrap_msg('1. read Giotto object \n')

  gobject_file = list.files(path_to_folder, pattern = 'gobject')

  if(grepl('.RDS', x = gobject_file)) {
    gobject = do.call('readRDS', c(file = paste0(path_to_folder,'/','gobject.RDS'), load_params))
  }

  if(grepl('.qs', x = gobject_file)) {
    package_check(pkg_name = 'qs', repository = 'CRAN')
    qread_fun = get("qread", asNamespace("qs"))
    gobject = do.call(qread_fun, c(file = paste0(path_to_folder,'/','gobject.qs'), load_params))
  }




  ## 2. read in features
  if(verbose) wrap_msg('2. read Giotto feature information \n')
  feat_files = list.files(path = paste0(path_to_folder, '/Features'), pattern = '.shp')
  if(length(feat_files) != 0) {
    feat_names = gsub(feat_files, pattern = '_feature_spatVector.shp', replacement = '')
    feat_paths = list.files(path = paste0(path_to_folder, '/Features'), pattern = '.shp', full.names = TRUE)
    for(feat_i in 1:length(feat_names)) {
      if(verbose) print(feat_paths[feat_i])
      spatVector = terra::vect(x = feat_paths[feat_i])
      feat_name = feat_names[feat_i]
      if(verbose) print(feat_name)
      gobject@feat_info[[feat_name]]@spatVector = spatVector
    }
  }

  ## 3. read in spatial polygons
  if(isTRUE(verbose)) wrap_msg('3. read Giotto spatial information \n')
  spat_paths = list.files(path = paste0(path_to_folder, '/SpatialInfo'), pattern = 'spatVector.shp', full.names = TRUE)
  spat_files = basename(spat_paths)
  if(length(spat_files) != 0) {

    ## 3.1. shapes
    if(isTRUE(verbose)) {
      wrap_msg('\n3.1 read Giotto spatial shape information \n')
      print(spat_files)
    }
    spat_names = gsub(spat_files, pattern = '_spatInfo_spatVector.shp', replacement = '')
    for(spat_i in 1:length(spat_names)) {
      spatVector = terra::vect(x = spat_paths[spat_i])
      spat_name = spat_names[spat_i]
      if(isTRUE(verbose)) message(spat_name)
      gobject@spatial_info[[spat_name]]@spatVector = spatVector
    }

    ## 3.2. centroids
    centroid_search_term = gsub(spat_files, pattern = '_spatInfo_spatVector.shp', replacement = '_spatInfo_spatVectorCentroids.shp')
    centroid_paths = sapply(centroid_search_term, function(gp_centroid) {
      list.files(path = paste0(path_to_folder, '/SpatialInfo'), pattern = gp_centroid, full.names = TRUE)
    }, USE.NAMES = FALSE)
    centroid_files = basename(centroid_paths)

    if(isTRUE(verbose)) {
      wrap_msg('\n3.2 read Giotto spatial centroid information \n')
      print(centroid_files)
    }
    if(length(centroid_files != 0)) {
      spat_names = gsub(centroid_files, pattern = '_spatInfo_spatVectorCentroids.shp', replacement = '')
      for(spat_i in 1:length(spat_names)) {
        spatVector = terra::vect(x = centroid_paths[spat_i])
        spat_name = spat_names[spat_i]
        if(isTRUE(verbose)) message(spat_name)
        gobject@spatial_info[[spat_name]]@spatVectorCentroids = spatVector
      }
    }


    ## 3.3. overlaps
    overlap_search_term = gsub(spat_files, pattern = '_spatInfo_spatVector.shp', replacement = '_spatInfo_spatVectorOverlaps.shp')
    overlap_files = list.files(path = paste0(path_to_folder, '/SpatialInfo'), pattern = 'spatVectorOverlaps.shp')

    if(isTRUE(verbose)) {
      wrap_msg('\n3.3 read Giotto spatial overlap information \n')
      print(overlap_files)
    }
    if(length(overlap_files != 0)) {
      # find overlaps per spatVector
      for(sv_i in seq_along(overlap_search_term)) {
        overlap_paths = list.files(path = paste0(path_to_folder, '/SpatialInfo'), pattern = overlap_search_term[sv_i], full.names = TRUE)
        overlap_filenames = basename(overlap_paths)

        for(spat_i in seq_along(overlap_filenames)) {
          spatVector = terra::vect(x = overlap_paths[spat_i])

          feat_name = gsub(overlap_filenames[spat_i], pattern = paste0('_', overlap_search_term[sv_i]), replacement = '')
          spat_name = gsub(overlap_filenames[spat_i], pattern = paste0(feat_name, '_'), replacement = '')
          spat_name = gsub(spat_name, pattern = '_spatInfo_spatVectorOverlaps.shp', replacement = '')

          if(isTRUE(verbose)) wrap_msg(spat_name, ' and ', feat_name)
          gobject@spatial_info[[spat_name]]@overlaps[[feat_name]] = spatVector
        }
      }
    }

  }




  ## 4. images
  if(verbose) wrap_msg('\n3. read Giotto image information \n')
  image_files = list.files(path = paste0(path_to_folder, '/Images'))
  if(length(image_files) != 0) {
    image_names = unique(gsub(image_files, pattern = '_spatRaster.*', replacement = ''))
    for(image_i in 1:length(image_names)) {
      image_name = image_names[image_i]
      if(verbose) image_name
      new_path = paste0(path_to_folder, '/Images','/', image_name,'_spatRaster')
      spatRaster = terra::rast(x = new_path)
      gobject@largeImages[[image_name]]@raster_object = spatRaster
    }
  }

  if(isTRUE(reconnect_giottoImage)) {
    if(list_images(gobject)[img_type == 'image', .N] > 0) {
      gobject = reconnectGiottoImage(gobject, reconnect_type = 'image')
    }
  }


  ## 5. Update python path
  identified_python_path = set_giotto_python_path(python_path = python_path)
  gobject = changeGiottoInstructions(gobject = gobject,
                                     params = c('python_path'),
                                     new_values = c(identified_python_path))

  ## 6. overallocate for data.tables
  # (data.tables when read from disk have a truelength of 0)
  gobject = giotto_alloc_dt_slots(gobject)


  return(gobject)

}



#' @title get10Xmatrix
#' @name get10Xmatrix
#' @description This function creates an expression matrix from a 10X structured folder
#' @param path_to_data path to the 10X folder
#' @param gene_column_index which column from the features or genes .tsv file to use for row ids
#' @param remove_zero_rows removes rows with sum equal to zero
#' @param split_by_type split into multiple matrices based on 3rd column of features.tsv(.gz)
#' @return sparse expression matrix from 10X
#' @details A typical 10X folder is named raw_feature_bc_matrix or filtered_feature_bc_matrix and it has 3 files:
#' \itemize{
#'   \item{barcodes.tsv(.gz)}
#'   \item{features.tsv(.gz) or genes.tsv(.gz)}
#'   \item{matrix.mtx(.gz)}
#' }
#' By default the first column of the features or genes .tsv file will be used, however if multiple
#' annotations are provided (e.g. ensembl gene ids and gene symbols) the user can select another column.
#' @export
get10Xmatrix = function(path_to_data,
                        gene_column_index = 1,
                        remove_zero_rows = TRUE,
                        split_by_type = TRUE) {

  # data.table variables
  total = gene_symbol = gene_id = gene_id_num = cell_id = cell_id_num = sort_gene_id_num = NULL

  # data directory
  files_10X = list.files(path_to_data)

  # get barcodes and create vector
  barcodes_file = grep(files_10X, pattern = 'barcodes', value = T)
  barcodesDT = data.table::fread(input = paste0(path_to_data,'/',barcodes_file), header = F)
  barcodes_vec = barcodesDT$V1
  names(barcodes_vec) = 1:nrow(barcodesDT)

  # get features and create vector
  features_file = grep(files_10X, pattern = 'features|genes', value = T)
  featuresDT = data.table::fread(input = paste0(path_to_data,'/',features_file), header = F)

  g_name = colnames(featuresDT)[gene_column_index]
  ## convert ensembl gene id to gene symbol ##
  ## TODO

  featuresDT[, total := .N, by = get(g_name)]
  featuresDT[, gene_symbol := ifelse(total > 1, paste0(get(g_name),'--',1:.N), get(g_name)), by = get(g_name)]
  features_vec = featuresDT$gene_symbol
  names(features_vec) = 1:nrow(featuresDT)

  # get matrix
  matrix_file = grep(files_10X, pattern = 'matrix', value = T)
  MMmatrix = Matrix::readMM(paste0(path_to_data,'/',matrix_file))
  rownames(MMmatrix) = features_vec
  colnames(MMmatrix) = barcodes_vec


  # Split by type of feature (features.tzv 3rd col)
  feat_classes = unique(featuresDT$V3)
  if(length(feat_classes) > 1) {
    result_list = list()

    for(fclass in feat_classes) {
      result_list[[fclass]] = MMmatrix[featuresDT$V3 == fclass, ]
    }

    if(isTRUE(remove_zero_rows)) {
      result_list = lapply(result_list, function(MMmatrix) {
        rowsums_result = rowSums_flex(MMmatrix)
        rowsums_bool = rowsums_result != 0
        MMmatrix = MMmatrix[rowsums_bool, ]
      })
    }

    return(result_list)
  } else {

    if(remove_zero_rows == TRUE) {
      rowsums_result = rowSums_flex(MMmatrix)
      rowsums_bool = rowsums_result != 0
      MMmatrix = MMmatrix[rowsums_bool, ]
    }

    return(MMmatrix)

  }

}



#' @title get10XmatrixOLD
#' @name get10XmatrixOLD
#' @description This function creates an expression matrix from a 10X structured folder
#' @param path_to_data path to the 10X folder
#' @param gene_column_index which column from the features or genes .tsv file to use for row ids
#' @return sparse expression matrix from 10X
#' @details A typical 10X folder is named raw_feature_bc_matrix or filtered_feature_bc_matrix and it has 3 files:
#' \itemize{
#'   \item{barcodes.tsv(.gz)}
#'   \item{features.tsv(.gz) or genes.tsv(.gz)}
#'   \item{matrix.mtx(.gz)}
#' }
#' By default the first column of the features or genes .tsv file will be used, however if multiple
#' annotations are provided (e.g. ensembl gene ids and gene symbols) the user can select another column.
get10XmatrixOLD = function(path_to_data, gene_column_index = 1) {

  # data.table variables
  total = gene_symbol = gene_id = gene_id_num = cell_id = cell_id_num = sort_gene_id_num = NULL

  # data directory
  files_10X = list.files(path_to_data)

  # get barcodes and create vector
  barcodes_file = grep(files_10X, pattern = 'barcodes', value = T)
  barcodesDT = data.table::fread(input = paste0(path_to_data,'/',barcodes_file), header = F)
  barcodes_vec = barcodesDT$V1
  names(barcodes_vec) = 1:nrow(barcodesDT)

  # get features and create vector
  features_file = grep(files_10X, pattern = 'features|genes', value = T)
  featuresDT = data.table::fread(input = paste0(path_to_data,'/',features_file), header = F)

  g_name = colnames(featuresDT)[gene_column_index]
  ## convert ensembl gene id to gene symbol ##
  ## TODO

  featuresDT[, total := .N, by = get(g_name)]
  featuresDT[, gene_symbol := ifelse(total > 1, paste0(get(g_name),'--',1:.N), get(g_name)), by = get(g_name)]
  features_vec = featuresDT$gene_symbol
  names(features_vec) = 1:nrow(featuresDT)

  # get matrix
  matrix_file = grep(files_10X, pattern = 'matrix', value = T)
  matrixDT = data.table::fread(input = paste0(path_to_data,'/',matrix_file), header = F, skip = 3)
  colnames(matrixDT) = c('gene_id_num', 'cell_id_num', 'umi')

  # extend matrixDT with missing cell IDs
  all_matrix_cell_ids = unique(matrixDT$cell_id_num)
  missing_barcodes_cell_ids = as.integer(names(barcodes_vec)[!names(barcodes_vec) %in% all_matrix_cell_ids])
  length_missing = length(missing_barcodes_cell_ids)

  if(length_missing > 0) {
    missing_matrixDT = data.table(gene_id_num = rep(1, length_missing),
                                  cell_id_num = missing_barcodes_cell_ids,
                                  umi = rep(0, length_missing))
    matrixDT = rbind(matrixDT, missing_matrixDT)
  }


  # convert barcodes and features
  matrixDT[, gene_id := features_vec[gene_id_num]]
  matrixDT[, cell_id := barcodes_vec[cell_id_num]]

  # make sure that gene id are consecutive
  sort_gene_id_vec = 1:length(unique(matrixDT$gene_id))
  names(sort_gene_id_vec) = unique(matrixDT$gene_id)
  matrixDT[, sort_gene_id_num := sort_gene_id_vec[gene_id]]

  sparsemat = Matrix::sparseMatrix(i = matrixDT$sort_gene_id_num,
                                   j = matrixDT$cell_id_num,
                                   x = matrixDT$umi,
                                   dimnames = list(unique(matrixDT$gene_id), unique(matrixDT$cell_id)))

  return(sparsemat)

}




#' @title get10Xmatrix_h5
#' @name get10Xmatrix_h5
#' @description This function creates an expression matrix from a 10X h5 file path
#' @param path_to_data path to the 10X .h5 file
#' @param gene_ids use gene symbols (default) or ensembl ids for the gene expression matrix
#' @inheritParams get10Xmatrix
#' @return (list of) sparse expression matrix from 10X
#' @details If the .h5 10x file has multiple classes of features (e.g. expression vs QC
#' probes) or modalities (e.g. RNA and protein), and \code{split_by_type} param is \code{TRUE},
#' multiple matrices will be returned
#' @export
get10Xmatrix_h5 = function(path_to_data,
                           gene_ids = c('symbols', 'ensembl'),
                           remove_zero_rows = TRUE,
                           split_by_type = TRUE) {

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

    barcodes = h5[[paste0(root, "/barcodes")]][]
    feature_tag_keys = h5[[paste0(root, "/features/_all_tag_keys")]][]
    feature_types = h5[[paste0(root, "/features/feature_type")]][]
    genome = unique(h5[[paste0(root, "/features/genome")]][])
    feature_id = h5[[paste0(root, "/features/id")]][]
    feature_names = h5[[paste0(root, "/features/name")]][]
    indices = h5[[paste0(root, "/indices")]][]
    indptr = h5[[paste0(root, "/indptr")]][]
    data_shape = h5[[paste0(root, "/shape")]][]

    # create a feature data.table
    features_dt = data.table::data.table(
      'id' = feature_id,
      'name' = feature_names,
      'feature_type' = feature_types,
      'genome' = genome
    )

  },
  finally = {
    h5$close_all()
  })

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

  # multiple feature classes (e.g. gene vs diagnostic or even modalities?)
  if(isTRUE(split_by_type)) {
    result_list = list()

    for(fclass in unique(feature_types)) {

      result_list[[fclass]] = sparsemat[features_dt$feature_type == fclass, ]

      # change names to gene symbols if it's expression
      if(fclass == 'Gene Expression' & gene_ids == 'symbols') {

        conv_vector = features_dt$uniq_name
        names(conv_vector) = features_dt$id

        current_names = rownames(result_list[[fclass]])
        new_names = conv_vector[current_names]
        rownames(result_list[[fclass]]) = new_names
      }
    }

    if(isTRUE(remove_zero_rows)) {
      result_list = lapply(result_list, function(sparsemat) {
        rowsums_result = rowSums_flex(sparsemat)
        rowsums_bool = rowsums_result != 0
        sparsemat = sparsemat[rowsums_bool, ]
      })
    }

    return(result_list)

  } else {

    if(isTRUE(remove_zero_rows)) {
      rowsums_result = rowSums_flex(sparsemat)
      rowsums_bool = rowsums_result != 0
      sparsemat = sparsemat[rowsums_bool, ]
    }

    return(list('Gene Expression' = sparsemat))
  }
}



#' @title convertEnsemblToGeneSymbol
#' @name convertEnsemblToGeneSymbol
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





#' @title readPolygonFilesVizgenHDF5
#' @name readPolygonFilesVizgenHDF5_old
#' @description Read and create polygons for all cells, or for only selected FOVs.
#' @param boundaries_path path to the cell_boundaries folder
#' @param fovs subset of fovs to use
#' @param custom_polygon_names a character vector to provide custom polygon names
#'   (optional)
#' @param polygon_feat_types a vector containing the polygon feature types
#' @param flip_x_axis flip x axis of polygon coordinates (multiply by -1)
#' @param flip_y_axis flip y axis of polygon coordinates (multiply by -1)
#' @param smooth_polygons smooth polygons (default = TRUE)
#' @param smooth_vertices number of vertices for smoothing
#' @param set_neg_to_zero set negative values to zero when smoothing
#' @param H5Fopen_flags see \code{\link[rhdf5]{H5Fopen}} for more details
#' @param cores cores to use
#' @param verbose be verbose
#' @seealso \code{\link{smoothGiottoPolygons}}
#' @details Set H5Fopen_flags to "H5F_ACC_RDONLY" if you encounter permission issues.
#' @export
readPolygonFilesVizgenHDF5_old = function(boundaries_path,
                                      fovs = NULL,
                                      polygon_feat_types = 0:6,
                                      custom_polygon_names = NULL,
                                      flip_x_axis = F,
                                      flip_y_axis = F,
                                      smooth_polygons = TRUE,
                                      smooth_vertices = 60,
                                      set_neg_to_zero = FALSE,
                                      H5Fopen_flags = "H5F_ACC_RDWR",
                                      cores = NA,
                                      verbose = TRUE) {

  # necessary pkgs
  package_check(pkg_name = 'rhdf5', repository = 'Bioc')

  cores = determine_cores(cores)

  # data.table vars
  x = y = cell_id = file_id = my_id = NULL

  # prepare poly feat names
  poly_feat_names = paste0('z', polygon_feat_types)
  poly_feat_indexes = paste0('zIndex_', polygon_feat_types)

  # provide your own custom names
  if(!is.null(custom_polygon_names)) {

    if(!is.character(custom_polygon_names)) {
      stop(wrap_txt('If custom_polygon_names are provided, it needs to be a character vector'))
    }

    if(length(custom_polygon_names) != length(poly_feat_names)) {
      stop(wrap_txt('length of custom names need to be same as polygon_feat_types'))
    } else {
      poly_feat_names = custom_polygon_names
    }
  }

  if(isTRUE(verbose)) wrap_msg('Reading from:', boundaries_path)
  # list all files in the folder
  hdf5_boundary_list = list.files(full.names = TRUE, boundaries_path)
  # only load subset of files if fov is given
  if(!is.null(fovs)) {

    selected_hdf5s = paste0('feature_data_', fovs, '.hdf5')
    selected_hdf5s_concatenated = paste0(selected_hdf5s, collapse = '|')
    hdf5_boundary_selected_list = grep(selected_hdf5s_concatenated, x = hdf5_boundary_list, value = TRUE)

  } else {
    hdf5_boundary_selected_list = hdf5_boundary_list
  }

  if(isTRUE(verbose)) wrap_msg('finished listing .hdf5 files
                               start extracting .hdf5 information')

  # open selected polygon files
  hdf5_list_length = length(hdf5_boundary_selected_list)

  # append data from all FOVs to single list
  init = proc.time()
  progressr::with_progress({
    pb = progressr::progressor(along = hdf5_boundary_selected_list)
    read_list = lapply_flex(seq_along(hdf5_boundary_selected_list),
                            cores = cores,
                            future.packages = c('rhdf5', 'Rhdf5lib'),
                            function(bound_i) {

      # get feature data
      read_file = rhdf5::H5Fopen(hdf5_boundary_selected_list[[bound_i]][[1]], flags = H5Fopen_flags)
      fov_info = read_file$featuredata

      # update progress
      print(basename(hdf5_boundary_selected_list[[bound_i]]))
      elapsed = (proc.time() - init)[[3L]]
      step_time = elapsed/bound_i
      est = (hdf5_list_length * step_time) - elapsed
      pb(message = c('// E:', time_format(elapsed), '| R:', time_format(est)))
      rhdf5::H5Fclose(read_file)
      return(fov_info)
    })
  })

  # # combine to FOV data single list
  read_list = Reduce('append', read_list)
  cell_names = names(read_list)


  # extract values for each z index and cell from read_list
  result_list = lapply_flex(seq_along(poly_feat_indexes), cores = cores, function(z_i) {
    # Reduce('append', lapply_flex(seq_along(read_list), cores = cores, function(bound_i) {
    #   cell_names = names(read_list[[bound_i]])
      # lapply_flex(seq_along(read_list[[bound_i]]), cores = cores, function(cell_i) {
      lapply_flex(seq_along(read_list), cores = cores, function(cell_i) {
        # singlearray = read_list[[bound_i]][[cell_i]][[poly_feat_indexes[z_i]]]$p_0$coordinates
        singlearray = read_list[[cell_i]][[poly_feat_indexes[z_i]]]$p_0$coordinates
        cell_name = cell_names[[cell_i]]
        if(!is.null(singlearray)) {
          singlearraydt = data.table::as.data.table(t_flex(as.matrix(singlearray[,,1])))
          data.table::setnames(singlearraydt, old = c('V1', 'V2'), new = c('x', 'y'))
          if(flip_x_axis) singlearraydt[, x := -1 * x]
          if(flip_y_axis) singlearraydt[, y := -1 * y]

          # singlearraydt[, file_id := paste0('file', bound_i)]
          singlearraydt[, cell_id := cell_name]
          # singlearraydt[, my_id := paste0('cell', cell_i)]
        }
      })
    # }))
  })
  result_list_rbind = lapply_flex(seq_along(result_list), cores = cores, function(z_i) {
    data.table::rbindlist(result_list[[z_i]])
  })



  if(isTRUE(verbose)) wrap_msg('finished extracting .hdf5 files
                               start creating polygons')


  # create Giotto polygons and add them to gobject
  # smooth_cell_polygons_list = list()
  progressr::with_progress({
    pb = progressr::progressor(along = result_list_rbind)
    smooth_cell_polygons_list = lapply_flex(seq_along(result_list_rbind), cores = cores, function(i) {
      dfr_subset = result_list_rbind[[i]][,.(x, y, cell_id)]
      cell_polygons = createGiottoPolygonsFromDfr(segmdfr = dfr_subset,
                                                  name = poly_feat_names[i],
                                                  verbose = verbose)

      pb(message = poly_feat_names[i])

      if(smooth_polygons == TRUE) {
        return(smoothGiottoPolygons(cell_polygons,
                                    vertices = smooth_vertices,
                                    set_neg_to_zero = set_neg_to_zero))
      } else {
        return(cell_polygons)
      }
    })
  })


  # TODO: add spatial centroids
  # needs to happen after smoothing to be correct

  return(smooth_cell_polygons_list)

}




#' @title readPolygonFilesVizgenHDF5
#' @name readPolygonFilesVizgenHDF5
#' @description Read and create polygons for all cells, or for only selected FOVs.
#' @param boundaries_path path to the cell_boundaries folder
#' @param fovs subset of fovs to use
#' @param z_indices z indices of polygons to use
#' @param segm_to_use segmentation results to use (usually = 1. Depends on if
#' alternative segmentations were generated)
#' @param custom_polygon_names a character vector to provide custom polygon names
#'   (optional)
#' @param polygon_feat_types deprecated. Use \code{z_indices}
#' @param flip_x_axis flip x axis of polygon coordinates (multiply by -1)
#' @param flip_y_axis flip y axis of polygon coordinates (multiply by -1)
#' @param smooth_polygons smooth polygons (default = TRUE)
#' @param smooth_vertices number of vertices for smoothing
#' @param set_neg_to_zero set negative values to zero when smoothing
#' @param calc_centroids calculate centroids (default = FALSE)
#' @param H5Fopen_flags see \code{\link[rhdf5]{H5Fopen}} for more details
#' @param cores cores to use
#' @param create_gpoly_parallel (default = TRUE) Whether to run gpoly creation in
#' parallel
#' @param create_gpoly_bin (Optional, default = FALSE) Parallelization option.
#' Accepts integer values as an binning size when generating giottoPolygon objects
#' @param verbose be verbose
#' @seealso \code{\link{smoothGiottoPolygons}}
#' @details Set H5Fopen_flags to "H5F_ACC_RDONLY" if you encounter permission issues.
#' @export
readPolygonFilesVizgenHDF5 = function(boundaries_path,
                                      fovs = NULL,
                                      z_indices = 1L:7L,
                                      segm_to_use = 1L,
                                      custom_polygon_names = NULL,
                                      flip_x_axis = FALSE,
                                      flip_y_axis = TRUE,
                                      calc_centroids = FALSE,
                                      smooth_polygons = TRUE,
                                      smooth_vertices = 60L,
                                      set_neg_to_zero = FALSE,
                                      H5Fopen_flags = "H5F_ACC_RDWR",
                                      cores = NA,
                                      create_gpoly_parallel = TRUE,
                                      create_gpoly_bin = FALSE,
                                      verbose = TRUE,
                                      polygon_feat_types = NULL) {

  # necessary pkgs
  package_check(pkg_name = 'rhdf5', repository = 'Bioc')

  cores = determine_cores(cores)

  # deprecation
  if(!is.null(polygon_feat_types)) {
    warning('polygon_feat_types is deprecated.\n Use z_indices instead')
    z_indices = polygon_feat_types + 1L
  }

  segm_to_use = paste0('p_', (segm_to_use - 1L))

  # data.table vars
  x = y = z = cell_id = poly_ID = file_id = my_id = NULL

  # provide your own custom names
  if(!is.null(custom_polygon_names)) {

    if(!is.character(custom_polygon_names)) {
      stop(wrap_txt('If custom_polygon_names are provided, it needs to be a character vector'))
    }

    if(length(custom_polygon_names) != length(z_indices)) {
      stop(wrap_txt('length of custom names need to be same as z_indices'))
    }
  }

  if(isTRUE(verbose)) wrap_msg('Reading from:', boundaries_path)
  # list all files in the folder
  hdf5_boundary_list = list.files(full.names = TRUE, boundaries_path)
  # only load subset of files if fov is given
  if(!is.null(fovs)) {

    selected_hdf5s = paste0('feature_data_', fovs, '.hdf5')
    selected_hdf5s_concatenated = paste0(selected_hdf5s, collapse = '|')
    hdf5_boundary_selected_list = grep(selected_hdf5s_concatenated, x = hdf5_boundary_list, value = TRUE)

  } else {
    hdf5_boundary_selected_list = hdf5_boundary_list
  }

  if(isTRUE(verbose)) wrap_msg('finished listing .hdf5 files
                               start extracting .hdf5 information')

  # open selected polygon files
  hdf5_list_length = length(hdf5_boundary_selected_list)

  # append data from all FOVs to single list
  init = Sys.time()
  progressr::with_progress({
    pb = progressr::progressor(length(hdf5_boundary_selected_list)/5)
    read_list = lapply_flex(seq_along(hdf5_boundary_selected_list),
                            future.packages = c('rhdf5', 'Rhdf5lib'),
                            function(init, z_indices, segm_to_use, bound_i) {
                              read_file = h5read_vizgen(h5File = hdf5_boundary_selected_list[[bound_i]][[1]],
                                                        z_indices = z_indices,
                                                        segm_to_use = segm_to_use,
                                                        H5Fopen_flags = H5Fopen_flags)

                              # update progress
                              print(basename(hdf5_boundary_selected_list[[bound_i]]))
                              if(bound_i %% 5 == 0) {
                                pb()
                              }

                              return(read_file)
                            },
                            cores = cores,
                            init = init,
                            z_indices = z_indices,
                            segm_to_use = segm_to_use)
  })

  # combine to FOV data single list
  read_DT = data.table::rbindlist(read_list)

  # perform any necessary flips
  if(flip_x_axis) read_DT[, x := -1 * x]
  if(flip_y_axis) read_DT[, y := -1 * y]

  # separate polygons by z index
  zvals = read_DT[, unique(z)]
  z_names = paste0('z', zvals)
  z_read_DT = lapply(seq_along(zvals), function(z_idx) {
    read_DT[z == zvals[z_idx],]
  })
  names(z_read_DT) = z_names
  if(!is.null(custom_polygon_names)) poly_names = custom_polygon_names
  else poly_names = z_names

  if(isTRUE(verbose)) wrap_msg('finished extracting .hdf5 files
                               start creating polygons')

  # create Giotto polygons and add them to gobject

  # **** sequential method ****
  if(!isTRUE(create_gpoly_parallel)) {
    progressr::with_progress({
      pb = progressr::progressor(along = z_read_DT)
      smooth_cell_polygons_list = lapply(seq_along(z_read_DT), function(i) {
        dfr_subset = z_read_DT[[i]][,.(x, y, cell_id)]
        data.table::setnames(dfr_subset, old = 'cell_id', new = 'poly_ID')
        cell_polygons = createGiottoPolygonsFromDfr(segmdfr = dfr_subset,
                                                    name = poly_names[i],
                                                    calc_centroids = FALSE,
                                                    skip_eval_dfr = TRUE,
                                                    copy_dt = FALSE,
                                                    verbose = verbose)
        if(isTRUE(smooth_polygons)) cell_polygons = smoothGiottoPolygons(gpolygon = cell_polygons,
                                                                         vertices = smooth_vertices,
                                                                         k = 3L,
                                                                         set_neg_to_zero = set_neg_to_zero)
        if(isTRUE(calc_centroids)) cell_polygons = calculate_centroids_polygons(gpolygon = cell_polygons,
                                                                                append_gpolygon = TRUE)
        pb(message = c(poly_names[i], ' (', i, '/', length(z_read_DT), ')'))
        return(cell_polygons)
      })
    })
    return(smooth_cell_polygons_list)
  }


  # **** parallel methods ****
  # no binning
  if(!is.numeric(create_gpoly_bin)) {

    progressr::with_progress({
      pb = progressr::progressor(along = z_read_DT)
      smooth_cell_polygons_list = lapply_flex(
        seq_along(z_read_DT),
        future.packages = c('terra', 'stats', 'data.table'),
        function(i) {
          dfr_subset = z_read_DT[[i]][,.(x, y, cell_id)]
          data.table::setnames(dfr_subset, old = 'cell_id', new = 'poly_ID')
          cell_polygons = gpoly_from_dfr_smoothed_wrapped(
            segmdfr = dfr_subset,
            name = poly_names[i],
            skip_eval_dfr = TRUE,
            copy_dt = FALSE,
            smooth_polygons = smooth_polygons,
            vertices = smooth_vertices,
            set_neg_to_zero = set_neg_to_zero,
            calc_centroids = calc_centroids,
            verbose = verbose
          )

          pb(message = c(poly_names[i], ' (', i, '/', length(z_read_DT), ')'))
          return(cell_polygons)
        }
      )
    })

    # unwrap results
    smooth_cell_polygons_list = lapply(smooth_cell_polygons_list, function(x) {
      slot(x, 'spatVector') = terra::vect(slot(x, 'spatVector'))
      if(isTRUE(calc_centroids)) {
        slot(x, 'spatVectorCentroids') = terra::vect(slot(x, 'spatVectorCentroids'))
      }
      return(x)
    })

  } else {
    # with binning

    dfr_subset = lapply(z_read_DT, function(bin, DT) {

      DT = DT[,.(x, y, cell_id)]
      data.table::setnames(DT, old = 'cell_id', new = 'poly_ID')
      pid = DT[, unique(poly_ID)]

      bin_pid = data.table::data.table(
        'poly_ID' = pid,
        'bin_ID' = as.numeric(
          cut(x = seq_along(pid),
              breaks = ceiling(length(pid)/bin))
        )
      )
      DT = data.table::merge.data.table(DT, bin_pid, by = 'poly_ID', all.x = TRUE)
      DT = split(DT, DT$bin_ID)

    }, bin = create_gpoly_bin)

    bin_steps = sum(unlist(lapply(dfr_subset, length)))

    progressr::with_progress({
      pb = progressr::progressor(steps = bin_steps)
      smooth_cell_polygons_list = lapply( # sequential across z index
        seq_along(dfr_subset),
        function(i) {
          lapply_flex(                    # parallelize across bins
            dfr_subset[[i]],
            future.packages = c('terra', 'stats', 'data.table'),
            function(bin_DT) {
              cell_polygons = gpoly_from_dfr_smoothed_wrapped(
                segmdfr = bin_DT,
                name = poly_names[i],
                skip_eval_dfr = TRUE,
                copy_dt = FALSE,
                smooth_polygons = smooth_polygons,
                vertices = smooth_vertices,
                set_neg_to_zero = set_neg_to_zero,
                calc_centroids = calc_centroids,
                verbose = verbose
              )

              pb(message = c(poly_names[i], ' (', i, '/', length(dfr_subset), ')'))
              return(cell_polygons)
            }
          )
        }
      )
    })

    # unwrap results
    smooth_cell_polygons_list = lapply(seq_along(smooth_cell_polygons_list), function(i) {
      p_list = lapply(smooth_cell_polygons_list[[i]], function(x) {
        slot(x, 'spatVector') = terra::vect(slot(x, 'spatVector'))
        if(isTRUE(calc_centroids)) {
          slot(x, 'spatVectorCentroids') = terra::vect(slot(x, 'spatVectorCentroids'))
        }
        return(x)
      })
      # rbind results
      names(p_list) = NULL
      return(do.call('rbind', p_list))
    })

  }


  return(smooth_cell_polygons_list)

}





#' @title readPolygonFilesVizgen
#' @name readPolygonFilesVizgen
#' @description Read selected polygon files for the FOVs present in the Giotto
#' object and add the smoothed polygons to the object
#' @param gobject giotto object
#' @param boundaries_path path to the cell_boundaries folder
#' @param fovs selected fovs, if NULL select all fovs within Giotto object
#' @param polygon_feat_types a vector containing the polygon feature types
#' @param flip_x_axis flip x axis of polygon coordinates (multiply by -1)
#' @param flip_y_axis flip y axis of polygon coordinates (multiply by -1)
#' @param smooth_polygons smooth polygons (default = TRUE)
#' @param smooth_vertices number of vertices for smoothing
#' @param set_neg_to_zero set negative values to zero when smoothing
#' @param return_gobject return giotto object
#' @param verbose be verbose
#' @seealso \code{\link{smoothGiottoPolygons}}
#' @export
readPolygonFilesVizgen = function(gobject,
                                  boundaries_path,
                                  fovs = NULL,
                                  polygon_feat_types = 0:6,
                                  flip_x_axis = F,
                                  flip_y_axis = F,
                                  smooth_polygons = TRUE,
                                  smooth_vertices = 60,
                                  set_neg_to_zero = FALSE,
                                  return_gobject = TRUE,
                                  verbose = TRUE) {

  # define names
  poly_feat_names = paste0('z', polygon_feat_types)
  poly_feat_indexes = paste0('zIndex_', polygon_feat_types)

  # select FOVs present in the subset
  if(is.null(fovs)) {
    subset_metadata = pDataDT(gobject)
    fovs = unique(subset_metadata$fov)
  }



  smooth_cell_polygons_list = readPolygonFilesVizgenHDF5(boundaries_path = boundaries_path,
                                                         fovs = fovs,
                                                         polygon_feat_types = polygon_feat_types,
                                                         flip_x_axis = flip_x_axis,
                                                         flip_y_axis = flip_y_axis,
                                                         smooth_polygons = smooth_polygons,
                                                         smooth_vertices = smooth_vertices,
                                                         set_neg_to_zero = set_neg_to_zero,
                                                         verbose = verbose)


  if(return_gobject) {
    # add cell polygons to Giotto object
    names(smooth_cell_polygons_list) = poly_feat_names
    gobject = addGiottoPolygons(gobject = gobject,
                                gpolygons = smooth_cell_polygons_list)
    return(gobject)
  } else {
    return(smooth_cell_polygons_list)
  }


}




#' @describeIn readPolygonFilesVizgen (internal) Optimized .hdf5 reading for vizgen
#' merscope output. Returns a data.table of xyz coords and cell_id
#' @keywords internal
h5read_vizgen = function(h5File,
                         z_indices = 1L:7L,
                         segm_to_use = 'p_0',
                         H5Fopen_flags = "H5F_ACC_RDWR") {

  # data.table vars
  group = name = cell = z_name = NULL

  h5_ls = data.table::setDT(rhdf5::h5ls(h5File, recursive = 5, datasetinfo = FALSE))
  cell_names = as.character(h5_ls[group == '/featuredata', name])
  z_names = h5_ls[grep('zIndex', name), unique(name)]

  dset_names = h5_ls[otype == 'H5I_DATASET' & name == 'coordinates',]
  # subset by segm_to_use
  dset_names = dset_names[grep(segm_to_use, group),]
  # tag cellnames
  dset_names[, cell := gsub(pattern = '/featuredata/|/zIndex.*$', replacement = '', x = group)]
  # tag z_names
  dset_names[, z_name := gsub(pattern = '^.*/(zIndex_\\d*).*$', replacement = '\\1', x = group)]
  # subset by z_indices
  dset_names = dset_names[z_name %in% z_names[z_indices],]
  # create full file location
  dset_names[, d_name := paste0(group, '/', name)]

  fid = rhdf5::H5Fopen(h5File, flags = H5Fopen_flags)
  dapl = rhdf5::H5Pcreate('H5P_DATASET_ACCESS')

  contents = lapply(cell_names, function(fid, dapl, cell_name) {

    zvals = .h5_read_bare(file = fid,
                          name = paste0(c('/featuredata', cell_name, 'z_coordinates'), collapse = '/'),
                          dapl = dapl)
    names(zvals) = z_names

    # subset to datasets related to cell
    cell_dsets = dset_names[cell == cell_name,]

    cell_data = lapply(seq(nrow(cell_dsets)), function(fid, dapl, zvals, d_i) {

      res = .h5_read_bare(file = fid, name = cell_dsets[d_i, d_name], dapl = dapl)
      res = t_flex(res[,,1L])
      res = cbind(res, zvals[cell_dsets[d_i, z_name]])
      colnames(res) = c('x', 'y', 'z')
      res

    }, fid = fid, dapl = dapl, zvals = zvals)
    cell_data = data.table::as.data.table(do.call('rbind', cell_data))
    cell_data[, cell_id := cell_name]
    cell_data

  }, fid = fid, dapl = dapl)

  rhdf5::H5Pclose(dapl)
  rhdf5::H5Fclose(fid)
  contents = data.table::rbindlist(contents)
  contents

}



#' @name .h5_read_bare
#' @title Read dataset from opened HDF5 with C functions
#' @param file opened HDF5 file
#' @param name dataset name within
#' @param dapl HDF5 property list (H5Pcreate('H5P_DATASET_ACCESS'))
#' @keywords internal
.h5_read_bare = function(file, name = "", dapl) {
  did = .Call("_H5Dopen", file@ID, name, dapl@ID, PACKAGE = "rhdf5")
  res = .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, FALSE,
              PACKAGE = "rhdf5")
  invisible(.Call("_H5Dclose", did, PACKAGE = "rhdf5"))
  res
}





#' @title getGEFtxCoords
#' @name getGEFtxCoords
#' @description Converts .gef file (output stereo-seq pipeline) into
#' transcript with coordinates
#' @param gef_file path to .gef file
#' @param bin_size bin size to select from .gef file
#' @export
getGEFtxCoords = function(gef_file,
                          bin_size = 'bin100') {

  # data.table vars
  genes = NULL

  # package check
  package_check(pkg_name = 'rhdf5', repository = 'Bioc')
  if(!file.exists(gef_file)) stop('File path to .gef file does not exist')

  # step 1: read expression and gene data from gef file
  geneExpData = rhdf5::h5read(file = gef_file, name = 'geneExp')
  exprDT = data.table::as.data.table(geneExpData[[bin_size]][['expression']])
  geneDT = data.table::as.data.table(geneExpData[[bin_size]][['gene']])

  # step 2: combine gene information from the geneDT to the exprDT
  exprDT[, genes := rep(x = geneDT$gene, geneDT$count)]

  return(exprDT)

}



#' @title Fread rows based on column matches
#' @name fread_colmatch
#' @param file path to file to load
#' @param col name of col to match from
#' @param sep grep term to match as column delimiters within the file
#' @param values_to_match values in \code{col} to match given as a vector
#' @param verbose whether to print the grep command
#' @param ... additional parameters to pass to \code{\link[data.table]{fread}}
#' @keywords internal
fread_colmatch = function(file,
                          col,
                          sep = NULL,
                          values_to_match,
                          verbose = FALSE,
                          ...) {

  # get colnames
  col_names = colnames(data.table::fread(file, nrows = 1L))
  col_num = which(col_names == col)

  # try to guess col separating char if not given
  if(is.null(sep)) {
    filename = basename(file)
    if(grepl(pattern = '.csv', x = filename)) {
      sep = '.*,'
    } else if(grepl(pattern = '.tsv', x = filename)) {
      sep = '.*\t'
    } else {
      stop('sep param cannot be guessed')
    }
  }

  # create grep search
  pattern = paste(values_to_match, collapse = '|')
  gpat = paste0('\'', strrep(x = sep, times = col_num - 1), '(', pattern, '),\' ')
  fread_cmd = paste0('grep -E ', gpat, file)
  if(isTRUE(verbose)) print(fread_cmd)

  file_DT = data.table::fread(cmd = fread_cmd, col.names = col_names, ...)
  return(file_DT)
}

