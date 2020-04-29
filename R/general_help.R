

#' @title getDistinctColors
#' @description Returns a number of distint colors based on the RGB scale
#' @param n number of colors wanted
#' @return number of distinct colors
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
    dist_mat <- as.matrix(dist(t(xxx)));
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


#' @title dt_to_matrix
#' @description converts data.table to matrix
#' @examples
#'     dt_to_matrix(x)
dt_to_matrix <- function(x) {
  rownames = as.character(x[[1]])
  mat = as(as.matrix(x[,-1]), 'Matrix')
  rownames(mat) = rownames
  return(mat)
}


#' @title mygini_fun
#' @description calculate gini coefficient
#' @return gini coefficient
mygini_fun <- function(x, weights = rep(1,length(x))) {

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
#' @return gini coefficient
extended_gini_fun <- function(x, weights = rep(1, length = length(x)), minimum_length = 16) {

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
#' @examples
#'     stitchFieldCoordinates(gobject)
stitchFieldCoordinates <- function(location_file,
                                   offset_file,
                                   cumulate_offset_x = F,
                                   cumulate_offset_y = F,
                                   field_col = 'Field of View',
                                   X_coord_col = 'X',
                                   Y_coord_col = 'Y',
                                   reverse_final_x = F,
                                   reverse_final_y = T) {


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
#' @details ...
#' @export
#' @examples
#'     stitchTileCoordinates(gobject)
stitchTileCoordinates <- function (location_file, Xtilespan, Ytilespan) {
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



#' @title get10Xmatrix
#' @description This function creates an expression matrix from a 10X structured folder
#' @param path_to_data path to the 10X folder
#' @param gene_column_index which column from the features or genes .tsv file to use for row ids
#' @return sparse expression matrix from 10X
#' @details A typical 10X folder is named raw_feature_bc_matrix or raw_feature_bc_matrix and tt has 3 files:
#' \itemize{
#'   \item{barcodes.tsv(.gz)}
#'   \item{features.tsv(.gz) or genes.tsv(.gz)}
#'   \item{matrix.mtx(.gz)}
#' }
#' By default the first column of the features or genes .tsv file will be used, however if multiple
#' annotations are provided (e.g. ensembl gene ids and gene symbols) the user can select another column.
#' @export
#' @examples
#'     get10Xmatrix(10Xmatrix)
get10Xmatrix = function(path_to_data, gene_column_index = 1) {
  
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
  #matrix_ab_mat = Giotto:::dt_to_matrix(matrix_ab)
  #matrix_ab_mat[is.na(matrix_ab_mat)] = 0

}


#' @title convertEnsemblToGeneSymbol
#' @description This function convert ensembl gene IDs from a matrix to official gene symbols
#' @param matrix an expression matrix with ensembl gene IDs as rownames
#' @param species species to use for gene symbol conversion
#' @return expression matrix with gene symbols as rownames
#' @details This function requires that the biomaRt library is installed
#' @export
#' @examples
#'     convertEnsemblToGeneSymbol(matrix)
convertEnsemblToGeneSymbol = function(matrix,
                                      species = c('mouse', 'human')) {


  if("biomaRt" %in% rownames(installed.packages()) == FALSE) {
    cat("\n package 'biomaRt' is not yet installed and is required for this function \n")
  }

  species = match.arg(species, choices = c('mouse', 'human'))

  if(species == 'mouse') {

    # ensembl IDs to change
    ensemblsIDS = rownames(matrix)

    # prepare ensembl database
    ensembl = useMart("ensembl",
                      dataset = "mmusculus_gene_ensembl")
    gene_names = getBM(attributes= c('mgi_symbol', 'ensembl_gene_id'),
                       filters = 'ensembl_gene_id',
                       values = ensemblsIDS,
                       mart = ensembl)
    gene_names_DT = as.data.table(gene_names)
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
    ensembl = useMart("ensembl",
                      dataset = "hsapiens_gene_ensembl")
    gene_names = getBM(attributes= c('hgnc_symbol', 'ensembl_gene_id'),
                       filters = 'ensembl_gene_id',
                       values = ensemblsIDS,
                       mart = ensembl)
    gene_names_DT = as.data.table(gene_names)
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



#' @title my_arowMeans
#' @description arithmic rowMeans that works for a single column
#' @examples
#'     my_arowMeans(x)
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
#' @examples
#'     my_growMeans(x)
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
#' @examples
#'     my_rowMeans(x)
my_rowMeans = function(x, method = c('arithmic', 'geometric'), offset = 0.1) {
  method = match.arg(method, c('arithmic', 'geometric'))
  if(method == 'arithmic') return(my_arowMeans(x))
  if(method == 'geometric') return(my_growMeans(x))
}

#' @title DT_removeNA
#' @name DT_removeNA
#' @description set NA values to 0 in a data.table object
DT_removeNA = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0]
  return(DT)
}

#' @title kmeans_binarize
#' @name kmeans_binarize
#' @description create binarized scores from a vector using kmeans
kmeans_binarize = function(x, nstart = 3, iter.max = 10) {

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

#' @title rank_binarize
#' @name rank_binarize
#' @description create binarized scores from a vector using arbitrary rank
rank_binarize = function(x, max_rank = 200) {

  sel_gene_rank = rank(-x, ties.method = 'average')

  sel_gene_bin = x
  sel_gene_bin[sel_gene_rank <= max_rank] = 1
  sel_gene_bin[sel_gene_rank > max_rank] = 0

  return(sel_gene_bin)

}


