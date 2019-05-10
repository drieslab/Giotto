

#' @title getDistinctColors
#' @description Returns a number of distint colors based on the RGB scale
#' @param n number of colors wanted
#' @return number of distinct colors
getDistinctColors <- function(n) {
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unique (base::unlist(base::mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, base::rownames(qual_col_pals))));
  stopifnot (n <= base::length(col_vector));
  xxx <- grDevices::col2rgb(col_vector);
  dist_mat <- base::as.matrix(dist(t(xxx)));
  base::diag(dist_mat) <- 1e10;
  while (base::length(col_vector) > n) {
    minv <- base::apply (dist_mat,1,function(x)min(x));
    idx <- base::which(minv==min(minv))[1];
    dist_mat <- dist_mat[-idx, -idx];
    col_vector <- col_vector[-idx]
  }
  col_vector
}



#' @title mygini_fun
#' @description calculate gini coefficient
#' @return gini coefficient
mygini_fun <- function(x, weights = rep(1, length = length(x))) {
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox]/sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu/nu[n]
  sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

#' @title extended_gini_fun
#' @description calculate gini coefficient on a minimum length vector
#' @return gini coefficient
extended_gini_fun <- function(x, weights = rep(1, length = length(x)), minimum_length = 16) {

  if(length(x) < minimum_length) {
    difference = minimum_length - length(x)
    x = c(x,rep(0,difference))
  }

  result <- mygini_fun(x = x, weights = weights)
  return(result)
}


#' @title stitchFieldCoordinates
#' @description Helper function to stitch field coordinates together to form one complete picture
#' @param location_file location dataframe with X and Y coordinates
#' @param offset_file dataframe that describes to offset for each field (see details)
#' @param cumulate_offset_x (boolean) Do the x-axis offset values need to be cumulated?
#' @param cumulate_offset_y (boolean) Do the y-axis offset values need to be cumulated?
#' @param field_col column that indicates the field within the location_file
#' @param X_coord_col column that indicates the x coordinates
#' @param Y_coord_col column that indicates the x coordinates
#' @param reverse_final_x (boolean) Do the final x coordinates need to be reversed?
#' @param reverse_final_y (boolean) Do the final y coordinates need to be reversed?
#' @return Updated location dataframe with new X ['X_final'] and Y ['Y_final'] coordinates
#' @details Describe how stitching works.
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

  copy_loc_file[, c('X_final', 'Y_final') := list(new_x_coord, new_y_coord)]

  return(copy_loc_file)
}

