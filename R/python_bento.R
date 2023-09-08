#' @title Create bento adata object from gobject
#' @name createBentoAdata
#' @description Create bento adata object from gobject
#' @param gobject Giotto object
#' @return bento_adata bento adata object
#' @export
createBentoAdata <- function(gobject){
  # Transcripts
  transcripts_df <- as.data.frame(sf::st_as_sf(gobject@feat_info$rna@spatVector))
  coordinates_df <- sf::st_coordinates(transcripts_df['geometry'])
  t_df <- cbind(coordinates_df, transcripts_df[c("feat_ID")])
  colnames(t_df) <- c('x','y','gene')
  t_df <- t_df[t_df$cell != '-1',]  # Remove extracellular points, required by Bento

  # Cell shapes
  # TODO: Add batch information based on?
  cell_poly <- spatVector_to_dt(gobject@spatial_info$cell@spatVector)
  cell_poly <- data.frame(cell_id = cell_poly$poly_ID, x = cell_poly$x, y = cell_poly$y, batch = 0L)

  # Nuclei shapes
  # TODO: Add batch information based on?
  nucleus_poly <- spatVector_to_dt(gobject@spatial_info$nucleus@spatVector)
  nucleus_poly <- data.frame(cell_id = nucleus_poly$poly_ID, x = nucleus_poly$x, y = nucleus_poly$y, batch = 0L)

  # Create AnnData object
  bento_adata <- Giotto::create_AnnData(trainscripts = t_df, cell_shape = cell_poly, nucleus_shape = nucleus_poly)
  
  return(bento_adata)
}
