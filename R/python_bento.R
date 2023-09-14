#' @title Create bento adata object from gobject
#' @name createBentoAdata
#' @description Create bento adata object from gobject
#' @param gobject Giotto object
#' @return bento_adata bento adata object
#' @export
createBentoAdata <- function(gobject){
  # Transcripts
  transcripts_df <- as.data.frame(sf::st_as_sf(gobject@feat_info$rna@spatVector))
  coordinates_df <- lapply(transcripts_df['geometry'], sf::st_coordinates)$geometry
  t_df <- as.data.frame(cbind(coordinates_df, transcripts_df[c("feat_ID")]))
  colnames(t_df) <- c('x','y','gene')

  # Cell shapes
  # TODO: Add batch information based on?
  cell_poly <- spatVector_to_dt(gobject@spatial_info$cell@spatVector)
  cell_poly <- data.frame(cell_id = cell_poly$poly_ID, x = cell_poly$x, y = cell_poly$y, batch = 0L)

  # Nuclei shapes
  # TODO: Add batch information based on?
  nucleus_poly <- spatVector_to_dt(gobject@spatial_info$nucleus@spatVector)
  nucleus_poly <- data.frame(cell_id = nucleus_poly$poly_ID, x = nucleus_poly$x, y = nucleus_poly$y, batch = 0L)

  # Create AnnData object
  g2bento_path <- system.file("python","g2bento.py",package="Giotto")
  reticulate::source_python(g2bento_path)
  bento_adata <- create_AnnData(trainscripts=t_df, cell_shape=cell_poly, nucleus_shape=nucleus_poly)
  
  return(bento_adata)
}
