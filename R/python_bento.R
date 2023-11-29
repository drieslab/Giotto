#' @title Create bento adata object from gobject
#' @name createBentoAdata
#' @description Create bento adata object from gobject
#' @param gobject Giotto object
#' @param env_to_use Python environment within which bento is installed.
#' If it is not already installed, the user
#' will be prompted to install `bento-tools`
#' DEFAULT: "giotto_env"
#' @return bento_adata bento adata object
#' @export
createBentoAdata <- function(gobject = NULL,
                             env_to_use = "giotto_env"){
  if(!c("giotto") %in% class(gobject)) stop(wrap_txt("Please provide a valid Giotto Object.", errWidth=TRUE))
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

  # Install bento-tools / Check python environment for bento-tools
  bento_installed = checkPythonPackage(github_package_url = "git+https://github.com/wwang-chcn/bento-tools.git",
                                       env_to_use = env_to_use)
  # Will crash downstream if installation unsuccessful/denied 
  # or if the package is not found.

  # Create AnnData object
  g2bento_path <- system.file("python","g2bento.py",package="Giotto")
  reticulate::source_python(g2bento_path)
  bento_adata <- create_AnnData(trainscripts=t_df, cell_shape=cell_poly, nucleus_shape=nucleus_poly)
  
  return(bento_adata)
}
