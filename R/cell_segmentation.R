
#' @title doCellSegmentation
#' @name doCellSegmentation
#' @description segment cells in Dapi image
#' @param raster_img raster image; Dapi image of cells
#' @param folder_path character; where to save the file
#' @param reduce_resolution numeric; the original Dapi image from Vizgen works
#' better in the Mesmer algorithm if its resolution is reduced 4 times.
#' @param overlapping_pixels numeric; the number of pixels to overlap when
#' calculating the rolling window
#' @param python_path specify specific path to python if required
#' @details
#' This function is a wrapper for the Mesmer algorithm implemented in python.
#' It segments the cells in tissue by applying a rolling window on the complete
#' image. It saves all the files in the specified folder with the coordinates
#' of the tile: sx (start x), ex (end x), sy, and ey.
#'
#' @export
doCellSegmentation <- function(raster_img,
                               folder_path,
                               reduce_resolution = 4,
                               overlapping_pixels = 50,
                               python_path = NULL){

  # python path (requires use of gobject)
  # if(is.null(python_path)) {
  #   python_path = readGiottoInstructions(gobject, param = "python_path")
  # }

  # prepare python path and segmentation script
  reticulate::use_python(required = T, python = python_path)
  python_segmentation_function = system.file("python",
                                             "python_segmentation.py",
                                             package = 'Giotto')
  reticulate::source_python(file = python_segmentation_function)

  # create mesmer app
  mesmer_app = python_create_mesmer_app()

  # get params for rolling window
  img_dim = dim(raster_img)
  xdim = img_dim[2]
  ydim = img_dim[1]

  mesmer_tile_dim = 512
  tile_dim = mesmer_tile_dim * reduce_resolution
  margin = overlapping_pixels
  window_size = tile_dim - margin

  nxwindow = xdim %/% window_size
  nywindow = ydim %/% window_size

  # sliding window
  start_x = 0
  end_x = start_x + tile_dim
  for (i in 1:nxwindow) {

    start_y = 0
    end_y = start_y + tile_dim
    for (j in 1:nywindow) {
      ext_crop = terra::ext(c(start_x, end_x, start_y, end_y))
      img_crop = terra::crop(raster_img, ext_crop, snap="in")
      img_crop_rescaled = terra::aggregate(img_crop, reduce_resolution)
      img_coordinates = as.character(paste('sx', start_x,
                                           'ex', end_x,
                                           'sy', start_y,
                                           'ey', end_y,
                                           sep='_'))
      file_name = file.path(folder_path, paste0(img_coordinates, '.png'))
      segmentation_result = python_segment_image(mesmer_app,
                                                 Matrix::as.matrix(img_crop_rescaled,
                                                                   wide=TRUE),
                                                 file_name)

      start_y = end_y - margin
      end_y = start_y + tile_dim
    }
    start_x = end_x - margin
    end_x = start_x + tile_dim
  }


  print(segmentation_result)
}
