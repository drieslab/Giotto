#' @title doCellSegmentation
#' @name doCellSegmentation
#' @description segment cells in Dapi image
#' @param image_files #name of image files - can be a single image or multiple images
#' @param img_channel cytoplasm, nucleus
#' @param image_diameter diameter chosen by user
#' @param model_type cytoplasm ("cyto") or nuclear ("nuclei")
#' @param python_path specify specific path to python where the cellpose package has been downloaded
#' @details
#' This function is based on the cellpose algorithm implemented in python.
#' It segments the cells and nuclei in tissue on the complete
#' image. It saves all the files in the specified folder.
#'
#' @export

doCellSegmentation <- function(image_files,
                               img_channels,
                               image_diameter,
                               model_type,
                               python_path){

  # prepare python path and segmentation script
  reticulate::use_python(required = T, python = python_path)
  python_segmentation_function = system.file("python",
                                             "python_segmentation.py",
                                             package = 'Giotto')
  reticulate::source_python(file = python_segmentation_function)

  # Running the cellpose model on user input image data
  results <- cellsegmentation_cellpose(image_files, img_channels, image_diameter)

  cellmasks <- results[[2]] # Matrix with polyIDs for each segmented cell (including boundaries and interior)
  cellmasks_outlines <- results[[1]] # Matrix with the segmented boundaries (1) and background (0)

  number <- max(cellmasks) # Determining the total number of segmented cells

  outline_matrix <- which(cellmasks_outlines == 1, arr.ind = TRUE) # Identifying polygon boundary coordinates

  # Identifying coordinates of the boundaries and polygon interior for each polygon ID
  cell_positions <- NULL
  for (i in 1:number) {
    cell_locations <- which(cellmasks == i, arr.ind = TRUE)
    positions <- as.matrix(cell_locations)
    positions <- cbind(positions, rep(i, nrow(positions)))
    cell_positions <- rbind(cell_positions, positions)
  }

  # Matching the segmented boundaries to their respective poly IDs
  outline_matrix_df <- data.frame(outline_matrix)
  colnames(outline_matrix_df) <- c("x", "y")
  cell_positions_df <- data.frame(cell_positions)
  colnames(cell_positions_df) <- c("x", "y", "poly ID")

  # Final result gives a dataframe with the x and y coordinates of the polygon boundaries and their respective poly IDs
  segmentationresults <- merge(outline_matrix_df, cell_positions_df[, c("x", "y", "poly ID")], how = "left", on = c("x", "y"))

  print(segmentationresults)

}
