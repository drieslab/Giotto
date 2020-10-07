

#' mini Giotto object for spatial single-cell resolution data
#'
#' Mini Giotto object created from the seqFISH+ data.
#'
#' @docType data
#'
#' @usage data(mini_giotto_single_cell)
#'
#' @format An object of class \code{"giotto"}; see \code{\link[Giotto]{createGiottoObject}}.
#'
#' @keywords datasets
#'
#' @references Eng et al. (2019) Nature
#' (\href{https://www.nature.com/articles/s41586-019-1049-y}{PubMed})
#'
#'
#' @examples
#' data(mini_giotto_single_cell)
#'
#' \dontrun{spatPlot2D(mini_giotto_single_cell,cell_color = 'cell_types', point_size = 5)}
"mini_giotto_single_cell"



#' mini Giotto object for spatial single-cell 3D data
#'
#' Mini Giotto object created from the STARmap data.
#'
#' @docType data
#'
#' @usage data(mini_giotto_3D)
#'
#' @format An object of class \code{"giotto"}; see \code{\link[Giotto]{createGiottoObject}}.
#'
#' @keywords datasets
#'
#' @references Wang et al. (2018) Science
#' (\href{https://pubmed.ncbi.nlm.nih.gov/29930089/}{PubMed})
#'
#'
#' @examples
#' data(mini_giotto_3D)
#'
#' \dontrun{spatPlot3D(mini_giotto_3D, cell_color = 'cell_types', point_size = 5)}
"mini_giotto_3D"





#' mini Giotto object for spatial multi-cell resolution data
#'
#' Mini Giotto object created from the Brain Visium 10X data.
#'
#' @docType data
#'
#' @usage data(mini_giotto_multi_cell)
#'
#' @format An object of class \code{"giotto"}; see \code{\link[Giotto]{createGiottoObject}}.
#'
#' @keywords datasets
#'
#' @references 10 Genomics Visium technology
#' (\href{https://www.10xgenomics.com/spatial-transcriptomics/}{10xgenomics})
#'
#'
#' @examples
#' data(mini_giotto_multi_cell)
#'
#' \dontrun{spatPlot(mini_giotto_multi_cell, cell_color = 'cell_types', point_size = 5)}
"mini_giotto_multi_cell"


