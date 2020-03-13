### plot voronoi tile and cell neighbors

plotVoronoiTile <- function(gobject,
                            color_by = NULL,
                            bound = NULL,
                            subset_vec = NULL,
                            show_plot = NA,
                            return_plot = NA,
                            save_plot = NA,
                            save_param = list(),
                            default_save_name = "voronoiPlot"
                            ){

  # create input data frame
  df = as.data.frame(cbind(gobject@spatial_locs,
                           gobject@cell_metadata[[color_by]]))
  names(df)[4] = "det_cell_type"
  if (is.null(subset_vec)){
    df_subset = df
  }else{
    df_subset = df[subset_vec,]
  }
  # calculate a proper rectangular bound
  if (is.null(bound)){
    bound_to_use = c(min(df$sdimx),
                     max(df$sdimx),
                     min(df$sdimy),
                     max(df$sdimy))
  }else{
    bound_to_use = bound
  }
  # generate veronoi plot

  # customize color pallette
  det_cell_type = df$det_cell_type
  factor_vec_master = as.factor(c(det_cell_type,"zzz_background_zzz"))
  color_vec_master = c(hue_pal()(length(levels(factor_vec_master))-1),"grey")
  names(color_vec_master) = levels(factor_vec_master)

  if (is.null(subset_vec)){
    factor_vec = as.factor(det_cell_type)
    color_vec = hue_pal()(length(levels(factor_vec)))
  }else{
    det_cell_type_subset = det_cell_type
    det_cell_type_subset[!subset_vec] = "zzz_background_zzz"
    factor_vec =  as.factor(det_cell_type_subset)
    color_vec = color_vec_master[levels(factor_vec)]
  }

  library(ggforce)
  pl = ggplot()
  pl = pl + theme_classic()
  pl = pl + geom_voronoi_tile(data = df,
                                  aes(x = sdimx, y = sdimy, group = -1L,fill = factor_vec),
                                  color="white",
                                  bound = bound_to_use) + scale_fill_manual(values = color_vec)
  pl = pl + geom_point(data = df_subset,
                           aes(x = sdimx, y = sdimy), size = 0.5)

  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject,
                                                              param = "show_plot"), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject,
                                                              param = "save_plot"), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject,
                                                                  param = "return_plot"), return_plot)
  if (show_plot == TRUE) {
    print(pl)
  }
  if (save_plot == TRUE) {
    do.call("all_plots_save_function", c(list(gobject = gobject,
                                              plot_object = pl, default_save_name = default_save_name),
                                         save_param))
  }
  if (return_plot == TRUE) {
    return(pl)
  }



}

plotDelaunayNeighbors <- function(gobject,
                                  color_by = NULL,
                                  bound = NULL,
                                  select_vec = NULL,
                                  show_plot = NA,
                                  return_plot = NA,
                                  save_plot = NA,
                                  save_param = list(),
                                  default_save_name = "delaunayNeighborsPlot"){

  # get delaunay neighbors of selected cell
  sp_network = annotateSpatialNetwork(gobject, spatial_network_name = "Delaunay_network", cluster_column = color_by,create_full_network = T)
  sp_network_select = sp_network[is.element(sp_network$from,gobject@cell_metadata$cell_ID[select_vec]),]
  sp_network_select_neighboring_cells = unique(sp_network_select$to)
  df_select = df[select_vec,]

  subset_vec = is.element(gobject@cell_metadata$cell_ID,sp_network_select_neighboring_cells)
  # create input data frame
  df = as.data.frame(cbind(gobject@spatial_locs,
                           gobject@cell_metadata[[color_by]]))
  names(df)[4] = "det_cell_type"
  if (is.null(subset_vec)){
    df_subset = df
  }else{
    df_subset = df[subset_vec,]
  }
  # calculate a proper rectangular bound
  if (is.null(bound)){
    bound_to_use = c(min(df$sdimx),
                     max(df$sdimx),
                     min(df$sdimy),
                     max(df$sdimy))
  }else{
    bound_to_use = bound
  }
  # generate veronoi plot

  # customize color pallette
  det_cell_type = df$det_cell_type
  factor_vec_master = as.factor(c(det_cell_type,"zzz_background_zzz"))
  color_vec_master = c(hue_pal()(length(levels(factor_vec_master))-1),"grey")
  names(color_vec_master) = levels(factor_vec_master)

  det_cell_type_subset = det_cell_type
  det_cell_type_subset[!subset_vec] = "zzz_background_zzz"
  factor_vec =  as.factor(det_cell_type_subset)
  color_vec = color_vec_master[levels(factor_vec)]

  library(ggforce)
  pl = ggplot()
  pl = pl + theme_classic()
  pl = pl + geom_voronoi_tile(data = df,
                                  aes(x = sdimx, y = sdimy, group = -1L,fill = factor_vec),
                                  color="white",
                                  bound = bound_to_use) + scale_fill_manual(values = color_vec)
  pl = pl + geom_point(data = df_select,
                           aes(x = sdimx, y = sdimy), size = 0.5)

  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject,
                                                              param = "show_plot"), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject,
                                                              param = "save_plot"), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject,
                                                                  param = "return_plot"), return_plot)
  if (show_plot == TRUE) {
    print(pl)
  }
  if (save_plot == TRUE) {
    do.call("all_plots_save_function", c(list(gobject = gobject,
                                              plot_object = pl, default_save_name = default_save_name),
                                         save_param))
  }
  if (return_plot == TRUE) {
    return(pl)
  }


}
