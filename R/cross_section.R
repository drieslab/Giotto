
# cross section helper functions ####

#' @title create_crossSection_object
#' @name create_crossSection_object
#' @description create a crossSection object
#' @param name name of cress section object. (default = cross_sectino)
#' @param method method to define the cross section plane.
#' @param thickness_unit unit of the virtual section thickness. If "cell", average size of the observed cells is used as length unit. If "natural", the unit of cell location coordinates is used.(default = cell)
#' @param slice_thickness thickness of slice
#' @param cell_distance_estimate_method method to estimate average distance between neighobring cells. (default = mean)
#' @param extend_ratio deciding the span of the cross section meshgrid, as a ratio of extension compared to the borders of the vitural tissue section. (default = 0.2)
#' @param plane_equation a numerical vector of length 4, in the form of c(A,B,C,D), which defines plane Ax+By+Cz=D.
#' @param mesh_grid_n numer of meshgrid lines to generate along both directions for the cross section plane.
#' @param mesh_obj object that stores the cross section meshgrid information.
#' @param cell_subset cells selected by the cross section
#' @param cell_subset_spatial_locations locations of cells selected by the cross section
#' @param cell_subset_projection_locations 3D projection coordinates of selected cells onto the cross section plane
#' @param cell_subset_projection_PCA pca of projection coordinates
#' @param cell_subset_projection_coords 2D PCA coordinates of selected cells in the cross section plane
create_crossSection_object <- function(name=NULL,
                                       method=NULL,
                                       thickness_unit=NULL,
                                       slice_thickness=NULL,
                                       cell_distance_estimate_method=NULL,
                                       extend_ratio=NULL,
                                       plane_equation=NULL,
                                       mesh_grid_n=NULL,
                                       mesh_obj=NULL,
                                       cell_subset=NULL,
                                       cell_subset_spatial_locations=NULL,
                                       cell_subset_projection_locations=NULL,
                                       cell_subset_projection_PCA=NULL,
                                       cell_subset_projection_coords=NULL){

  crossSection_obj = list("method"=method,
                          "thickness_unit"=thickness_unit,
                          "slice_thickness" = slice_thickness,
                          "plane_equation"=plane_equation,
                          "mesh_grid_n"=mesh_grid_n,
                          "mesh_obj"=mesh_obj,
                          "cell_subset"=cell_subset,
                          "cell_subset_spatial_locations"=cell_subset_spatial_locations,
                          "cell_subset_projection_locations"=cell_subset_projection_locations,
                          "cell_subset_projection_PCA"=cell_subset_projection_PCA,
                          "cell_subset_projection_coords"=cell_subset_projection_coords)
}

#' @title read_crossSection
#' @name read_crossSection
#' @description read a cross section object from a giotto object
#' @param gobject gobject
#' @param name name
#' @param spatial_network_name spatial_network_name
#' @keywords internal
read_crossSection <- function(gobject,
                              name=NULL,
                              spatial_network_name=NULL){
  if(is.null(spatial_network_name)){
    stop("spatial_network_name is not specified.")
  }else if (!is.element(spatial_network_name,names(gobject@spatial_network))){
    stop(paste0(spatial_network_name, " has not been created."))
  }else {
    sp_network_obj = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = TRUE)
    if (length(sp_network_obj$crossSectionObjects)==0){
      stop("No cross section object has been created.")
    }else if (is.null(name)){
      sprintf("cross section object is not specified, reading the last one %s from the existing list",
              names(sp_network_obj$crossSectionObjects)[length(sp_network_obj$crossSectionObjects)])
      crossSection_obj = sp_network_obj$crossSectionObjects[[length(sp_network_obj$crossSectionObjects)]]
    }else if(!is.element(name,names(sp_network_obj$crossSectionObjects))){
      stop(paste0(name, " has not been created."))
    }
    else{
      crossSection_obj = sp_network_obj$crossSectionObjects[[name]]
    }
  }
  return(crossSection_obj)
}

#' @title get_distance
#' @name get_distance
#' @description estimate average distance between neighboring cells with network table as input
#' @param networkDT networkDT
#' @param method method
#' @keywords internal
get_distance <- function(networkDT,
                         method=c("mean","median")
                         ){

  if (method=="median"){
    distance = stats::median(networkDT$distance)
  }else if(method=="mean"){
    distance = mean(networkDT$distance)
  }
  return(distance)
}

#' @title estimateCellCellDistance
#' @name estimateCellCellDistance
#' @description estimate average distance between neighboring cells
#' @param gobject gobject
#' @param spatial_network_name spatial_network_name
#' @param method method
#' @keywords internal
estimateCellCellDistance <- function(gobject,
                                     spatial_network_name="Delaunay_network",
                                     method=c("mean","median")
                                     ){

  delaunay_network_DT = gobject@spatial_network[[spatial_network_name]]$networkDT

  CellCellDistance = get_distance(networkDT= delaunay_network_DT,
                                              method=method)
  return(CellCellDistance)

}
#' @title get_sectionThickness
#' @name get_sectionThickness
#' @description get section thickness
#' @param gobject gobject
#' @param thickness_unit thickness_unit
#' @param spatial_network_name spatial_network_name
#' @param cell_distance_estimate_method cell_distance_estimate_method
#' @param plane_equation plane_equation
#' @keywords internal
get_sectionThickness <- function(gobject,thickness_unit=c("cell","natural"),
                                 slice_thickness = 2,
                                 spatial_network_name="Delaunay_network",
                                 cell_distance_estimate_method = c("mean","median"),
                                 plane_equation=NULL){

  thickness_unit = match.arg(thickness_unit, c("cell", "natural"))

  if (thickness_unit == "cell"){
    CellCellDistance = estimateCellCellDistance(gobject,
                                                 method = cell_distance_estimate_method,
                                            spatial_network_name = spatial_network_name)
    sectionThickness = CellCellDistance*slice_thickness
  }else if (thickness_unit=="natural"){
    sectionThickness = slice_thickness
  }
  return(sectionThickness)
}

#' @title projection_fun
#' @name projection_fun
#' @description project a point onto a plane
#' @param point_to_project point_to_project
#' @param plane_point plane_point
#' @param plane_norm plane_norm
#' @keywords internal
projection_fun <- function(point_to_project,plane_point,plane_norm){

  a = plane_norm[1]
  b = plane_norm[2]
  c = plane_norm[3]
  x = point_to_project[1]
  y = point_to_project[2]
  z = point_to_project[3]
  d = plane_point[1]
  e = plane_point[2]
  f = plane_point[3]
  t = (a*d - a*x + b*e - b*y + c*f - c*z)/(a^2+b^2+c^2)
  xp = x + t*a
  yp = y + t*b
  zp = z + t*c
  projection = c(xp,yp,zp)
  return(projection)
}

#' @title adapt_aspect_ratio
#' @name adapt_aspect_ratio
#' @description adapt the aspact ratio after inserting cross section mesh grid lines
#' @param current_ratio current_ratio
#' @param cell_locations cell_locations
#' @param sdimx sdimx
#' @param sdimy sdimy
#' @param sdimz sdimz
#' @param mesh_obj mesh_obj
#' @keywords internal
adapt_aspect_ratio <-function(current_ratio,cell_locations,
                              sdimx = NULL,sdimy = NULL,sdimz = NULL,
                              mesh_obj=NULL){
  x_range = max(cell_locations[[sdimx]]) - min(cell_locations[[sdimx]])
  y_range = max(cell_locations[[sdimy]]) - min(cell_locations[[sdimy]])
  z_range = max(cell_locations[[sdimz]]) - min(cell_locations[[sdimz]])

  x_mesh_range = max(mesh_obj$mesh_grid_lines$mesh_grid_lines_X) - min(mesh_obj$mesh_grid_lines$mesh_grid_lines_X)
  y_mesh_range = max(mesh_obj$mesh_grid_lines$mesh_grid_lines_Y) - min(mesh_obj$mesh_grid_lines$mesh_grid_lines_Y)
  z_mesh_range = max(mesh_obj$mesh_grid_lines$mesh_grid_lines_Z) - min(mesh_obj$mesh_grid_lines$mesh_grid_lines_Z)

  if (x_mesh_range>x_range){
    x_adapt =  x_mesh_range/x_range
  }else{
    x_adapt = 1
  }
  if (y_mesh_range>y_range){
    y_adapt =  y_mesh_range/y_range
  }else{
    y_adapt = 1
  }
  if (z_mesh_range>z_range){
    z_adapt =  z_mesh_range/z_range
  }else{
    z_adapt = 1
  }

  new_ratio = as.numeric(current_ratio)*c(as.numeric(x_adapt),as.numeric(y_adapt),as.numeric(z_adapt))
  new_ratio = new_ratio/min(new_ratio)
  return(new_ratio)
}

# mesh grid line helper functions ####

#' @title extend_vector
#' @name extend_vector
#' @description extend the range of a vector by a given ratio
#' @param x x
#' @param extend_ratio extend_ratio
#' @keywords internal
extend_vector <- function(x,extend_ratio){

  x_center = (max(x)+min(x))/2
  y = (x-x_center)*(extend_ratio+1)+x_center

  return(y)
}

#' @title find_x_y_ranges
#' @name find_x_y_ranges
#' @description get the extended ranges of x and y
#' @param data data
#' @param extend_ratio extend_ratio
#' @keywords internal
find_x_y_ranges <- function(data,extend_ratio){

  x_extend = extend_vector(data[,1],extend_ratio)
  y_extend = extend_vector(data[,2],extend_ratio)

  x_min = min(x_extend)
  x_max = max(x_extend)
  y_min = min(y_extend)
  y_max = max(y_extend)

  out = list("x_min"=x_min,
             "x_max"=x_max,
             "y_min"=y_min,
             "y_max"=y_max
  )
}

#' @title create_2d_mesh_grid_line_obj
#' @name create_2d_mesh_grid_line_obj
#' @description create 2d mesh grid line object
#' @param x_min x_min
#' @param x_max x_max
#' @param y_min y_min
#' @param y_max y_max
#' @param mesh_grid_n mesh_grid_n
#' @keywords internal
create_2d_mesh_grid_line_obj <- function(x_min,x_max,y_min,y_max,mesh_grid_n){

  x_grid = seq(x_min,x_max,length.out = mesh_grid_n)
  y_grid = seq(y_min,y_max,length.out = mesh_grid_n)

  mesh_grid_lines_X = cbind(matrix(rep(x_grid,mesh_grid_n),nrow = mesh_grid_n,byrow = T),
                            matrix(rep(x_grid,mesh_grid_n),nrow = mesh_grid_n,byrow = F))

  mesh_grid_lines_Y = cbind(matrix(rep(y_grid,mesh_grid_n),nrow = mesh_grid_n,byrow = F),
                            matrix(rep(y_grid,mesh_grid_n),nrow = mesh_grid_n,byrow = T))


  mesh_grid_line_obj_2d = list("mesh_grid_lines_X"=mesh_grid_lines_X,
                               "mesh_grid_lines_Y"=mesh_grid_lines_Y)
  return(mesh_grid_line_obj_2d)
}

#' @title reshape_to_data_point
#' @name reshape_to_data_point
#' @description reshape a mesh grid line object to data point matrix
#' @param mesh_grid_obj mesh_grid_obj
#' @keywords internal
reshape_to_data_point <- function(mesh_grid_obj){

  if (length(mesh_grid_obj)==3){
    data_points = cbind(as.vector(mesh_grid_obj[[1]]),
                        as.vector(mesh_grid_obj[[2]]),
                        as.vector(mesh_grid_obj[[3]]))
  }else if (length(mesh_grid_obj)==2){
    data_points = cbind(as.vector(mesh_grid_obj[[1]]),
                        as.vector(mesh_grid_obj[[2]])
    )
  }
  return(data_points)
}

#' @title reshape_to_mesh_grid_obj
#' @name reshape_to_mesh_grid_obj
#' @description reshape a data point matrix to a mesh grid line object
#' @param data_points data_points
#' @param mesh_grid_n mesh_grid_n
#' @keywords internal
reshape_to_mesh_grid_obj <- function(data_points,mesh_grid_n){

  if (dim(data_points)[2]==2){

    mesh_grid_lines_X = matrix(data_points[,1],nrow=mesh_grid_n,byrow=F)
    mesh_grid_lines_Y = matrix(data_points[,2],nrow=mesh_grid_n,byrow=F)

    mesh_grid_obj = list("mesh_grid_lines_X"=mesh_grid_lines_X,
                         "mesh_grid_lines_Y"=mesh_grid_lines_Y)

  }else if (dim(data_points)[2]==3){
    mesh_grid_lines_X = matrix(data_points[,1],nrow=mesh_grid_n,byrow=F)
    mesh_grid_lines_Y = matrix(data_points[,2],nrow=mesh_grid_n,byrow=F)
    mesh_grid_lines_Z = matrix(data_points[,3],nrow=mesh_grid_n,byrow=F)
    mesh_grid_obj = list("mesh_grid_lines_X"=mesh_grid_lines_X,
                         "mesh_grid_lines_Y"=mesh_grid_lines_Y,
                         "mesh_grid_lines_Z"=mesh_grid_lines_Z)
  }
  return(mesh_grid_obj)
}


#' @title transform_2d_mesh_to_3d_mesh
#' @name transform_2d_mesh_to_3d_mesh
#' @description transform 2d mesh to 3d mesh by reversing PCA
#' @param mesh_line_obj_2d mesh_line_obj_2d
#' @param pca_out pca_out
#' @param center_vec center_vec
#' @param mesh_grid_n mesh_grid_n
#' @keywords internal
transform_2d_mesh_to_3d_mesh <- function(mesh_line_obj_2d,pca_out,center_vec,mesh_grid_n){

  data_point_2d = reshape_to_data_point(mesh_line_obj_2d)
  center_mat = matrix(rep(center_vec,dim(data_point_2d)[1]),nrow=dim(data_point_2d)[1],byrow=T)
  data_point_3d = cbind(data_point_2d,rep(0,dim(data_point_2d)[1])) %*% t((pca_out$rotation))+center_mat
  mesh_grid_line_obj_3d = reshape_to_mesh_grid_obj(data_point_3d,mesh_grid_n)

  return(mesh_grid_line_obj_3d)
}

#' @title get_cross_section_coordinates
#' @name get_cross_section_coordinates
#' @description get local coordinates within cross section plane
#' @param cell_subset_projection_locations cell_subset_projection_locations
#' @keywords internal
get_cross_section_coordinates <- function(cell_subset_projection_locations){

  cell_subset_projection_PCA = stats::prcomp(cell_subset_projection_locations)

  cell_subset_projection_coords = cell_subset_projection_PCA$x[,c("PC1","PC2")]

  return(cell_subset_projection_coords)
}

#' @title create_mesh_grid_lines
#' @name create_mesh_grid_lines
#' @description create mesh grid lines for cross section
#' @param cell_subset_projection_locations cell_subset_projection_locations
#' @param extend_ratio extend_ratio
#' @param mesh_grid_n mesh_grid_n
#' @keywords internal
create_mesh_grid_lines <- function(cell_subset_projection_locations,extend_ratio,mesh_grid_n){

  cell_subset_projection_PCA = stats::prcomp(cell_subset_projection_locations)

  cell_subset_projection_coords = cell_subset_projection_PCA$x[,c("PC1","PC2")]

  x_y_ranges = find_x_y_ranges(cell_subset_projection_coords,extend_ratio)

  mesh_line_obj_2d = create_2d_mesh_grid_line_obj(x_y_ranges$x_min,
                                                  x_y_ranges$x_max,
                                                  x_y_ranges$y_min,
                                                  x_y_ranges$y_max,
                                                  mesh_grid_n)
  center_vec = apply(cell_subset_projection_locations,2,function(x) mean(x))
  mesh_grid_line_obj_3d = transform_2d_mesh_to_3d_mesh(mesh_line_obj_2d,
                                                       cell_subset_projection_PCA,
                                                       center_vec,
                                                       mesh_grid_n)
  return(mesh_grid_line_obj_3d)
}


# cross section creation function ####

#' @title createCrossSection
#' @description Create a virtual 2D cross section.
#' @param gobject giotto object
#' @param name name of cress section object. (default = cross_sectino)
#' @param spatial_network_name name of spatial network object. (default = Delaunay_network)
#' @param thickness_unit unit of the virtual section thickness. If "cell", average size of the observed cells is used as length unit. If "natural", the unit of cell location coordinates is used.(default = cell)
#' @param slice_thickness thickness of slice. default = 2
#' @param cell_distance_estimate_method method to estimate average distance between neighobring cells. (default = mean)
#' @param extend_ratio deciding the span of the cross section meshgrid, as a ratio of extension compared to the borders of the vitural tissue section. (default = 0.2)
#' @param method method to define the cross section plane.
#' If equation, the plane is defined by a four element numerical vector (equation) in the form of c(A,B,C,D), corresponding to a plane with equation Ax+By+Cz=D.
#' If 3 points, the plane is define by the coordinates of 3 points, as given by point1, point2, and point3.
#' If point and norm vector, the plane is defined by the coordinates of one point (point1) in the plane and the coordinates of one norm vector (normVector) to the plane.
#' If point and two plane vector, the plane is defined by the coordinates of one point (point1) in the plane and the coordinates of two vectors (planeVector1, planeVector2) in the plane.
#' (default = equation)
#' @param equation equation required by method "equation".equations needs to be a numerical vector of length 4, in the form of c(A,B,C,D), which defines plane Ax+By+Cz=D.
#' @param point1 coordinates of the first point required by method "3 points","point and norm vector", and "point and two plane vectors".
#' @param point2 coordinates of the second point required by method "3 points"
#' @param point3 coordinates of the third point required by method "3 points"
#' @param normVector coordinates of the norm vector required by method "point and norm vector"
#' @param planeVector1 coordinates of the first plane vector required by method "point and two plane vectors"
#' @param planeVector2 coordinates of the second plane vector required by method "point and two plane vectors"
#' @param mesh_grid_n numer of meshgrid lines to generate along both directions for the cross section plane.
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object with updated spatial network slot
#' @details Creates a virtual 2D cross section object for a given spatial network object. The users need to provide the definition of the cross section plane (see method).
#' @export
createCrossSection <- function(gobject,
                               name="cross_section",
                               spatial_network_name = "Delaunay_network",
                               thickness_unit = c("cell","natural"),
                               slice_thickness = 2,
                               cell_distance_estimate_method = "mean",
                               extend_ratio = 0.2,
                               method=c("equation","3 points","point and norm vector","point and two plane vectors"),
                               equation=NULL,
                               point1=NULL,point2=NULL,point3=NULL,
                               normVector=NULL,
                               planeVector1=NULL,planeVector2=NULL,
                               mesh_grid_n = 20,
                               return_gobject = TRUE
){

  # read spatial locations
  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, grepl("sdim", colnames(spatial_locations)),
                                        with = F]
  spatial_locations = as.matrix(spatial_locations)
  rownames(spatial_locations) = gobject@cell_ID
  cell_ID_vec = c(1:nrow(spatial_locations))
  names(cell_ID_vec) = rownames(spatial_locations)

  # generate section plane equation

  method = match.arg(method, c("equation","3 points","point and norm vector","point and two plane vectors"))

  if (method == "equation"){
    if (is.null(equation)){
      print("equation was not provided.")
    }else{
      plane_equation = equation
      plane_equation[4] = -equation[4]
    }
  }else if (method == "point and norm vector"){
    if (is.null(point1)|is.null(normVector)){
      print("either point or norm vector was not provided.")
    }else{
      plane_equation = c()
      plane_equation[1:3] = normVector
      plane_equation[4] = -point1 %*% normVector
    }
  }else if (method == "point and two plane vectors"){
    if(is.null(point1)|is.null(planeVector1)|is.null(planeVector2)){
      print("either point or any of the two plane vectors was not provided.")
    }else{
      normVector = crossprod(planeVector1,planeVector2)
      plane_equation[1:3] = normVector
      plane_equation[4] = -point1 %*% normVector
    }
  }else if (method == "3 points"){
    if (is.null(point1)|is.null(point2)|is.null(point3)){
      print("not all three points were provided.")
    }else{
      planeVector1 = point2-point1;
      planeVector2 = point3-point1;
      normVector = crossprod(planeVector1,planeVector2)
      plane_equation[1:3] = normVector
      plane_equation[4] = -point1 %*% normVector
    }
  }
  names(plane_equation)=c("A","B","C","D")

  # determine section thickness
  thickness_unit = match.arg(thickness_unit, c("cell", "natural"))
  sectionThickness = get_sectionThickness(gobject,thickness_unit=thickness_unit,
                                          slice_thickness = slice_thickness,
                                          spatial_network_name=spatial_network_name,
                                          cell_distance_estimate_method = cell_distance_estimate_method,
                                          plane_equation=plane_equation)

  max_distance_to_section_plane = sectionThickness/2

  # calculate distances to cross section
  spatial_locations_mat = cbind(spatial_locations,as.matrix(rep(1,dim(spatial_locations)[1])))
  norm_vec <- function(x) sqrt(sum(x^2))
  distance_to_plane_vector = abs(spatial_locations_mat %*% as.matrix(plane_equation)/norm_vec(plane_equation[1:3]))

  # select cells within section ###
  cell_subset = distance_to_plane_vector<=max_distance_to_section_plane

  # project the selected cells onto the section plane ###
  cell_subset_spatial_locations = spatial_locations[cell_subset,]

  ## find a point on the section plane ##
  if (plane_equation["A"]!=0){
    plane_point = c(-plane_equation["D"]/plane_equation["A"],0,0)
  }else if (plane_equation["B"]!=0){
    plane_point = c(0,-plane_equation["D"]/plane_equation["B"],0)
  }else if (plane_equation["C"]!=0){
    plane_point = c(0,0,-plane_equation["D"]/plane_equation["C"])
  }
  ## find the projection Xp,Yp,Zp coordinates ##
  cell_subset_projection_locations = t(apply(cell_subset_spatial_locations,1,function(x) projection_fun(x,plane_point = plane_point, plane_norm = plane_equation[1:3])))

  # get the local coordinates of selected cells on the section plane
  cell_subset_projection_PCA = stats::prcomp(cell_subset_projection_locations)
  cell_subset_projection_coords = get_cross_section_coordinates(cell_subset_projection_locations)

  # create mesh grid lines for the cross section ###
  mesh_grid_lines = create_mesh_grid_lines(cell_subset_projection_locations,extend_ratio,mesh_grid_n)
  mesh_obj = list("mesh_grid_lines" = mesh_grid_lines)

  ### save and update the spatial object ###

  crossSection_obj <- create_crossSection_object(method=method,
                                                 thickness_unit=thickness_unit,
                                                 slice_thickness=slice_thickness,
                                                 cell_distance_estimate_method=cell_distance_estimate_method,
                                                 extend_ratio=extend_ratio,
                                                 plane_equation=plane_equation,mesh_grid_n=mesh_grid_n,
                                                 mesh_obj=mesh_obj,cell_subset=cell_subset,
                                                 cell_subset_spatial_locations=cell_subset_spatial_locations,
                                                 cell_subset_projection_locations=cell_subset_projection_locations,
                                                 cell_subset_projection_PCA=cell_subset_projection_PCA,
                                                 cell_subset_projection_coords=cell_subset_projection_coords)


  if (return_gobject == TRUE) {

    cs_names = names(gobject@spatial_network[[spatial_network_name]]$crossSectionObjects)
    if (name %in% cs_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    gobject@spatial_network[[spatial_network_name]]$crossSectionObjects[[name]] = crossSection_obj

    return(gobject)
  }
  else {
    return(crossSection_obj)
  }


}


# cross section visual functions ####

####
#' @title crossSectionGenePlot
#' @name crossSectionGenePlot
#' @description Visualize cells and gene expression in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj crossSection object
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatGenePlot2D
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatGenePlot3D}} and \code{\link{spatGenePlot2D}}
#'
crossSectionGenePlot <-function(gobject=NULL,
                                crossSection_obj=NULL,
                                name=NULL,
                                spatial_network_name = "Delaunay_network",
                                default_save_name = "crossSectionGenePlot",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }

  cell_subset = crossSection_obj$cell_subset
  cell_subset_projection_coords = crossSection_obj$cell_subset_projection_coords
  # modify gobject based on crossSection object
  subset_cell_IDs = gobject@cell_metadata$cell_ID[cell_subset]
  temp_gobject = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)
  temp_gobject@spatial_locs$sdimx=cell_subset_projection_coords[,1]
  temp_gobject@spatial_locs$sdimy=cell_subset_projection_coords[,2]
  temp_gobject@spatial_locs$sdimz=rep(0,dim(cell_subset_projection_coords)[1])
  # call spatGenePlot2D to generate the plots
  spatGenePlot2D(gobject = temp_gobject,
                 spatial_network_name = spatial_network_name,
                 default_save_name = default_save_name,
                 ...)
}
####

#' @title crossSectionPlot
#' @name crossSectionPlot
#' @description Visualize cells in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatPlot2D
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{crossSectionPlot}}
crossSectionPlot <-function(gobject,
                            crossSection_obj = NULL,
                            name=NULL,
                            spatial_network_name = "Delaunay_network",
                            default_save_name = "crossSectionPlot",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }


  cell_subset = crossSection_obj$cell_subset
  cell_subset_projection_coords = crossSection_obj$cell_subset_projection_coords
  # modify gobject based on crossSection object
  subset_cell_IDs = gobject@cell_metadata$cell_ID[cell_subset]
  temp_gobject = subsetGiotto(gobject, cell_ids = subset_cell_IDs)
  temp_gobject@spatial_locs$sdimx=cell_subset_projection_coords[,1]
  temp_gobject@spatial_locs$sdimy=cell_subset_projection_coords[,2]
  temp_gobject@spatial_locs$sdimz=rep(0,dim(cell_subset_projection_coords)[1])
  # call spatGenePlot2D to generate the plots
  spatPlot2D(gobject = temp_gobject,
             spatial_network_name = spatial_network_name,
             default_save_name = default_save_name,...)


}

####
#' @title crossSectionGenePlot3D
#' @name crossSectionGenePlot3D
#' @description Visualize cells and gene expression in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param other_cell_color color of cells outside the cross section. default = transparent.
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatGenePlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
crossSectionGenePlot3D <-function(gobject,
                                  crossSection_obj = NULL,
                                  name=NULL,
                                  spatial_network_name = "Delaunay_network",
                                  other_cell_color = alpha("lightgrey", 0),
                                  default_save_name = "crossSectionGenePlot3D",...){


  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }


  cell_subset = crossSection_obj$cell_subset
  cell_subset_projection_coords = crossSection_obj$cell_subset_projection_coords
  # modify gobject based on crossSection object
  subset_cell_IDs = gobject@cell_metadata$cell_ID[cell_subset]
  # call spatGenePlot3D to generate the plots
  spatGenePlot3D(gobject,
                 select_cells = subset_cell_IDs,
                 other_cell_color = other_cell_color,
                 default_save_name = default_save_name,...)
}
####
#' @title crossSectionPlot3D
#' @name crossSectionPlot3D
#' @description Visualize cells in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of cells outside the cross section. default = transparent.
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatPlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
crossSectionPlot3D <-function(gobject,
                              crossSection_obj = NULL,
                              name=NULL,
                              spatial_network_name = "Delaunay_network",
                              show_other_cells = T,
                              other_cell_color = alpha("lightgrey", 0),
                              default_save_name = "crossSection3D",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }

  cell_subset = crossSection_obj$cell_subset
  cell_subset_projection_coords = crossSection_obj$cell_subset_projection_coords
  # modify gobject based on crossSection object
  subset_cell_IDs = gobject@cell_metadata$cell_ID[cell_subset]
  # temp_gobject = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)
  # temp_gobject@spatial_locs$sdimx=cell_subset_projection_coords[,1]
  # temp_gobject@spatial_locs$sdimy=cell_subset_projection_coords[,2]
  # temp_gobject@spatial_locs$sdimz=rep(0,dim(cell_subset_projection_coords)[1])
  #
  # call spatPlot3D to generate the plots
  spatPlot3D(gobject=gobject,
             ##
             select_cells = subset_cell_IDs,
             ##
             show_other_cells = show_other_cells,
             other_cell_color = other_cell_color,
             default_save_name = default_save_name,...)
}


####
#' @title insertCrossSectionSpatPlot3D
#' @name insertCrossSectionSpatPlot3D
#' @description Visualize the meshgrid lines of cross section together with cells
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param mesh_grid_color color for the meshgrid lines
#' @param mesh_grid_width width for the meshgrid lines
#' @param mesh_grid_style style for the meshgrid lines
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimy')
#' @param show_other_cells display not selected cells
#' @param axis_scale axis_scale
#' @param custom_ratio custom_ratio
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatPlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
insertCrossSectionSpatPlot3D <- function(gobject,
                                         crossSection_obj=NULL,
                                         name=NULL,
                                         spatial_network_name = "Delaunay_network",
                                         mesh_grid_color = "#1f77b4",
                                         mesh_grid_width = 3,
                                         mesh_grid_style = "dot",
                                         sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
                                         show_other_cells = F,
                                         axis_scale = c("cube", "real", "custom"),
                                         custom_ratio = NULL,
                                         default_save_name = "spat3D_with_cross_section",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }



  pl = spatPlot3D(gobject,
                  sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
                  show_other_cells = show_other_cells,
                  show_plot = FALSE,
                  return_plot = TRUE,
                  save_plot = FALSE,
                  default_save_name = default_save_name,...)

  for (i in 1:dim(crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X)[2]){

    pl = pl %>% plotly::add_trace(x = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X[,i],
                                  y = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Y[,i],
                                  z = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Z[,i],
                                  mode = 'lines',type = 'scatter3d',
                                  line = list(color = mesh_grid_color, width = mesh_grid_width,dash = mesh_grid_style))
  }

  current_ratio = plotly_axis_scale_3D(gobject@spatial_locs,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                       mode = axis_scale,custom_ratio = custom_ratio)

  new_ratio = adapt_aspect_ratio(current_ratio,gobject@spatial_locs,
                                 sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mesh_obj=crossSection_obj$mesh_obj)

  pl = pl %>% plotly::layout(showlegend = FALSE,
                             scene = list(
                               aspectmode='manual',
                               aspectratio = list(x=new_ratio[[1]],
                                                  y=new_ratio[[2]],
                                                  z=new_ratio[[3]])))

  return(pl)


}
####
#' @title insertCrossSectionGenePlot3D
#' @name insertCrossSectionGenePlot3D
#' @description Visualize cells and gene expression in a virtual cross section according to spatial coordinates
#' @param gobject giotto object
#' @param crossSection_obj cross section object as alternative input. default = NULL.
#' @param name name of virtual cross section to use
#' @param spatial_network_name name of spatial network to use
#' @param mesh_grid_color color for the meshgrid lines
#' @param mesh_grid_width width for the meshgrid lines
#' @param mesh_grid_style style for the meshgrid lines
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimy')
#' @param show_other_cells display not selected cells
#' @param axis_scale axis_scale
#' @param custom_ratio custom_ratio
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for spatGenePlot3D
#' @return ggplot
#' @details Description of parameters.
#' @export
insertCrossSectionGenePlot3D <- function(gobject,
                                         crossSection_obj=NULL,
                                         name=NULL,
                                         spatial_network_name = "Delaunay_network",
                                         mesh_grid_color = "#1f77b4",
                                         mesh_grid_width = 3,
                                         mesh_grid_style = "dot",
                                         sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
                                         show_other_cells = F,
                                         axis_scale = c("cube", "real", "custom"),
                                         custom_ratio = NULL,
                                         show_plot = NA, return_plot = NA, save_plot = NA,
                                         save_param = list(),
                                         default_save_name = "spatGenePlot3D_with_cross_section",...){

  # load cross section object
  if (!is.null(crossSection_obj)){
    crossSection_obj = crossSection_obj
  }else{
    crossSection_obj = read_crossSection(gobject,name=name,spatial_network_name = spatial_network_name)
  }

  pl = spatGenePlot3D(gobject,
                      show_other_cells = F,
                      axis_scale = axis_scale,
                      custom_ratio = custom_ratio,
                      show_plot = FALSE,
                      return_plot = TRUE,
                      save_plot = FALSE,
                      default_save_name = default_save_name,...)
  for (i in 1:dim(crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X)[2]){

    pl = pl %>% plotly::add_trace(x = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_X[,i],
                                  y = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Y[,i],
                                  z = crossSection_obj$mesh_obj$mesh_grid_lines$mesh_grid_lines_Z[,i],
                                  mode = 'lines+markers',type = 'scatter3d',color = mesh_grid_color,
                                  marker = list(color=alpha(mesh_grid_color,0)),
                                  line = list(color = mesh_grid_color, width = mesh_grid_width,dash = mesh_grid_style))
  }

  current_ratio = plotly_axis_scale_3D(gobject@spatial_locs,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                       mode = axis_scale,custom_ratio = custom_ratio)

  new_ratio = adapt_aspect_ratio(current_ratio,gobject@spatial_locs,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mesh_obj = crossSection_obj$mesh_obj)

  pl = pl %>% plotly::layout(showlegend = FALSE,
                             scene = list(
                               aspectmode='manual',
                               aspectratio = list(x=new_ratio[[1]],
                                                  y=new_ratio[[2]],
                                                  z=new_ratio[[3]])))

  cowplot = pl
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject,
                                                              param = "show_plot"), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject,
                                                              param = "save_plot"), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject,
                                                                  param = "return_plot"), return_plot)
  if (show_plot == TRUE) {
    print(cowplot)
  }
  if (save_plot == TRUE) {
    do.call("all_plots_save_function", c(list(gobject = gobject,
                                              plot_object = cowplot, default_save_name = default_save_name),
                                         save_param))
  }
  if (return_plot == TRUE) {
    return(cowplot)
  }

}
