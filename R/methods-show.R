
#' @include classes.R
#' @importFrom methods show
#' @importFrom methods setMethod
NULL


# ------------------------------------------------------ #







# Giotto ####

#' show method for giotto class
#' @param object giotto object
#' @aliases show,giotto-method
#' @docType methods
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giotto",
  definition = function(object) {
    spat_unit = feat_type = prints = name = img_type = name = NULL

    cat("An object of class",  class(object), "\n")


    # active spat_unit and feat_type
    active_su = try(instructions(object, 'active_spat_unit'), silent = TRUE)
    active_ft = try(instructions(object, 'active_feat_type'), silent = TRUE)
    if(!inherits(active_su, 'try-error')) {
      cat('>Active spat_unit: ', active_su, '\n')
    }
    if(!inherits(active_ft, 'try-error')) {
      cat('>Active feat_type: ', active_ft, '\n')
    }


    cat('[SUBCELLULAR INFO]\n')
    if(!is.null(object@spatial_info)) cat('polygons      :', wrap_txt(list_spatial_info_names(object)), '\n')
    if(!is.null(object@feat_info)) cat('features      :', wrap_txt(list_feature_info_names(object)), '\n')


    mini_avail_print = function(avail_dt) {
      if(!'spat_unit' %in% colnames(avail_dt)) {
        avail_dt[, spat_unit := '']
      } else avail_dt[, spat_unit := paste0('[', spat_unit, ']')]
      if(!'feat_type' %in% colnames(avail_dt)) {
        avail_dt[, feat_type := '']
      } else avail_dt[, feat_type := paste0('[', feat_type, ']')]
      avail_dt[, prints := paste0(spat_unit, feat_type)]

      unique_entry = avail_dt[, unique(prints)]
      for(entry in unique_entry) {
        cat('  ', entry, paste0(' ', wrap_txt(avail_dt[prints == entry, name])), '\n', sep = '')
      }
    }


    cat('[AGGREGATE INFO]\n')
    avail_expr = list_expression(object)
    if(!is.null(avail_expr)) {
      cat('expression -----------------------\n')
      mini_avail_print(avail_expr)
    }

    avail_sl = list_spatial_locations(object)
    if(!is.null(avail_sl)) {
      cat('spatial locations ----------------\n')
      mini_avail_print(avail_sl)
    }

    avail_sn = list_spatial_networks(object)
    if(!is.null(avail_sn)) {
      cat('spatial networks -----------------\n')
      mini_avail_print(avail_sn)
    }

    avail_se = list_spatial_enrichments(object)
    if(!is.null(avail_se)) {
      cat('spatial enrichments --------------\n')
      mini_avail_print(avail_se)
    }

    avail_dim = list_dim_reductions(object)
    if(!is.null(avail_dim)) {
      cat('dim reduction --------------------\n')
      mini_avail_print(avail_dim)
    }

    avail_nn = list_nearest_networks(object)
    if(!is.null(avail_nn)) {
      cat('nearest neighbor networks --------\n')
      mini_avail_print(avail_nn)
    }

    avail_im = list_images(object)
    if(!is.null(avail_im)) {
      cat('attached images ------------------\n')
      if('image' %in% avail_im$img_type) {
        if(sum(avail_im$img_type == 'image') > 3) {
          cat(  'giottoLargeImage :', sum(avail_im$img_type == 'image'), 'items...\n')
        } else {
          cat('giottoImage      :', wrap_txt(avail_im[img_type == 'image', name]), '\n')
        }
      }
      if('largeImage' %in% avail_im$img_type) {
        if(sum(avail_im$img_type == 'largeImage') > 3) {
          cat(  'giottoLargeImage :', sum(avail_im$img_type == 'largeImage'), 'items...\n')
        } else {
          cat(  'giottoLargeImage :', wrap_txt(avail_im[img_type == 'largeImage', name]), '\n')
        }
      }
    }

    cat(wrap_txt('\n\nUse objHistory() to see steps and params used\n\n'))
    invisible(x = NULL)
  }
)




# packedGiotto ####
setMethod("show", signature(object='packedGiotto'),
          function(object) {
            print(paste("This is a", class(object), "object. Use 'Giotto::unwrap()' to unpack it"))
          }
)








# subobjects ####

## exprObj ####

#' show method for exprObj class
#' @param object expression object
#' @aliases show,exprObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('exprObj'), function(object) {

    show_class_and_name(object)

    # print spat/feat and provenance info
    show_spat_and_feat(object)
    show_prov(object)

    cat('\ncontains:\n')
    # preview matrix

    # * Matrix sparseMatrix specific *
    if(inherits(slot(object, 'exprMat'), 'sparseMatrix')) {
      print_cap = capture.output(Matrix::printSpMatrix2(x = slot(object, 'exprMat'),
                                                        zero.print = ".",
                                                        col.names = FALSE,
                                                        note.dropping.colnames = FALSE,
                                                        suppRows = TRUE,
                                                        suppCols = TRUE,
                                                        width = 40,
                                                        maxp = 80))
      print_cap = print_cap[-which(print_cap == ' ..............................')]
      writeLines(gsub(pattern = "in show(.*?))'", replacement = '', x = print_cap))
      cat('\n First four colnames:')
      cat('\n', wrap_txt(head(colnames(slot(object, 'exprMat')), 4), strWidth = 40), '\n')
    } else if(inherits(slot(object, 'exprMat'), 'denseMatrix')) {
      abb_mat(object, nrows = 10, ncols = 10, header = FALSE)
    } else {
      # * other matrices *
      print(slot(object, 'exprMat'))
      cat('\n')
    }

  })











## cellMetaObj ####
setMethod('show', signature('cellMetaObj'), function(object) {

  cat('An object of class', class(object), '\n')
  show_spat_and_feat(object)
  show_prov(object)
  cat('\n')
  if(!is.null(object[])) print(head(object[], 3L))

})




## featMetaObj ####
setMethod('show', signature('featMetaObj'), function(object) {

  cat('An object of class', class(object), '\n')
  show_spat_and_feat(object)
  show_prov(object)
  cat('\n')
  if(!is.null(object[])) print(head(object[], 3L))

})











## dimObj ####

#' Show method for dimObj class
#' @param object dimension reduction object
#' @aliases show,dimObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('dimObj'), function(object) {

    show_class_and_name(object)
    if(!is.null(object@reduction_method)) cat('--| Contains dimension reduction generated with:', object@reduction_method, '\n')
    if(!is.null(object@feat_type) & !is.null(object@spat_unit)) {
      cat('----| for feat_type:', object@feat_type, '\n')
      cat('----|     spat_unit:', object@spat_unit, '\n\n')
    }

    if(!is.null(object@coordinates)) cat('  ', ncol(object@coordinates), 'dimensions for', nrow(object@coordinates),'data points\n\n')

    if(!is.null(object@misc)) {
      cat('Additional included info:\n')
      print(names(object@misc))
      cat('\n')
    }

  })










## nnNetObj ####
#' Show method for nnNetObj class
#' @param object nearest neigbor network object
#' @aliases show,nnNetObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('nnNetObj'), function(object) {

    show_class_and_name(object)
    if(!is.null(object@nn_type)) cat('--| Contains nearest neighbor network generated with:', object@nn_type, '\n')
    if(!is.null(object@feat_type) & !is.null(object@spat_unit)) {
      cat('----| for feat_type:', object@feat_type, '\n')
      cat('----|     spat_unit:', object@spat_unit, '\n')
    }
    if(!is.null(object@provenance)) cat('----|     provenance:', object@provenance, '\n\n')

    if(!is.null(object@igraph)) {
      print(object@igraph)
      cat('\n\n')
    }

    if(!is.null(object@misc)) {
      cat('Additional included info:\n')
      print(names(object@misc))
      cat('\n')
    }

  })








## spatLocsObj ####

#' show method for spatLocsObj class
#' @param object spatial locations object
#' @aliases show,spatLocsObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('spatLocsObj'), function(object) {

    sdimx = sdimy = NULL

    show_class_and_name(object)
    show_spat(object)
    show_prov(object)

    cat('   ------------------------\n\npreview:\n')
    if(!is.null(slot(object, 'coordinates'))) show(head(slot(object, 'coordinates'), 3L))

    cat('\nranges:\n')
    try(expr = print(sapply(slot(object, 'coordinates')[,.(sdimx,sdimy)], range)),
        silent = TRUE)

    cat('\n\n')

  })










## spatialNetworkObj ####

#' show method for spatialNetworkObj class
#' @param object spatial network object
#' @aliases show,spatialNetworkObj-method
#' @docType methods
#' @importFrom methods show
#' @importFrom graphics segments

#' @rdname show-methods
setMethod(
  f = "show", signature('spatialNetworkObj'), function(object) {

    show_class_and_name(object)
    if(!is.na(object@method)) cat('Contains spatial network generated with:', object@method, '\n')
    show_spat(object)
    show_prov(object)

    if(!is.null(object@networkDT)) cat('  ', nrow(object@networkDT), 'connections (filtered)\n')
    if(!is.null(object@networkDT_before_filter)) cat('  ', nrow(object@networkDT_before_filter), 'connections (before filter)\n\n')

  })













## spatialGridObj ####

#' show method for spatialGridObj class
#' @param object spatial grid object
#' @aliases show,spatialGridObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('spatialGridObj'), function(object) {

    # define for data.table
    x_start = x_end = y_start = y_end = z_start = z_end = NULL

    show_class_and_name(object)
    cat('Contains annotations for spatial unit: "', slot(object, 'spat_unit'),'"', sep = '')
    if(!is.na(slot(object, 'feat_type'))) {
      cat(' and feature type: "', slot(object, 'feat_type'), '"\n', sep = '')
    } else cat('\n')

    # find grid spatial extent
    gridNames = colnames(slot(object, 'gridDT'))
    sdimx_max = slot(object, 'gridDT')[, max(x_start, x_end)]
    sdimx_min = slot(object, 'gridDT')[, min(x_start, x_end)]
    sdimy_max = slot(object, 'gridDT')[, max(y_start, y_end)]
    sdimy_min = slot(object, 'gridDT')[, min(y_start, y_end)]
    sdimx_uniques = slot(object, 'gridDT')[, length(unique(x_start))]
    sdimy_uniques = slot(object, 'gridDT')[, length(unique(y_start))]
    cat('Contains spatial grid defined for:\n  ', sdimx_uniques,'intervals from x range:',
        sdimx_min, 'to', sdimx_max, '\n  ' ,sdimy_uniques,'intervals from y range:',
        sdimy_min, 'to', sdimy_max)
    if('z_start' %in% gridNames & 'z_end' %in% gridNames) {
      sdimz_max = slot(object, 'gridDT')[, max(z_start, z_end)]
      sdimz_min = slot(object, 'gridDT')[, min(z_start, z_end)]
      sdimz_uniques = slot(object, 'gridDT')[, length(unique(z_start))]
      cat('\n  ', sdimz_uniques, 'intervals from z range:', sdimz_min, 'to', sdimz_max, '\n\n')
    } else {
      cat('\n\n')
    }

    if(!is.null(slot(object, 'method'))) cat('Contains spatial grid generated with:', slot(object, 'method'), '\n\n')

    if(!is.null(slot(object, 'parameters'))) {
      cat('Parameters used:\n')
      for(param in names(slot(object, 'parameters'))) {
        cat(paste0('  ',param, ': ', slot(object, 'parameters')[[param]], '\n'))
      }
      cat('\n')
    }

    if(!is.null(slot(object, 'misc'))) {
      cat('Additional included info:\n')
      print(names(slot(object, 'misc')))
      cat('\n')
    }
  })












## spatEnrObj ####

#' show method for spatEnrObj class
#' @param object spatial locations object
#' @aliases show,spatEnrObj-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods
setMethod(
  f = "show", signature('spatEnrObj'), function(object) {

    show_class_and_name(object)
    show_spat_and_feat(object)
    show_prov(object)

    cat('   ------------------------\n\npreview:\n')
    if(!is.null(slot(object, 'enrichDT'))) {
      enr_cols = ncol(slot(object, 'enrichDT'))
      if(enr_cols > 10L) {
        show(slot(object, 'enrichDT')[1:3, 1:10])
        cat(rep(' ', times = getOption('width')/2.5 - 10L), rep('.', 20L),'\n', sep = '')
        show(slot(object, 'enrichDT')[1:3, 'cell_ID'])
        cat('...', enr_cols - 11L, ' cols omitted\n', sep = '')
      } else {
        show(slot(object, 'enrichDT')[1:3])
      }
    }

    cat('\n...first 20 remaining colnames:\n')
    cat('\n', wrap_txt(head(colnames(slot(object, 'enrichDT'))[-c(1:10)], 20L), strWidth = 40L), '\n')

    cat('\n\n')

  })












## giottoPolygon ####
setMethod('show', signature = 'giottoPolygon', function(object) {

  cat('An object of class giottoPolygon\n')
  show_spat(object)
  cat('Spatial Information:\n')
  print(object@spatVector)

  if(!is.null(object@spatVectorCentroids)) {
    cat(' centroids   : calculated\n')
  } else {
    cat(' centroids   : NULL\n')
  }

  if(!is.null(object@overlaps)) {
    cat(' overlaps    : calculated')
  } else {
    cat(' overlaps    : NULL')
  }

})






## packedGiottoPolygon ####
setMethod("show", signature(object='packedGiottoPolygon'),
          function(object) {
            print(paste("This is a", class(object), "object. Use 'Giotto::unwrap()' to unpack it"))
          }
)









## giottoPoints ####
setMethod('show', signature = 'giottoPoints', function(object) {

  cat('An object of class giottoPoints\n')
  show_feat(object)
  cat('Feature Information:\n')
  print(object@spatVector)

  if(!is.null(object@networks)) {
    cat(' feat. net.  :')
    print(object@networks)
  }

})







## packedGiottoPoints ####
setMethod("show", signature(object='packedGiottoPoints'),
          function(object) {
            print(paste("This is a", class(object), "object. Use 'Giotto::unwrap()' to unpack it"))
          }
)











## giottoImage ####

#' show method for giottoImage class
#' @param object giottoImage object
#' @aliases show,giottoImage-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giottoImage",
  definition = function(object) {

    cat("An object of class '",  class(object), "' with name ", object@name, "\n \n")

    cat("Min and max values are: \n",
        "Max on x-axis: ", object@minmax[['xmax_sloc']], "\n",
        "Min on x-axis: ", object@minmax[['xmin_sloc']], "\n",
        "Max on y-axis: ", object@minmax[['ymax_sloc']], "\n",
        "Min on y-axis: ", object@minmax[['ymin_sloc']], "\n",
        "\n")

    cat("Boundary adjustment are: \n",
        "Max adjustment on x-axis: ", object@boundaries[['xmax_adj']], "\n",
        "Min adjustment on x-axis: ", object@boundaries[['xmin_adj']], "\n",
        "Max adjustment on y-axis: ", object@boundaries[['ymax_adj']], "\n",
        "Min adjustment on y-axis: ", object@boundaries[['ymin_adj']], "\n",
        "\n")

    cat("Boundaries are: \n",
        "Image x-axis max boundary: ", object@minmax[['xmax_sloc']] + object@boundaries[['xmax_adj']], "\n",
        "Image x-axis min boundary: ", object@minmax[['xmin_sloc']] - object@boundaries[['xmin_adj']], "\n",
        "Image y-axis max boundary: ", object@minmax[['ymax_sloc']] + object@boundaries[['ymax_adj']], "\n",
        "Image y-axis min boundary: ", object@minmax[['ymin_sloc']] - object@boundaries[['ymin_adj']], "\n",
        "\n")

    cat("Scale factor: \n")
    print(object@scale_factor)

    cat("\n Resolution: \n")
    print(object@resolution)

    cat("\n File Path: \n")
    print(object@file_path)

    # print(object@mg_object)

  }
)










# giottoLargeImage ####

#' show method for giottoLargeImage class
#' @param object giottoLargeImage object
#' @aliases show,giottoLargeImage-method
#' @docType methods
#' @importFrom methods show
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giottoLargeImage",
  definition = function(object) {

    e = ext(object)
    img_dim = dim(object)[c(2,1,3)] # x, y, layers
    x_scalefactor = diff(e[c(1,2)]) / img_dim[1]
    y_scalefactor = diff(e[c(3,4)]) / img_dim[2]

    show_class_and_name(object)
    cat("Image extent            :", show_ext(object))
    cat("Original image extent   :", show_ext(object@overall_extent))
    cat('Scale factor            :', paste(x_scalefactor, y_scalefactor, sep = ', '), '(x, y)\n')
    cat("Resolution              :", paste(1/x_scalefactor, 1/y_scalefactor, sep = ', '), '(x, y)\n')
    cat('Layers                  :' , img_dim[3], '\n')
    cat('Estimated max intensity :', object@max_intensity, '\n')
    cat('Estimated min intensity :', object@min_intensity, '\n')
    if(object@is_int == TRUE) cat('Values                  : integers\n')
    if(object@is_int == FALSE) cat('Values                  : floating point\n')
    cat(paste0('File path               : \'', object@file_path, '\'\n'))

  }
)








# show helpers ####

#' @noRd
show_class_and_name = function(object) {
  cat("An object of class ",  class(object), ' : \"', objName(object), '\"\n', sep = '')
}

#' @noRd
show_spat_and_feat = function(object) {
  show_spat(object)
  show_feat(object)
  # cat(paste0('for spatial unit: "', spatUnit(object), '" and feature type: "', featType(object),'" \n'))
}

#' @noRd
show_spat = function(object) {
  cat(paste0('spat_unit : "', spatUnit(object), '\"\n'))
}

#' @noRd
show_feat = function(object) {
  cat(paste0('feat_type : "', featType(object), '\"\n'))
}

#' @noRd
show_prov = function(object) {
  if(!is.null(object@provenance)) cat('provenance:', object@provenance, '\n')
}

#' @noRd
show_ext = function(object) {
  paste0(paste0(ext(object)[], collapse = (', ')), ' (xmin, xmax, ymin, ymax)\n')
}




