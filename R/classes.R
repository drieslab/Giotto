

# common internals ####
abbrev_path <- function(path, head = 15, tail = 35L) {
    nch <- nchar(path)
    if (nch > 60L) {
        p1 <- substring(path, first = 0L, last = head)
        p2 <- substring(path, first = nch - tail, last = nch)
        path <- paste0(p1, "[...]", p2)
    }
    return(path)
}

.reader_fun_prints <- function(x, pre) {
    nfun <- length(x@calls)
    funs <- names(x@calls)
    if (nfun > 0L) {
        pre_funs <- format(c(pre, rep("", nfun - 1L)))
        for (i in seq_len(nfun)) {
            cat(pre_funs[i], " ", funs[i], "()\n", sep = "")
        }
    }
}

.filetype_prints <- function(x, pre) {
    nftype <- length(x@filetype)
    datatype <- format(names(x@filetype))
    pre_ftypes <- format(c(pre, rep("", nftype - 1L)))
    cat(sprintf("%s %s -- %s\n",
                pre_ftypes,
                datatype,
                x@filetype),
        sep = "")
}

# pattern - list.files pattern to use to search for specific files/dirs
# warn - whether to warn when a pattern does not find any files
# first - whether to only return the first match
.detect_in_dir <- function(
        path, pattern, platform, warn = TRUE, first = TRUE
) {
    f <- list.files(path, pattern = pattern, full.names = TRUE)
    lenf <- length(f)
    if (lenf == 1L) return(f) # one match
    else if (lenf == 0L) { # no matches
        if (warn) {
            warning(sprintf(
                "%s not detected in %s directory",
                pattern,
                platform
            ),
            call. = FALSE)
        }
        return(NULL)
    }

    # more than one match
    if (first) {
        return(f[[1L]])
    } else {
        return(f)
    }
}










# Xenium ####

setClass(
    "XeniumReader",
    slots = list(
        xenium_dir = "character",
        filetype = "list",
        qv = "ANY",
        calls = "list"
    ),
    prototype = list(
        filetype = list(
            transcripts = "parquet",
            boundaries = "parquet",
            expression = "h5",
            cell_meta = "parquet"
        ),
        qv = 20,
        calls = list()
    )
)

# * show ####
setMethod("show", signature("XeniumReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "XeniumReader"))
    print_slots <- c("dir", "filetype", "qv_cutoff", "funs")
    pre <- sprintf(
        "%s :", format(print_slots)
    )
    names(pre) <- print_slots

    # dir
    d <- object@xenium_dir
    if (length(d) > 0L) {
        d <- abbrev_path(d)
        cat(pre["dir"], d, "\n")
    } else {
        cat(pre["dir"], "\n")
    }

    # qv
    qv <- object@qv
    cat(pre["qv_cutoff"], paste(qv, collapse = ", "), "\n")

    # filetype
    .filetype_prints(x = object, pre = pre["filetype"])

    # funs
    .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("XeniumReader"), function(x, ...) show(x))

# * init ####
setMethod(
    "initialize", signature("XeniumReader"),
    function(
        .Object,
        xenium_dir,
        filetype,
        qv_cutoff
    ) {
        .Object <- callNextMethod(.Object)

        # provided params (if any)
        if (!missing(xenium_dir)) {
            checkmate::assert_directory_exists(xenium_dir)
            .Object@xenium_dir <- xenium_dir
        }
        if (!missing(filetype)) {
            .Object@filetype <- filetype
        }
        if (!missing(qv_cutoff)) {
            .Object@qv <- qv_cutoff
        }


        # check filetype
        ftype_data <- c("transcripts", "boundaries", "expression", "cell_meta")
        if (!all(ftype_data %in% names(.Object@filetype))) {
            stop(wrap_txt("`$filetype` must have entries for each of:\n",
                          paste(ftype_data, collapse = ", ")))
        }

        ftype <- .Object@filetype
        ft_tab <- c("csv", "parquet")
        ft_exp <- c("h5", "mtx", "zarr")
        if (!ftype$transcripts %in% ft_tab) {
            stop(wrap_txt("`$filetype$transcripts` must be one of",
                          paste(ft_tab, collapse = ", ")),
                 call. = FALSE)
        }
        if (!ftype$boundaries %in% ft_tab) {
            stop(wrap_txt("`$filetype$boundaries` must be one of",
                          paste(ft_tab, collapse = ", ")),
                 call. = FALSE)
        }
        if (!ftype$cell_meta %in% ft_tab) {
            stop(wrap_txt("`$filetype$cell_meta` must be one of",
                          paste(ft_tab, collapse = ", ")),
                 call. = FALSE)
        }
        if (!ftype$expression %in% ft_exp) {
            stop(wrap_txt("`$filetype$expression` must be one of",
                          paste(ft_tab, collapse = ", ")),
                 call. = FALSE)
        }


        # detect paths and subdirs
        p <- .Object@xenium_dir
        .xenium_detect <- function(pattern, ...) {
            .detect_in_dir(
                pattern = pattern, ...,
                path = p, platform = "Xenium",
            )
        }

        cell_meta_path <- .xenium_detect("cells", first = FALSE)
        panel_meta_path <- .xenium_detect("panel") # json
        experiment_info_path <- .xenium_detect(".xenium") # json

        # 3D stack - DAPI
        img_path <- .xenium_detect("morphology.", warn = FALSE)
        # 2D fusion images
        # - DAPI
        # - stainings for multimodal segmentation
        img_focus_path <- .xenium_detect("morphology_focus", warn = FALSE)
        # Maximum intensity projection (MIP) of the morphology image.
        # (Xenium Outputs v1.0 - 1.9. only)
        img_mip_path <- .xenium_detect("morphology_mip", warn = FALSE)

        tx_path <- .xenium_detect("transcripts", first = FALSE)
        cell_bound_path <- .xenium_detect("cell_bound", first = FALSE)
        nuc_bound_path <- .xenium_detect("nucleus_bound", first = FALSE)

        expr_path <- .xenium_detect("cell_feature_matrix", first = FALSE)

        .xenium_ftype <- function(paths, ftype) {
            paths[grepl(pattern = paste0(".", ftype), x = paths)]
        }


        # select file formats based on reader settings
        tx_path <- .xenium_ftype(tx_path, ftype$transcripts)
        cell_bound_path <- .xenium_ftype(cell_bound_path, ftype$boundaries)
        nuc_bound_path <- .xenium_ftype(nuc_bound_path, ftype$boundaries)
        expr_path <- .xenium_ftype(expr_path, ftype$expression)
        cell_meta_path <- .xenium_ftype(cell_meta_path, ftype$cell_meta)


        # transcripts load call
        tx_fun <- function(
        path = tx_path,
        feat_type = c(
            "rna",
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        split_keyword = list(
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        dropcols = c(),
        qv_threshold = .Object@qv,
        cores = determine_cores(),
        verbose = NULL
        ) {
            .xenium_transcript(
                path = path,
                feat_type = feat_type,
                split_keyword = split_keyword,
                dropcols = dropcols,
                qv_threshold = qv_threshold,
                cores = cores,
                verbose = verbose
            )
        }
        .Object@calls$load_transcripts <- tx_fun

        # load polys call
        poly_fun <- function(
        path = cell_bound_path,
        name = "cell",
        calc_centroids = TRUE,
        cores = determine_cores(),
        verbose = NULL
        ) {
            .xenium_poly(
                path = path,
                name = name,
                calc_centroids = calc_centroids,
                cores = cores,
                verbose = verbose
            )
        }
        .Object@calls$load_polys <- poly_fun

        # load cellmeta
        cmeta_fun <- function(
        path = cell_meta_path,
        dropcols = c(),
        cores = determine_cores(),
        verbose = NULL
        ) {
            .xenium_cellmeta(
                path = path,
                dropcols = dropcols,
                cores = cores,
                verbose = verbose
            )
        }
        .Object@calls$load_cellmeta <- cmeta_fun

        # load featmeta
        fmeta_fun <- function(
        path = panel_meta_path,
        dropcols = c(),
        cores = determine_cores(),
        verbose = NULL
        ) {
            .xenium_featmeta(
                path = path,
                gene_ids,
                dropcols = dropcols,
                verbose = verbose
            )
        }
        .Object@calls$load_featmeta <- fmeta_fun

        # load expression call
        expr_fun <- function(
        path,
        gene_ids = "symbols",
        remove_zero_rows = TRUE,
        split_by_type = TRUE,
        verbose = NULL
        ) {
            .xenium_expression(
                path = path,
                gene_ids = gene_ids,
                remove_zero_rows = remove_zero_rows,
                split_by_type = split_by_type,
                verbose = verbose
            )
        }
        .Object@calls$load_expression <- expr_fun

        # load image call




        # create giotto object call
        gobject_fun <- function(
        transcript_path = tx_path,
        load_bounds = list(
            cell = "cell",
            nucleus = "nucleus"
        ),
        expression_path = expr_path,
        metadata_path = meta_path,
        feat_type = c(
            "rna",
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        split_keyword = list(
            "NegControlProbe",
            "UnassignedCodeword",
            "NegControlCodeword"
        ),
        load_images = list(
            morphology = "focus",
        ),
        load_expression = FALSE,
        load_cellmeta = FALSE
        ) {
            load_expression <- as.logical(load_expression)
            load_cellmeta <- as.logical(load_cellmeta)

            if (!is.null(load_images)) {
                checkmate::assert_list(load_images)
                if (is.null(names(load_images))) {
                    stop("Images paths provided to 'load_images' must be named")
                }
            }
            if (!is.null(load_bounds)) {
                checkmate::assert_list(load_bounds)
                if (is.null(names(load_bounds))) {
                    stop("bounds paths provided to 'load_bounds' must be named")
                }
            }



            funs <- .Object@calls

            # init gobject
            g <- giotto()


            # transcripts
            tx_list <- funs$load_transcripts(
                path = transcript_path,
                feat_type = feat_type,
                split_keyword = split_keyword
            )
            for (tx in tx_list) {
                g <- setGiotto(g, tx)
            }


            # polys
            if (!is.null(load_bounds)) {
                # replace convenient shortnames
                load_bounds[load_bounds == "cell"] <- cell_bound_path
                load_bounds[load_bounds == "nucleus"] <- nuc_bound_path

                blist <- list()
                bnames <- names(load_bounds)
                for (b_i in seq_along(load_bounds)) {
                    b <- funs$load_polys(
                        path = load_bounds[[b_i]],
                        name = bnames[[b_i]]
                    )
                    blist <- c(blist, b)
                }
                for (gpoly_i in seq_along(blist)) {
                    g <- setGiotto(g, blist[[gpoly_i]])
                }
            }


            # feat metadata
            fx <- funs$load_featmeta(
                path =
            )


            # expression
            if (load_expression) {

            }


            # cell metadata
            if (load_cellmeta) {

            }


            # images
            if (!is.null(load_images)) {
                # replace convenient shortnames
                load_images[load_images == "focus"] <- img_focus_path
            }




        }
        .Object@calls$create_gobject <- gobject_fun


        return(.Object)
    }
)




# access ####

#' @export
setMethod("$", signature("XeniumReader"), function(x, name) {
    basic_info <- c("xenium_dir", "filetype", "qv")
    if (name %in% basic_info) return(methods::slot(x, name))

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("XeniumReader"), function(x, name, value) {
    basic_info <- c("xenium_dir", "filetype", "qv")
    if (name %in% basic_info) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }

    stop(sprintf("Only items in '%s' can be set",
                 paste0(basic_info, collapse = "', '")))
})

#' @export
`.DollarNames.XeniumReader` <- function(x, pattern) {
    dn <- c("xenium_dir", "filetype", "qv")
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}






# CosMx ####

setClass(
    "CosmxReader",
    slots = list(
        cosmx_dir = "character",
        slide = "numeric",
        fovs = "numeric",
        micron = "logical",
        px2mm = "numeric",
        offsets = "ANY",
        calls = "list"
    ),
    prototype = list(
        slide = 1,
        micron = FALSE,
        px2mm = 0.12028, # from cosmx output help files
        offsets = NULL,
        calls = list()
    )
)

# * show ####
setMethod("show", signature("CosmxReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "CosmxReader"))
    print_slots <- c("dir", "slide", "fovs", "micron", "offsets", "funs")
    pre <- sprintf(
        "%s :", format(print_slots)
    )
    names(pre) <- print_slots

    # dir
    d <- object@cosmx_dir
    if (length(d) > 0L) {
        nch <- nchar(d)
        d <- abbrev_path(d)
        cat(pre["dir"], d, "\n")
    } else {
        cat(pre["dir"], "\n")
    }

    # slide
    slide <- object@slide
    cat(pre["slide"], slide, "\n")

    # fovs
    fovs <- object@fovs %none% "all"
    cat(pre["fovs"], paste(fovs, collapse = ", "), "\n")

    # micron scaling
    micron <- ifelse(object@micron, object@px2mm / 1000, FALSE)
    cat(pre["micron"], micron, "\n")

    # offsets
    offs_status <- ifelse(nrow(object@offsets) > 0L, "found", "none")
    cat(pre["offsets"], offs_status, "\n")

    # funs
    .fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("CosmxReader"), function(x, ...) show(x))

# * plot ####
setMethod(
    "plot", signature(x = "CosmxReader", y = "missing"),
    function(x, cex = 0.8, ...) {
        a <- list(...)
        dat <- x@offsets

        if (is.null(dat)) { # don't run if no offsets
            cat("no offsets to plot\n")
            return(invisible(NULL))
        }

        plot(y ~ x, data = dat, asp = 1L, type = "n", ...)
        text(y ~ x, data = dat, labels = dat$fov, cex = cex, ...)
    })




#' @title Import a Nanostring CosMx Assay
#' @name importCosMx
#' @description
#' Giotto import functionalities for CosMx datasets. This function generates
#' a `CosmxReader` instance that has convenient reader functions for converting
#' individual pieces of CosMx data into Giotto-compatible representations when
#' the params `cosmx_dir` and `fovs` (if only a subset is desired) are provided.
#' A function that creates the full `giotto` object is also available.
#' These functions should have all param values provided as defaults, but
#' can be flexibly modified to do things such as look in alternative
#' directories or paths.
#' @param cosmx_dir CosMx output directory
#' @param slide numeric. Slide number. Defaults to 1
#' @param fovs numeric. (optional) If provided, will load specific fovs.
#' Otherwise, all FOVs will be loaded
#' @param micron logical. Whether to scale spatial information as micron
#' instead of the default pixels
#' @param px2mm numeric. Scalefactor from pixels to mm. Defaults to 0.12028
#' based on `CosMx-ReadMe.html` info
#' @details
#' Loading functions are generated after the `cosmx_dir` is added.
#' Transcripts, expression, and metadata loading are all expected to be done
#' from the top level of the directory. Loading of polys, and any image sets
#' are expected to be from specific subdirectories containing only those
#' images for the set of FOVs.
#' @returns CosmxReader object
#' @examples
#' # Create a `CosmxReader` object
#' reader <- importCosMx()
#'
#' \dontrun{
#' # Set the cosmx_dir and fov parameters
#' reader$cosmx_dir <- "path to cosmx dir"
#' reader$fov <- c(1, 4)
#'
#' plot(reader) # displays FOVs (top left corner) in px scale.
#'
#' # Load polygons, transcripts, and images
#' polys <- reader$load_polys()
#' tx <- reader$load_transcripts()
#' imgs <- reader$load_images()
#'
#' # Create a `giotto` object and add the loaded data
#' g <- giotto()
#' g <- setGiotto(g, tx[["rna"]])
#' g <- setGiotto(g, polys)
#' g <- addGiottoLargeImage(g, largeImages = imgs)
#' force(g)
#' }
#' @export
importCosMx <- function(
        cosmx_dir = NULL, slide = 1, fovs = NULL, micron = FALSE, px2mm = 0.12028
) {
    # get params
    a <- list(Class = "CosmxReader")
    if (!is.null(cosmx_dir)) {
        a$cosmx_dir <- cosmx_dir
    }
    if (!is.null(fovs)) {
        a$fovs <- fovs
    }
    a$slide <- slide
    a$micron <- micron
    a$px2mm <- px2mm

    do.call(new, args = a)
}

# * init ####
setMethod("initialize", signature("CosmxReader"), function(
        .Object, cosmx_dir, slide, fovs, micron, px2mm
) {
    # provided params (if any)
    if (!missing(cosmx_dir)) {
        checkmate::assert_directory_exists(cosmx_dir)
        .Object@cosmx_dir <- cosmx_dir
    }
    if (!missing(slide)) {
        .Object@slide <- slide
    }
    if (!missing(fovs)) {
        .Object@fovs <- fovs
    }
    if (!missing(micron)) {
        .Object@micron <- micron
    }
    if (!missing(px2mm)) {
        .Object@px2mm <- px2mm
    }

    # NULL case
    if (length(.Object@cosmx_dir) == 0) {
        return(.Object) # return early if no path given
    }


    # detect paths and subdirs
    p <- .Object@cosmx_dir
    .cosmx_detect <- function(pattern) {
        .detect_in_dir(pattern = pattern, path = p, platform = "CosMx")
    }

    shifts_path <- .cosmx_detect("fov_positions_file")
    meta_path <- .cosmx_detect("metadata_file")
    tx_path <- .cosmx_detect("tx_file")
    mask_dir <- .cosmx_detect("CellLabels")
    expr_path <- .cosmx_detect("exprMat_file")
    composite_img_dir <- .cosmx_detect("CellComposite")
    overlay_img_dir <- .cosmx_detect("CellOverlay")
    compart_img_dir <- .cosmx_detect("CompartmentLabels")


    # load fov offsets through one of several methods
    if (is.null(.Object@offsets)) { # only run if not already existing
        pos <- NULL

        if (!is.null(shifts_path)) {
            fov_shifts <- data.table::fread(shifts_path)
            if (!"X_mm" %in% colnames(fov_shifts)) {
                # older version has fov, x, y (all numeric) in px shifts
                data.table::setnames(fov_shifts, new = c("fov", "x", "y"))
                pos <- fov_shifts
            }
        }

        # proceed with other possible methods of inferring shifts if present
        if (!is.null(meta_path) && is.null(pos)) {
            pos <- .cosmx_infer_fov_shifts(
                meta_dt = data.table::fread(meta_path),
                flip_loc_y = TRUE
            )
        } else if (!is.null(tx_path) && is.null(pos)) {
            warning(wrap_txt(
                "metadata_file not found:
                Detecting fov shifts from tx_file. (This is slower)"
            ), call. = FALSE)
            pos <- .cosmx_infer_fov_shifts(
                tx_dt = data.table::fread(tx_path),
                flip_loc_y = TRUE
            )
        }
        else {
            pos <- data.table::data.table()
            warning(wrap_txt(
                "fov_positions_file, tx_file, and metadata_file not auto detected.
                One of these must be provided to infer FOV shifts"
            ))
        }

        .Object@offsets <- pos
    }



    # transcripts load call
    tx_fun <- function(
        path = tx_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb"),
        dropcols = c(
            "x_local_px",
            "y_local_px",
            "cell_ID",
            "cell"
        ),
        verbose = NULL
    ) {
        .cosmx_transcript(
            path = path,
            fovs = .Object@fovs %none% NULL,
            feat_type = feat_type,
            split_keyword = split_keyword,
            dropcols = dropcols,
            micron = .Object@micron,
            px2mm = .Object@px2mm,
            cores = determine_cores(),
            verbose = verbose
        )
    }
    .Object@calls$load_transcripts <- tx_fun



    # mask load call
    mask_fun <- function(
        path = mask_dir,
        # VERTICAL FLIP + NO VERTICAL SHIFT
        flip_vertical = TRUE,
        flip_horizontal = FALSE,
        shift_vertical_step = FALSE,
        shift_horizontal_step = FALSE,
        remove_background_polygon = TRUE,
        verbose = NULL
    ) {
        .cosmx_poly(
            path = path,
            fovs = .Object@fovs %none% NULL,
            flip_vertical = flip_vertical,
            flip_horizontal = flip_horizontal,
            shift_vertical_step = shift_vertical_step,
            shift_horizontal_step = shift_horizontal_step,
            remove_background_polygon = remove_background_polygon,
            micron = .Object@micron,
            px2mm = .Object@px2mm,
            offsets = .Object@offsets,
            verbose = verbose
        )
    }
    .Object@calls$load_polys <- mask_fun


    # expression load call
    expr_fun <- function(
        path = expr_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list("NegPrb")
    ) {
        .cosmx_expression(
            path = path,
            fovs = .Object@fovs %none% NULL,
            feat_type = feat_type,
            split_keyword = split_keyword
        )
    }
    .Object@calls$load_expression <- expr_fun


    # images load call
    img_fun <- function(
        path = composite_img_dir,
        img_type = "composite",
        img_name_fmt = paste0(img_type, "_fov%03d"),
        negative_y = TRUE,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        verbose = NULL
    ) {
        .cosmx_image(
            path = path,
            fovs = .Object@fovs %none% NULL,
            img_type = img_type,
            img_name_fmt = img_name_fmt,
            negative_y = negative_y,
            flip_vertical = flip_vertical,
            flip_horizontal = flip_horizontal,
            micron = .Object@micron,
            px2mm = .Object@px2mm,
            offsets = .Object@offsets,
            verbose = verbose
        )
    }
    .Object@calls$load_images <- img_fun


    # meta load call
    meta_fun <- function(
        path = meta_path,
        dropcols = c(
            "CenterX_local_px",
            "CenterY_local_px",
            "CenterX_global_px",
            "CenterY_global_px",
            "cell_id"
        ),
        verbose = NULL
    ) {
        .cosmx_cellmeta(
            path = path,
            fovs = .Object@fovs %none% NULL,
            dropcols = dropcols,
            cores = determine_cores(),
            verbose = verbose
        )
    }
    .Object@calls$load_cellmeta <- meta_fun


    # build gobject call
    gobject_fun <- function(
        transcript_path = tx_path,
        cell_labels_dir = mask_dir,
        expression_path = expr_path,
        metadata_path = meta_path,
        feat_type = c("rna", "negprobes"),
        split_keyword = list(
            "NegPrb"
        ),
        load_images = list(
            composite = "composite",
            overlay = "overlay"
        ),
        load_expression = FALSE,
        load_cellmeta = FALSE
    ) {
        load_expression <- as.logical(load_expression)
        load_cellmeta <- as.logical(load_cellmeta)

        if (!is.null(load_images)) {
            checkmate::assert_list(load_images)
            if (is.null(names(load_images))) {
                stop("Images directories provided to 'load_images' must be named")
            }
        }

        funs <- .Object@calls

        # init gobject
        g <- giotto()

        # transcripts
        tx_list <- funs$load_transcripts(
            path = transcript_path,
            feat_type = feat_type,
            split_keyword = split_keyword
        )
        for (tx in tx_list) {
            g <- setGiotto(g, tx)
        }

        # polys
        polys <- funs$load_polys(
            path = cell_labels_dir,
            verbose = FALSE
        )
        g <- setGiotto(g, polys)

        # images
        if (!is.null(load_images)) {
            # replace convenient shortnames
            load_images[load_images == "composite"] <- composite_img_dir
            load_images[load_images == "overlay"] <- overlay_img_dir

            imglist <- list()
            dirnames <- names(load_images)
            for (imdir_i in seq_along(load_images)) {
                dir_imgs <- funs$load_images(
                    path = load_images[[imdir_i]],
                    img_type = dirnames[[imdir_i]],
                )
                imglist <- c(imglist, dir_imgs)
            }
            g <- addGiottoLargeImage(g, largeImages = imglist)
        }

        # expression & meta
        # Need to check that names agree for poly/expr/meta
        allowed_ids <- spatIDs(polys)

        if (load_expression) {
            exlist <- funs$load_expression(
                path = expression_path,
                feat_type = feat_type,
                split_keyword = split_keyword
            )

            # only keep allowed cells and set into gobject
            for (ex in exlist) {
                bool <- colnames(ex[]) %in% allowed_ids
                ex[] <- ex[][, bool]
                g <- setGiotto(g, ex)
            }
        }

        if (load_cellmeta) {
            cx <- funs$load_cellmeta(
                path = metadata_path
            )

            cx[] <- cx[][cell_ID %in% allowed_ids,]
            g <- setGiotto(g, cx)
        }

        return(g)
    }
    .Object@calls$create_gobject <- gobject_fun

    return(.Object)
})





# * access ####

#' @export
setMethod("$", signature("CosmxReader"), function(x, name) {
    basic_info <- c("cosmx_dir", "slide", "fovs", "micron", "px2mm", "offsets")
    if (name %in% basic_info) return(methods::slot(x, name))

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("CosmxReader"), function(x, name, value) {
    basic_info <- c("cosmx_dir", "slide", "fovs", "micron", "px2mm")
    if (name %in% basic_info) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }

    if (name == "offsets") {
        methods::slot(x, name) <- data.table::setDT(value)
        return(initialize(x))
    }

    stop(sprintf("Only items in '%s' can be set",
                 paste0(basic_info, collapse = "', '")))
})

#' @export
`.DollarNames.CosmxReader` <- function(x, pattern) {
    dn <- c("cosmx_dir", "slide", "fovs", "micron", "px2mm", "offsets")
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}




