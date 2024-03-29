

setClass(
    "CosmxReader",
    slots = list(
        cosmx_dir = "character",
        slide = "numeric",
        fovs = "numeric",
        mm = "logical",
        px2mm = "numeric",
        offsets = "ANY",
        calls = "list"
    ),
    prototype = list(
        slide = 1,
        mm = FALSE,
        px2mm = 0.12028, # from cosmx output help files
        offsets = NULL,
        calls = list()
    )
)

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
#' @param mm logical. Whether to scale spatial information as millimeters
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
        cosmx_dir = NULL, slide = 1, fovs = NULL, mm = FALSE, px2mm = 0.12028
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
    a$mm <- mm
    a$px2mm <- px2mm

    do.call(new, args = a)
}

setMethod("initialize", signature("CosmxReader"), function(
        .Object, cosmx_dir, slide, fovs, mm, px2mm
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
    if (!missing(mm)) {
        .Object@mm <- mm
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
    .detect_in_dir <- function(pattern) {
        f <- list.files(p, pattern = pattern, full.names = TRUE)
        lenf <- length(f)
        if (lenf == 1L) return(f)
        else if (lenf == 0L) {
            warning(pattern, " not detected in CosMx directory", call. = FALSE)
            return(NULL)
        }
        return(f[[1L]]) # more than one match
    }

    shifts_path <- .detect_in_dir("fov_positions_file")
    meta_path <- .detect_in_dir("metadata_file")
    tx_path <- .detect_in_dir("tx_file")
    mask_dir <- .detect_in_dir("CellLabels")
    expr_path <- .detect_in_dir("exprMat_file")
    composite_img_dir <- .detect_in_dir("CellComposite")
    overlay_img_dir <- .detect_in_dir("CellOverlay")
    compart_img_dir <- .detect_in_dir("CompartmentLabels")


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
                flip_loc_y = FALSE
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
            mm = .Object@mm,
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
        verbose = NULL
    ) {
        .cosmx_poly(
            path = path,
            fovs = .Object@fovs %none% NULL,
            flip_vertical = flip_vertical,
            flip_horizontal = flip_horizontal,
            shift_vertical_step = shift_vertical_step,
            shift_horizontal_step = shift_horizontal_step,
            mm = .Object@mm,
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
            mm = .Object@mm,
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
            "CenterY_global_px"
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

            cx[] <- c[][cell_ID %in% allowed_ids,]
            g <- setGiotto(g, cx)
        }

        return(g)
    }
    .Object@calls$create_gobject <- gobject_fun

    return(.Object)
})





# access ####

#' @export
setMethod("$", signature("CosmxReader"), function(x, name) {
    basic_info <- c("cosmx_dir", "slide", "fovs", "mm", "px2mm", "offsets")
    if (name %in% basic_info) return(methods::slot(x, name))

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("CosmxReader"), function(x, name, value) {
    basic_info <- c("cosmx_dir", "slide", "fovs", "mm", "px2mm")
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
    dn <- c("cosmx_dir", "slide", "fovs", "mm", "px2mm", "offsets")
    if (length(methods::slot(x, "calls")) > 0) {
        dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
    }
    return(dn)
}


# show ####
setMethod("show", signature("CosmxReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "CosmxReader"))
    print_slots <- c("dir", "slide", "fovs", "mm", "offsets", "funs")
    pre <- sprintf(
        "%s :", format(print_slots)
    )
    names(pre) <- print_slots

    # dir
    d <- object@cosmx_dir
    if (length(d) > 0L) {
        nch <- nchar(d)
        if (nch > 60L) {
            d1 <- substring(d, first = 0L, last = 15L)
            d2 <- substring(d, first = nch - 35L, last = nch)
            d <- paste0(d1, "[...]", d2)
        }
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

    # mm scaling
    mm <- ifelse(object@mm, object@px2mm, FALSE)
    cat(pre["mm"], mm, "\n")

    # offsets
    offs_status <- ifelse(nrow(object@offsets) > 0L, "found", "none")
    cat(pre["offsets"], offs_status, "\n")

    # funs
    nfun <- length(object@calls)
    funs <- names(object@calls)
    if (nfun > 0L) {
        pre_funs <- format(c(pre["funs"], rep("", nfun - 1L)))
        for (i in seq_len(nfun)) {
            cat(pre_funs[i], " ", funs[i], "()\n", sep = "")
        }
    }
})

setMethod("print", signature("CosmxReader"), function(x, ...) show(x))



