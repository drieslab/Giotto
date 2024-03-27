

setClass(
    "CosmxReader",
    slots = list(
        cosmx_dir = "character",
        fovs = "numeric",
        offsets = "data.frame",
        calls = "list"
    ),
    prototype = list(
        calls = list()
    )
)

#' @title Import a CosMx Assay
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
#' @param fovs numeric. (optional) If provided, will load specific fovs.
#' Otherwise, all FOVs will be loaded
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
importCosMx <- function(cosmx_dir = NULL, fovs = NULL) {
    # get params
    a <- list(Class = "CosmxReader")
    if (!is.null(cosmx_dir)) {
        a$cosmx_dir <- cosmx_dir
    }
    if (!is.null(fovs)) {
        a$fovs <- fovs
    }

    do.call(new, args = a)
}

setMethod("initialize", signature("CosmxReader"), function(.Object, cosmx_dir, fovs) {

    if (!missing(cosmx_dir)) {
        checkmate::assert_directory_exists(cosmx_dir)
        .Object@cosmx_dir <- cosmx_dir
    }
    if (!missing(fovs)) {
        checkmate::assert_numeric(fovs)
        .Object@fovs <- fovs
    }

    if (length(.Object@cosmx_dir) == 0) {
        return(.Object) # return early if no path given
    }

    p <- .Object@cosmx_dir
    .detect_in_dir <- function(pattern) {
        list.files(p, pattern = pattern, full.names = TRUE)
    }[[1L]]

    # detect paths and dirs
    pos_path <- .detect_in_dir("fov_positions_file")
    meta_path <- .detect_in_dir("metadata_file")
    tx_path <- .detect_in_dir("tx_file")
    mask_dir <- .detect_in_dir("CellLabels")
    expr_path <- .detect_in_dir("exprMat_file")
    composite_img_path <- .detect_in_dir("CellComposite")
    overlay_img_path <- .detect_in_dir("CellOverlay")
    compart_img_path <- .detect_in_dir("CompartmentLabels")


    # load fov offsets through one of several methods if not already existing
    if (nrow(.Object@offsets) == 0L) {
        if (!is.null(pos_path)) {
            pos <- data.table::fread(pos_path)
            data.table::setnames(pos, new = c("fov", "x", "y"))
        }
        else if (!is.null(meta_path)) {
            pos <- .cosmx_infer_fov_shifts(
                meta_dt = data.table::fread(meta_path),
                flip_loc_y = FALSE
            )
        } else if (!is.null(tx_path)) {
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
        gpoints_params = list(
            feat_type = c("rna", "negprobes"),
            split_keyword = list("NegPrb")
        ),
        verbose = NULL
    ) {
        .cosmx_transcript(
            path = path,
            fovs = .Object@fovs %none% NULL,
            gpoints_params = gpoints_params,
            cores = determine_cores(),
            verbose = verbose
        )
    }
    .Object@calls$load_transcripts <- tx_fun



    # mask load call
    mask_fun <- function(
        path = mask_dir,
        mask_params = list(
            # VERTICAL FLIP + NO VERTICAL SHIFT
            flip_vertical = TRUE,
            flip_horizontal = FALSE,
            shift_vertical_step = FALSE,
            shift_horizontal_step = FALSE,
            ID_fmt = NULL
        ),
        verbose = NULL
    ) {
        .cosmx_poly(
            path = path,
            fovs = .Object@fovs %none% NULL,
            mask_params = mask_params,
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
        path = composite_img_path,
        img_name_fmt = "composite_fov%03d",
        negative_y = FALSE,
        flip_vertical = FALSE,
        flip_horizontal = FALSE,
        verbose = NULL
    ) {
        .cosmx_image(
            path = path,
            fovs = .Object@fovs %none% NULL,
            img_name_fmt = img_name_fmt,
            negative_y = negative_y,
            flip_vertical = flip_vertical,
            flip_horizontal = flip_horizontal,
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
        load_images = list(
            composite = "composite",
            overlay = "overlay"
        ),
        load_expression = FALSE,
        load_cellmeta = FALSE
    ) {
        if (!is.null(load_images)) {
            checkmate::assert_list(load_images)
            if (is.null(names(load_images))) {
                stop("Images directories provided to 'load_images' must be named")
            }
        }
        g <- giotto()

        tx_list <- .Object@calls$load_transcripts()
        polys <- .Object@calls$load_polys()

        if (!is.null(load_images)) {
            # convenient shortnames
            load_images[load_images == "composite"] <- composite_img_path
            load_images[load_images == "overlay"] <- overlay_img_path

            imglist <- list()
            dirnames <- names(load_images)
            for (imdir_i in seq_along(load_images)) {
                dir_imgs <- .Object@calls$load_images(
                    path = load_images[[imdir_i]],
                    img_name_fmt = paste0(dirnames, "_fov%03d")
                )
                imglist <- c(imglist, dir_imgs)
            }
        }

        g <- setGiotto(g, gpoly)
        for (tx in tx_list) {
            g <- setGiotto(g, tx)
        }
        g@largeImages <- imglist

        # TODO expression & meta
        # Will need to check that names agree for poly/expr/meta

        return(g)
    }
    .Object@calls$create_gobject <- gobject_fun

    return(.Object)
})

#' @export
setMethod("$", signature("CosmxReader"), function(x, name) {
    basic_info <- c("offsets", "fovs", "cosmx_dir")
    if (name %in% basic_info) return(methods::slot(x, name))

    return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("CosmxReader"), function(x, name, value) {
    basic_info <- c("offsets", "fovs", "cosmx_dir")
    if (name %in% basic_info) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }

    stop(sprintf("Only items in '%s' can be set",
                 paste0(basic_info, collapse = "', '")))
})

#' @export
`.DollarNames.CosmxReader` <- function(x, pattern) {
    basic_info <- c("offsets", "fovs", "cosmx_dir")
    c(basic_info, paste0(names(methods::slot(x, "calls")), "()"))
}

setMethod("show", signature("CosmxReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", "CosmxReader"))
    pre <- sprintf("%s :", format(c("dir", "fovs", "offsets", "funs")))
    d <- object@cosmx_dir
    nch <- nchar(d)
    if (nch > 60) {
        d1 <- substring(d, first = 0L, last = 10L)
        d2 <- substring(d, first = nch - 40, last = nch)
        d <- paste0(d1, "[...]", d2)
    }
    cat(pre[1], d, "\n")
    fovs <- object@fovs %none% "all"
    cat(pre[2], paste(fovs, collapse = ", "), "\n")
    offs_status <- ifelse(nrow(object@offsets) > 0L, "found", "none")
    cat(pre[3], offs_status, "\n")

    nfun <- length(object@calls)
    funs <- names(object@calls)
    if (nfun > 0L) {
        pre_funs <- format(c(pre[4], rep("", nfun - 1L)))
        for (i in seq_len(nfun)) {
            cat(pre_funs[i], " ", funs[i], "()\n", sep = "")
        }
    }
})

setMethod("print", signature("CosmxReader"), function(x, ...) show(x))



