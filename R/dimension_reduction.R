# * Dimension Reduction Object Creation #
# ! Moved to classes.R



## * PCA  ####
# ---------- #



#' @title Run PCA - factominer
#' @name .run_pca_factominer
#' @description performs PCA based on the factominer package
#' @param x matrix or object that can be converted to matrix
#' @param ncp number of principal components to calculate
#' @param scale scale features
#' @param rev reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @keywords internal
#' @returns list of eigenvalues, loadings and pca coordinates
.run_pca_factominer <- function(x,
    ncp = 100,
    scale = TRUE,
    rev = FALSE,
    set_seed = TRUE,
    seed_number = 1234,
    ...) {
    # verify if optional package is installed
    package_check(pkg_name = "FactoMineR", repository = "CRAN")

    if (!methods::is(x, "matrix")) {
        x <- as.matrix(x)
    }

    if (isTRUE(rev)) {
        x <- t_flex(x)

        if (ncp > nrow(x)) {
            warning("ncp > nrow(x), will be set to nrow(x)")
            ncp <- nrow(x)
        }

        # start seed
        if (isTRUE(set_seed)) {
            set.seed(seed = seed_number)
        }

        pca_res <- FactoMineR::PCA(
            X = x, ncp = ncp, scale.unit = scale, graph = FALSE, ...)

        # exit seed
        if (isTRUE(set_seed)) {
            set.seed(seed = Sys.time())
        }

        # eigenvalues
        eigenvalues <- pca_res$eig[, 1]

        # PC loading
        loadings <- pca_res$ind$coord
        rownames(loadings) <- rownames(x)
        colnames(loadings) <- paste0("Dim.", 1:ncol(loadings))

        # coordinates
        coords <- sweep(pca_res$var$coord, 
                        2, sqrt(eigenvalues[1:ncp]), FUN = "/")
        rownames(coords) <- colnames(x)
        colnames(coords) <- paste0("Dim.", 1:ncol(coords))

        result <- list(
            eigenvalues = eigenvalues, loadings = loadings, coords = coords)
    } else {
        if (ncp > ncol(x)) {
            warning("ncp > ncol(x), will be set to ncol(x)")
            ncp <- ncol(x)
        }

        # start seed
        if (isTRUE(set_seed)) {
            set.seed(seed = seed_number)
        }

        pca_res <- FactoMineR::PCA(
            X = x, ncp = ncp, scale.unit = scale, graph = FALSE, ...)

        # exit seed
        if (isTRUE(set_seed)) {
            set.seed(seed = Sys.time())
        }

        # eigenvalues
        eigenvalues <- pca_res$eig[, 1]

        # PC loading
        loadings <- sweep(
            pca_res$var$coord, 2, sqrt(eigenvalues[1:ncp]), FUN = "/")
        rownames(loadings) <- colnames(x)
        colnames(loadings) <- paste0("Dim.", 1:ncol(loadings))

        # coordinates
        coords <- pca_res$ind$coord
        rownames(coords) <- rownames(x)
        colnames(coords) <- paste0("Dim.", 1:ncol(coords))

        result <- list(
            eigenvalues = eigenvalues, loadings = loadings, coords = coords)
    }

    vmsg(.is_debug = TRUE, 
        "finished .run_pca_factominer, method == factominer")

    return(result)
}


#' @title Run PCA - BiocSingular
#' @name .run_pca_biocsingular
#' @description Performs PCA based on the biocSingular package
#' @param x matrix or object that can be converted to matrix
#' @param ncp number of principal components to calculate
#' @param center center the matrix before pca
#' @param scale scale features
#' @param rev reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param BSPARAM method to use
#' @param BPPARAM BiocParallelParam object
#' @keywords internal
#' @returns list of eigenvalues, loadings and pca coordinates
.run_pca_biocsingular <- function(x,
    ncp = 100,
    center = TRUE,
    scale = TRUE,
    rev = FALSE,
    set_seed = TRUE,
    seed_number = 1234,
    BSPARAM = c("irlba", "exact", "random"),
    BPPARAM = BiocParallel::SerialParam(),
    ...) {
    BSPARAM <- match.arg(BSPARAM, choices = c("irlba", "exact", "random"))

    min_ncp <- min(dim(x))

    if (ncp >= min_ncp) {
        warning("ncp >= minimum dimension of x, will be set to 
                minimum dimension of x - 1")
        ncp <- min_ncp - 1
    }

    # start seed
    if (isTRUE(set_seed)) {
        set.seed(seed = seed_number)
    }

    if (isTRUE(rev)) {
        x <- t_flex(x)

        if (BSPARAM == "irlba") {
            pca_res <- BiocSingular::runPCA(
                x = x, rank = ncp,
                center = center, scale = scale,
                BSPARAM = BiocSingular::IrlbaParam(),
                BPPARAM = BPPARAM,
                ...
            )
        } else if (BSPARAM == "exact") {
            pca_res <- BiocSingular::runPCA(
                x = x, rank = ncp,
                center = center, scale = scale,
                BSPARAM = BiocSingular::ExactParam(),
                BPPARAM = BPPARAM,
                ...
            )
        } else if (BSPARAM == "random") {
            pca_res <- BiocSingular::runPCA(
                x = x, rank = ncp,
                center = center, scale = scale,
                BSPARAM = BiocSingular::RandomParam(),
                BPPARAM = BPPARAM,
                ...
            )
        }



        # eigenvalues
        eigenvalues <- pca_res$sdev^2
        # PC loading
        loadings <- pca_res$x
        rownames(loadings) <- rownames(x)
        colnames(loadings) <- paste0("Dim.", 1:ncol(loadings))
        # coordinates
        coords <- pca_res$rotation
        rownames(coords) <- colnames(x)
        colnames(coords) <- paste0("Dim.", 1:ncol(coords))
        result <- list(
            eigenvalues = eigenvalues, loadings = loadings, coords = coords)
    } else {
        if (BSPARAM == "irlba") {
            pca_res <- BiocSingular::runPCA(
                x = x, rank = ncp,
                center = center, scale = scale,
                BSPARAM = BiocSingular::IrlbaParam(),
                BPPARAM = BPPARAM,
                ...
            )
        } else if (BSPARAM == "exact") {
            pca_res <- BiocSingular::runPCA(
                x = x, rank = ncp,
                center = center, scale = scale,
                BSPARAM = BiocSingular::ExactParam(),
                BPPARAM = BPPARAM,
                ...
            )
        } else if (BSPARAM == "random") {
            pca_res <- BiocSingular::runPCA(
                x = x, rank = ncp,
                center = center, scale = scale,
                BSPARAM = BiocSingular::RandomParam(),
                BPPARAM = BPPARAM,
                ...
            )
        }

        # eigenvalues
        eigenvalues <- pca_res$sdev^2
        # PC loading
        loadings <- pca_res$rotation
        rownames(loadings) <- colnames(x)
        colnames(loadings) <- paste0("Dim.", 1:ncol(loadings))
        # coordinates
        coords <- pca_res$x
        rownames(coords) <- rownames(x)
        colnames(coords) <- paste0("Dim.", 1:ncol(coords))
        result <- list(
            eigenvalues = eigenvalues, loadings = loadings, coords = coords)
    }

    # exit seed
    if (isTRUE(set_seed)) {
        set.seed(seed = Sys.time())
    }

    vmsg(.is_debug = TRUE, "finished .run_pca_biocsingular, method ==", BSPARAM)

    return(result)
}






#' @title Create features to use matrix
#' @name .create_feats_to_use_matrix
#' @description subsets matrix based on vector of genes or hvf column
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param sel_matrix selected expression matrix
#' @param feats_to_use feats to use, character or vector of features
#' @param verbose verbosity
#' @keywords internal
#' @returns subsetted matrix based on selected features
.create_feats_to_use_matrix <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    sel_matrix,
    feats_to_use,
    verbose = FALSE) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # cell metadata
    feat_metadata <- fDataDT(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # for hvf features
    if (is.character(feats_to_use) && length(feats_to_use) == 1) {
        if (feats_to_use %in% colnames(feat_metadata)) {
            vmsg(
                .v = verbose, str_double_quote(feats_to_use),
                "column was found in the feats metadata information and will be 
                used to select highly variable features"
            )
            feats_to_use <- feat_metadata[
                get(feats_to_use) == "yes"][["feat_ID"]]
            sel_matrix <- sel_matrix[rownames(sel_matrix) %in% feats_to_use, ]
        } else {
            vmsg(
                .v = verbose, str_double_quote(feats_to_use),
                "was not found in the gene metadata information.
                all genes will be used."
            )
        }
    } else {
        vmsg(.v = verbose, 
            "a custom vector of genes will be used to subset the matrix")
        sel_matrix <- sel_matrix[rownames(sel_matrix) %in% feats_to_use, ]
    }

    vmsg(.v = verbose, .is_debug = TRUE, 
        "class of selected matrix: ", class(sel_matrix))

    return(sel_matrix)
}



#' @title runPCA
#' @name runPCA
#' @description runs a Principal Component Analysis
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param name arbitrary name for PCA run
#' @param feats_to_use subset of features to use for PCA
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param center center data first (default = TRUE)
#' @param scale_unit scale features before PCA (default = TRUE)
#' @param ncp number of principal components to calculate (see details)
#' @param method which implementation to use
#' @param method_params BiocParallelParam object
#' @param rev do a reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param verbose verbosity of the function
#' @param ... additional parameters for PCA (see details)
#' @returns giotto object with updated PCA dimension recuction
#' @details See \code{\link[BiocSingular]{runPCA}} and
#' \code{\link[FactoMineR]{PCA}} for more information about other parameters.
#' With the feats_to_use param, you can control which features are used to
#' calculate your PCA, which can be useful for making sure that your downstream
#' dimension reduction and clusterings are based on your features of interest.
#' \itemize{
#'   \item feats_to_use = NULL: will use all features from the selected matrix
#'   \item feats_to_use = <hvg name>: can be used to select a column name of
#'   highly variable features, created by (see \code{\link{calculateHVF}})
#'   \item feats_to_use = c('geneA', 'geneB', ...): will use all manually 
#'   provided features
#' }
#' By default the number of principle components that we calculate is 100, which
#' may not encompass all the variation within the dataset. Setting ncp to NULL
#' will calculate all the principle components.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' runPCA(g)
#' @export
runPCA <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    name = NULL,
    feats_to_use = "hvf",
    return_gobject = TRUE,
    center = TRUE,
    scale_unit = TRUE,
    ncp = 100,
    method = c("irlba", "exact", "random", "factominer"),
    method_params = BiocParallel::SerialParam(),
    rev = FALSE,
    set_seed = TRUE,
    seed_number = 1234,
    verbose = TRUE,
    ...) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # specify name to use for pca
    if (is.null(name)) {
        if (feat_type == "rna") {
            name <- "pca"
        } else {
            name <- paste0(feat_type, ".", "pca")
        }
    }

    # expression values to be used
    values <- match.arg(
        expression_values, 
        unique(c("normalized", "scaled", "custom", expression_values)))
    expr_values <- getExpression(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        values = values,
        output = "exprObj"
    )

    provenance <- prov(expr_values)

    expr_values <- expr_values[] # extract matrix

    # set max ncp if NULL was provided
    if (is.null(ncp)) {
        ncp <- nrow(expr_values)
    }

    ## subset matrix
    if (!is.null(feats_to_use)) {
        expr_values <- .create_feats_to_use_matrix(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sel_matrix = expr_values,
            feats_to_use = feats_to_use,
            verbose = verbose
        )
    }


    # do PCA dimension reduction
    reduction <- match.arg(reduction, c("cells", "feats"))

    # PCA implementation
    method <- match.arg(method, c("irlba", "exact", "random", "factominer"))

    if (reduction == "cells") {
        # PCA on cells
        if (method %in% c("irlba", "exact", "random")) {
            pca_object <- .run_pca_biocsingular(
                x = t_flex(expr_values),
                center = center,
                scale = scale_unit,
                ncp = ncp,
                rev = rev,
                set_seed = set_seed,
                seed_number = seed_number,
                BSPARAM = method,
                BPPARAM = method_params,
                ...
            )
        } else if (method == "factominer") {
            pca_object <- .run_pca_factominer(
                x = t_flex(expr_values),
                scale = scale_unit,
                ncp = ncp, rev = rev,
                set_seed = set_seed,
                seed_number = seed_number,
                ...
            )
        } else {
            stop("only PCA methods from the BiocSingular and factominer 
                package have been implemented")
        }
    } else {
        # PCA on genes
        if (method %in% c("irlba", "exact", "random")) {
            pca_object <- .run_pca_biocsingular(
                x = expr_values,
                center = center,
                scale = scale_unit,
                ncp = ncp,
                rev = rev,
                set_seed = set_seed,
                seed_number = seed_number,
                BSPARAM = method,
                BPPARAM = method_params,
                ...
            )
        } else if (method == "factominer") {
            pca_object <- .run_pca_factominer(
                x = expr_values,
                scale = scale_unit, ncp = ncp, rev = rev,
                set_seed = set_seed, seed_number = seed_number, ...
            )
        } else {
            stop("only PCA methods from the irlba and factominer package have 
                been implemented")
        }
    }




    if (isTRUE(return_gobject)) {
        pca_names <- list_dim_reductions_names(
            gobject = gobject,
            data_type = reduction,
            spat_unit = spat_unit,
            feat_type = feat_type,
            dim_type = "pca"
        )

        if (name %in% pca_names) {
            cat(name, " has already been used, will be overwritten")
        }

        if (reduction == "cells") {
            my_row_names <- colnames(expr_values)
        } else {
            my_row_names <- rownames(expr_values)
        }

        dimObject <- create_dim_obj(
            name = name,
            feat_type = feat_type,
            spat_unit = spat_unit,
            provenance = provenance,
            reduction = reduction,
            reduction_method = "pca",
            coordinates = pca_object$coords,
            misc = list(
                eigenvalues = pca_object$eigenvalues,
                loadings = pca_object$loadings
            ),
            my_rownames = my_row_names
        )

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_dimReduction(gobject = gobject, dimObject = dimObject)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


        ## update parameters used ##
        gobject <- update_giotto_params(gobject, description = "_pca")
        return(gobject)
    } else {
        return(pca_object)
    }
}







## * PCA projections  ####
# ---------------------- #


#' @title Run PCA - BiocSingular irlba projection
#' @name .run_pca_biocsingular_irlba_projection
#' @description Performs PCA based on the biocSingular package on a
#' subset of the matrix. It uses the obtained loadings to predicted coordinates
#' for the remaining matrix.
#' @param x matrix or object that can be converted to matrix
#' @param ncp number of principal components to calculate
#' @param center center the matrix before pca
#' @param scale scale features
#' @param rev reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param BPPARAM BiocParallelParam object
#' @param random_subset random subset to perform PCA on
#' @param verbose verbosity level
#' @keywords internal
#' @returns list of eigenvalues, loadings and pca coordinates
.run_pca_biocsingular_irlba_projection <- function(x,
    ncp = 100,
    center = TRUE,
    scale = TRUE,
    rev = FALSE,
    set_seed = TRUE,
    seed_number = 1234,
    BPPARAM = BiocParallel::SerialParam(),
    random_subset = 500,
    verbose = TRUE,
    ...) {
    x <- scale(x, center = center, scale = scale)

    min_ncp <- min(dim(x))

    if (ncp >= min_ncp) {
        warning("ncp >= minimum dimension of x, will be set to minimum 
                dimension of x - 1")
        ncp <- min_ncp - 1
    }

    # seed
    if (isTRUE(set_seed)) {
        set.seed(seed = seed_number)
    }
    on.exit(random_seed(), add = TRUE)




    if (isTRUE(rev)) {
        x <- t_flex(x)

        # store cell ID order information
        cell_ID_order <- rownames(x)

        # create random selection
        random_selection <- sort(sample(1:nrow(x), random_subset))
        subsample_matrix <- x[random_selection, ]


        if (verbose) message("pca random subset: start")

        # pca on random selection
        pca_res <- BiocSingular::runPCA(
            x = subsample_matrix,
            rank = ncp,
            center = FALSE, scale = FALSE,
            BSPARAM = BiocSingular::IrlbaParam(),
            BPPARAM = BPPARAM,
            ...
        )

        if (verbose) message("pca random subset: done")
        if (verbose) message("pca prediction: start")

        # create leftover matrix
        leftover_matrix <- x[-random_selection, ]

        # predict on leftover matrix
        class(pca_res) <- "prcomp"
        pca_res$center <- FALSE
        pca_res$scale <- FALSE
        project_results <- predict(pca_res, leftover_matrix)

        if (verbose) message("pca prediction: done")

        # combine subsample + predicted coordinates
        coords <- rbind(pca_res$x, project_results)
        coords <- coords[match(cell_ID_order, rownames(coords)), ]

        # eigenvalues
        eigenvalues <- pca_res$sdev^2

        # PC loading
        loadings <- coords
        rownames(loadings) <- rownames(x)
        colnames(loadings) <- paste0("Dim.", 1:ncol(loadings))

        # coordinates
        coords <- pca_res$rotation
        rownames(coords) <- colnames(x)
        colnames(coords) <- paste0("Dim.", 1:ncol(coords))

        result <- list(
            eigenvalues = eigenvalues, loadings = loadings, coords = coords)
    } else {
        # store cell ID order information
        cell_ID_order <- rownames(x)

        # create random selection
        random_selection <- sort(sample(1:nrow(x), random_subset))
        subsample_matrix <- x[random_selection, ]

        if (verbose) message("pca random subset: start")

        pca_res <- BiocSingular::runPCA(
            x = subsample_matrix,
            rank = ncp,
            center = FALSE, scale = FALSE,
            BSPARAM = BiocSingular::IrlbaParam(),
            BPPARAM = BPPARAM,
            ...
        )

        if (verbose) message("pca random subset: done")
        if (verbose) message("pca prediction: start")

        # create leftover matrix
        leftover_matrix <- x[-random_selection, ]

        # predict on leftover matrix
        class(pca_res) <- "prcomp"
        pca_res$center <- FALSE
        pca_res$scale <- FALSE
        project_results <- predict(pca_res, leftover_matrix)

        if (verbose) message("pca prediction: done")

        # combine subsample + predicted coordinates
        coords <- rbind(pca_res$x, project_results)
        coords <- coords[match(cell_ID_order, rownames(coords)), ]

        # eigenvalues
        eigenvalues <- pca_res$sdev^2

        # PC loading
        loadings <- pca_res$rotation
        rownames(loadings) <- colnames(x)
        colnames(loadings) <- paste0("Dim.", 1:ncol(loadings))

        # coordinates
        colnames(coords) <- paste0("Dim.", 1:ncol(coords))

        result <- list(
            eigenvalues = eigenvalues, loadings = loadings, coords = coords)
    }

    return(result)
}









#' @title runPCAprojection
#' @name runPCAprojection
#' @description runs a Principal Component Analysis on a random 
#' subset + projection
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param random_subset numeric. random subset to perform PCA on
#' @param name arbitrary name for PCA run
#' @param feats_to_use subset of features to use for PCA
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param center center data first (default = TRUE)
#' @param scale_unit scale features before PCA (default = TRUE)
#' @param ncp number of principal components to calculate
#' @param method which implementation to use
#' @param method_params BiocParallelParam object
#' @param rev do a reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param verbose verbosity of the function
#' @param ... additional parameters for PCA (see details)
#' @returns giotto object with updated PCA dimension recuction
#' @details See \code{\link[BiocSingular]{runPCA}} and 
#' \code{\link[FactoMineR]{PCA}} for more information about other parameters.
#' This PCA implementation is similar to  \code{\link{runPCA}}, except that it
#' performs PCA on a subset of the cells or features, and predict on the others.
#' This can significantly increase speed without sacrificing accuracy too much.
#' \itemize{
#'   \item feats_to_use = NULL: will use all features from the selected matrix
#'   \item feats_to_use = <hvg name>: can be used to select a column name of
#'   highly variable features, created by (see \code{\link{calculateHVF}})
#'   \item feats_to_use = c('geneA', 'geneB', ...): will use all manually 
#'   provided features
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' runPCAprojection(g)
#' @export
runPCAprojection <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    random_subset = 500,
    name = "pca.projection",
    feats_to_use = "hvf",
    return_gobject = TRUE,
    center = TRUE,
    scale_unit = TRUE,
    ncp = 100,
    method = c("irlba"),
    method_params = BiocParallel::SerialParam(),
    rev = FALSE,
    set_seed = TRUE,
    seed_number = 1234,
    verbose = TRUE,
    ...) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # specify name to use for pca
    if (is.null(name)) {
        if (feat_type == "rna") {
            name <- "pca"
        } else {
            name <- paste0(feat_type, ".", "pca")
        }
    }

    # expression values to be used
    values <- match.arg(
        expression_values, 
        unique(c("normalized", "scaled", "custom", expression_values)))
    expr_values <- getExpression(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        values = values,
        output = "exprObj"
    )

    provenance <- prov(expr_values)

    if (!is.null(slot(gobject, "h5_file"))) {
        expr_path <- slot(expr_values, "exprMat")

        expr_values <- HDF5Array::h5mread(
            filepath = slot(gobject, "h5_file"),
            name = paste0(
                "expression/",
                feat_type, "/",
                values
            )
        )

        expr_dimnames <- HDF5Array::h5readDimnames(
            filepath = slot(gobject, "h5_file"),
            name = paste0(
                "expression/",
                feat_type, "/",
                values
            )
        )

        rownames(expr_values) <- expr_dimnames[[1]]
        colnames(expr_values) <- expr_dimnames[[2]]
    } else {
        expr_values <- expr_values[] # extract matrix
    }



    ## subset matrix
    if (!is.null(feats_to_use)) {
        expr_values <- .create_feats_to_use_matrix(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sel_matrix = expr_values,
            feats_to_use = feats_to_use,
            verbose = verbose
        )
    }


    # do PCA dimension reduction
    reduction <- match.arg(reduction, c("cells", "feats"))

    # PCA implementation
    method <- match.arg(method, c("irlba"))

    if (reduction == "cells") {
        # PCA on cells
        pca_object <- .run_pca_biocsingular_irlba_projection(
            x = t_flex(expr_values),
            ncp = ncp,
            center = center,
            scale = scale_unit,
            rev = rev,
            set_seed = set_seed,
            seed_number = seed_number,
            BPPARAM = method_params,
            random_subset = random_subset,
            verbose = verbose,
            ...
        )
    } else {
        # PCA on features
        pca_object <- .run_pca_biocsingular_irlba_projection(
            x = expr_values,
            ncp = ncp,
            center = center,
            scale = scale_unit,
            rev = rev,
            set_seed = set_seed,
            seed_number = seed_number,
            BPPARAM = method_params,
            random_subset = random_subset,
            verbose = verbose,
            ...
        )
    }

    if (isTRUE(return_gobject)) {
        pca_names <- list_dim_reductions_names(
            gobject = gobject,
            data_type = reduction,
            spat_unit = spat_unit,
            feat_type = feat_type,
            dim_type = "pca"
        )

        if (name %in% pca_names) {
            cat(name, " has already been used, will be overwritten")
        }

        if (reduction == "cells") {
            my_row_names <- colnames(expr_values)
        } else {
            my_row_names <- rownames(expr_values)
        }

        dimObject <- create_dim_obj(
            name = name,
            feat_type = feat_type,
            spat_unit = spat_unit,
            provenance = provenance,
            reduction = reduction,
            reduction_method = "pca",
            coordinates = pca_object$coords,
            misc = list(
                eigenvalues = pca_object$eigenvalues,
                loadings = pca_object$loadings
            ),
            my_rownames = my_row_names
        )

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_dimReduction(gobject = gobject, dimObject = dimObject)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


        ## update parameters used ##
        gobject <- update_giotto_params(gobject, description = "_pca")
        return(gobject)
    } else {
        return(pca_object)
    }
}





#' @title runPCAprojectionBatch
#' @name runPCAprojectionBatch
#' @description runs a Principal Component Analysis on multiple random 
#' batches + projection
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param random_subset random subset to perform PCA on
#' @param batch_number number of random batches to run
#' @param name arbitrary name for PCA run
#' @param feats_to_use subset of features to use for PCA
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param center center data first (default = TRUE)
#' @param scale_unit scale features before PCA (default = TRUE)
#' @param ncp number of principal components to calculate
#' @param method which implementation to use
#' @param method_params BiocParallelParam object
#' @param rev do a reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param verbose verbosity of the function
#' @param ... additional parameters for PCA (see details)
#' @returns giotto object with updated PCA dimension reduction
#' @details See \code{\link[BiocSingular]{runPCA}} and 
#' \code{\link[FactoMineR]{PCA}} for more information about other parameters.
#' This PCA implementation is similar to  \code{\link{runPCA}} and  
#' \code{\link{runPCAprojection}},
#' except that it performs PCA on multiple subsets (batches) of the cells or 
#' features,
#' and predict on the others. This can significantly increase speed without 
#' sacrificing accuracy too much.
#' \itemize{
#'   \item feats_to_use = NULL: will use all features from the selected matrix
#'   \item feats_to_use = <hvg name>: can be used to select a column name of
#'   highly variable features, created by (see \code{\link{calculateHVF}})
#'   \item feats_to_use = c('geneA', 'geneB', ...): will use all manually 
#'   provided features
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' runPCAprojectionBatch(g)
#' @export
runPCAprojectionBatch <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    random_subset = 500,
    batch_number = 5,
    name = "pca.projection.batch",
    feats_to_use = "hvf",
    return_gobject = TRUE,
    center = TRUE,
    scale_unit = TRUE,
    ncp = 100,
    method = c("irlba"),
    method_params = BiocParallel::SerialParam(),
    rev = FALSE,
    set_seed = TRUE,
    seed_number = 1234,
    verbose = TRUE,
    ...) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # specify name to use for pca
    if (is.null(name)) {
        if (feat_type == "rna") {
            name <- "pca"
        } else {
            name <- paste0(feat_type, ".", "pca")
        }
    }

    # expression values to be used
    values <- match.arg(
        expression_values, 
        unique(c("normalized", "scaled", "custom", expression_values)))
    expr_values <- getExpression(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        values = values,
        output = "exprObj"
    )

    provenance <- prov(expr_values)

    if (!is.null(slot(gobject, "h5_file"))) {
        expr_path <- slot(expr_values, "exprMat")

        expr_values <- HDF5Array::h5mread(
            filepath = slot(gobject, "h5_file"),
            name = paste0(
                "expression/",
                feat_type, "/",
                values
            )
        )

        expr_dimnames <- HDF5Array::h5readDimnames(
            filepath = slot(gobject, "h5_file"),
            name = paste0(
                "expression/",
                feat_type, "/",
                values
            )
        )

        rownames(expr_values) <- expr_dimnames[[1]]
        colnames(expr_values) <- expr_dimnames[[2]]
    } else {
        expr_values <- expr_values[] # extract matrix
    }



    ## subset matrix
    if (!is.null(feats_to_use)) {
        expr_values <- .create_feats_to_use_matrix(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sel_matrix = expr_values,
            feats_to_use = feats_to_use,
            verbose = verbose
        )
    }


    # do PCA dimension reduction
    reduction <- match.arg(reduction, c("cells", "feats"))

    # PCA implementation
    method <- match.arg(method, c("irlba"))

    if (reduction == "cells") {
        pca_batch_results <- list()

        for (batch in 1:batch_number) {
            if (verbose) wrap_msg("start batch ", batch)

            if (isTRUE(set_seed)) {
                seed_batch <- seed_number + batch
            } else {
                seed_batch <- NULL
            }


            # PCA on cells
            pca_object <- .run_pca_biocsingular_irlba_projection(
                x = t_flex(expr_values),
                ncp = ncp,
                center = center,
                scale = scale_unit,
                rev = rev,
                set_seed = set_seed,
                seed_number = seed_batch,
                BPPARAM = method_params,
                random_subset = random_subset,
                verbose = verbose,
                ...
            )

            # adjust the sign of the coordinates and loadings vector relative 
            # to first batch this is necessary for the next averaging step
            if (batch == 1) {
                pca_batch_results[[batch]] <- pca_object
            } else {
                for (dimension in 1:ncol(pca_object[["coords"]])) {
                    sum_evaluation <- sum(sign(pca_batch_results[[1]][[
                        "coords"]][1:20, dimension]) *
                        sign(pca_object[["coords"]][1:20, dimension]))
                    if (sum_evaluation < 0) {
                        pca_object$coords[, dimension] <- -1 * pca_object$coords[, dimension]
                        pca_object$loadings[, dimension] <- -1 * pca_object$loadings[, dimension]
                    }
                }
                pca_batch_results[[batch]] <- pca_object
            }
        }

        if (verbose) message("start averaging pca results of batches")

        # calculate average eigenvalues, coordinates and loadings
        # TODO: test out DT approach, might be faster and more efficient for
        # large datasets

        # eigenvalues
        eigenvalues_list <- lapply(
            pca_batch_results, FUN = function(x) x$eigenvalues)
        eigenvalues_matrix <- do.call("cbind", eigenvalues_list)
        eigenvalues_mean <- rowMeans_flex(eigenvalues_matrix)

        # coordinates
        coords_list <- lapply(pca_batch_results, FUN = function(x) x$coords)
        coords_vector <- do.call("c", coords_list)
        coords_array <- array(
            data = coords_vector, 
            dim = c(ncol(expr_values), ncp, length(pca_batch_results)))
        coords_all <- apply(coords_array, MARGIN = c(1:2), function(arr) {
            mean(arr, na.rm = TRUE)
        })
        rownames(coords_all) <- rownames(pca_batch_results[[1]][["coords"]])
        colnames(coords_all) <- colnames(pca_batch_results[[1]][["coords"]])

        # loadings
        loadings_list <- lapply(pca_batch_results, FUN = function(x) x$loadings)
        loadings_vector <- do.call("c", loadings_list)
        loadings_array <- array(
            data = loadings_vector, 
            dim = c(nrow(expr_values), ncp, length(pca_batch_results)))
        loadings_all <- apply(loadings_array, MARGIN = c(1:2), function(arr) {
            mean(arr, na.rm = TRUE)
        })
        rownames(loadings_all) <- rownames(pca_batch_results[[1]][["loadings"]])
        colnames(loadings_all) <- colnames(pca_batch_results[[1]][["loadings"]])


        pca_object <- list(
            eigenvalues = eigenvalues_mean, 
            loadings = loadings_all, coords = coords_all)
    } else {
        pca_batch_results <- list()

        for (batch in 1:batch_number) {
            if (verbose) wrap_msg("start batch ", batch)

            if (isTRUE(set_seed)) {
                seed_batch <- seed_number + batch
            } else {
                seed_batch <- NULL
            }


            # PCA on features
            pca_object <- .run_pca_biocsingular_irlba_projection(
                x = expr_values,
                ncp = ncp,
                center = center,
                scale = scale_unit,
                rev = rev,
                set_seed = set_seed,
                seed_number = seed_number,
                BPPARAM = method_params,
                random_subset = random_subset,
                verbose = verbose,
                ...
            )


            # adjust the sign of the coordinates and loadings vector relative 
            # to first batch this is necessary for the next averaging step
            if (batch == 1) {
                pca_batch_results[[batch]] <- pca_object
            } else {
                for (dimension in 1:ncol(pca_object[["coords"]])) {
                    sum_evaluation <- sum(sign(pca_batch_results[[1]][[
                        "coords"]][1:20, dimension]) *
                        sign(pca_object[["coords"]][1:20, dimension]))
                    if (sum_evaluation < 0) {
                        pca_object$coords[, dimension] <- -1 * pca_object$coords[, dimension]
                        pca_object$loadings[, dimension] <- -1 * pca_object$loadings[, dimension]
                    }
                }
                pca_batch_results[[batch]] <- pca_object
            }
        }

        if (verbose) wrap_msg("start averaging pca results of batches")

        # calculate average eigenvalues, coordinates and loadings
        # TODO: test out DT approach, might be faster and more efficient for
        # large datasets

        # eigenvalues
        eigenvalues_list <- lapply(
            pca_batch_results, FUN = function(x) x$eigenvalues)
        eigenvalues_matrix <- do.call("cbind", eigenvalues_list)
        eigenvalues_mean <- rowMeans_flex(eigenvalues_matrix)

        # coordinates
        coords_list <- lapply(pca_batch_results, FUN = function(x) x$coords)
        coords_vector <- do.call("c", coords_list)
        coords_array <- array(
            data = coords_vector, 
            dim = c(ncol(expr_values), ncp, length(pca_batch_results)))
        coords_all <- apply(coords_array, MARGIN = c(1:2), function(arr) {
            mean(arr, na.rm = TRUE)
        })
        rownames(coords_all) <- rownames(pca_batch_results[[1]][["coords"]])
        colnames(coords_all) <- colnames(pca_batch_results[[1]][["coords"]])

        # loadings
        loadings_list <- lapply(pca_batch_results, FUN = function(x) x$loadings)
        loadings_vector <- do.call("c", loadings_list)
        loadings_array <- array(
            data = loadings_vector, 
            dim = c(nrow(expr_values), ncp, length(pca_batch_results)))
        loadings_all <- apply(loadings_array, MARGIN = c(1:2), function(arr) {
            mean(arr, na.rm = TRUE)
        })
        rownames(loadings_all) <- rownames(pca_batch_results[[1]][["loadings"]])
        colnames(loadings_all) <- colnames(pca_batch_results[[1]][["loadings"]])


        pca_object <- list(
            eigenvalues = eigenvalues_mean, 
            loadings = loadings_all, coords = coords_all)
    }


    if (return_gobject == TRUE) {
        pca_names <- list_dim_reductions_names(
            gobject = gobject,
            data_type = reduction,
            spat_unit = spat_unit,
            feat_type = feat_type,
            dim_type = "pca"
        )

        if (name %in% pca_names) {
            cat(name, " has already been used, will be overwritten")
        }

        if (reduction == "cells") {
            my_row_names <- colnames(expr_values)
        } else {
            my_row_names <- rownames(expr_values)
        }

        dimObject <- create_dim_obj(
            name = name,
            feat_type = feat_type,
            spat_unit = spat_unit,
            provenance = provenance,
            reduction = reduction,
            reduction_method = "pca",
            coordinates = pca_object$coords,
            misc = list(
                eigenvalues = pca_object$eigenvalues,
                loadings = pca_object$loadings
            ),
            my_rownames = my_row_names
        )

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_dimReduction(gobject = gobject, dimObject = dimObject)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


        ## update parameters used ##
        gobject <- update_giotto_params(gobject, description = "_pca")
        return(gobject)
    } else {
        return(pca_object)
    }
}





## * PC estimates ####
# ------------------ #


#' @title screePlot
#' @name screePlot
#' @description identify significant principal components (PCs) using an 
#' screeplot (a.k.a. elbowplot)
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams create_screeplot
#' @param name name of PCA object if available
#' @param expression_values expression values to use
#' @param reduction cells or features
#' @param method which implementation to use
#' @param rev do a reverse PCA
#' @param feats_to_use subset of features to use for PCA
#' @param center center data before PCA
#' @param scale_unit scale features before PCA
#' @param verbose be verbose
#' @param ... additional arguments to pca function, see \code{\link{runPCA}}
#' @returns ggplot object for scree method
#' @details
#'  Screeplot works by plotting the explained variance of each
#'  individual PC in a barplot allowing you to identify which PC provides a 
#'  significant contribution (a.k.a 'elbow method'). \cr
#'  Screeplot will use an available pca object, based on the parameter 'name', 
#'  or it will create it if it's not available (see \code{\link{runPCA}})
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' screePlot(g)
#' @export
screePlot <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    name = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    method = c("irlba", "exact", "random", "factominer"),
    rev = FALSE,
    feats_to_use = NULL,
    center = FALSE,
    scale_unit = FALSE,
    ncp = 100,
    ylim = c(0, 20),
    verbose = TRUE,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "screePlot",
    ...) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # specify name to use for screeplot
    if (is.null(name)) {
        if (feat_type == "rna") {
            name <- "pca"
        } else {
            name <- paste0(feat_type, ".", "pca")
        }
    }

    # select direction of reduction
    reduction <- match.arg(reduction, c("cells", "feats"))
    pca_obj <- getDimReduction(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        reduction = reduction,
        reduction_method = "pca",
        name = name,
        output = "dimObj"
    )

    # if pca already exists plot
    if (!is.null(pca_obj)) {
        if (isTRUE(verbose)) 
            wrap_msg("PCA with name: ", name, 
                    " already exists and will be used for the screeplot")

        screeplot <- create_screeplot(
            eigs = slot(pca_obj, "misc")$eigenvalues, ncp = ncp, ylim = ylim)
    } else {
        # if pca doesn't exists, then create pca and then plot
        if (isTRUE(verbose)) 
            wrap_msg("PCA with name: ", name, 
                    " does NOT exist, PCA will be done first")

        # expression values to be used
        values <- match.arg(
            expression_values, 
            unique(c("normalized", "scaled", "custom", expression_values)))
        expr_values <- getExpression(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            values = values,
            output = "exprObj"
        )

        provenance <- prov(expr_values)
        expr_values <- expr_values[] # extract matrix

        # PCA implementation
        biocsingular_methods <- c("irlba", "exact", "random", "factominer")
        method <- match.arg(method, choices = biocsingular_methods)

        ## subset matrix
        if (!is.null(feats_to_use)) {
            expr_values <- .create_feats_to_use_matrix(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                sel_matrix = expr_values,
                feats_to_use = feats_to_use,
                verbose = verbose
            )
        }

        # reduction of cells
        if (reduction == "cells") {
            # PCA on cells
            if (method %in% biocsingular_methods) {
                pca_object <- .run_pca_biocsingular(
                    x = t_flex(expr_values),
                    center = center,
                    scale = scale_unit,
                    ncp = ncp,
                    rev = rev,
                    BSPARAM = method,
                    BPPARAM = BiocParallel::SerialParam(),
                    ...
                )
            } else if (method == "factominer") {
                pca_object <- .run_pca_factominer(
                    x = t_flex(expr_values), 
                    scale = scale_unit, ncp = ncp, rev = rev, ...)
            } else {
                stop("only PCA methods from the irlba and factominer package 
                    have been implemented")
            }

            dimObject <- create_dim_obj(
                name = name,
                feat_type = feat_type,
                spat_unit = spat_unit,
                provenance = provenance,
                reduction = reduction,
                reduction_method = "pca",
                coordinates = pca_object$coords,
                misc = list(
                    eigenvalues = pca_object$eigenvalues,
                    loadings = pca_object$loadings
                ),
                my_rownames = colnames(expr_values)
            )

            screeplot <- create_screeplot(
                eigs = slot(dimObject, "misc")$eigenvalues, 
                ncp = ncp, ylim = ylim)
        }
    }

    return(plot_output_handler(
        gobject = gobject,
        plot_object = screeplot,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}





#' @title create_screeplot
#' @name create_screeplot
#' @description create screeplot with ggplot
#' @param eigs numeric. Vector of pca eigenvalues
#' @param ncp numeric. max number of principal components to plot
#' @param ylim numeric. y-axis limits on scree plot
#' @returns ggplot
#' @examples
#' dr <- GiottoData::loadSubObjectMini("dimObj")
#' scree <- create_screeplot(methods::slot(dr, "misc")$eigenvalues)
#' scree
#' @export
create_screeplot <- function(eigs, ncp = 20, ylim = c(0, 20)) {
    checkmate::assert_numeric(eigs)
    checkmate::assert_numeric(ncp, len = 1L)
    checkmate::assert_numeric(ylim, len = 2L)

    # DT vars
    PC <- NULL

    eigs <- sort(eigs, decreasing = TRUE)

    # variance explained
    var_expl <- eigs / sum(eigs) * 100
    var_expl_cum <- cumsum(eigs) / sum(eigs) * 100

    # create data.table
    screeDT <- data.table::data.table(
        "PC" = paste0("PC.", seq_along(var_expl)),
        "var_expl" = var_expl,
        "var_expl_cum" = var_expl_cum
    )
    screeDT[, "PC" := factor(PC, levels = PC)]

    max_ncp <- length(eigs)
    ncp <- ifelse(ncp > max_ncp, max_ncp, ncp)

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()
    pl <- pl + ggplot2::geom_bar(
        data = screeDT[1:ncp], 
        ggplot2::aes(x = PC, y = var_expl), stat = "identity")
    pl <- pl + ggplot2::coord_cartesian(ylim = ylim)
    pl <- pl + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
    pl <- pl + ggplot2::labs(x = "", y = "% of variance explained per PC")

    cpl <- ggplot2::ggplot()
    cpl <- cpl + ggplot2::theme_bw()
    cpl <- cpl + ggplot2::geom_bar(
        data = screeDT[1:ncp], 
        ggplot2::aes(x = PC, y = var_expl_cum), stat = "identity")
    cpl <- cpl + ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1))
    cpl <- cpl + ggplot2::labs(x = "", y = "cumulative % of variance explained")

    savelist <- list(pl, cpl)

    ## combine plots with cowplot
    combo_plot <- cowplot::plot_grid(
        plotlist = savelist,
        ncol = 1,
        rel_heights = c(1),
        rel_widths = c(1),
        align = "v"
    )

    return(combo_plot)
}








#' @title jackstrawPlot
#' @name jackstrawPlot
#' @description identify significant prinicipal components (PCs)
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param feats_to_use subset of features to use for PCA
#' @param center center data before PCA
#' @param scale_unit scale features before PCA
#' @param ncp number of principal components to calculate
#' @param ylim y-axis limits on jackstraw plot
#' @param iter number of interations for jackstraw
#' @param threshold p-value threshold to call a PC significant
#' @param verbose show progress of jackstraw method
#' @returns ggplot object for jackstraw method
#' @details
#'  The Jackstraw method uses the \code{\link[jackstraw]{permutationPA}} 
#'  function. By systematically permuting genes it identifies robust, and thus 
#'  significant, PCs.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' jackstrawPlot(gobject = g)
#' @export
jackstrawPlot <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    feats_to_use = NULL,
    center = FALSE,
    scale_unit = FALSE,
    ncp = 20,
    ylim = c(0, 1),
    iter = 10,
    threshold = 0.01,
    verbose = TRUE,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "jackstrawPlot") {
    package_check(pkg_name = "jackstraw", repository = "CRAN")

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # print message with information #
    if (verbose) 
        message("using 'jackstraw' to identify significant PCs If used in 
        published research, please cite: 
        Neo Christopher Chung and John D. Storey (2014).
        'Statistical significance of variables driving systematic variation in 
        high-dimensional data. Bioinformatics")

    # select direction of reduction
    reduction <- match.arg(reduction, c("cells", "feats"))

    # expression values to be used
    values <- match.arg(
        expression_values, 
        unique(c("normalized", "scaled", "custom", expression_values)))
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )


    ## subset matrix
    if (!is.null(feats_to_use)) {
        expr_values <- .create_feats_to_use_matrix(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sel_matrix = expr_values,
            feats_to_use = feats_to_use,
            verbose = verbose
        )
    }

    # reduction of cells
    if (reduction == "cells") {
        if (scale_unit == TRUE | center == TRUE) {
            expr_values <- t_flex(scale(
                t_flex(expr_values), center = center, scale = scale_unit))
        }

        jtest <- jackstraw::permutationPA(
            dat = as.matrix(expr_values), 
            B = iter, threshold = threshold, verbose = verbose)

        ## results ##
        nr_sign_components <- jtest$r
        if (verbose) 
            cat("number of estimated significant components: ", 
                nr_sign_components)
        final_results <- jtest$p
        jackplot <- create_jackstrawplot(
            jackstraw_data = final_results, 
            ncp = ncp, ylim = ylim, threshold = threshold)
    }

    return(plot_output_handler(
        gobject = gobject,
        plot_object = jackplot,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}



#' @title create_jackstrawplot
#' @name create_jackstrawplot
#' @description create jackstrawplot with ggplot
#' @param jackstraw_data result from jackstraw function (`testresult$p`)
#' @param ncp number of principal components to calculate
#' @param ylim y-axis limits on jackstraw plot
#' @param threshold p.value threshold to call a PC significant
#' @keywords internal
#' @returns ggplot
#' @export
create_jackstrawplot <- function(jackstraw_data,
    ncp = 20,
    ylim = c(0, 1),
    threshold = 0.01) {
    checkmate::assert_numeric(ncp, len = 1L)
    checkmate::assert_numeric(ylim, len = 2L)
    checkmate::assert_numeric(threshold, len = 1L)

    # data.table variables
    PC <- p.val <- NULL

    testDT <- data.table::data.table(
        PC = paste0("PC.", seq_along(jackstraw_data)),
        p.val = jackstraw_data
    )
    testDT[, "PC" := factor(PC, levels = PC)]
    testDT[, "sign" := ifelse(p.val <= threshold, "sign", "n.s.")]

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()
    pl <- pl + ggplot2::geom_point(
        data = testDT[1:ncp], 
        ggplot2::aes(x = PC, y = p.val, fill = sign), shape = 21)
    pl <- pl + ggplot2::scale_fill_manual(
        values = c("n.s." = "lightgrey", "sign" = "darkorange"))
    pl <- pl + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
    pl <- pl + ggplot2::coord_cartesian(ylim = ylim)
    pl <- pl + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
    pl <- pl + ggplot2::labs(x = "", y = "p-value per PC")

    return(pl)
}








#' @title signPCA
#' @name signPCA
#' @description identify significant prinicipal components (PCs)
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param name name of PCA object if available
#' @param method method to use to identify significant PCs
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param pca_method which implementation to use
#' @param rev do a reverse PCA
#' @param feats_to_use subset of features to use for PCA
#' @param center center data before PCA
#' @param scale_unit scale features before PCA
#' @param ncp number of principal components to calculate
#' @param scree_ylim y-axis limits on scree plot
#' @param jack_iter number of interations for jackstraw
#' @param jack_threshold p-value threshold to call a PC significant
#' @param jack_ylim y-axis limits on jackstraw plot
#' @param verbose be verbose
#' @returns ggplot object for scree method and maxtrix of p-values for jackstraw
#' @details Two different methods can be used to assess the number of relevant 
#' or significant prinicipal components (PC's). \cr
#'  1. Screeplot works by plotting the explained variance of each
#'  individual PC in a barplot allowing you to identify which PC provides a 
#'  significant
#'  contribution  (a.k.a. 'elbow method'). \cr
#'  2. The Jackstraw method uses the \code{\link[jackstraw]{permutationPA}} 
#'  function. By systematically permuting genes it identifies robust, and thus 
#'  significant, PCs.
#'  \cr
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' signPCA(g)
#' @export
signPCA <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    name = NULL,
    method = c("screeplot", "jackstraw"),
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    pca_method = c("irlba", "factominer"),
    rev = FALSE,
    feats_to_use = NULL,
    center = TRUE,
    scale_unit = TRUE,
    ncp = 50,
    scree_ylim = c(0, 10),
    jack_iter = 10,
    jack_threshold = 0.01,
    jack_ylim = c(0, 1),
    verbose = TRUE,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "signPCA") {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # specify name to use
    if (!is.null(name)) {
        if (feat_type == "rna") {
            name <- "pca"
        } else {
            name <- paste0(feat_type, ".", "pca")
        }
    }

    # select method
    method <- match.arg(method, choices = c("screeplot", "jackstraw"))

    # select PCA method
    pca_method <- match.arg(pca_method, choices = c("irlba", "factominer"))

    # select direction of reduction
    reduction <- match.arg(reduction, c("cells", "feats"))

    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )

    ## subset matrix
    if (!is.null(feats_to_use)) {
        expr_values <- .create_feats_to_use_matrix(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sel_matrix = expr_values,
            feats_to_use = feats_to_use,
            verbose = verbose
        )
    }

    # reduction of cells
    if (reduction == "cells") {
        if (method == "screeplot") {
            screeplot <- screePlot(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                name = name,
                expression_values = values,
                reduction = reduction,
                feats_to_use = feats_to_use,
                center = center,
                scale_unit = scale_unit,
                ncp = ncp,
                rev = rev,
                method = pca_method,
                ylim = scree_ylim,
                verbose = verbose,
                show_plot = FALSE,
                return_plot = TRUE,
                save_plot = FALSE,
                save_param = list(),
                default_save_name = "screePlot"
            )

            return(GiottoVisuals::plot_output_handler(
                gobject = gobject,
                plot_object = screeplot,
                save_plot = save_plot,
                return_plot = return_plot,
                show_plot = show_plot,
                default_save_name = default_save_name,
                save_param = save_param,
                else_return = NULL
            ))
        } else if (method == "jackstraw") {
            jackplot <- jackstrawPlot(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                expression_values = values,
                reduction = reduction,
                feats_to_use = feats_to_use,
                center = center,
                scale_unit = scale_unit,
                ncp = ncp,
                ylim = jack_ylim,
                iter = jack_iter,
                threshold = jack_threshold,
                verbose = verbose,
                show_plot = FALSE,
                return_plot = TRUE,
                save_plot = FALSE,
                save_param = list(),
                default_save_name = "jackstrawPlot"
            )

            return(GiottoVisuals::plot_output_handler(
                gobject = gobject,
                plot_object = jackplot,
                save_plot = save_plot,
                return_plot = return_plot,
                show_plot = show_plot,
                default_save_name = default_save_name,
                save_param = save_param,
                else_return = jackplot 
                # TODO potentially return all results instead
            ))
        }
    } else {
        message("gene reduction not yet implemented")
    }
}








## * Dim reduction algos ####
# ------------------------- #

#' @title Run UMAP dimension reduction
#' @name runUMAP
#' @description run UMAP
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name arbitrary name for UMAP run
#' @param feats_to_use if dim_reduction_to_use = NULL, which features to use
#' @param return_gobject logical: return giotto object (default = TRUE)
#' @param n_neighbors UMAP param: number of neighbors
#' @param n_components UMAP param: number of components
#' @param n_epochs UMAP param: number of epochs
#' @param min_dist UMAP param: minimum distance
#' @param n_threads UMAP param: threads/cores to use
#' @param spread UMAP param: spread
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param verbose verbosity of function
#' @param toplevel_params parameters to extract
#' @inheritDotParams uwot::umap -X -n_neighbors -n_components -n_epochs 
#' -min_dist -n_threads -spread -seed -scale -pca -pca_center -pca_method
#' @returns giotto object with updated UMAP dimension reduction
#' @details See \code{\link[uwot]{umap}} for more information about these and
#' other parameters.
#' \itemize{
#'   \item Input for UMAP dimension reduction can be another dimension reduction
#'   (default = 'pca')
#'   \item To use gene expression as input set dim_reduction_to_use = NULL
#'   \item If dim_reduction_to_use = NULL, feats_to_use can be used to select a
#'   column name of highly variable features (see \code{\link{calculateHVF}}) or
#'   simply provide a vector of features
#'   \item multiple UMAP results can be stored by changing the \emph{name} of
#'   the analysis
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' runUMAP(g)
#' @export
runUMAP <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    dim_reduction_to_use = "pca",
    dim_reduction_name = NULL,
    dimensions_to_use = 1:10,
    name = NULL,
    feats_to_use = NULL,
    return_gobject = TRUE,
    n_neighbors = 40,
    n_components = 2,
    n_epochs = 400,
    min_dist = 0.01,
    n_threads = NA,
    spread = 5,
    set_seed = TRUE,
    seed_number = 1234L,
    verbose = TRUE,
    toplevel_params = 2L,
    ...) {
    # NSE vars
    cell_ID <- NULL

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    reduction <- match.arg(reduction, choices = c("cells", "feats"))


    # specify dim_reduction_name to use for pca input for umap
    if (!is.null(dim_reduction_to_use)) {
        if (dim_reduction_to_use == "pca") {
            if (is.null(dim_reduction_name)) {
                if (feat_type == "rna") {
                    dim_reduction_name <- "pca"
                } else {
                    dim_reduction_name <- paste0(feat_type, ".", "pca")
                }
            }
        }
    }



    # specify name to use for umap
    if (is.null(name)) {
        if (feat_type == "rna") {
            name <- "umap"
        } else {
            name <- paste0(feat_type, ".", "umap")
        }
    }


    # set cores to use
    n_threads <- determine_cores(cores = n_threads)

    ## umap on cells ##
    if (reduction == "cells") {
        ## using dimension reduction ##
        if (!is.null(dim_reduction_to_use)) {
            ## TODO: check if reduction exists
            dimObj_to_use <- get_dimReduction(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                reduction = reduction,
                reduction_method = dim_reduction_to_use,
                name = dim_reduction_name,
                output = "dimObj"
            )

            provenance <- prov(dimObj_to_use)
            matrix_to_use <- dimObj_to_use[]

            if (any(!dimensions_to_use %in% seq(ncol(matrix_to_use)))) {
                warning(wrap_txt(
                    "dimensions_to_use requested is outside what is available.",
                    paste0("(1 to", ncol(matrix_to_use), ")"),
                    "Ignoring dimensions_to_use that are outside the range."
                ))
                dimensions_to_use <- dimensions_to_use[
                    dimensions_to_use %in% seq(ncol(matrix_to_use))]
            }

            matrix_to_use <- matrix_to_use[, dimensions_to_use]
        } else {
            ## using original matrix ##
            # expression values to be used
            values <- match.arg(
                expression_values, 
                unique(c("normalized", "scaled", "custom", expression_values)))

            expr_values <- getExpression(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                values = values,
                output = "exprObj"
            )

            provenance <- prov(expr_values)
            expr_values <- expr_values[] # extract matrix


            ## subset matrix
            if (!is.null(feats_to_use)) {
                expr_values <- .create_feats_to_use_matrix(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    sel_matrix = expr_values,
                    feats_to_use = feats_to_use,
                    verbose = verbose
                )
            }

            matrix_to_use <- t_flex(expr_values)
        }

        # start seed
        if (isTRUE(set_seed)) {
            set.seed(seed = seed_number)
            on.exit(
                {
                    GiottoUtils::random_seed(set.seed = TRUE)
                },
                add = TRUE
            )
        }

        ## run umap ##
        uwot_clus <- uwot::umap(
            X = matrix_to_use, # as.matrix(matrix_to_use) necessary?
            n_neighbors = n_neighbors,
            n_components = n_components,
            n_epochs = n_epochs,
            min_dist = min_dist,
            n_threads = n_threads,
            spread = spread,
            ...
        )

        uwot_clus_pos_DT <- data.table::as.data.table(uwot_clus)

        uwot_clus_pos_DT[, cell_ID := rownames(matrix_to_use)]

        # exit seed
        if (isTRUE(set_seed)) {
            set.seed(seed = Sys.time())
        }


        if (return_gobject == TRUE) {
            umap_names <- list_dim_reductions_names(
                gobject = gobject,
                data_type = reduction,
                spat_unit = spat_unit,
                feat_type = feat_type,
                dim_type = "umap"
            )

            if (name %in% umap_names) {
                message(name, " has already been used, will be overwritten")
            }


            coordinates <- uwot_clus
            rownames(coordinates) <- rownames(matrix_to_use)

            dimObject <- create_dim_obj(
                name = name,
                feat_type = feat_type,
                spat_unit = spat_unit,
                reduction = reduction,
                provenance = provenance,
                reduction_method = "umap",
                coordinates = coordinates,
                misc = NULL
            )


            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- set_dimReduction(gobject = gobject, 
                                        dimObject = dimObject)
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###



            ## update parameters used ##
            gobject <- update_giotto_params(gobject,
                description = "_umap",
                return_gobject = TRUE,
                toplevel = toplevel_params
            )
            return(gobject)
        } else {
            return(uwot_clus_pos_DT)
        }
    } else if (reduction == "feats") {
        message("Feats reduction is not yet implemented")
    }
}







#' @title Run UMAP dimension reduction
#' @name runUMAPprojection
#' @description run UMAP on subset and project on the rest
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param reduction 'cells' or 'feats'
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param random_subset random subset to perform UMAP on
#' @param name arbitrary name for UMAP run
#' @param feats_to_use if dim_reduction_to_use = NULL, which features to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param n_neighbors UMAP param: number of neighbors
#' @param n_components UMAP param: number of components
#' @param n_epochs UMAP param: number of epochs
#' @param min_dist UMAP param: minimum distance
#' @param n_threads UMAP param: threads/cores to use
#' @param spread UMAP param: spread
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param verbose verbosity of function
#' @param toplevel_params parameters to extract
#' @param ... additional UMAP parameters
#' @returns giotto object with updated UMAP dimension reduction
#' @details See \code{\link[uwot]{umap}} for more information about these and 
#' other parameters.
#' \itemize{
#'   \item Input for UMAP dimension reduction can be another dimension reduction (default = 'pca')
#'   \item To use gene expression as input set dim_reduction_to_use = NULL
#'   \item If dim_reduction_to_use = NULL, feats_to_use can be used to select a column name of
#'   highly variable genes (see \code{\link{calculateHVF}}) or simply provide a vector of genes
#'   \item multiple UMAP results can be stored by changing the \emph{name} of the analysis
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' runUMAPprojection(g)
#' @export
runUMAPprojection <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    dim_reduction_to_use = "pca",
    dim_reduction_name = NULL,
    dimensions_to_use = 1:10,
    random_subset = 500,
    name = NULL,
    feats_to_use = NULL,
    return_gobject = TRUE,
    n_neighbors = 40,
    n_components = 2,
    n_epochs = 400,
    min_dist = 0.01,
    n_threads = NA,
    spread = 5,
    set_seed = TRUE,
    seed_number = 1234,
    verbose = TRUE,
    toplevel_params = 2,
    ...) {
    # NSE vars
    cell_ID <- NULL

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    reduction <- match.arg(reduction, choices = c("cells", "feats"))


    # specify dim_reduction_name to use for pca input for umap
    if (!is.null(dim_reduction_to_use)) {
        if (dim_reduction_to_use == "pca") {
            if (is.null(dim_reduction_name)) {
                if (feat_type == "rna") {
                    dim_reduction_name <- "pca"
                } else {
                    dim_reduction_name <- paste0(feat_type, ".", "pca")
                }
            }
        }
    }



    # specify name to use for umap
    if (is.null(name)) {
        if (feat_type == "rna") {
            name <- "umap.projection"
        } else {
            name <- paste0(feat_type, ".", "umap.projection")
        }
    }


    # set cores to use
    n_threads <- determine_cores(cores = n_threads)

    ## umap on cells ##
    if (reduction == "cells") {
        ## using dimension reduction ##
        if (!is.null(dim_reduction_to_use)) {
            ## TODO: check if reduction exists
            dimObj_to_use <- getDimReduction(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                reduction = reduction,
                reduction_method = dim_reduction_to_use,
                name = dim_reduction_name,
                output = "dimObj"
            )

            provenance <- prov(dimObj_to_use)
            matrix_to_use <- dimObj_to_use[]

            matrix_to_use <- matrix_to_use[, dimensions_to_use]
        } else {
            ## using original matrix ##
            # expression values to be used
            values <- match.arg(
                expression_values,
                unique(c("normalized", "scaled", "custom", expression_values))
            )

            expr_values <- getExpression(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                values = values,
                output = "exprObj"
            )

            provenance <- prov(expr_values)
            expr_values <- expr_values[] # extract matrix

            ## subset matrix
            if (!is.null(feats_to_use)) {
                expr_values <- .create_feats_to_use_matrix(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    sel_matrix = expr_values,
                    feats_to_use = feats_to_use,
                    verbose = verbose
                )
            }

            matrix_to_use <- t_flex(expr_values)
        }

        # start seed
        if (isTRUE(set_seed)) {
            set.seed(seed = seed_number)
        }


        ## run umap ##
        cell_ID_order <- rownames(matrix_to_use)

        # create random selection
        random_selection <- sort(sample(1:nrow(matrix_to_use), random_subset))
        subsample_matrix <- matrix_to_use[random_selection, ]

        uwot_clus_subset <- uwot::umap(
            X = subsample_matrix,
            n_neighbors = n_neighbors,
            n_components = n_components,
            n_epochs = n_epochs,
            min_dist = min_dist,
            n_threads = n_threads,
            spread = spread,
            ret_model = TRUE,
            ...
        )

        # create leftover matrix
        leftover_matrix <- matrix_to_use[-random_selection, ]

        # make prediction on leftover matrix
        uwot_clus_pred <- uwot::umap_transform(
            X = leftover_matrix, model = uwot_clus_subset
        )

        # combine subset and prediction
        coords_umap <- rbind(uwot_clus_subset$embedding, uwot_clus_pred)
        coords_umap <- coords_umap[
            match(cell_ID_order, rownames(coords_umap)), ]

        coords_umap_DT <- data.table::as.data.table(coords_umap)
        coords_umap_DT[, cell_ID := rownames(coords_umap)]

        # exit seed
        if (isTRUE(set_seed)) {
            set.seed(seed = Sys.time())
        }
    } else if (reduction == "feats") {
        message("Feats reduction is not yet implemented")
    }



    if (isTRUE(return_gobject)) {
        umap_names <- list_dim_reductions_names(
            gobject = gobject,
            data_type = reduction,
            spat_unit = spat_unit,
            feat_type = feat_type,
            dim_type = "umap"
        )

        if (name %in% umap_names) {
            message(name, " has already been used, will be overwritten")
        }


        coordinates <- coords_umap

        dimObject <- create_dim_obj(
            name = name,
            feat_type = feat_type,
            spat_unit = spat_unit,
            reduction = reduction,
            provenance = provenance,
            reduction_method = "umap",
            coordinates = coordinates,
            misc = NULL
        )


        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_dimReduction(gobject = gobject, dimObject = dimObject)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        ## update parameters used ##
        gobject <- update_giotto_params(
            gobject,
            description = "_umap",
            return_gobject = TRUE,
            toplevel = toplevel_params
        )
        return(gobject)
    } else {
        return(coords_umap_DT)
    }
}









#' @title Run tSNE dimensional reduction
#' @name runtSNE
#' @description run tSNE
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param reduction cells or feats
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name arbitrary name for tSNE run
#' @param feats_to_use if dim_reduction_to_use = NULL, which features to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param dims tSNE param: number of dimensions to return
#' @param perplexity tSNE param: perplexity
#' @param theta tSNE param: theta
#' @param do_PCA_first tSNE param: do PCA before tSNE (default = FALSE)
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param verbose verbosity of the function
#' @param ... additional tSNE parameters
#' @returns giotto object with updated tSNE dimension recuction
#' @details See \code{\link[Rtsne]{Rtsne}} for more information about these and 
#' other parameters. \cr
#' \itemize{
#'   \item Input for tSNE dimension reduction can be another dimension reduction (default = 'pca')
#'   \item To use gene expression as input set dim_reduction_to_use = NULL
#'   \item If dim_reduction_to_use = NULL, feats_to_use can be used to select a column name of
#'   highly variable genes (see \code{\link{calculateHVF}}) or simply provide a vector of genes
#'   \item multiple tSNE results can be stored by changing the \emph{name} of the analysis
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' runtSNE(g)
#' @export
runtSNE <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = c("cells", "feats"),
    dim_reduction_to_use = "pca",
    dim_reduction_name = NULL,
    dimensions_to_use = 1:10,
    name = NULL,
    feats_to_use = NULL,
    return_gobject = TRUE,
    dims = 2,
    perplexity = 30,
    theta = 0.5,
    do_PCA_first = FALSE,
    set_seed = TRUE,
    seed_number = 1234,
    verbose = TRUE,
    ...) {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    reduction <- match.arg(reduction, choices = c("cells", "feats"))


    # specify dim_reduction_name to use for pca input for tsne
    if (!is.null(dim_reduction_to_use)) {
        if (dim_reduction_to_use == "pca") {
            if (is.null(dim_reduction_name)) {
                if (feat_type == "rna") {
                    dim_reduction_name <- "pca"
                } else {
                    dim_reduction_name <- paste0(feat_type, ".", "pca")
                }
            }
        }
    }


    # specify name to use for umap
    if (is.null(name)) {
        if (feat_type == "rna") {
            name <- "tsne"
        } else {
            name <- paste0(feat_type, ".", "tsne")
        }
    }




    ## tsne on cells ##
    if (reduction == "cells") {
        ## using dimension reduction ##
        if (!is.null(dim_reduction_to_use)) {
            ## TODO: check if reduction exists
            dimObj_to_use <- getDimReduction(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                reduction = reduction,
                reduction_method = dim_reduction_to_use,
                name = dim_reduction_name,
                output = "dimObj"
            )

            provenance <- prov(dimObj_to_use)
            matrix_to_use <- dimObj_to_use[]
            matrix_to_use <- matrix_to_use[, dimensions_to_use]
        } else {
            ## using original matrix ##
            # expression values to be used
            values <- match.arg(
                expression_values, 
                unique(c("normalized", "scaled", "custom", expression_values)))
            expr_values <- getExpression(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                values = values,
                output = "exprObj"
            )

            provenance <- prov(expr_values)
            expr_values <- expr_values[] # extract matrix

            ## subset matrix
            if (!is.null(feats_to_use)) {
                expr_values <- .create_feats_to_use_matrix(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    sel_matrix = expr_values,
                    feats_to_use = feats_to_use,
                    verbose = verbose
                )
            }

            matrix_to_use <- t_flex(expr_values)
        }

        # start seed
        if (isTRUE(set_seed)) {
            set.seed(seed = seed_number)
        }

        ## run tSNE ##
        tsne_clus <- Rtsne::Rtsne(
            X = matrix_to_use,
            dims = dims,
            perplexity = perplexity,
            theta = theta,
            pca = do_PCA_first, ...
        )

        tsne_clus_pos_DT <- data.table::as.data.table(tsne_clus$Y)

        # data.table variables
        cell_ID <- NULL
        tsne_clus_pos_DT[, cell_ID := rownames(matrix_to_use)]

        # exit seed
        if (isTRUE(set_seed)) {
            set.seed(Sys.time())
        }


        if (isTRUE(return_gobject)) {
            tsne_names <- list_dim_reductions_names(
                gobject = gobject, data_type = reduction,
                spat_unit = spat_unit, feat_type = feat_type,
                dim_type = "tsne"
            )

            if (name %in% tsne_names) {
                cat(name, " has already been used, will be overwritten")
            }


            coordinates <- tsne_clus$Y
            rownames(coordinates) <- rownames(matrix_to_use)

            dimObject <- create_dim_obj(
                name = name,
                feat_type = feat_type,
                spat_unit = spat_unit,
                provenance = provenance,
                reduction = reduction,
                reduction_method = "tsne",
                coordinates = coordinates,
                misc = tsne_clus
            )

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- set_dimReduction(gobject = gobject, 
                                        dimObject = dimObject)
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

            ## update parameters used ##
            gobject <- update_giotto_params(gobject, description = "_tsne")
            return(gobject)
        } else {
            return(tsne_clus_pos_DT)
        }
    } else if (reduction == "feats") {
        message("Not yet implemented")
    }
}






## * Data Integration ####
# ---------------------- #


#' @title runGiottoHarmony
#' @name runGiottoHarmony
#' @description run UMAP
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param vars_use If meta_data is dataframe, this defines which variable(s) to
#'   remove (character vector).
#' @param do_pca Whether to perform PCA on input matrix.
#' @param expression_values expression values to use
#' @param reduction reduction on cells or features
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name arbitrary name for Harmony run
#' @param feats_to_use if dim_reduction_to_use = NULL, which feats to use
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param toplevel_params parameters to extract
#' @param verbose be verbose
#' @param ... additional \code{\link[harmony]{HarmonyMatrix}} parameters
#' @returns giotto object with updated Harmony dimension reduction
#' @details This is a simple wrapper for the HarmonyMatrix function in the 
#' Harmony package \doi{10.1038/s41592-019-0619-0}.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' 
#' runGiottoHarmony(g, vars_use = "leiden_clus")
#' @export
runGiottoHarmony <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    vars_use = "list_ID",
    do_pca = FALSE,
    expression_values = c("normalized", "scaled", "custom"),
    reduction = "cells",
    dim_reduction_to_use = "pca",
    dim_reduction_name = NULL,
    dimensions_to_use = 1:10,
    name = NULL,
    feats_to_use = NULL,
    set_seed = TRUE,
    seed_number = 1234,
    toplevel_params = 2,
    return_gobject = TRUE,
    verbose = NULL,
    ...) {
    # verify if optional package is installed
    package_check(pkg_name = "harmony", repository = "CRAN")


    # print message with information #
    message("using 'Harmony' to integrate different datasets. If used in 
            published research, please cite:")

    wrap_msg("Korsunsky, I., Millard, N., Fan, J. et al.
    Fast, sensitive and accurate integration of single-cell data with Harmony.
    Nat Methods 16, 1289-1296 (2019).
    https://doi.org/10.1038/s41592-019-0619-0",
        .initial = "  ",
        .prefix = "    "
    )


    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )


    # specify dim_reduction_name to use for pca input for umap
    if (!is.null(dim_reduction_to_use)) {
        if (dim_reduction_to_use == "pca") {
            if (is.null(dim_reduction_name)) {
                if (feat_type == "rna") {
                    dim_reduction_name <- "pca"
                } else {
                    dim_reduction_name <- paste0(feat_type, ".", "pca")
                }
            }
        }
    }


    # specify name to use for harmony
    if (is.null(name)) {
        if (feat_type == "rna") {
            name <- "harmony"
        } else {
            name <- paste0(feat_type, ".", "harmony")
        }
    }




    # set cores to use
    # n_threads = determine_cores(cores = n_threads)


    ## using dimension reduction ##
    if (!is.null(dim_reduction_to_use)) {
        ## TODO: check if reduction exists
        matrix_to_use <- get_dimReduction(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = reduction, # set to spat_unit?
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "dimObj"
        )
        provenance <- prov(matrix_to_use)
        matrix_to_use <- matrix_to_use[]

        matrix_to_use <- matrix_to_use[, dimensions_to_use]
    } else {
        ## using original matrix ##
        # expression values to be used
        values <- match.arg(
            expression_values, 
            unique(c("normalized", "scaled", "custom", expression_values)))
        expr_values <- getExpression(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            values = values,
            output = "exprObj"
        )

        provenance <- prov(expr_values)
        expr_values <- expr_values[] # extract matrix


        ## subset matrix
        if (!is.null(feats_to_use)) {
            expr_values <- .create_feats_to_use_matrix(
                gobject = gobject,
                feat_type = feat_type,
                spat_unit = spat_unit,
                sel_matrix = expr_values,
                feats_to_use = feats_to_use,
                verbose = verbose
            )
        }

        matrix_to_use <- t_flex(expr_values)
    }

    # get metadata
    metadata <- pDataDT(gobject, feat_type = feat_type, spat_unit = spat_unit)

    # start seed
    if (isTRUE(set_seed)) {
        set.seed(seed = seed_number)
        on.exit(GiottoUtils::random_seed())
    }

    # run harmony
    harmony_results <- harmony::RunHarmony(
        data_mat = matrix_to_use,
        meta_data = metadata,
        vars_use = vars_use,
        do_pca = do_pca,
        ...
    )


    colnames(harmony_results) <- paste0("Dim.", 1:ncol(harmony_results))
    rownames(harmony_results) <- rownames(matrix_to_use)

    harmdimObject <- create_dim_obj(
        name = name,
        spat_unit = spat_unit,
        feat_type = feat_type,
        provenance = provenance,
        reduction = "cells", # set to spat_unit?
        reduction_method = "harmony",
        coordinates = harmony_results,
        misc = NULL
    )

    # return giotto object or harmony results
    if (isTRUE(return_gobject)) {

        harmony_names <- list_dim_reductions_names(
            gobject = gobject,
            data_type = reduction,
            spat_unit = spat_unit,
            feat_type = feat_type,
            dim_type = "harmony"
        )

        if (name %in% harmony_names) {
            cat(name, 
                " has already been used with harmony, will be overwritten")
        }

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_dimReduction(gobject = gobject, 
                                    dimObject = harmdimObject)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


        ## update parameters used ##
        gobject <- update_giotto_params(gobject,
            description = "_harmony",
            return_gobject = TRUE,
            toplevel = toplevel_params
        )
        return(gobject)
    } else {
        return(harmdimObject)
    }
}
