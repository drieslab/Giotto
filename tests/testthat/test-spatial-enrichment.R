# Characterization tests for the sign-matrix enrichment methods.
#
# These pin what PAGE, rank and hypergeometric produce *today*, before the
# analyzeParam refactor moves the arithmetic behind a verb generic. Nothing
# here asserts that the numbers are right -- only that they do not change. The
# reference values were produced by the pre-refactor implementations on the
# fixture below; a diff in any of them means the refactor changed behaviour.
#
# Where a test pins a bug rather than a behaviour it says so, and is updated
# by the commit that fixes it.

skip_if_no_mini <- function() {
    skip_if_not_installed("GiottoData")
}

# 634 genes x 624 spots. `python_path = NA` skips the conda probe, which this
# object does not need -- nothing below touches a python module.
.enrich_fixture <- function(n_markers = 150L, seed = 1L) {
    withr::local_options(giotto.use_conda = FALSE, giotto.verbose = FALSE)
    g <- suppressMessages(
        GiottoData::loadGiottoMini("visium", python_path = NA)
    )
    genes <- rownames(getExpression(g, values = "normalized",
                                    output = "matrix"))
    set.seed(seed)
    types <- c("typeA", "typeB", "typeC")
    sm <- matrix(0L, nrow = length(genes), ncol = length(types),
                 dimnames = list(genes, types))
    for (j in seq_along(types)) sm[sample(length(genes), n_markers), j] <- 1L
    list(g = g, sm = sm, genes = genes)
}


test_that("runPAGEEnrich output is unchanged", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    res <- runPAGEEnrich(f$g, sign_matrix = f$sm,
                         return_gobject = FALSE, verbose = FALSE)

    # PAGE is the odd one out: it returns list(DT=, matrix=) where the three
    # siblings return a bare spatEnrObj. Pinned so the refactor has to make a
    # deliberate choice about it rather than an accidental one.
    expect_type(res, "list")
    expect_named(res, c("DT", "matrix"))
    expect_s4_class(res$matrix, "spatEnrObj")

    dt <- res$matrix[]
    expect_identical(nrow(dt), 624L)
    expect_setequal(names(dt), c("cell_ID", "typeA", "typeB", "typeC"))
    expect_equal(head(dt$typeA, 3),
                 c(1.660300, -0.065675, -0.416552), tolerance = 1e-5)
    expect_equal(head(dt$typeB, 3),
                 c(0.411777, 1.100740, -0.997041), tolerance = 1e-5)
    expect_equal(head(dt$typeC, 3),
                 c(0.167189, 0.625430, 0.233294), tolerance = 1e-5)

    expect_identical(res$matrix@method, "PAGE")
    expect_identical(res$matrix@name, "PAGE")
    expect_setequal(names(res$matrix@misc), c(
        "expr_values_used", "reverse_log_scale", "logbase",
        "p_values_calculated", "output_enrichment_scores",
        "include_depletion", "nr_permutations"
    ))
})


test_that("runRankEnrich output is unchanged", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # `expression_values` must be given explicitly: the default vector cannot
    # survive this function's own match.arg call (see the bug test below).
    res <- runRankEnrich(f$g, sign_matrix = f$sm,
                         expression_values = "normalized",
                         return_gobject = FALSE)

    expect_s4_class(res, "spatEnrObj")
    dt <- res[]
    expect_identical(nrow(dt), 624L)
    expect_setequal(names(dt), c("cell_ID", "typeA", "typeB", "typeC"))
    expect_equal(head(dt$typeA, 3),
                 c(0.887582, 0.880433, 0.877802), tolerance = 1e-5)
    expect_equal(head(dt$typeB, 3),
                 c(0.882421, 0.883287, 0.875983), tolerance = 1e-5)
    expect_equal(head(dt$typeC, 3),
                 c(0.883911, 0.888095, 0.882159), tolerance = 1e-5)

    expect_identical(res@method, "rank")
    expect_identical(res@name, "rank")
})


test_that("runHyperGeometricEnrich output is unchanged", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    res <- runHyperGeometricEnrich(f$g, sign_matrix = f$sm,
                                   return_gobject = FALSE)

    expect_s4_class(res, "spatEnrObj")
    dt <- res[]
    expect_identical(nrow(dt), 624L)
    expect_setequal(names(dt), c("cell_ID", "typeA", "typeB", "typeC"))
    expect_equal(head(dt$typeA, 3),
                 c(1.288880, 0.326247, 0.427120), tolerance = 1e-5)
    expect_equal(head(dt$typeB, 3),
                 c(2.091817, 1.166187, 0.886146), tolerance = 1e-5)
    expect_equal(head(dt$typeC, 3),
                 c(0.694707, 0.819049, 0.632128), tolerance = 1e-5)

    expect_identical(res@method, "hypergeometric")
    expect_identical(res@name, "hypergeometric")
})


test_that("the enrichment result lands in the gobject", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    g2 <- runPAGEEnrich(f$g, sign_matrix = f$sm, verbose = FALSE)
    expect_s4_class(g2, "giotto")
    expect_true("PAGE" %in% list_spatial_enrichments_names(
        g2, spat_unit = "cell", feat_type = "rna"
    ))

    enr <- getSpatialEnrichment(g2, name = "PAGE", output = "spatEnrObj")
    direct <- runPAGEEnrich(f$g, sign_matrix = f$sm,
                            return_gobject = FALSE, verbose = FALSE)$matrix
    expect_equal(enr[], direct[])

    # every run appends one parameter-history entry
    expect_gt(length(g2@parameters), length(f$g@parameters))
    last <- g2@parameters[[length(g2@parameters)]]
    expect_identical(unname(last[["method used"]]), "PAGE")
})


test_that("runSpatialEnrich routes to the same result as the direct call", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    for (m in c("PAGE", "rank", "hypergeometric")) {
        via <- runSpatialEnrich(f$g, enrich_method = m, sign_matrix = f$sm,
            expression_values = "normalized",
            return_gobject = FALSE, verbose = FALSE)
        direct <- switch(m,
            PAGE = runPAGEEnrich(f$g, sign_matrix = f$sm,
                expression_values = "normalized",
                return_gobject = FALSE, verbose = FALSE)$matrix,
            rank = runRankEnrich(f$g, sign_matrix = f$sm,
                expression_values = "normalized", return_gobject = FALSE),
            hypergeometric = runHyperGeometricEnrich(f$g, sign_matrix = f$sm,
                expression_values = "normalized", return_gobject = FALSE)
        )
        if (m == "PAGE") via <- via$matrix
        expect_equal(via[], direct[], info = m)
    }
})


# --- bug fixes ---------------------------------------------------------------

test_that("runRankEnrich works on its own default expression_values", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # `expression_values` defaults to c("normalized", "raw", "scaled",
    # "custom") but the choices were built as unique(c("normalized", "scaled",
    # "custom", expression_values)) -- the same four in a different order. Not
    # identical to the arg, so match.arg refused a length-4 arg and the
    # function could not be called without naming a value explicitly.
    res <- runRankEnrich(f$g, sign_matrix = f$sm, return_gobject = FALSE)
    expect_s4_class(res, "spatEnrObj")

    explicit <- runRankEnrich(f$g, sign_matrix = f$sm,
        expression_values = "normalized", return_gobject = FALSE)
    expect_equal(res[], explicit[])
})


test_that("runPAGEEnrich honours output_enrichment", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # The wrapper passed the literal c("original", "zscore") down to
    # .page_dt_method(), which match.arg'd it back to "original". The user's
    # choice was discarded, so PAGE always returned unscaled scores.
    orig <- runPAGEEnrich(f$g, sign_matrix = f$sm, output_enrichment = "original",
        return_gobject = FALSE, verbose = FALSE)$matrix[]
    zsc <- runPAGEEnrich(f$g, sign_matrix = f$sm, output_enrichment = "zscore",
        return_gobject = FALSE, verbose = FALSE)$matrix[]

    expect_false(isTRUE(all.equal(orig$typeA, zsc$typeA)))
    # "zscore" standardizes within cell type
    expect_equal(mean(zsc$typeA), 0, tolerance = 1e-8)
    expect_equal(stats::sd(zsc$typeA), 1, tolerance = 1e-8)
    # the default is still "original"
    expect_equal(
        runPAGEEnrich(f$g, sign_matrix = f$sm,
            return_gobject = FALSE, verbose = FALSE)$matrix[],
        orig
    )
})


test_that("runRankEnrich(p_value = TRUE) returns p-values", {
    skip_if_no_mini()
    f <- .enrich_fixture()

    # The permutation branch recursed into runRankEnrich() without
    # return_gobject = FALSE, so it got a giotto object back and then subset it
    # as a table. p_value = TRUE errored in fitdistrplus every time.
    res <- runRankEnrich(f$g, sign_matrix = f$sm, p_value = TRUE,
                         n_times = 20, return_gobject = FALSE)
    expect_s4_class(res, "spatEnrObj")

    dt <- res[]
    score_cols <- c("typeA", "typeB", "typeC")
    for (cl in score_cols) {
        expect_true(all(dt[[cl]] >= 0 & dt[[cl]] <= 1), info = cl)
    }
    # cell_ID must survive as an ID, not be swept into the gamma transform
    expect_type(dt$cell_ID, "character")
    expect_identical(nrow(dt), 624L)

    # and p_value = FALSE still gives the scores, unchanged
    plain <- runRankEnrich(f$g, sign_matrix = f$sm, return_gobject = FALSE)
    expect_equal(head(plain[]$typeA, 3),
                 c(0.887582, 0.880433, 0.877802), tolerance = 1e-5)
})
