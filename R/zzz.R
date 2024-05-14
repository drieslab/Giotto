# Run on library loading


.onAttach <- function(libname, pkgname) {
    ## print version number ##
    packageStartupMessage("Giotto Suite ", utils::packageVersion("Giotto"))

    check_ver <- getOption("giotto.check_version", TRUE)
    if (isTRUE(check_ver)) {
        check_github_suite_ver()
        options("giotto.check_version" = FALSE)
    }


    ## cores detection ##
    check_core <- getOption("giotto.check_core", TRUE)
    if (isTRUE(check_core)) {
        cores <- determine_cores(cores = NA)
        data.table::setDTthreads(threads = cores)
        options("giotto.check_core" = FALSE)
    }
}
