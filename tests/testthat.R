# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(Giotto)

# additional needed packages
suppressWarnings({
  if(!require(remotes)){
    install.packages('R.utils', repos = 'http://cran.us.r-project.org')
    install.packages('remotes', repos = 'http://cran.us.r-project.org')
  }

  if(!require(GiottoData)){
    remotes::install_github('drieslab/GiottoData', build = FALSE)
  }
})

# install giotto environment
if (!checkGiottoEnvironment()) {
  installGiottoEnvironment()
}

test_check("Giotto")
