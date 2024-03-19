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
  # needed for reading .gz files.
  if(!require(R.utils)) {
    install.packages('R.utils', repos = 'http://cran.us.r-project.org')
  }
})


test_check("Giotto")
