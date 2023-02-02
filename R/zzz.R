# Run on library loading

# print version number
.onAttach = function(libname, pkgname) {
  packageStartupMessage('Giotto Suite ', utils::packageVersion('Giotto'))

  check_ver = getOption('giotto.check_version', TRUE)
  if(!isTRUE(check_ver)) {
    check_github_suite_ver()
    options('giotto.check_version' = FALSE)
  }

}
