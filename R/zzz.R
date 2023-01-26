# Run on library loading

# print version number
.onAttach = function(libname, pkgname) {
  packageStartupMessage('Giotto Suite ', utils::packageVersion('Giotto'))
}
