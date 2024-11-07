.onLoad <- function(libname, pkgname) {
  msg <- package_loader()
  setup_logger()
  invisible()
}

