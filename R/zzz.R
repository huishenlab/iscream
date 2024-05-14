.onAttach <- function(libname, pkgname) {
  opt_threads <- getOption("iscream.threads")
  threads <- get_threads(verbose = FALSE)

  if (is.null(opt_threads)) {
    options(iscream.threads = threads["use threads"])
  }

  msg <- paste0("iscream using ", threads["use threads"], " threads of ", threads["available threads"], ". See ?get_threads")
  packageStartupMessage(msg)
  invisible()
}

