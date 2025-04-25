package_loader <- function() {
  threads <- get_threads()
  use_threads <- threads["use_threads"]
  opt_set <- threads["opt_set"]
  avail_threads <- threads["avail_threads"]
  options("iscream.threads" = unname(use_threads))
  msg <- "iscream using "
  help_msg <- " See `?set_threads` for information on multithreading"

  if (use_threads == avail_threads) {
    msg <- paste0(
      msg,
      use_threads,
      " threads of detected ",
      avail_threads,
      " threads.",
      help_msg,
      "."
    )
  } else {
    msg <- paste0(
      msg,
      use_threads,
      ifelse(use_threads > 1, " threads", " thread"),
      ifelse(opt_set, " based on 'options(iscream.threads)'", " by default"),
      " but parallelly::availableCores() detects ",
      avail_threads,
      " possibly available threads.",
      help_msg,
      " before trying to use more."
    )
  }

  if (Sys.which("tabix") == "") {
    options("tabix.method" = "htslib")
    msg <- paste0(
      msg,
      "\n'tabix' executable not found in $PATH.",
      " tabix() will use htslib to make queries instead which can be slower.",
      " See ?tabix for details."
    )
  } else {
    options("tabix.method" = "shell")
  }

  return(msg)
}

.onAttach <- function(libname, pkgname) {
  msg <- package_loader()
  packageStartupMessage(msg)
  setup_logger()
  invisible()
}
