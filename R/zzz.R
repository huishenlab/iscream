.onLoad <- function(libname, pkgname) {
  threads <- get_threads()
  use_threads <- threads["use_threads"]
  opt_set <- threads["opt_set"]
  avail_threads <- threads["avail_threads"]
  options("iscream.threads" = unname(use_threads))
  msg <- "iscream using "
  help_msg <- " See `?set_threads` for information on multithreading"

  if (use_threads == avail_threads) {
    msg <- paste0(
      msg, use_threads, " threads of detected ",
      avail_threads, " threads.", help_msg, "."
    )
  } else {
    msg <- paste0(
      msg, use_threads,
      ifelse(use_threads > 1, " threads", " thread"),
      ifelse(opt_set, " based on 'options(iscream.threads)'", " by default"),
      " but parallelly::availableCores() detects ", avail_threads,
      " possibly available threads.",
      help_msg, " before trying to use more.")
  }

  packageStartupMessage(msg)
  setup_logger()
  invisible()
}
