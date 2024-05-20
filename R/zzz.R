.onAttach <- function(libname, pkgname) {
  use_threads <- get_threads()

  msg <- paste0(
    "iscream using ",
    ifelse(use_threads == 1, "only ", ""),
    use_threads, " ",
    ifelse(use_threads > 1, "threads", "thread"),
    ". See `?set_threads` for information on multithreading."
  )
  packageStartupMessage(msg)
  invisible()
}
