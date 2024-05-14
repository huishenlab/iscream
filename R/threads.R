#' Get the numebr of available threads
#'
#' This function returns the half the number of available OpenMP threads. To
#' set the number of threads use setDTthreads sets it to half the available
#' number.
#' @param verbose Whether to be verbose on available omp threads
#'
#' @export
get_threads <- function(verbose = FALSE) {
  avail_threads <- get_omp_threads(verbose = FALSE)
  use_threads <- ifelse(avail_threads < 1, 1, round(avail_threads / 2))
  return(c("available threads" = avail_threads, "use threads" = use_threads))
}

#' Set the numebr of available threads
#'
#' This function sets the 'iscream.thread' option to `n_threads`. To
#' see how many threads are available, run `get_omp_threads`.
#' @param n_threads The number of threads to use
#'
#' @export
set_threads <- function(n_threads) {
  avail_threads <- get_omp_threads(verbose = FALSE)
  if (avail_threads < n_threads) {
    stop(paste0("Cannot use ", n_threads, ". Your system has only ", avail_threads, " threads."))
  }
  options(iscream.threads = n_threads)
}

