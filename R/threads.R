#' Get the number of available threads
#'
#' This function sets the number of threads iscream uses to 1 if the
#' `iscream.threads` option is not set. To set the number of threads use
#' setDTthreads or set the `iscream.threads` option in your `~/.Rprofile`. See
#' `?set_threads` for more information.
#' @param verbose Whether to be verbose on available omp threads
#' @importFrom parallelly availableCores
#'
#' @export
#' @examples
#' get_threads()
get_threads <- function() {
  opt_threads <- getOption("iscream.threads")
  avail_threads <- availableCores()
  if (is.null(opt_threads)) {
    return(1)
  }

  if (avail_threads < opt_threads) {
    stop(
      "The `iscream.threads` option is set to ", opt_threads,
      " but your system has only ", avail_threads,
      " threads. See parallelly::availableCores(which = 'all') for more information on available resources."
    )
  }
  opt_threads
}

#' Set the number of available threads
#'
#' This function sets the 'iscream.thread' option to `n_threads`. To
#' see how many threads you have available see `parallelly::availableCores`.
#' @param n_threads The number of threads to use
#' @importFrom parallelly availableCores
#'
#' @details iscream uses OpenMP to parallelize certain functions. You can use
#' as many threads as are available to you on your system to varying degrees of
#' performance improvements. The `parallelly::availableCores()` function will
#' report the number of available threads. To get more information about OpenMP
#' threads, run `get_omp_threads(verbose = TRUE)`. However, on high preformance
#' computers (HPCs) with resource allocating job schedulers like SLURM, OpenMP
#' may detect all available threads across the HPC and not limit itself to the
#' cores that were allocated to you by the scheduler. If your system
#' administrator has not set up any limits, this may result in your job taking
#' resources from other jobs. If there are limits, trying to use more threads
#' that those available will reduce iscream's performance. Job schedulers will
#' typically have an environment variable (e.g. SLURM_CPUS_ON_NODE with SLURM)
#' that gives you the actual number of available cores. Further, on
#' hyperthreaded systems, this count may be double that of the available
#' processors. Using hyperthreading does not guarantee any performance
#' improvement - it may be better to set the number of threads to half the
#' reported number. `parallelly::availableCores()` takes HPC
#' scheduler/CRAN/Bioconductor limits into account when reporting the number of
#' available threads but it may not reliably report hyperthreading ('system' or
#' 'nproc'). To set the number of threads without having to call
#' `set_threads()` in every session, put
#'
#' ```
#' options(iscream.threads = [n_threads])
#' ```
#' in your .Rprofile. See help('Rprofile') for information on startup options.
#'
#' Functions currently using OpenMP:
#' - `region_map()`
#' @export
#' @examples
#' \dontrun{
#' set_threads(4)
#' }
set_threads <- function(n_threads) {
  avail_threads <- availableCores()
  if (avail_threads < n_threads) {
    stop(paste0(
      "Cannot use ", n_threads,
      " threads. Your system has only ", avail_threads,
      " threads. See parallelly::availableCores(which = 'all') for more informaion on available resources"
    ))
  }
  options(iscream.threads = n_threads)
  message(paste0("iscream now using ", n_threads, " of ", avail_threads, " available threads."))
}
