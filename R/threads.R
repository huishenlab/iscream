#' Get the number of available threads
#'
#' Gets the number of threads iscream is currently set to use, whether the
#' `"iscream.threads"` option is set and how many threads are available for
#' use. To set the number of threads use `set_threads()` or set the
#' `iscream.threads` option in your `~/.Rprofile`. See `?set_threads` for more
#' information.
#' @importFrom parallelly availableCores
#' @returns A named vector:
#' - `use_threads` = the number of threads iscream will use
#' - `opt_set` = whether the option was set by the user
#' - `avail_threads` = The number of available threads as reported by
#' `parallelly::availableCores`
#'
#' @export
#'
#' @examples
#' get_threads()
get_threads <- function() {
  opt_threads <- unname(getOption("iscream.threads"))
  avail_threads <- unname(availableCores())
  if (is.null(opt_threads)) {
    return(c("use_threads" = 1, "opt_set" = FALSE, "avail_threads" = avail_threads))
  }
  check_thread_count(opt_threads, avail_threads, TRUE)
  c("use_threads" = opt_threads, "opt_set" = TRUE, "avail_threads" = avail_threads)
}

# Get the number of threads to use in functions from param or option
.get_threads <- function(nthreads) {
  ifelse(
    is.null(nthreads),
    getOption("iscream.threads"),
    check_thread_count(nthreads)
  )
}

#' Set the number of available threads
#'
#' Sets the `"iscream.threads"` option to `n_threads`. To see how many threads
#' you have available see `?get_threads()`.

#' @param n_threads The number of threads to use
#' @importFrom parallelly availableCores
#' @returns None. Sets the `iscream.threads` option to the requested number of
#' threads if available
#'
#' @details iscream uses OpenMP to parallelize certain functions. You can use
#' as many threads as are available to you on your system to varying degrees of
#' performance improvements. The `get_threads()` function uses
#' `parallelly::availableCores()` to report the number of available threads.
#' Although OpenMP can detect the number of available cores, on high
#' preformance computers (HPCs) with resource allocating job schedulers like
#' SLURM, OpenMP may detect all available threads across the HPC and not limit
#' itself to the cores that were allocated to you by the scheduler. If your
#' system administrator has not set up any limits, this may result in your job
#' taking resources from other jobs. If there are limits, trying to use more
#' threads that those available will reduce iscream's performance. Job
#' schedulers will typically have an environment variable (e.g.
#' `SLURM_CPUS_ON_NODE` with SLURM) that gives you the actual number of
#' available cores. Further, on hyperthreaded systems, this count may be double
#' that of the available processors. Using hyperthreading does not guarantee
#' any performance improvement - it may be better to set the number of threads
#' to half the reported number. `parallelly::availableCores()` takes HPC
#' scheduler/CRAN/Bioconductor limits into account when reporting the number of
#' available threads but it may not reliably report hyperthreading ('system' or
#' 'nproc'). To set the number of threads without having to call
#' `set_threads()` in every session, put
#'
#' ```
#' options(iscream.threads = [n_threads])
#' ```
#' in your `.Rprofile` See `help('Rprofile')` for information on startup options.
#'
#' Functions currently using multithreading:
#' - `tabix()`, `tabix_gr(), 'tabix_raw()`
#' - `query_chroms()`
#' - `make_mat()`, `make_mat_se()`, `make_mat_gr()`, `make_mat_bsseq()`
#' - `summarize_regions()`, `summarize_meth_regions()`
#'
#' @export
#'
#' @examples
#' (ncores <- parallelly::availableCores())
#' \dontrun{
#' set_threads(ncores)
#' }
set_threads <- function(n_threads) {
  check_thread_count(n_threads)
  options(iscream.threads = n_threads)
  message("iscream now using ", n_threads, " of ", availableCores(), " available threads.")
}

#' Check that the required threads are available
#'
#' @param n_threads The number of threads to check availability for
#' @param avail_threads The number of threads that are available on the system.
#' Defaults to `parallelly::availableCores()`
#' @param opt_set Whether the `iscream.threads` options is set
#'
#' @importFrom parallelly availableCores
#' @returns `n_threads` if the requested number of threads are available and
#' stops if not
#'
#' @keywords internal
check_thread_count <- function(
  n_threads,
  avail_threads = availableCores(),
  opt_set = FALSE
) {
  if (n_threads <= avail_threads) {
    return(n_threads)
  }
  if (opt_set) {
    msg <- paste0("The `iscream.threads` option is set to ", n_threads, " threads but your")
  } else {
    msg <- paste0("Cannot use ", n_threads, " threads. Your")
  }
  stop(
    msg,
    " system has only ",
    avail_threads,
    " threads. See parallelly::availableCores(which = 'all') for more informaion on available resources"
  )
}
