#' Query the chromosomes or seqnames from a vector of BED files
#' @param bedfiles The vector of bedfile paths
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#' @return A vector of seqnames
#'
#' @export
query_chroms <- function(bedfiles, nthreads = NULL) {
  verify_files_or_stop(bedfiles)

  n_threads <- .get_threads(nthreads)

  Cpp_query_chroms(bedfiles, n_threads)
}
