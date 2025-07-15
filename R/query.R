#' Query the chromosomes or seqnames from a vector of BED files
#' @param bedfiles The vector of bedfile paths
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#' @returns A vector of seqnames
#'
#' @export
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#' query_chroms(bedfiles)
query_chroms <- function(bedfiles, nthreads = NULL) {
  verify_files_or_stop(bedfiles)

  n_threads <- .get_threads(nthreads)

  Cpp_query_chroms(bedfiles, n_threads)
}
