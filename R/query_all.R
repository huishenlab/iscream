#' Make a bsseq object of CpGs from input genomic regions
#'
#' Queries the provided regions and produces M and Coverage matrices that are
#' wrapped in a BSSeq object
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector of genomic regions strings
#' @param nthreads Set number of threads to use. Should not be necessary as this
#' is set by `option('iscream.threads')`
#' @importFrom fs file_exists
#' @return A data.frame
#'
#' @export
query_all <- function(bedfiles, regions, nthreads = NULL) {

  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  verify_regions_or_stop(regions)

  n_threads <- ifelse(
    is.null(nthreads),
    getOption("iscream.threads"),
    check_thread_count(nthreads)
  )
  Cpp_query_all(bedfiles, regions, nthreads = n_threads)
}

