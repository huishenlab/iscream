#' Make a bsseq object of CpGs from input genomic regions
#'
#' Queries the provided regions and produces M and Coverage matrices that are
#' wrapped in a BSSeq object. Parallelized across files using threads from the
#' `"iscream.threads"` option.
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector of genomic regions strings
#' @param bismark Whether the input is a bismark coverage file
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#' @importFrom fs file_exists
#' @return A data.frame
#'
#' @export
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#' # examine the bedfiles
#' colnames <- c("chr", "start", "end", "beta", "coverage")
#' lapply(bedfiles, function(i) knitr::kable(data.table::fread(i, col.names = colnames)))
#'
#' # make a vector of regions
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' query_all(bedfiles, regions)
query_all <- function(bedfiles, regions, bismark = FALSE, nthreads = NULL) {

  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  verify_regions_or_stop(regions)

  n_threads <- ifelse(
    is.null(nthreads),
    getOption("iscream.threads"),
    check_thread_count(nthreads)
  )
  Cpp_query_all(bedfiles, regions, bismark, nthreads = n_threads)
}

