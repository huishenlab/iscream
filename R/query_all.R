#' Query all CpGs from input genomic regions
#'
#' Queries the provided regions and produces M and Coverage matrices and their
#' genomic positions. Parallelized across files using threads from the
#' `"iscream.threads"` option. The output of `query_all` may be used to create
#' a BSseq object: `do.call(BSseq, query_all(...))`.
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector of genomic regions strings
#' @param bismark Whether the input is a bismark coverage file
#' @param merged Whether the input strands have been merged/collapsed
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#' @importFrom fs file_exists
#' @return A named list of
#' - M and coverage matrices
#' - a character vector of chromosomes and numeric vector of corresponding CpG
#' base positions
#' - a character vector of the input sample names
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
#' # for BSseq object run
#' \dontrun{
#' library(bsseq)
#' do.call(BSseq, query_all(bedfiles, regions))
#' }
query_all <- function(bedfiles, regions, bismark = FALSE, merged = TRUE, nthreads = NULL) {

  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  verify_regions_or_stop(regions)

  n_threads <- ifelse(
    is.null(nthreads),
    getOption("iscream.threads"),
    check_thread_count(nthreads)
  )
  Cpp_query_all(bedfiles, regions, bismark, merged, nthreads = n_threads)
}

