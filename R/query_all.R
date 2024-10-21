#' Query all CpGs from input genomic regions
#'
#' Queries the provided regions and produces M and Coverage matrices and their
#' genomic positions. Parallelized across files using threads from the
#' `"iscream.threads"` option. The output of `query_all` may be used to create
#' a BSseq object: `do.call(BSseq, query_all(...))`.
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector of genomic regions strings
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt".
#' @param merged Whether the input strands have been merged/collapsed
#' @param sparse Whether to return M and coverage matrices as sparse matrices
#' ("dgCMatrix"). Set this `TRUE` only for scWGBS data
#' @param prealloc The number of rows to initialize the matrices with. If the
#' number of methyltion loci are approximately known, this can reduce runtime
#' as fewer resizes need to be made.
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
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
#' lapply(bedfiles, function(i) knitr::kable(read.table(i, col.names = colnames)))
#'
#' # make a vector of regions
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' query_all(bedfiles, regions)
#' # for BSseq object run
#' \dontrun{
#' library(bsseq)
#' do.call(BSseq, query_all(bedfiles, regions))
#' }
query_all <- function(
  bedfiles,
  regions,
  aligner = "biscuit",
  merged = TRUE,
  sparse = FALSE,
  prealloc = 10000,
  nthreads = NULL
) {

  verify_aligner_or_stop(aligner)
  verify_filetype(bedfiles, aligner, stop_on_error = TRUE)
  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  verify_regions_or_stop(regions)

  n_threads <- .get_threads(nthreads)
  validate_log_level(n_threads = n_threads)

  Cpp_query_all(
    bedfiles = bedfiles,
    regions = regions,
    bismark = aligner != "biscuit",
    merged = merged,
    sparse = sparse,
    prealloc = prealloc,
    nthreads = n_threads
  )
}

