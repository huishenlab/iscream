#' Query all CpGs from input genomic regions
#'
#' Queries the provided regions and produces M/beta and Coverage matrices and
#' their genomic positions. Parallelized across files using threads from the
#' `"iscream.threads"` option. The output of `query_all` may be used to create
#' a BSseq object: `do.call(BSseq, query_all(...))`.
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector, data frame or GenomicRanges of genomic regions. See
#' details.
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt".
#' @param mval Whether to return M-values or beta-values with the coverage
#' matrix. Defaults to M-value. Set `mval=FALSE` to get beta value matrix.
#' @param merged Whether the input strands have been merged/collapsed
#' @param sparse Whether to return M and coverage matrices as sparse matrices
#' ("dgCMatrix"). Set this `TRUE` only for scWGBS data
#' @param prealloc The number of rows to initialize the matrices with. If the
#' number of methyltion loci are approximately known, this can reduce runtime
#' as fewer resizes need to be made.
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#' @return A named list of
#' - coverage and either a beta- or M-value matrix
#' - a character vector of chromosomes and numeric vector of corresponding CpG
#' base positions
#' - a character vector of the input sample names
#'
#' @details
#' The input regions may be string vector in the form "chr:start-end"
#' or a GRanges object. If a data frame is provided, they must have "chr",
#' "start", and "end" columns.
#'
#' @importFrom Matrix drop0
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
  mval = TRUE,
  merged = TRUE,
  sparse = FALSE,
  prealloc = 10000,
  nthreads = NULL
) {

  verify_aligner_or_stop(aligner)
  verify_filetype(bedfiles, aligner, stop_on_error = TRUE)
  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  if (class(regions)[1] == "GRanges"){
    regions <- get_granges_string(regions)
  } else if ("data.frame" %in% class(regions)) {
    regions <- get_df_string(regions)
  }

  n_threads <- .get_threads(nthreads)
  validate_log_level(n_threads = n_threads)

  b <- Cpp_query_all(
    bedfiles = bedfiles,
    regions = regions,
    bismark = aligner != "biscuit",
    merged = merged,
    sparse = sparse,
    prealloc = prealloc,
    nthreads = n_threads
  )

  if (sparse) {
    get_beta_m <- ifelse(mval, get_m_sparse, get_beta_sparse)
  } else {
    get_beta_m <- ifelse(mval, get_m, get_beta)
  }

  if (sparse) {
    get_beta_m(b$M)
    b$M <- drop0(b$M)
    get_cov_sparse(b$Cov)
  } else {
    get_beta_m(b$M, n_threads)
    get_cov(b$Cov, n_threads)
  }
  if (!mval) {
    names(b)[which(names(b) == "M")] = "beta"
  }
  b
}

