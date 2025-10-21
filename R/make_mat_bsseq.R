#' Make M/beta and coverage matrices from WGBS BED files
#'
#' Queries the CpG/CpH loci from provided regions and produces M/beta and
#' coverage matrices with their genomic positions. Parallelized across files
#' using threads from the `"iscream.threads"` option. The output of
#' `make_mat_bsseq` may be used to create a BSseq object: `do.call(BSseq,
#' make_mat_bsseq(...))`.
#'
#' @inheritParams make_mat
#'
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt".
#' @param mval Whether to return M-values or beta-values with the coverage
#' matrix. Defaults to M-value. Set `mval=FALSE` to get beta value matrix.
#' @param merged Whether the input strands have been merged/collapsed
#'
#' @returns A named list of
#' - coverage and either a beta- or M-value matrix
#' - a character vector of chromosomes and numeric vector of corresponding CpG
#' base positions
#' - a character vector of the input sample names
#'
#' @inherit make_mat details
#'
#' @section Bitpacking limits:
#' `make_mat_bsseq()` makes two matrices: M-value (or beta-value) and coverage.
#' For speed and memory efficiency these two values are bitpacked during matrix
#' creation so that only one matrix needs to be populated and resized. This
#' matrix is unpacked into the two required matrices only after the matrix
#' dimensions are known after querying all input files. The two values are
#' packed using the INT16 type, which has an upper limit of 32,767, into one
#' INT32. If the coverage values exceed 32,767, the upper limit of a 16-bit
#' signed integer, it will be capped at the limit. Beta values will also be
#' capped similarly, but any such beta values would indicate a bug in the
#' aligner that produced the data.
#'
#' @importFrom Matrix drop0
#' @importFrom methods is
#'
#' @export
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#' # examine the BED files
#' colnames <- c("chr", "start", "end", "beta", "coverage")
#' lapply(bedfiles, function(i) knitr::kable(read.table(i, col.names = colnames)))
#'
#' # make a vector of regions
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' mat <- make_mat_bsseq(bedfiles, regions)
#' # for BSseq object run
#' if (requireNamespace("bsseq", quietly = TRUE)) {
#'   do.call(bsseq::BSseq, mat)
#' }
make_mat_bsseq <- function(
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
  if (is(regions, "GRanges")) {
    regions <- get_granges_string(regions)
  } else if ("data.frame" %in% class(regions)) {
    regions <- get_df_string(regions)
  }

  n_threads <- .get_threads(nthreads)
  validate_log_level(n_threads = n_threads)

  b <- Cpp_query_all(
    bedfiles = bedfiles,
    regions = regions,
    aligner = aligner,
    valInd = 0,
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
    names(b)[which(names(b) == "M")] <- "beta"
  }
  b
}
