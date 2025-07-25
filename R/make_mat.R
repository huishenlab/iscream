#' Make a matrix from a numeric column of BED files
#'
#' Queries the provided regions and produces a matrix along with genomic
#' positions as a named list (`make_mat()`), a `RangedSummarizedExperiment`
#' (`make_mat_se()`), `GRanges` (`make_mat_gr()`). Parallelized across files
#' using threads from the `"iscream.threads"` option.
#'
#' @param bedfiles A vector of BED file paths
#' @param regions A vector, data frame or GenomicRanges of genomic regions. See details.
#' @param column The index of the data column needed for the matrix
#' @param mat_name What to name the matrix in the returned object
#' @param sparse Whether to return a sparse matrix
#' @param prealloc The number of rows to initialize the matrices with. If the
#' number of loci are approximately known, this can reduce runtime as fewer
#' resizes need to be made.
#' @param nthreads Set the number of threads to use. Overrides the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @returns
#' - `make_mat()`: A named list of
#'   - the matrix with the value of interest
#'   - a character vector of chromosomes and numeric vector of base positions
#'   - a character vector of the input sample BED file names
#' - `make_mat_gr()`: if `GenomicRanges` is available, a `GRanges`
#' - `make_mat_se()`: if `SummarizedExperiment` is available, a `RangedSummarizedExperiment`
#'
#' @details
#' The input regions may be string vector in the form "chr:start-end"
#' or a GRanges object. If a data frame is provided, they must have "chr",
#' "start", and "end" columns.
#'
#' @importFrom Matrix drop0
#' @importFrom methods is
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
#' # make matrix of beta values
#' make_mat(bedfiles, regions, column = 4)
make_mat <- function(
  bedfiles,
  regions,
  column,
  mat_name = "value",
  sparse = FALSE,
  prealloc = 10000,
  nthreads = NULL
) {
  if (column < 4) {
    stop("`col` < 3 - must be a the index of a numeric data column not any of chr, start or end ")
  }
  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  if (is(regions, "GRanges")) {
    regions <- get_granges_string(regions)
  } else if ("data.frame" %in% class(regions)) {
    regions <- get_df_string(regions)
  }

  n_threads <- .get_threads(nthreads)
  validate_log_level(n_threads = n_threads)

  mat <- Cpp_query_all(
    bedfiles = bedfiles,
    regions = regions,
    aligner = "none",
    valInd = column,
    merged = FALSE,
    sparse = sparse,
    prealloc = prealloc,
    nthreads = n_threads
  )

  names(mat)[which(names(mat) == "M")] <- mat_name
  mat
}

#' @rdname make_mat
#' @export
make_mat_se <- function(
  bedfiles,
  regions,
  column,
  mat_name = "value",
  sparse = FALSE,
  prealloc = 10000,
  nthreads = NULL
) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("The 'SummarizedExperiment' package must be installed for this functionality")
  }
  mat <- make_mat(bedfiles, regions, column, mat_name, sparse, prealloc, nthreads)
  matlist <- list(mat[[mat_name]])
  names(matlist) <- mat_name
  gr <- getGR(mat$chr, mat$pos)
  SummarizedExperiment::SummarizedExperiment(assays = matlist, rowRanges = gr)
}

#' @rdname make_mat
#' @export
make_mat_gr <- function(
  bedfiles,
  regions,
  column,
  mat_name = "value",
  prealloc = 10000,
  nthreads = NULL
) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The 'GenomicRanges' package must be installed for this functionality")
  }
  mat <- make_mat(
    bedfiles,
    regions,
    column,
    mat_name,
    sparse = FALSE,
    prealloc,
    nthreads
  )
  gr <- getGR(mat$chr, mat$pos)
  GenomicRanges::mcols(gr) <- mat[[mat_name]]
  gr
}
