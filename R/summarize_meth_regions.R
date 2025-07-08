#' Summarize CpGs  methylation information over genomic regions
#'
#' Run summarizing functions on the CpGs in bedfiles across genomic regions.
#' Parallelized across files using threads from the `"iscream.threads"` option.
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector, data frame or GenomicRanges of genomic regions. See
#' details.
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt".
#' @param feature_col If the input is a dataframe, the column to use as the
#' feature label instead of the genomic region string
#' @param fun Function(s) to apply over the region. See details.
#' @param mval Whether to calculate the M value (coverage \eqn{\times \beta})
#' or use the beta value when applying the function.
#' @param set_region_rownames Use the region strings as the returned data
#' frame's rownames. Can be useful if you have a named regions and want both
#' the regions strings rownames and the feature names.
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @details
#' The input regions may be string vector in the form "chr:start-end"
#' or a GRanges object. If a data frame is provided, they must have "chr",
#' "start", and "end" columns. If the string vector and GenomicRanges inputs
#' are named, the names will be used to describe each feature in the output
#' dataframe. If input dataframes have a feature column, set `feature_col` to
#' that column name to populate the output's feature column.
#'
#' Supported `fun` arguments are given below. For each of these functions,
#' setting `mval = FALSE` will use the beta values instead of the M value:
#' - Sum: `"sum"`
#' - Mean: `"mean"`
#' - Median: `"median"`
#' - Standard deviation: `"stddev"`
#' - Variance: `"variance"`
#' - Minimum: `"min"`
#' - Maximum: `"max"`
#' - Range: `"range"`
#' - No. of CpGs in the region: `"cpg_count"`
#'
#' The summarizing computations are backed by the Armadillo library. See
#' <https://arma.sourceforge.net/docs.html#stats_fns> for futher details on the
#' supported functions
#'
#' @importFrom methods is
#'
#' @return A data.frame
#'
#' @export
#'
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#' # examine the bedfiles
#' colnames <- c("chr", "start", "end", "beta", "coverage")
#' lapply(bedfiles, function(i) knitr::kable(read.table(i, col.names = colnames)))
#'
#' # make a vector of regions
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' summarize_meth_regions(bedfiles, regions)
#' names(regions) <- c("A", "B", "C")
#' summarize_meth_regions(bedfiles, regions, fun = c("mean", "stddev"), mval = FALSE)
#' summarize_meth_regions(bedfiles, regions, fun = "sum")
summarize_meth_regions <- function(
  bedfiles,
  regions,
  fun = "all",
  aligner = "biscuit",
  feature_col = NULL,
  mval = TRUE,
  set_region_rownames = FALSE,
  nthreads = NULL
) {
  supported_funcs <- c("sum", "mean", "median", "stddev", "variance", "min", "max", "range", "count")

  if (length(fun) > 1) {
    if ("all" %in% fun) {
      stop("'all' can't be used with other summary funcions")
    }
    stopifnot("Selected function not supported" = all(fun %in% supported_funcs))
    fun_to_use <- fun
  } else {
    stopifnot("Selected function not supported" = fun %in% c(supported_funcs, "all"))
    fun_to_use <- supported_funcs
    if (fun != "all") {
      fun_to_use <- fun
    }
  }

  if (aligner != "general") {
    col_names <- c("coverage", ifelse(mval, "M", "beta"))
  }

  stopifnot("'mval' must be TRUE or FALSE" = mval %in% c(TRUE, FALSE))

  verify_aligner_or_stop(aligner)
  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  verify_filetype(bedfiles, aligner, stop_on_error = TRUE)
  if (is(regions, "GRanges")) {
    regions <- get_granges_string(regions)
  } else if ("data.frame" %in% class(regions)) {
    regions <- get_df_string(regions, feature_col)
  }

  n_threads <- .get_threads(nthreads)
  validate_log_level(n_threads = n_threads)

  df <- Cpp_summarize_regions(
    bedfiles = bedfiles,
    regions = regions,
    col_indices = c(4, 5),
    col_names = col_names,
    fun_vec = fun_to_use,
    aligner = aligner,
    mval = mval,
    region_rownames = set_region_rownames,
    nthreads = n_threads
  )
  df[df == -99] <- NA

  mval_count <- paste0(ifelse(mval, "M", "beta"), ".count")
  if (mval_count %in% colnames(df)) {
    df <- df[, !(names(df) %in% "coverage.count")]
  }

  colnames(df)[which(colnames(df) == mval_count)] <- "cpg_count"
  df
}
