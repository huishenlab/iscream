#' Summarize information over genomic regions from any BED file
#'
#' Run summarizing functions on BED file records across genomic regions.
#' Parallelized across files using threads from the `"iscream.threads"` option.
#' @param bedfiles A vector of BED file paths
#' @param regions A vector, data frame or GenomicRanges of genomic regions. See
#' details.
#' @param columns A vector of indices of the numeric columns to be summarized
#' @param col_names A vector of names to use for `columns` in the output
#' @param fun Function(s) to apply over the region. See details.
#' @param feature_col Column name of the input `regions` data frame or GRanges
#' `mcols` containing a name for each genomic region. Set only if the using a
#' data frame-like object as the input regions format, not for string vectors.
#' See details.
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @details
#'
#' # Supported functions
#'
#' - Sum: `"sum"`
#' - Mean: `"mean"`
#' - Median: `"median"`
#' - Mode: `"mode"`
#' - Anti-mode: `"antimode"`
#' - Sample standard deviation: `"stddev"` (`sstdev` in `bedtools map`)
#' - Population standard deviation: `"pstddev"` (`stdev` in `bedtools map`)
#' - Variance: `"variance"`
#' - Minimum: `"min"`
#' - Maximum: `"max"`
#' - Minimum of absolute values: `"absmin"`
#' - Maximum of absolute values: `"absmax"`
#' - Range: `"range"`
#' - First element: `"first"`
#' - Last element: `"last"`
#' - No. of records in the region: `"count"`
#' - No. of records in the region with unique data values: `"count_distinct"`
#'
#' Most summarizing computations are backed by the Armadillo library. See
#' <https://arma.sourceforge.net/docs.html#stats_fns> for futher details on the
#' supported functions
#'
#' # Using feature identifiers
#'
#' `regions` may be a string vector in the form "chr:start-end", a GRanges
#' object or a data frame with "chr", "start", and "end" columns. If the input
#' data.frame or GRanges has a column with feature identifiers, like gene names
#' for a set of gene regions, pass that column's name to `feature_col`. If
#' `regions` is a vector, set its `names()` to those identifiers. These will be
#' used to populate a 'feature' column in the summary. See examples.
#'
#' @returns A data.table
#'
#' @importFrom methods is
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
#' summarize_regions(bedfiles, regions, columns = c(4, 5), col_names = c("beta", "cov"))
#'
#' # select functions
#' summarize_regions(
#'   bedfiles,
#'   regions,
#'   fun = c("mean", "stddev"),
#'   columns = c(4, 5),
#'   col_names = c("beta", "cov")
#' )
#'
#' # add names to the regions
#' names(regions) <- c("gene1", "gene2", "gene3")
#' summarize_regions(
#'   bedfiles,
#'   regions,
#'   fun = "sum",
#'   columns = 5,
#'   col_names = "coverage"
#' )
#'
#' # using `feature_col`
#' library(data.table)
#'
#' # convert string vector to a data.table
#' regions_df <- data.table::as.data.table(regions) |>
#' _[, tstrsplit(regions, ":|-", fixed = FALSE, names = c("chr", "start", "end"))] |>
#' _[, feature := names(regions)][]
#' regions_df
#'
#' summarize_regions(
#'   bedfiles,
#'   regions_df,
#'   fun = "sum",
#'   columns = 5,
#'   col_names = "coverage",
#'   feature_col = "feature"
#' )
summarize_regions <- function(
  bedfiles,
  regions,
  columns,
  col_names = NULL,
  fun = "all",
  feature_col = NULL,
  nthreads = NULL
) {
  n_threads <- .get_threads(nthreads)
  validate_log_level(n_threads = n_threads)
  verify_files_or_stop(bedfiles, verify_tabix = TRUE)

  supported_funcs <- c(
    "sum",
    "mean",
    "median",
    "mode",
    "antimode",
    "stddev",
    "pstddev",
    "variance",
    "min",
    "max",
    "absmin",
    "absmax",
    "range",
    "first",
    "last",
    "count_distinct",
    "count"
  )

  fun_to_use <- validate_summary_function(fun, supported_funcs)

  regions_str <- get_string_input_regions(regions, feature_col)
  regions_df <- get_regions_as_df(regions, regions_str)

  col_names <- col_names %||% paste0("V", seq_len(length(columns)))

  df <- Cpp_summarize_regions(
    bedfiles = bedfiles,
    regions = regions_str,
    col_indices = columns,
    col_names = col_names,
    fun_vec = fun_to_use,
    regions_df = regions_df,
    aligner = "general",
    mval = FALSE,
    nthreads = n_threads
  )
  setDT(df)
  df[df == -99] <- NA

  if ("count" %in% fun_to_use) {
    count_colnames <- paste0(col_names, ".count")
    df[, eval(count_colnames[-1]) := NULL]
    colnames(df)[which(colnames(df) == count_colnames[1])] <- "count"
  }
  df
}

#' Return functions to use if the input summarizing functions are valid
#' @keywords internal
validate_summary_function <- function(fun, supported_funcs) {
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
  fun_to_use
}

#' Return the input regions as a data frame so the output can have positions as columns
#' @keywords internal
get_regions_as_df <- function(regions, regions_str) {
  start <- end <- NULL
  if (!is(regions, "data.frame")) {
    return(get_df_from_string(regions_str))
  } else if (is(regions, "data.table")) {
    regions_df <- regions
  } else {
    regions_df <- as.data.table(regions)
  }
  regions_df[, `:=`(start = as.integer(start), end = as.integer(end))]
}
