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
#' @param feature_col Column name of the input `regions` data frame containing
#' a name for each genomic region. Set only if the using a data frame as the
#' input regions format. See details.
#' @param set_region_rownames Use the region strings as the returned data
#' frame's rownames. Can be useful if you have a named regions and want both
#' the regions strings rownames and the feature names. See details.
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @details
#'
#' ## Supported `fun` functions
#'
#' - Sum: `"sum"`
#' - Mean: `"mean"`
#' - Median: `"median"`
#' - Standard deviation: `"stddev"`
#' - Variance: `"variance"`
#' - Minimum: `"min"`
#' - Maximum: `"max"`
#' - Range: `"range"`
#' - No. of records in the region: `"count"`
#'
#' The summarizing computations are backed by the Armadillo library. See
#' <https://arma.sourceforge.net/docs.html#stats_fns> for futher details on the
#' supported functions
#'
#' ## Using feature identifiers
#'
#' `regions` may be string vector in the form "chr:start-end", a GRanges
#' object or a data frame with "chr", "start", and "end" columns. The `feature`
#' column of the output will contain a "chr:start-end" identifier for each
#' summarized region. To use other identifiers, like a gene name for a region
#' instead of the coordinates, set the names of the vector or GRanges to those
#' identifiers. These names will be used instead of the genomic region string
#' to describe each feature in the output dataframe. If `regions` is a data
#' frame make an additional column with the identifiers and pass that column
#' name to `feature_col`. See examples.
#'
#' @returns A data.frame
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
#' names(regions) <- c("A", "B", "C")
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
#' _[, start := as.integer(start)] |>
#' _[, feature := LETTERS[.I]][]
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

  col_names <- col_names %||% paste0("V", seq_len(length(columns)))

  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
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
    col_indices = columns,
    col_names = col_names,
    fun_vec = fun_to_use,
    region_rownames = set_region_rownames,
    aligner = "general",
    mval = FALSE,
    nthreads = n_threads
  )
  df[df == -99] <- NA
  df

  count_colnames <- paste0(col_names, ".count")
  if (any(count_colnames %in% colnames(df))) {
    df <- df[, !(names(df) %in% count_colnames[-1])]
  }

  colnames(df)[which(colnames(df) == count_colnames[1])] <- "count"
  df
}
