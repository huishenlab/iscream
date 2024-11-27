#' Summarize CpGs  methylation information over genomic regions
#'
#' Run summarizing functions on the CpGs in bedfiles across genomic regions.
#' Parallelized across files using threads from the `"iscream.threads"` option.
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector, data frame or GenomicRanges of genomic regions. See
#' details.
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt".
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
#' dataframe. If input dataframes have a 'name' column, it will be used to
#' populate the output's feature column.
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
#' summarize_regions(bedfiles, regions)
#' summarize_regions(bedfiles, regions, fun = c("mean", "stddev"), mval = FALSE)
#' summarize_regions(bedfiles, regions, fun = "sum")
summarize_regions <- function(
  bedfiles,
  regions,
  fun = "all",
  aligner = "biscuit",
  mval = TRUE,
  set_region_rownames = FALSE,
  nthreads = NULL
) {

  supported_funcs <- c("sum", "mean", "median", "stddev", "variance", "min", "max", "range", "cpg_count")

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

  stopifnot("'mval' must be TRUE or FALSE" = mval %in% c(TRUE, FALSE))

  verify_aligner_or_stop(aligner)
  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  verify_filetype(bedfiles, aligner, stop_on_error = TRUE)
  if (class(regions)[1] == "GRanges"){
    regions <- get_granges_string(regions)
  } else if ("data.frame" %in% class(regions)) {
    regions <- get_df_string(regions)
  }

  n_threads <- .get_threads(nthreads)
  validate_log_level(n_threads = n_threads)

  df <- Cpp_summarize_regions(
    bedfiles = bedfiles,
    regions = regions,
    fun_vec = fun_to_use,
    bismark = aligner != "biscuit",
    mval = mval,
    region_rownames = set_region_rownames,
    nthreads = n_threads
  )
  df[df == -99] <- NA
  df

  mval_cpg_count <- paste0(ifelse(mval, "M", "beta"), ".cpg_count")
  if (mval_cpg_count %in% colnames(df)) {
    df <- df[, !(names(df) %in% "coverage.cpg_count")]
  }

  colnames(df)[which(colnames(df) == mval_cpg_count)] <- "cpg_count"
  df
}
#' Read WGBS BED file and get summary statistics
#'
#' @param bedfile Path to the BED file
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' @param merged Whether the input strands have been merged/collapsed
#' "bismark", "bsbolt".
#'
#' @return A data.table
#'
#' @importFrom data.table fread :=
#' @importFrom parallel mclapply
#' @importFrom stats var sd median
#'
#' @keywords internal
fsummarize <- function(bedfile, aligner, merged) {

  beta <- coverage <- U <- M <- NULL
  `.` <- list()

  if (aligner != "biscuit") {
    col_names <- c("chr", "start", "end", "beta", "U", "M")
  } else {
    col_names <- c("chr", "start", "end", "beta", "coverage")
  }

  if (merged) {
    dt <- fread(bedfile, col.names = col_names, drop = 6)
  } else {
    dt <- fread(bedfile, col.names = col_names)
  }

  if (aligner != "biscuit") {
      dt[, coverage := U + M][, beta := beta / 100]
  }

  dt[, .(
    Sample = basename(tools::file_path_sans_ext(bedfile, compression = T)),
    beta.mean = mean(beta),
    coverage.mean = mean(coverage),
    beta.median = median(beta),
    coverage.median = median(coverage),
    beta.min = min(beta),
    coverage.min = min(coverage),
    beta.max = max(beta),
    coverage.max = max(coverage),
    beta.range = max(beta) - min(beta),
    coverage.range = max(coverage) - min(coverage),
    beta.variance = var(beta),
    coverage.variance = var(coverage),
    beta.stddev = sd(beta),
    coverage.stddev = sd(coverage)
  )]
}

#' Summarize CpGs methylation information over full files
#'
#' @param bedfiles A `summarize_regions()` table with chromosomes as regions
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt".
#' @param merged Whether the input strands have been merged/collapsed
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @return A data.table
#'
#' @importFrom data.table rbindlist
#' @importFrom parallel mclapply
#'
#' @export
#'
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#'
#' summary <- summarize_samples(bedfiles, merged = FALSE, aligner = "biscuit")
summarize_samples <- function(bedfiles, aligner, merged, nthreads = NULL) {
  verify_filetype(bedfiles, aligner, stop_on_error = TRUE)
  mclapply(
    bedfiles,
    function(bedfile) {
      fsummarize(bedfile, aligner, merged)
    },
    mc.cores = .get_threads(nthreads)
  ) |> rbindlist()
}
