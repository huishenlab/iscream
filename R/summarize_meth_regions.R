#' Summarize methylation information over genomic regions
#'
#' Run summarizing functions on the CpG/CpH loci in BED files across genomic
#' regions. Parallelized across files using threads from the
#' `"iscream.threads"` option.
#'
#' @inheritParams summarize_regions
#'
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt".
#' @param mval Whether to calculate the M value (coverage \eqn{\times \beta})
#' or use the beta value when applying the function.
#'
#' @inheritSection summarize_regions Supported functions
#'
#' @inheritSection summarize_regions Using feature identifiers
#'
#' @importFrom methods is
#'
#' @returns A data.table
#'
#' @export
#'
#' @examples
#' # also see examples from ?summarize_regions
#'
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#'
#' # make a vector of regions
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' summarize_meth_regions(bedfiles, regions)
#'
#' # add names to the regions to populate the 'feature' column
#' names(regions) <- c("gene1", "gene2", "gene3")
#' summarize_meth_regions(bedfiles, regions, fun = c("mean", "stddev"), mval = FALSE)
#' summarize_meth_regions(bedfiles, regions, fun = "sum")
summarize_meth_regions <- function(
  bedfiles,
  regions,
  fun = "all",
  aligner = "biscuit",
  feature_col = NULL,
  mval = TRUE,
  nthreads = NULL
) {
  n_threads <- .get_threads(nthreads)
  validate_log_level(n_threads = n_threads)
  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  verify_aligner_or_stop(aligner)
  verify_filetype(bedfiles, aligner, stop_on_error = TRUE)

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

  if (aligner != "general") {
    col_names <- c("coverage", ifelse(mval, "M", "beta"))
  }
  stopifnot("'mval' must be TRUE or FALSE" = mval %in% c(TRUE, FALSE))

  df <- Cpp_summarize_regions(
    bedfiles = bedfiles,
    regions = regions_str,
    col_indices = c(4, 5),
    col_names = col_names,
    fun_vec = fun_to_use,
    regions_df = regions_df,
    aligner = aligner,
    mval = mval,
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
