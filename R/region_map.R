#' Run a defined function over genomic regions
#'
#' Run a function on the CpGs in bedfiles across genomic regions. Currently
#' supported functions are aggregate and average. Parallelized across files
#' using threads from the `"iscream.threads"` option.
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector of genomic regions strings
#' @param bismark Whether the input is a bismark coverage file
#' @param fun Function to apply over the region. See details.
#' @param mval Whether to calculate the M value (coverage \eqn{\times \beta})
#' or use the beta value
#' when applying the function.
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#' @details
#' Available functions:
#' - `"aggregate"` sums the values in the region with aggregated beta values if
#' `mval =` FALSE and aggregated M values if `mval =` TRUE
#' - `"average"` averages the values in the region with average beta values if
#' `mval =` FALSE and average M values if `mval =` TRUE
#' @importFrom fs file_exists
#' @return A data.frame
#'
#' @export
#'
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#' # examine the bedfiles
#' colnames <- c("chr", "start", "end", "beta", "coverage")
#' lapply(bedfiles, function(i) knitr::kable(data.table::fread(i, col.names = colnames)))
#'
#' # make a vector of regions
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' region_map(bedfiles, regions)
#' region_map(bedfiles, regions, mval = FALSE)
#' region_map(bedfiles, regions, fun = "average")
region_map <- function(bedfiles, regions, fun = "aggregate", bismark = FALSE, mval = TRUE, nthreads = NULL) {

  supported_funcs <- c("aggregate", "average")
  stopifnot("Selected function not supported" = fun %in% supported_funcs)

  verify_files_or_stop(bedfiles, verify_tabix = TRUE)
  verify_regions_or_stop(regions)

  stopifnot("'mval' must be TRUE or FALSE" = mval %in% c(TRUE, FALSE))

  n_threads <- ifelse(
    is.null(nthreads),
    getOption("iscream.threads"),
    check_thread_count(nthreads)
  )
  Cpp_region_map(bedfiles, regions, fun, mval, bismark, nthreads = n_threads)
}
