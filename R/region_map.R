#' Run a defined function over genomic regions
#'
#' Run a function on the CpGs in bedfiles across genomic regions. Currently
#' supported functions are aggregate, average, and ...
#' @param befiles A vector of bedfile paths
#' @param regions A vector of genomic regions strings
#' @param fun Function to apply over the region. See details.
#' @param mval Whether to calculate the M value (coverage \eqn{\times \beta})
#' or use the beta value
#' when applying the function.
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
region_map <- function(bedfiles, regions, fun = "aggregate", mval = TRUE) {

  supported_funcs <- c("aggregate", "average")
  stopifnot("Selected function not supported" = fun %in% supported_funcs)

  valid_files <- sapply(bedfiles, fs::file_exists)
  missing_files <- bedfiles[!valid_files]
  if (length(missing_files != 0)) {
    stop(paste0("File: ", missing_files, " could not be found"))
  }

  valid_regions <- sapply(regions, function(i) grepl("^chr", i))
  invalid_regions <- regions[!valid_regions]
  if (length(invalid_regions != 0)) {
    stop(paste0(invalid_regions, " are invalid regions\n"))
  }

  stopifnot("'mval' must be TRUE or FALSE" = mval %in% c(T, F))

  n_threads <- getOption("iscream.threads")
  Cpp_region_map(bedfiles, regions, fun, mval, nthreads = n_threads)
}
