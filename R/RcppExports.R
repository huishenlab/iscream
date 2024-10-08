# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Make a Bsseq object
#' @param bedfiles A vector of bedfiles
#' @param regions A vector of regions
#' @param bismark Whether the input is a bismark coverage file
#' @param prealloc The number of rows to initialize the matrices with
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @keywords internal
#' @export
Cpp_query_all <- function(bedfiles, regions, bismark, merged, sparse, prealloc, nthreads) {
    .Call(`_iscream_Cpp_query_all`, bedfiles, regions, bismark, merged, sparse, prealloc, nthreads)
}

#' Get the number of available threads from OpenMP.
#'
#' This queries the number of available threads usign OpenMP, but will not
#' reliably provide an accurate available thread count. To get a more reliable
#' count that accounts for environment variables and HPC schedulers, use
#' get_threads()`. This function was pulled from
#' github.com/rdatatable/data.table
#' @param verbose Whether to be verbose on available omp threads
#'
#' @keywords internal
#' @export
#'
#' @examples
#' get_omp_threads(verbose = TRUE)
get_omp_threads <- function(verbose) {
    .Call(`_iscream_get_omp_threads`, verbose)
}

setup_logger <- function(logname = "iscream") {
    invisible(.Call(`_iscream_setup_logger`, logname))
}

#' spdlog Logging Lever Setter
#'
#' A helper function to turn a logging level given as string
#' into the current logging level
#'
#' @param name A string with the logging level. Value understood are,
#' in decreasing verbosity \sQuote{trace}, \sQuote{debug}, \sQuote{info},
#' \sQuote{warning}, \sQuote{error}, \sQuote{critical}, and \sQuote{off}.
#' Unrecognised names are equivalent to \sQuote{off}.
#' @return Nothing is returned.
#' @keywords internal
Cpp_set_log_level <- function(name) {
    invisible(.Call(`_iscream_Cpp_set_log_level`, name))
}

#' Get the current log level
#'
#' Can handle all of spdlogs levels, but iscream functions only supports
#' "info" and "debug"
#' @return The current logging level as a string
#' @export
get_log_level <- function() {
    .Call(`_iscream_get_log_level`)
}

#' Query a genomic interval from a opened htsFile and return the reads in it
#'
#' @param region Genomic region string in the form "chr:start-end"
#' @param bedFile The opened htslib bed file stream
#' @param tbx The bed file's tab-index
#' @returns A vector of strings from the matching region of the bed file
NULL

#' Get reads from multiple genomic regions from a tabixed bed file
#'
#' @param bedfile The name of the bed file - must have a corresponding tabix
#' file with the same name and .tbi extension
#' @param regions A vector of region strings in the form "chr:start-end"
NULL

#' Get reads from single genomic regions from multiple tabixed bed file.
#'
#' @param bedfiles A vector of bedfile names - must have corresponding tabix
#' files with the same name and .tbi extension
#' @param region A vector regions string in the form "chr:start-end"
NULL

#' Get reads from a single genomic region from one tabixed bed file.
#'
#' @param bedfile The name of the bed file - must have a corresponding tabix
#' file with the same name and .tbi extension
#' @param region The region string in the form "chr:start-end"
NULL

#' Query the chromosomes or seqnames from a vector of files
#' @param bedfile_vec The vector of bedfile paths
#' @return A vector of seqnames
#'
#' @keywords internal
#' @export
Cpp_query_chroms <- function(bedfile_vec, nthreads) {
    .Call(`_iscream_Cpp_query_chroms`, bedfile_vec, nthreads)
}

#' Get reads from a single genomic region from one tabixed bed file to return as CharacterVector
#'
#' @param bedfile The name of the bed file - must have a corresponding tabix
#' file with the same name and .tbi extension
#' @param regions A vector of region strings in the form "chr:start-end"
#'
#' @keywords internal
#' @export
Cpp_query_interval <- function(bedfile, regions) {
    .Call(`_iscream_Cpp_query_interval`, bedfile, regions)
}

#' Get namde list of reads from a single genomic region from one tabixed bed file
#'
#' @param bedfile The name of the bed file - must have a corresponding tabix
#' file with the same name and .tbi extension
#' @param regions A vector of region strings in the form "chr:start-end"
#'
#' @keywords internal
#' @export
scan_tabix <- function(bedfile, regions) {
    .Call(`_iscream_scan_tabix`, bedfile, regions)
}

#' Apply a function over CpGs within features
#'
#' This function should be called from `summarize_regions()` since there are few
#' sanity checks on the C++ side.
#' @param bedfiles A vector of bedfile paths
#' @param regions A vector of genomic regions
#' @param fun_vec Vector of the armadillo-supported stats functions to apply over the
#' CpGs in the ' regions: `"sum"`, `"mean"`, `"median"`, `"stddev"`,
#' `"variance"`, `"range"`.
#' @param mval Calculates M values when TRUE, use beta values when FALSE
#' @param bismark If the input is in the bismark column format instead of BISCUIT
#' @param region_rownames Whether to set rownames to the regions strings. Not
#' necessary if your regions vector is unnamed. If its names, then the "Feature"
#' column is set to the names and the rownames are set to the regions string
#' @param nthreads Number of cores to use. See details.
#'
#' @details
#' The optimal number of threads depends on the number of bedfiles, but is set
#' to half the available OpenMP cores. See `?get_threads` for more details. It
#' can be manaully set with `set_threads()`.
#'
#' @keywords internal
#' @export
Cpp_summarize_regions <- function(bedfiles, regions, fun_vec, mval, bismark, region_rownames = FALSE, nthreads = 1L) {
    .Call(`_iscream_Cpp_summarize_regions`, bedfiles, regions, fun_vec, mval, bismark, region_rownames, nthreads)
}

