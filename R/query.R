#' Query the chromosomes or seqnames from a vector of files
#' @param bedfiles The vector of bedfile paths
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#' @return A vector of seqnames
#'
#' @export
query_chroms <- function(bedfiles, nthreads = NULL) {
  verify_files_or_stop(bedfiles)

  n_threads <- .get_threads(nthreads)

  Cpp_query_chroms(bedfiles, n_threads)
}

#' Query lines from a tabixed bedfile
#' @param bedfile The bedfile to be queried
#' @param regions A vector of genomic region strings
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt". Will set the result data.table's column names based on
#' this argument.
#' @param colnames A vector of column names for the result data.table. Set if
#' you your bedfile is not from the supported aligners or is a general bedfile.
#' @param raw Set true to give a named list of raw strings from the regions in
#' the style of `Rsamtools::scanTabix` instead of a data.table
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @importFrom data.table as.data.table tstrsplit set
#' @return A data.table
#'
#' @export
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' tabix(bedfiles[1], regions, colnames = c("chr", "start", "end", "beta", "coverage"))
tabix <- function(bedfile, regions, aligner = "biscuit", colnames = NULL, raw = FALSE, nthreads = NULL) {
  verify_files_or_stop(bedfile)
  verify_regions_or_stop(regions)
  verify_aligner_or_stop(aligner)
  verify_filetype(bedfile, aligner)
  base_colnames <- c("chr", "start", "end")
  biscuit_colnames <- c("beta", "coverage")

  if (length(bedfile) > 1) {
    stop("Cannot tabix multiple files - only single-file queries are currently supported")
  }

  mergecg <- FALSE
  if (!is.null(colnames)) {
    result_colnames <- colnames
  } else if (aligner == "biscuit") {
    result_colnames <- c(base_colnames, biscuit_colnames)
    if (grepl("mergecg", bedfile)) {
      result_colnames <- c(result_colnames, "mergecg")
      mergecg <- TRUE
    }
  } else {
    result_colnames <- c(base_colnames, bismark_colnames)
  }



  n_threads <- .get_threads(nthreads)
  if (raw) return(scan_tabix(bedfile, regions))

  lines <- Cpp_query_interval(bedfile, regions)
  if (length(lines) == 0) {
    warning("No records found")
    return(NULL)
  }
  lines_dt <- as.data.table(lines)

  lines_dt <- lines_dt[, tstrsplit(lines, "\t", fixed = TRUE)]
  n_col <- ncol(lines_dt)
  if (length(result_colnames) < n_col) {
    warning(paste(
        "Did not use input 'colnames' - only",
        length(colnames), "names provided for", n_col, "data.table"
      ))
    return(lines_dt)
  } else if (length(result_colnames) > n_col) {
    warning("Fewer columns in data than provided colnames")
  }

  colnames(lines_dt) <- result_colnames[1:n_col]
  end_col <- ifelse(mergecg, ncol(lines_dt) - 1, ncol(lines_dt))
  for (i in 2:end_col) set(lines_dt, j = i, value = as.numeric(lines_dt[[i]]))

  return(lines_dt)
}
