#' Query lines from a tabixed bedfile
#' @param bedfiles The bedfiles to be queried
#' @param regions A vector, data frame or GenomicRanges of genomic regions. See
#' details.
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt". Will set the result data.table's column names based on
#' this argument.
#' @param colnames A vector of column names for the result data.table. Set if
#' your bedfile is not from the supported aligners or is a general bedfile.
#' @param raw Set true to give a named list of raw strings from the regions in
#' the style of `Rsamtools::scanTabix` instead of a data.table
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @details
#' The input regions may be string vector in the form "chr:start-end", a
#' dataframe with "chr", "start" and "end" columns or a GRanges object. If the
#' input is a GRanges, the output will also be GRanges with any associated
#' metadata columns (joined onto the result using
#' `GenomicRanges::findOverlaps()`)
#'
#' @importFrom data.table as.data.table tstrsplit set := rbindlist
#' @importFrom parallel mclapply
#' @importFrom tools file_path_sans_ext
#' @importFrom stats setNames
#' @return A data.table
#'
#' @export
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' tabix(bedfiles[1], regions, colnames = c("chr", "start", "end", "beta", "coverage"))
tabix <- function(bedfiles, regions, aligner = "biscuit", colnames = NULL, raw = FALSE, nthreads = NULL) {
  verify_files_or_stop(bedfiles)

  if (class(regions)[1] == "GRanges") {
    input_regions <- get_granges_string(regions)
  } else if ("data.frame" %in% class(regions)) {
    input_regions <- get_df_string(regions)
  } else {
    input_regions <- regions
  }

  verify_aligner_or_stop(aligner)
  verify_filetype(bedfiles, aligner)

  if (raw) {
    if (length(bedfiles) == 1) {
      return (scan_tabix(bedfiles, input_regions))
    } else {
      bedline_list <- mclapply(bedfiles, function(file) {
        scan_tabix(file, input_regions)
      }, mc.cores = .get_threads(nthreads)) |>
      setNames(
        nm = file_path_sans_ext(basename(bedfiles), compression = TRUE),
        object = _
      )
      return(bedline_list)
    }
  }

  base_colnames <- c("chr", "start", "end")
  biscuit_colnames <- c("beta", "coverage")
  bismark_colnames <- c("methylation.percentage", "count.methylated", "count.unmethylated")

  mergecg <- FALSE
  if (!is.null(colnames)) {
    result_colnames <- colnames
  } else if (aligner == "biscuit") {
    result_colnames <- c(base_colnames, biscuit_colnames)
    if (grepl("mergecg", bedfiles[1])) {
      result_colnames <- c(result_colnames, "mergecg")
      mergecg <- TRUE
    }
  } else {
    result_colnames <- c(base_colnames, bismark_colnames)
  }

  if (length(bedfiles) == 1) {
      result <- single_tabix(
        bedfile = bedfiles,
        regions = input_regions,
        result_colnames = result_colnames,
        mergecg = mergecg
      )
  } else {
    dt_list <- mclapply(bedfiles, function(file) {
      tbx_query <- single_tabix(
        bedfile = file,
        regions = input_regions,
        result_colnames = result_colnames,
        mergecg = mergecg
      )
      if (!is.null(tbx_query)) {
        tbx_query[, sample := file_path_sans_ext(basename(file), compression = TRUE)]
      }
      return(tbx_query)
    }, mc.cores = .get_threads(nthreads))
    result <- rbindlist(dt_list)
  }

  if (class(regions)[1] == "GRanges") {
    result.gr <- GenomicRanges::GRanges(result)
    overlaps <- GenomicRanges::findOverlaps(result.gr, regions)

    if (dim(GenomicRanges::mcols(regions))[2] > 0) {
      mcols.colnames <- colnames(GenomicRanges::mcols(regions))
      mcols.subjectHits <- GenomicRanges::mcols(regions)[S4Vectors::subjectHits(overlaps), mcols.colnames]
      GenomicRanges::mcols(result.gr)[S4Vectors::queryHits(overlaps), mcols.colnames] <-
        GenomicRanges::mcols(regions)[S4Vectors::subjectHits(overlaps), mcols.colnames]
    }

    return(result.gr)
  }

  return(result)
}

single_tabix <- function(bedfile, regions, result_colnames, mergecg) {
  lines <- Cpp_query_interval(bedfile, regions)
  if (length(lines) == 0) {
    warning(paste("No records found in", bedfile, "- if this is unexpected check that your region format matches your bedfiles"))
    return(NULL)
  }
  lines_dt <- as.data.table(lines)

  lines_dt <- lines_dt[, tstrsplit(lines, "\t", fixed = TRUE)]
  n_col <- ncol(lines_dt)
  if (length(result_colnames) < n_col) {
    warning(paste(
        "Did not use input 'colnames' - only",
        length(result_colnames), "names provided for", n_col, "column data.table"
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
