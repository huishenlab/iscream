#' Query lines from a tabixed BED file
#' @param bedfiles The bedfiles to be queried
#' @param regions A vector, data frame or GenomicRanges of genomic regions. See
#' details.
#' @param aligner The aligner used to produce the BED files - one of "biscuit",
#' "bismark", "bsbolt". Will set the result data.table's column names based on
#' this argument.
#' @param col.names A vector of column names for the data columns of the
#' result.table, not including "chr", "start", and "end". Set if your BED file
#' is not from the supported aligners or is a general BED file.
#' @param zero_based Whether the input BED file has a zero-based start column -
#' used when coverting the result data frame to GenomicRanges.
#' @param raw Set true to give a named list of raw strings from the regions in
#' the style of `Rsamtools::scanTabix` instead of a data.table
#' @param grlist Whether to return a GRangesList or a single GRanges object
#' when querying multiple files
#' @param nthreads Set number of threads to use overriding the
#' `"iscream.threads"` option. See `?set_threads` for more information.
#'
#' @details
#'
#' ## Query method
#' '*iscream* has two methods to query records from BED files:
#' - the *tabix* shell executable: fast since its output can be redirected to a
#' file (which `data.table::fread()` can then read) instead of having to
#' allocate memory and store it during the query
#'
#' - *iscream's* tabix implementation, based on the *tabix* executable using
#' *htslib*, but slower on large queries since it stores the records as they
#' are found instead of writing to a file. However it's able to store each
#' regions records independently instead of in a single file and is used in
#' `make_mat()`, `make_bsseq_mat()`, and `summarize_regions()`.
#'
#' When *iscream* is attached, it checks that the *tabix* executable is
#' available with `Sys.which()` and, if available, sets `options("tabix.method"
#' = "shell")`. `tabix()` then uses the *tabix* executable to make
#' queries, except if `raw = TRUE`. If *tabix* is not found, *iscream* uses its
#' tabix implementation. To use only *iscream's* tabix implementation, set
#' `options("tabix.method" = "htslib")`.
#'
#' ## Input region formats
#' The input regions may be string vector in the form "chr:start-end", a
#' dataframe with "chr", "start" and "end" columns or a GRanges object. Input
#' regions must be 1-based.  If the input is a GRanges, the output will also be
#' GRanges with any associated metadata columns (joined onto the result using
#' `GenomicRanges::findOverlaps()`). When making `GRanges`, the 0-based records
#' from BED-files will be converted to 1-based with
#' `GenomicRanges::makeGRangesFromDataFrame()`. Bismark's coverage files will
#' not be converted as they are already 1-based and the `ranges` slot will be
#' only one position.
#'
#' @importFrom data.table as.data.table tstrsplit set := rbindlist fread fwrite setnames
#' @importFrom parallel mclapply
#' @importFrom tools file_path_sans_ext
#' @importFrom stats setNames
#' @importFrom methods is
#' @return A data frame or `GRanges` if the input was `GRanges`.
#'
#' @export
#' @examples
#' bedfiles <- system.file("extdata", package = "iscream") |>
#'   list.files(pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
#' regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
#' tabix(bedfiles[1], regions, col.names = c("beta", "coverage"))
tabix <- function(
  bedfiles,
  regions,
  aligner = NULL,
  col.names = NULL,
  zero_based = TRUE,
  raw = FALSE,
  grlist = TRUE,
  nthreads = NULL
) {
  verify_files_or_stop(bedfiles)
  if (!is.null(aligner)) {
    verify_aligner_or_stop(aligner)
    verify_filetype(bedfiles, aligner)
  }

  if (raw) {
    input_regions <- get_string_input_regions(regions)
    return(run_scan_tabix(bedfiles, input_regions, nthreads))
  }

  # make the query
  if (getOption("tabix.method") == "htslib") {
    input_regions <- get_string_input_regions(regions)
    result <- tabix.htslib(bedfiles, input_regions, nthreads)
  } else {
    regions_df <- get_df_input_regions(regions)
    result <- tabix.shell(bedfiles, regions_df, nthreads)
  }

  # single-file empties returns null, multi file returns a non-null empty data.table
  if (is.null(result)) {
    return(NULL)
  } else if (nrow(result) == 0) {
    warning("No records found in any file - if this is unexpected check that your region format matches your bedfiles")
    return(NULL)
  }

  base_colnames <- c("chr", "start", "end")
  if (is.null(aligner)) {
    result_colnames <- get_bed_colnames(col.names, bedfiles, result, base_colnames)
  } else {
    result_colnames <- get_meth_colnames(aligner, bedfiles, result, base_colnames)
  }

  if (length(bedfiles) > 1) {
    result_colnames <- c(result_colnames, "file")
  }

  setnames(result, result_colnames)

  # get GRanges
  if (is(regions, "GRanges")) {
    result.gr <- GenomicRanges::makeGRangesFromDataFrame(
      result,
      starts.in.df.are.0based = zero_based,
      keep.extra.columns = TRUE
    )
    overlaps <- GenomicRanges::findOverlaps(result.gr, regions)

    if (dim(GenomicRanges::mcols(regions))[2] > 0) {
      mcols.colnames <- colnames(GenomicRanges::mcols(regions))
      mcols.subjectHits <- GenomicRanges::mcols(regions)[S4Vectors::subjectHits(overlaps), mcols.colnames]
      GenomicRanges::mcols(result.gr)[S4Vectors::queryHits(overlaps), mcols.colnames] <-
        GenomicRanges::mcols(regions)[S4Vectors::subjectHits(overlaps), mcols.colnames]
    }
    if ("file" %in% colnames(GenomicRanges::mcols(result.gr)) & grlist) {
      return(GenomicRanges::split(result.gr, as.factor(result.gr$file)))
    }
    return(result.gr)
  }
  return(result)
}

tabix.shell <- function(bedfiles, regions_df, nthreads) {
  if (length(bedfiles) == 1) {
    return(tabix.shell.single(bedfiles, regions_df))
  }

  mclapply(
    bedfiles,
    function(bedfile) {
      tbx_query <- tabix.shell.single(bedfile, regions_df)
      if (!is.null(tbx_query)) {
        tbx_query[, file := file_path_sans_ext(basename(bedfile), compression = TRUE)]
      }
    },
    mc.cores = .get_threads(nthreads)
  ) |>
    rbindlist()
}

tabix.shell.single <- function(bedfile, regions_df) {
  query.tmpfile <- tempfile(pattern = "regions", fileext = ".tsv")
  if (!is.null(regions_df)) {
    write_bed(regions_df, query.tmpfile)
  }

  cmd <- paste("tabix", bedfile, "-R", query.tmpfile)
  result <- suppressWarnings(fread(cmd = cmd))

  if (is_empty(result, bedfile)) {
    return(NULL)
  }
  result
}

tabix.htslib <- function(bedfiles, input_regions, nthreads) {
  if (length(bedfiles) == 1) {
    result <- tabix.htslib.single(
      bedfile = bedfiles,
      regions = input_regions
    )
  } else {
    dt_list <- mclapply(
      bedfiles,
      function(file) {
        tbx_query <- tabix.htslib.single(
          bedfile = file,
          regions = input_regions
        )
        if (!is.null(tbx_query)) {
          tbx_query[, file := file_path_sans_ext(basename(file), compression = TRUE)]
        }
        return(tbx_query)
      },
      mc.cores = .get_threads(nthreads)
    )
    result <- rbindlist(dt_list)
  }
}

tabix.htslib.single <- function(bedfile, regions) {
  lines <- Cpp_query_interval(bedfile, regions)
  lines_dt <- as.data.table(lines)
  if (is_empty(lines_dt, bedfile)) {
    return(NULL)
  }

  lines_dt[, tstrsplit(lines, "\t", fixed = TRUE, type.convert = TRUE)]
}

run_scan_tabix <- function(bedfiles, input_regions, nthreads) {
  if (length(bedfiles) == 1) {
    return(scan_tabix(bedfiles, input_regions))
  } else {
    bedline_list <- mclapply(bedfiles, scan_tabix, input_regions, mc.cores = .get_threads(nthreads)) |>
      setNames(
        nm = file_path_sans_ext(basename(bedfiles), compression = TRUE),
        object = _
      )
    return(bedline_list)
  }
}

# helpers

get_string_input_regions <- function(regions) {
  if (is(regions, "GRanges")) {
    get_granges_string(regions)
  } else if ("data.frame" %in% class(regions)) {
    get_df_string(regions)
  } else {
    regions
  }
}

get_df_input_regions <- function(regions) {
  if (is(regions, "GRanges")) {
    regions_df <- as.data.table(regions)[, 1:3]
    colnames(regions_df)[1] <- "chr"
    return(regions_df)
  } else if ("data.frame" %in% class(regions)) {
    regions
  } else {
    get_df_from_string(regions)
  }
}

get_meth_colnames <- function(aligner, bedfiles, result, base_colnames) {
  base_colnames <- c("chr", "start", "end")
  biscuit_colnames <- c("beta", "coverage")
  bismark_colnames <- c("methylation.percentage", "count.methylated", "count.unmethylated")

  if (aligner == "biscuit") {
    result_colnames <- c(base_colnames, biscuit_colnames)
    if (grepl("mergecg", bedfiles[1])) {
      result_colnames <- c(result_colnames, "mergecg")
    }
  } else {
    result_colnames <- c(base_colnames, bismark_colnames)
  }
  result_colnames
}

get_bed_colnames <- function(col.names, bedfiles, result, base_colnames) {
  data_cols <- ifelse(length(bedfiles) > 1, ncol(result) - 4, ncol(result) - 3)
  input_length <- length(col.names)

  if (input_length < data_cols) {
    data_colnames <- c(base_colnames, paste0("V", seq_len(data_cols)))
    if (input_length == 0) {
      return(data_colnames)
    }
    warning(
      "Did not use input data 'colnames' - only ",
      length(col.names),
      " names provided for ",
      data_cols,
      " data columns"
    )
    return(data_colnames)
  } else if (length(col.names) > data_cols) {
    warning("Fewer columns in data than provided colnames")
  }

  c(base_colnames, col.names[seq_len(data_cols)])
}


write_bed <- function(regions_df, outfile) {
  fwrite(
    regions_df[, 1:3],
    outfile,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE
  )
}

is_empty <- function(result, bedfile) {
  if (nrow(result) == 0) {
    warning(
      "No records found in ",
      bedfile,
      " - if this is unexpected check that your region format matches your bedfiles"
    )
    return(TRUE)
  }
  return(FALSE)
}
