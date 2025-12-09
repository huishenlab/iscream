#' GRanges to region strings
#'
#' Coerces GenomicRanges to `chr:start-end` strings with `as.character`. If any
#' regions have the same start and end, `as.character` returns `chr:start`
#' strings which are invalid for the htslib API. These are corrected to
#' `chr:start-start`.
#'
#' @param gr A GRanges object
#' @param feature_col The `mcols` column to use as the names of the output string vector
#' @returns A character vector
#'
#' @export
#' @examples
#' if (requireNamespace("GenomicRanges", quietly = TRUE)) {
#'   get_granges_string(GenomicRanges::GRanges(c("chr1:1-10", "chr2:15-20")))
#' }
get_granges_string <- function(gr, feature_col = NULL) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The 'GenomicRanges' package must be installed for this functionality")
  }
  region_str <- as.character(gr)
  singles <- which(GenomicRanges::width(gr) == 1)
  single_count <- length(singles)
  if (single_count > 0) {
    region_str[singles] <- gsub(":([0-9]+)", ":\\1-\\1", region_str[singles])
    message("Corrected ", single_count, " invalid 'chr:start' region strings to 'chr:start-start'")
  }

  if (!is.null(feature_col) && feature_col %in% colnames(GenomicRanges::mcols(gr))) {
    names(region_str) <- GenomicRanges::mcols(gr)[, feature_col]
  } else {
    names(region_str) <- names(gr)
  }
  region_str
}

#' DataFrame to region strings
#'
#' Convert DataFrame to a vector of strings. Set feature names in a "name" column
#'
#' @param regions_df A data frame with "chr", "start" and "end" columns
#' @param feature_col The data frame column to use as the names of the output string vector
#'
#' @importFrom data.table setDT
#'
#' @returns A character vector
#'
#' @export
#' @examples
#' (df <- data.frame(chr = c("chr1", "chr2"), start = c(1, 5), end = c(4, 10)))
#' get_df_string(df)
get_df_string <- function(regions_df, feature_col = NULL) {
  colnames.check <- colnames(regions_df)[seq_len(3)]
  stopifnot(
    "colnames must be 'chr', 'start' and 'end'" = colnames.check == c("chr", "start", "end")
  )
  chr <- start <- end <- NULL
  regions.dt <- setDT(regions_df)
  regions <- regions.dt[, paste0(chr, ":", start, "-", end)]
  if (!is.null(feature_col) && feature_col %in% colnames(regions_df)) {
    names(regions) <- regions_df[[feature_col]]
  }
  return(regions)
}

get_df_from_string <- function(regions) {
  start <- NULL
  as.data.table(regions)[, tstrsplit(regions, ":|-", fixed = FALSE, names = c("chr", "start", "end"))][,
    start := as.integer(start)
  ]
}

# Get GRanges from chr and pos vector
getGR <- function(chr, pos) {
  if (requireNamespace("GenomicRanges", quietly = TRUE)) {
    GenomicRanges::GRanges(chr, IRanges::IRanges(pos, pos, width = 1))
  }
}

#' Return region strings from GRanges or data.frame region inputs
#' @keywords internal
get_string_input_regions <- function(regions, feature_col = NULL) {
  if (is(regions, "GRanges")) {
    get_granges_string(regions, feature_col)
  } else if (is(regions, "data.frame")) {
    get_df_string(regions, feature_col)
  } else {
    regions
  }
}

#' Return data.frame from string region inputs to write to disk
#' @keywords internal
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
