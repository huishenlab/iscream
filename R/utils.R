#' Validate provided aligner
#'
#' Only "biscuit", "bismark", and "bsbolt" are currently supported
#' @param aligner The input alinger
#' @returns true; quits if the input is not among supported_aligners
#'
#' @keywords internal
verify_aligner_or_stop <- function(aligner) {
  supported_aligners <- c("biscuit", "bismark", "bsbolt")
  if (!(aligner %in% supported_aligners)) {
    stop(paste(
      "aligner =",
      aligner,
      "not supported. Use one of the supported aligners:",
      paste(supported_aligners, collapse = ", ")
    ))
  }
}

#' GRanges to region strings
#'
#' @param gr A GRanges object
#' @returns A character vector
#'
#' @export
#' @examples
#' if (requireNamespace("GenomicRanges", quietly = TRUE)) {
#'   get_granges_string(GenomicRanges::GRanges(c("chr1:1-10", "chr2:15-20")))
#' }
get_granges_string <- function(gr) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The 'GenomicRanges' package must be installed for this functionality")
  }
  region_str <- as.character(gr)
  names(region_str) <- names(gr)
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
