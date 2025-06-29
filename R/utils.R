#' Validate provided aligner
#'
#' Only "biscuit", "bismark", and "bsbolt" are currently supported
#' @param aligner The input alinger
#' @return true; quits if the input is not among supported_aligners
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
#' Convert GRanges object to a vector of strings
#'
#' This function was adapted from
#' [stuart-lab/signac](https://github.com/stuart-lab/signac/blob/0b33b1154f9a610897d1efad2c0065081c4e7132/R/utilities.R#L732C1-L756C2).
#'
#' @param gr A GRanges object
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @return A character vector
#'
#' @export
get_granges_string <- function(gr, sep = c(":", "-")) {
  if (requireNamespace("GenomicRanges", quietly = TRUE)) {
    region_str <- paste0(
      as.character(x = GenomicRanges::seqnames(x = gr)),
      sep[[1]],
      GenomicRanges::start(x = gr),
      sep[[2]],
      GenomicRanges::end(x = gr)
    )
    names(region_str) <- names(gr)
    return(region_str)
  } else {
    stop("The 'GenomicRanges' package must be installed for this functionality")
  }
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
#' @return A character vector
#'
#' @export
get_df_string <- function(regions_df, feature_col = NULL) {
  stopifnot("colnames must be 'chr', 'start' and 'end'" = colnames(regions_df)[1:3] == c("chr", "start", "end"))
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
