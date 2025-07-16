#' @rdname tabix
#' @param zero_based Whether the input BED file has a zero-based start column -
#' used when coverting the result data frame to GenomicRanges.
#' @export
tabix_gr <- function(
  bedfiles,
  regions,
  aligner = NULL,
  col.names = NULL,
  zero_based = TRUE,
  nthreads = NULL
) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The 'GenomicRanges' package must be installed for this functionality")
  }
  result <- tabix(
    bedfiles,
    regions,
    aligner,
    col.names,
    nthreads
  )
  if (is.null(result)) {
    stop("No records found")
  }
  result.gr <- GenomicRanges::makeGRangesFromDataFrame(
    result,
    starts.in.df.are.0based = zero_based,
    keep.extra.columns = TRUE
  )
  if (is(regions, "character")) {
    regions <- GenomicRanges::GRanges(regions)
  } else if ("data.frame" %in% class(regions)) {
    regions <- GenomicRanges::makeGRangesFromDataFrame(regions)
  }
  overlaps <- GenomicRanges::findOverlaps(result.gr, regions)

  if (dim(GenomicRanges::mcols(regions))[2] > 0) {
    mcols.colnames <- colnames(GenomicRanges::mcols(regions))
    mcols.subjectHits <- GenomicRanges::mcols(regions)[S4Vectors::subjectHits(overlaps), mcols.colnames]
    GenomicRanges::mcols(result.gr)[S4Vectors::queryHits(overlaps), mcols.colnames] <-
      GenomicRanges::mcols(regions)[S4Vectors::subjectHits(overlaps), mcols.colnames]
  }
  if ("file" %in% colnames(GenomicRanges::mcols(result.gr))) {
    return(GenomicRanges::split(result.gr, as.factor(result.gr$file)))
  }
  result.gr
}
