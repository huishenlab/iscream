#' @rdname tabix
#' @export
tabix_raw <- function(
  bedfiles,
  regions,
  col.names = NULL,
  zero_based = TRUE,
  raw = FALSE,
  nthreads = NULL
) {
  input_regions <- get_string_input_regions(regions)
  run_scan_tabix(bedfiles, input_regions, .get_threads(nthreads))
}
