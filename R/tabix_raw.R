#' @rdname tabix
#' @export
tabix_raw <- function(
  bedfiles,
  regions,
  nthreads = NULL,
  BPPARAM = NULL
) {
  input_regions <- get_string_input_regions(regions)
  run_scan_tabix(bedfiles, input_regions, nthreads, BPPARAM)
}

run_scan_tabix <- function(bedfiles, input_regions, nthreads, BPPARAM) {
  nthreads <- .get_threads(nthreads)
  bpparam <- BPPARAM %||% MulticoreParam(workers = nthreads, progressbar = TRUE)
  nworkers <- bpworkers(bpparam)
  check_thread_count(nworkers)

  if (length(bedfiles) == 1) {
    return(scan_tabix(bedfiles, input_regions))
  } else {
    bedline_list <- bplapply(bedfiles, scan_tabix, input_regions, BPPARAM = bpparam) |>
      setNames(
        nm = file_path_sans_ext(basename(bedfiles), compression = TRUE),
        object = _
      )
    return(bedline_list)
  }
}
