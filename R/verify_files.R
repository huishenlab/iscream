#' Check that files exist
#'
#' @param files_vec A vector of file paths
#' @param error_file_prefix Error message prefix for 'Bedfile' vs 'Tabix file'
#' @importFrom fs file_exists
#' @return TRUE if all input bedfiles have an associated tabix index file.
#' FALSE if not
check_files_exist <- function(files_vec, error_file_prefix = "Bedfile") {
  valid_files <- file_exists(files_vec)
  missing_files <- files_vec[!valid_files]
  if (length(missing_files != 0)) {
    stop(paste0(error_file_prefix, ": ", missing_files, " could not be found\n"))
  }
}

#' Verify that bedfiles are tabixed
#'
#' @param bedfiles A vector of bedfile paths
#' @param verify_tabix Whether to verify the presence of tabix files
#' @importFrom fs file_exists
#' @return TRUE if all input bedfiles have an associated tabix index file.
#' FALSE if not
verify_files_or_stop <- function(bedfiles, verify_tabix = TRUE) {
  check_files_exist(bedfiles)
  if (verify_tabix) {
    tbi_files <- paste0(bedfiles, ".tbi")
    check_files_exist(tbi_files, "Tabix file")
  }
}

#' Verify that regions are valid
#'
#' @param regions A vector of genomic regions
#' @return TRUE if all input regions start with 'chr' FALSE if not
verify_regions_or_stop <- function(regions) {
  valid_regions <- sapply(regions, function(i) grepl("^chr", i))
  invalid_regions <- regions[!valid_regions]
  if (length(invalid_regions != 0)) {
    stop(paste0(invalid_regions, " are invalid regions\n"))
  }
}
