check_files_exist <- function(files_vec, error_file_prefix = "File") {
  valid_files <- file_exists(files_vec)
  missing_files <- files_vec[!valid_files]
  if (length(missing_files != 0)) {
    stop(paste0(error_file_prefix, ": ", missing_files, " could not be found"))
  }
}

#' Verify that bedfiles are tabixed
#'
#' @param befiles A vector of bedfile paths
#' @importFrom fs file_exists
#' @return TRUE if all input bedfiles have an associated tabix index file.
#' FALSE if not
#'
#' @export
verify_files_or_stop <- function(bedfiles, verify_tabix = T) {
  check_files_exist(bedfiles)
  if (verify_tabix) {
    tbi_files <- paste0(bedfiles, ".tbi")
    check_files_exist(tbi_files, "Tabix file")
  }
}
