#' Check that files exist
#'
#' @param files_vec A vector of file paths
#' @param error_file_prefix Error message prefix for 'Bedfile' vs 'Tabix file'
#' @return TRUE if all input bedfiles have an associated tabix index file.
#' FALSE if not
#'
#' @keywords internal
check_files_exist <- function(files_vec, error_file_prefix = "Bedfile") {
  valid_files <- file.exists(files_vec)
  missing_files <- files_vec[!valid_files]
  if (length(missing_files != 0)) {
    stop(error_file_prefix, ": ", missing_files, " could not be found\n")
  }
}

#' Verify that bedfiles are tabixed
#'
#' @param bedfiles A vector of bedfile paths
#' @param verify_tabix Whether to verify the presence of tabix files
#' @return TRUE if all input bedfiles have an associated tabix index file.
#' FALSE if not
#'
#' @keywords internal
verify_files_or_stop <- function(bedfiles, verify_tabix = TRUE) {
  check_files_exist(bedfiles)
  if (verify_tabix) {
    tbi_files <- paste0(bedfiles, ".tbi")
    check_files_exist(tbi_files, "Tabix file")
  }
}

#' Verify that the input bedfiles are of the type specified by the input aligner
#'
#' @param bedfiles A vector of bedfile paths
#' @param aligner The aligner chosen
#' @param stop_on_error Whether to warn or stop on aligner-filename mismatch
#'
#' @importFrom stringfish sf_grepl
#'
#' @return TRUE if all input bedfiles have an associated tabix index file.
#' FALSE if not
#'
#' @keywords internal
verify_filetype <- function(bedfiles, aligner, stop_on_error = FALSE) {
  warn_stop <- ifelse(stop_on_error, stop, warning)
  check_name_warning <- "- verify aligner"
  warning_msg <- ifelse(stop_on_error, check_name_warning, paste(check_name_warning, "and data frame colnames"))

  if (aligner == "biscuit" & any(sf_grepl(bedfiles, pattern = ".cov"))) {
    warn_stop("'aligner' set to 'biscuit' but files found with '.cov', extension ", warning_msg)
  }
  if (aligner != "biscuit" & any(!sf_grepl(bedfiles, pattern = ".cov"))) {
    warn_stop("'aligner' set to ", aligner, " but no files found with '.cov', extension ", warning_msg)
  }
}
