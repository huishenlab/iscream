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
    stop(paste0(error_file_prefix, ": ", missing_files, " could not be found\n"))
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
    warn_stop(paste("'aligner' set to 'biscuit' but files found with '.cov', extension", warning_msg))
  }
  if (aligner != "biscuit" & any(!sf_grepl(bedfiles, pattern = ".cov"))) {
    warn_stop(paste("'aligner' set to", aligner, "but no files found with '.cov', extension", warning_msg))
  }
}

validate_chrom_field <- function(chrom_no) {
  alpha_chroms <- c("X", "Y", "M")
  if (!grepl("^chr", chrom_no) & !chrom_no %in% alpha_chroms) {
    stopifnot("seqname is not numeric" = !is.na(as.numeric(chrom_no)))
  } else {
    chrom_no <- gsub("chr", "", chrom_no)
    if (! chrom_no %in% c("X", "Y", "M")) {
      stopifnot("chr is not in the format 'chr[1-22/X/Y]'" = !is.na(as.numeric(chrom_no)))
    }
  }
}

#' Validate a region string
#'
#' CAUTION: very basic validation is done here, checking only that there are
#' one or three elements and that the last two are integers and the first is
#' a valid integer or among "X", "Y", "M", optionally prefixed with 'chr'
#'
#' @param region A vector of genomic regions
#' @return TRUE if region is valid
#'
#' @importFrom stringfish sf_split
#'
#' @keywords internal
validate_region <- function(region) {
  split_region <- unlist(sf_split(region, ":|-", fixed = FALSE))
  n_elems <- length(split_region)
  chr <- split_region[1]
  validate_chrom_field(chr)
  if (n_elems == 3) {
    start <- as.numeric(split_region[2])
    end <- as.numeric(split_region[3])
    stopifnot("start is not numeric" = !is.na(start))
    stopifnot("end is not numeric" = !is.na(end))
    stopifnot("end is smaller than start" = end >= start)
  }
}

#' Verify that regions are valid
#'
#' @param regions A vector of genomic regions
#' @return TRUE if all input regions start with 'chr' FALSE if not
#'
#' @importFrom stringfish sf_split
#'
#' @keywords internal
verify_regions_or_stop <- function(regions) {
  sapply(regions, validate_region)
}
