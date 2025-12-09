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

# Check if package is loaded
# https://github.com/HenrikBengtsson/R.utils/blob/74def095eaa244e355d05fdf790ee6393dad1d99/R/isPackageLoaded.R#L33-L43
is_package_loaded <- function(package, caller, fail) {
  loaded_packages <- gsub("package:", "", search())
  if (package %in% loaded_packages) {
    return(TRUE)
  }
  if (fail) {
    msg <- sprintf("Please load the '%s' package to use '%s' output", package, caller)
    stop(msg)
  }
  FALSE
}
