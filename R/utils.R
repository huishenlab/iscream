#' Validate provided aligner
#'
#' Only "biscuit", "bismark", and "bsbolt" are currently supported
#' @param aligner The input alinger
#' @return true; quits if the input is not among supported_aligners
#'
#' @keywords internal
verify_aligner_or_stop <- function(aligner) {
  supported_aligners <- c("biscuit", "bismark", "bsbolt")
  if (!aligner %in% supported_aligners) {
    stop(paste("aligner =", aligner, "not supported. Use one of the supported aligners:", paste(supported_aligners, collapse = ", ")))
  }
}

