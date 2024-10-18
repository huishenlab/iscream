#' Set logging level
#'
#' @param level The logging verbosity level to use
#' - `"info"`: the default that gives provides basic information about the
#' number of files and regions used in a function
#' - `"debug"`: more verbose about row allocations, how many CpGs were found in
#' a region, filename parsing etc. This mode cannot be used on more than one
#' thread as R cannot output messages from multiple threads without crashing.
#' @return None; sets the log level to the provided level
#'
#' @export
#'
#' @examples
#' set_log_level("info")
set_log_level <- function(level = "info") {
  n_threads <- get_threads()["use_threads"]
  validate_log_level(level, n_threads)
  Cpp_set_log_level(level)
}

#' Validate provided logging level
#'
#' Only "info" and "debug" are currently supported, with "debug" only supported
#' when using 1 thread
#' @param level The logging level to validate
#' @param n_threads The number of threads that the next iscream function call
#' will use
#' @return None; sets the log level to the provide level
#'
#' @keywords internal
#'
#' @examples
#' set_log_level("info")
validate_log_level <- function(level = get_log_level(), n_threads) {
  quiet_logging <- c("error", "info", "off")
  supported <- c(quiet_logging, "trace", "debug", "warn")
  unsupported <- c("critical")

  if (!(level %in% supported)) {
    if (level %in% unsupported) {
      stop(paste(level, "logging is currently not supported"))
    } else {
      stop(paste(level, "is not a supported spdlog log level"))
    }
  }

  if (n_threads > 1 && !(level %in% quiet_logging)) {
    stop(paste(
      n_threads, "threads were requested",
      "but", level, "logging cannot use more than 1 thread.",
      "Run `set_log_level(\"info\")` to use multiple threads.",
      "See `?set_log_level` for more details."
    ))
  }
}
