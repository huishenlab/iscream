#' Sample reader
#'
#' Reads scWGBS data from a sample bed file and generates a data.table
#' @param sample_name The name of the sample in the bed file
#' @param merged Whether the bed file has CGs merged or unmerged
#' @return The sample as a data.table
#'
#' @importFrom data.table fread := .SD
#' @importFrom fs path
#'
#' @keywords experimental
#' @export
#'
#' @examples
# TODO: write working examples
read_sample <- function(sample_name, merged) {
  sample_col_names <- c("chr", "start", "end", "beta", "cov")
  ifelse(merged,
    sample_data <- fread(sample_name, col.names = sample_col_names, drop = 6),
    sample_data <- fread(sample_name, col.names = sample_col_names)
  )

  # suppress 'no visible binding' error
  encoded <- NULL
  beta <- NULL
  cov <- NULL

  message("Found ", nrow(sample_data), " methylation loci\n")
  vencoder_cols <- c("beta", "cov")
  sample_data[, encoded := vencoder(beta, cov)]
}

#' Cpg and sample joiner
#'
#' Takes a sample data.table and joins it to the CpG index table on the
#' chromosome and start columns
#' @param cpg_table The table with CpG loci and any existing samples
#' @param sample_dt The sample data table from read_sample
#' @param sample_name name of the sample getting joined
#' @return The cpg_table with the sample_table tacked on
#'
#' @keywords experimental
#' @importFrom data.table :=
#'
joiner <- function(cpg_table, sample_dt, sample_name) {

  # suppress 'no visible binding' error
  chr <- NULL
  start <- NULL
  i.encoded <- NULL
  . <- NULL

  cpg_table[sample_dt, on = .(chr, start), paste0(sample_name) := i.encoded][]
}

#' BSSeq Sample reader
#'
#' Reads scWGBS data from a sample bed file and generates a data.table
#' @param sample_name The name of the sample in the bed file
#' @param merged Whether the bed file has CGs merged or unmerged
#' @return The sample as a data.table
#'
#' @importFrom data.table fread
#' @importFrom data.table :=
#' @importFrom data.table .SD
#' @importFrom data.table as.data.table
#'
#' @keywords experimental
#' @export
#'
#' @examples
# TODO: write working examples
read_sample_bs <- function(sample_name, merged) {
  sample_col_names <- c("chr", "start", "end", "beta", "cov")
  ifelse(merged,
    sample_data <- fread(sample_name, col.names = sample_col_names, drop = 6),
    sample_data <- fread(sample_name, col.names = sample_col_names)
  )
  message("Found ", nrow(sample_data), " methylation loci\n")

  # suppress 'no visible binding' error
  beta <- NULL
  cov <- NULL

  vencoder_cols <- c("beta", "cov")
  sample_data[, "M" := round(beta * cov)][, beta := NULL]
}

