#' Sample reader
#'
#' Reads scWGBS data from a sample bed file and generates a data.table
#' @param sample_name The name of the sample in the bed file
#' @param merged Whether the bed file has CGs merged or unmerged
#' @return The sample as a data.table
#'
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' sample1 <- read_sample("sample1", TRUE)
read_sample <- function(sample_name, merged) {
  sample_col_names <- c("chr", "start", "end", "beta", "cov")
  ifelse(merged,
    sample_data <- fread(sample_name, col.names = sample_col_names, drop = 6),
    sample_data <- fread(sample_name, col.names = sample_col_names)
  )
  message("Found ", nrow(sample_data), " methylation loci\n")
  sample_data[, encoded := mapply(encoder, beta, cov)][]
}

.joiner <- function(cpg_table, sample_dt, sample_name) {
  cpg_table[sample_dt, on = .(chr, start), paste0(sample_name) := i.encoded][]
}

#' Methylation matrix
#'
#' Make a methylation matrix with given scWGBS sample bed files
#' @param sample_list A list of sample names for the bed files
#' @param sample_path Path to the directory containing sample bed files
#' @param cpg_bed_file The gzipped bed file with CpG loci. This is used to make the matrix index
#' @param merged Whether the bed files have CPs merged or unmerged
#'
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' sample_list <- c("sample1", "sample2", "sample3")
#' sample_matrix <- make_meth_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
make_meth_mat <- function(
  sample_list,
  sample_path,
  cpg_bed_file,
  merged = FALSE
) {
  cpgs <- fread(cpg_bed_file, col.names = c("chr", "start", "end"))
  mapply(function(i) {
    sample_data <- read_sample(
      paste0(sample_path, sample_list[i], ".bed.gz"),
      merged
    )
    .joiner(cpgs, sample_data, sample_list[i])
  }, seq_along(sample_list)
  )
  cpgs
}
