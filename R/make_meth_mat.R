#' Methylation matrix
#'
#' Make a methylation matrix with given scWGBS sample bed files
#' @param sample_list A list of sample names for the bed files
#' @param sample_path Path to the directory containing sample bed files
#' @param cpg_bed_file The gzipped bed file with CpG loci. This is used to make the matrix index
#' @param file_ext The extension of the cpg index and sample files
#' @param merged Whether the bed files have CPs merged or unmerged
#'
#' @importFrom data.table fread
#' @importFrom data.table .SD
#' @export
#'
#' @examples
#' sample_list <- c("sample1", "sample2", "sample3")
#' sample_matrix <- make_meth_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
make_meth_mat <- function(
  sample_list,
  sample_path,
  cpg_bed_file,
  file_ext = ".bed.gz",
  merged = FALSE
) {
  cpgs <- fread(cpg_bed_file, col.names = c("chr", "start"), drop = 3)
  mapply(function(i) {
    sample_data <- read_sample(
      path(sample_path, paste0(sample_list[i], file_ext)),
      merged
    )
    joiner(cpgs, sample_data, sample_list[i])
  }, seq_along(sample_list)
  )
  cpgs
}
