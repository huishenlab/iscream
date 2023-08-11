#' Sample reader
#'
#' Reads scWGBS data from a sample bed file and generates a data.table
#' @param sample_name The name of the sample in the bed file
#' @param merged Whether the bed file has CGs merged or unmerged
#' @return The sample as a data.table
#'
#' @importFrom data.table fread
#' @importFrom data.table :=
#' @importFrom data.table .SD
#' @importFrom fs path
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
  vencoder_cols <- c("beta", "cov")
  sample_data[, encoded := do.call(vencoder, setNames(.SD, names(vencoder_cols))),
  .SDcols = vencoder_cols][]
}

#' Cpg and sample joiner
#'
#' Takes a sample data.table and joins it to the CpG index table on the
#' chromosome and start columns
#' @param cpg_table The table with CpG loci and any existing samples
#' @param sample_dt The sample data table from read_sample
#' @return The cpg_table with the sample_table tacked on
#'
#' @importFrom data.table :=
#'
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
    .joiner(cpgs, sample_data, sample_list[i])
  }, seq_along(sample_list)
  )
  cpgs
}

#' Beta-value matrix
#'
#' Make a beta matrix with methylation matrix
#' @param encoded_matrix A matrix CpG loci and samples with encoded beta values
#' @importFrom fs path
#' @export
#'
#' @examples
#' sample_matrix <- make_meth_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
#' beta_mat <- make_beta_mat(sample_matrix)
make_beta_mat <- function(encoded_matrix) {
  encoded_matrix[complete.cases(encoded_matrix)][
  , lapply(.SD, decode_beta), .SDcols = 4:ncol(encoded_matrix)][]
}

make_cpg_index <- function(cpg_bed_file) {
  cpgs <- fread(cpg_bed_file, col.names = c("chr", "start"), drop = 3)
  cpgs[, "cpg_id" := seq.int(nrow(cpgs))][]
}

#' Sparse Methylation matrix
#'
#' Make a methylation matrix with given scWGBS sample bed files
#' @param cpg_bed_file The gzipped bed file with CpG loci. This is used to make the matrix index
#' @param sample_list A list of sample names for the bed files
#' @param sample_path Path to the directory containing sample bed files
#'
#' @importFrom data.table := setcolorder rbindlist
#' @export
#'
#' @examples
#' sample_list <- c("sample1", "sample2", "sample3")
#' sample_matrix <- make_meth_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
make_sparse_mat <- function(cpg_bed_file, sample_list, sample_path, file_ext = ".bed.gz") {

  cpg_index <- make_cpg_index(cpg_bed_file)

  meth_mat_list <- vector("list", length(sample_list))
  mapply(function(i) {
      meth_mat_list[[i]] <<- read_sample(path(sample_path, paste0(sample_list[i], file_ext)), TRUE)[
        cpg_index, on = .(chr, start), .(cpg_id, encoded), nomatch = 0L
      ][, "sample_id" := i]
    }, seq_along(sample_list))

  result <- rbindlist(meth_mat_list)
  setcolorder(result, c("cpg_id", "sample_id", "encoded"))
  list(cpg_index = cpg_index, sample_index = sample_list, sparse_mat = result)
}
