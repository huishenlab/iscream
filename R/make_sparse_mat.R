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

make_cpg_index <- function(cpg_bed_file) {
  cpgs <- fread(cpg_bed_file, col.names = c("chr", "start"), drop = 3)
  cpgs[, "cpg_id" := seq.int(nrow(cpgs))][]
}

