#' Sparse Methylation matrix
#'
#' Make a sparse methylation matrix with given scWGBS sample bed files. The
#' matrix is sparse in that it stores the cpg locus id, the sample id and the
#' associated endcoded beta and coverage value. It is not a true sparse matrix
#' as you would find in the Matrix package and has no functions associated with
#' it.
#' @param cpg_bed_file The gzipped bed file with CpG loci. This is used to make the matrix index
#' @param sample_list A list of sample names for the bed files
#' @param sample_path Path to the directory containing sample bed files
#' @param to_mm Whether to write the matrix to a matrix market file
#' @param file_ext The file extension of the sample files
#'
#' @importFrom data.table := setcolorder rbindlist fwrite
#' @export
#'
#' @examples
#' sample_list <- c("sample1", "sample2", "sample3")
#' sparse_mat <- make_sparse_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
make_sparse_mat <- function(cpg_bed_file, sample_list, sample_path, to_mm = TRUE, file_ext = ".bed.gz") {

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

  # TODO: option if chunking is necessary

  if (to_mm) {
    message("Writing cpg_index")
    fwrite(cpg_index, "cpg_index.gz", sep = "\t", compress = "gzip")
    message("Writing sample_list")
    write(sample_list, "sample_list.txt")
    message("Writing methylation matrix")
    mm_writer(
      DT = result,
      row_count = nrow(cpg_index),
      col_count = length(sample_list),
      filename = "result.mtx.gz"
    )
    message("Done")
  } else {
    list(cpg_index = cpg_index, sample_index = sample_list, sparse_mat = result)
  }
}

make_cpg_index <- function(cpg_bed_file) {
  cpgs <- fread(cpg_bed_file, col.names = c("chr", "start"), drop = 3)
  cpgs[, "cpg_id" := seq.int(nrow(cpgs))][]
}

#' BSSeq Sparse Methylation matrix
#'
#' Make a sparse methylation matrix with given scWGBS sample bed files. The
#' matrix is sparse in that it stores the cpg locus id, the sample id and the
#' associated endcoded beta and coverage value. It is not a true sparse matrix
#' as you would find in the Matrix package and has no functions associated with
#' it.
#' @param cpg_bed_file The gzipped bed file with CpG loci. This is used to make the matrix index
#' @param sample_list A list of sample names for the bed files
#' @param sample_path Path to the directory containing sample bed files
#' @param file_ext The file extension of the sample files
#'
#' @importFrom data.table rbindlist
#' @export
#'
#' @examples
#' sample_list <- c("sample1", "sample2", "sample3")
#' sample_matrix <- make_sparse_mat(sample_list, "./pileup", "./data/cpgs.bed.gz", to_mm = FALSE)
#' make_sparse_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
make_sparse_mat_bs <- function(cpg_bed_file, sample_list, sample_path, file_ext = ".bed.gz") {

  cpg_index <- make_cpg_index(cpg_bed_file)

  meth_mat_list <- vector("list", length(sample_list))
  mapply(function(i) {
      meth_mat_list[[i]] <<- read_sample_bs(path(sample_path, paste0(sample_list[i], file_ext)), TRUE)[
        cpg_index, on = .(chr, start), .(cpg_id, cov, M), nomatch = 0L
      ][, "sample_id" := i]
    }, seq_along(sample_list))

  result <- rbindlist(meth_mat_list)
  setcolorder(result, c("cpg_id", "sample_id", "cov", "M"))
  list(cpg_index = cpg_index, sparse_mat = result)
}
