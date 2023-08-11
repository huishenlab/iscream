#' Make a BSSeq object
#'
#' Create a BSSeq object using a encoded sparse methylation matrix
#' @param encoded_matrix A sparse matrix of CpG loci and samples with encoded beta values
#' @importFrom data.table dcast setnafill
#'
#' @return A BSSeq object
#'
#' @export
#'
make_bsseq_from_sparse <- function(sp, sample_list) {

  message("Making Cov matrix")
  C <- DelayedArray::DelayedArray(
    setnafill(dcast(
      sp$sparse_mat[, c(1, 2, 3)],
      cpg_id ~ sample_id,
      value.var = "cov"), fill = 0L)
  )
  gc()

  message("Making GRanges")
  gr <- GRanges(sp$cpg_index[C[, 1]][, end := start + 2][, c(1, 2, 4)]); gc()
  C <- C[, -1]
  gc()

  message("Making M matrix")
  M <- DelayedArray::DelayedArray(
    setnafill(dcast(
      sp$sparse_mat[, c(1, 2, 4)],
      cpg_id ~ sample_id, value.var = "M"),
      fill = 0L)
  )[, -1]
  gc()

  rownames(C) <- NULL
  colnames(C) <- NULL
  rownames(M) <- NULL
  colnames(M) <- NULL

  message("Making BSseq object")
  bsseq::BSseq(
    gr = gr,
    M = M,
    Cov = C,
    sampleNames = sample_list
  )
}

