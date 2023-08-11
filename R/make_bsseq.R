#' Make a BSSeq object
#'
#' Create a BSSeq object using a encoded sparse methylation matrix. This method
#' uses more memory than `make_bsseq_lm` since it uses dense coverage and
#' methylation matrices when creating the BSSeq object.
#' @param encoded_matrix A sparse matrix of CpG loci and samples with encoded beta values
#' @importFrom data.table dcast setnafill
#'
#' @return A BSSeq object
#'
#' @export
#'
make_bsseq_hm <- function(sp, sample_list) {

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

#' Make a BSSeq object using less memory
#'
#' Create a BSSeq object using a encoded sparse methylation matrix. This method
#' uses less memory than `make_bsseq_hm` since it uses sparse coverage and
#' methylation matrices when creating the BSseq object. It is however slower in
#' the BSseq object creation.
#' @param encoded_matrix A sparse matrix of CpG loci and samples with encoded beta values
#' @importFrom data.table dcast setnafill
#' @importFrom Matrix readMM
#'
#' @return A BSSeq object
#'
#' @export
#'
make_bsseq_lm <- function(
    sparse_mat,
    sample_list,
    cov_filename = "cov.mtx.gz",
    m_filename = "M.mtx.gz"
) {
  row_count <- nrow(sparse_mat$cpg_index)
  col_count <- length(sample_list)

  message("Writing coverage matrix market file")
  mm_writer(sparse_mat$sparse_mat[, c(1, 2, 3)], row_count, col_count, cov_filename)

  message("Writing M matrix market file")
  mm_writer(sparse_mat$sparse_mat[, c(1, 2, 4)], row_count, col_count, m_filename)

  message("Making GRanges")
  gr <- GRanges(sparse_mat$cpg_index[, end := start + 2][, cpg_id := NULL])
  gc()

  cov_mat <- DelayedArray(readMM(cov_filename))
  M_mat <- DelayedArray(readMM(m_filename))

  message("Making BSseq object")
  BSseq(
    gr = gr,
    Cov = cov_mat,
    M = M_mat,
    sampleNames = sample_list
  )
}

