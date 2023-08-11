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

  C <- make_coverage_mat(sp$sparse_mat)
  M <- make_m_mat(sp$sparse_mat)

  C <- DelayedArray::DelayedArray(
    setnafill(
      dcast(cbind(
        sp$sparse_mat[, 1:2], M),
        cpg_id ~ sample_id,
        value.var = "encoded"
      )[],
      fill = 0L)
  )

  gr <- GenomicRanges::GRanges(
    sp$cpg_index[
    list(cpg_id = C[, 1]),
    on = "cpg_id", .(chr, start)
    ][, end := start + 2][]
  )
  C <- C[, -1]

  M <- DelayedArray::DelayedArray(
    setnafill(
      dcast(cbind(
        sp$sparse_mat[, 1:2], M),
        cpg_id ~ sample_id,
        value.var = "encoded"
      )[, cpg_id := NULL][],
      fill = 0L
    )
  )

  rownames(C) <- NULL
  colnames(C) <- NULL
  rownames(M) <- NULL
  colnames(M) <- NULL

  bsseq::BSseq(
    gr = gr,
    M = M,
    Cov = C,
    sampleNames = sample_list
  )
}

