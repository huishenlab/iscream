#' Beta matrix
#'
#' Make a cov matrix with methylation matrix. Wrapper for `decoder`.
#' @param encoded_matrix A matrix CpG loci and samples with encoded beta values
#' @return A data.table of coverage values for each sample
#'
#' @export
#'
#' @examples
#' sample_matrix <- make_meth_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
#' beta_mat <- make_beta_mat(sample_matrix)
make_beta_mat <- function(encoded_matrix) {
  decoder(encoded_matrix, 1)
}

#' Cov matrix
#'
#' Make a cov matrix with methylation matrix. Wrapper for `decoder`.
#' @param encoded_matrix A matrix CpG loci and samples with encoded beta values
#' @return A data.table of coverage values for each sample
#'
#' @importFrom data.table copy
#'
#' @export
#'
#' @examples
#' sample_matrix <- make_meth_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
#' beta_mat <- make_beta_mat(sample_matrix)
make_coverage_mat <- function(encoded_matrix) {
  decoder(encoded_matrix, 2)
}

#' M value matrix
#'
#' Make a beta matrix with methylation matrix. Wrapper for `decoder`.
#' @param encoded_matrix A matrix CpG loci and samples with encoded beta values
#' @return A data.table of M values for each sample
#'
#' @export
#'
#' @examples
#' sample_matrix <- make_meth_mat(sample_list, "./pileup", "./data/cpgs.bed.gz")
#' beta_mat <- make_beta_mat(sample_matrix)
make_m_mat <- function(encoded_matrix) {
  decoder(encoded_matrix, 3)
}

#' Get decoded values
#'
#' Get the decoded values from an encoded data.table methylation matrix
#' @param encoded_matrix A matrix CpG loci and samples with encoded beta values
#' @importFrom data.table copy
#' @param encoded_matrix The matrix with bit-packed beta and cov values
#' @param measure The required decoded value: 1 for beta, 2 for coverage and 3 for M
#' @param sparse Whether the matrix is a sparse matrix object
#' @return A data.table of the beta, coverage or M value by sample
#'
#' @export
#'
decoder <- function(encoded_matrix, measure, sparse = TRUE) {

  stopifnot(
    "Unknown measure: Enter a measure of 1 for beta, 2 for coverage, and 3 for M value" = measure %in% 1:3
  )

  encoded <- NULL # suppress 'no visible binding' error
  if (sparse) {
    measures <- c(
      "1" = "beta",
      "2" = "coverage",
      "3" = "M"
    )
    encoded_matrix[, c(1, 2, 3)][, paste(measures[measure]) := vdecoder(encoded, measure)][, 3][]
  } else {
    beta_mat <- copy(encoded_matrix)
    beta_mat[, lapply(.SD, vdecoder, measure), .SDcols = 3:ncol(encoded_matrix)][]
  }
}
