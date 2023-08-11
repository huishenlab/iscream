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
