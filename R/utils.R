#' Matrix Market writer
#'
#' Write a sparse matrix data.table to disk in the Matrix Market format.
#'
#' @details
#' Each row of the data.table has to be of the form:
#'      row_num  col_num value
#' This function requires `row_count` and `col_count` because the number of
#' entries can be calculated from the input data.table.
#' @param DT A data.table with the entries of the matrix to write
#' @param row_count The number of rows in the matrix - this is not the same as `nrow(DT)`
#' @param col_count The number of columns in the matrix - this is not the same as `ncol(DT)`
#' @param filename The name of the file to write to
#' @param compress: Whether to write a gzip compressed file or not
#' @importFrom data.table fwrite
#'
#' @export
#'
mm_writer <- function(
    DT,
    row_count,
    col_count,
    filename,
    compress = TRUE) {

  entry_count <- nrow(DT)

  stopifnot("invalid MatrixMarket format: DT contains more than 3 columns" = ncol(DT) == 3)

  fwrite(list("%%MatrixMarket matrix coordinate integer general"), filename, compress = "gzip")
  fwrite(
    list(row_count, col_count, entry_count),
    filename,
    sep = " ", compress = "gzip", append = TRUE
  )
  fwrite(DT, filename, sep = " ", append = TRUE, compress = "gzip")
}

