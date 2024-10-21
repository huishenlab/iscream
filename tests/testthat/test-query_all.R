library(data.table)
library(Matrix)
extdata <- system.file("extdata", package = "iscream")
biscuit_bedfiles <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = T)
bismark_bedfiles <- list.files(extdata, pattern = "[a|b|c|d].cov.gz$", full.names = T)
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")

M <- matrix(c(1, 0, 0, 1,
              1, 0, 2, 2,
              0, 2, 0, 0,
              0, 1, 0, 0,
              1, 0, 1, 1,
              2, 0, 0, 0,
              3, 1, 0, 1),
            nrow = 7,
            byrow = TRUE)
Cov <- matrix(c(1, 2, 0, 1,
                 1, 0, 2, 2,
                 2, 2, 0, 0,
                 1, 1, 2, 1,
                 2, 0, 1, 2,
                 2, 2, 0, 0,
                 3, 1, 0, 1),
               nrow = 7,
               byrow = TRUE)
M_sp <- Matrix(M, sparse = T)
Cov_sp <- Matrix(Cov, sparse = T)
row_names <- c("chr1:2", "chr1:4", "chr1:6", "chr1:8", "chr1:10", "chr1:12", "chr1:14")
col_names <- c("a", "b", "c", "d")
rownames(M) <- row_names
colnames(M) <- col_names
rownames(Cov) <- row_names
colnames(Cov) <- col_names
rownames(M_sp) <- row_names
colnames(M_sp) <- col_names
rownames(Cov_sp) <- row_names
colnames(Cov_sp) <- col_names

exp_pos <- seq(1, 13, by = 2)
exp_seq <- rep("chr1", 7)

set_threads(1)
biscuit_test <- query_all(biscuit_bedfiles, regions, prealloc = 2)
bismark_test <- query_all(bismark_bedfiles, regions, aligner = "bismark")
biscuit_sparse_test <- query_all(biscuit_bedfiles, regions, sparse = T)
bismark_sparse_test <- query_all(bismark_bedfiles, regions, sparse = T, aligner = "bismark")
results <- list(biscuit_test, bismark_test, biscuit_sparse_test, bismark_sparse_test)

dense_class <- "matrix"
sparse_class <- "dgCMatrix"

test_types <- function(result_obj, sparse) {

  if (sparse) {
    cl <- sparse_class
  } else {
    cl <- dense_class
  }

  expect_equal(class(result_obj), "list")
  expect_equal(class(result_obj$M)[1], cl)
  expect_equal(class(result_obj$Cov)[1], cl)

}

test_that("matrix metadata", {
  test_types(biscuit_test, sparse = F)
  test_types(bismark_test, sparse = F)
  test_types(biscuit_sparse_test, sparse = T)
  test_types(bismark_sparse_test, sparse = T)
})

test_content <- function(result_obj, sparse) {

  if (sparse) {
    exp_m <- M_sp
    exp_c <- Cov_sp
  } else {
    exp_m <- M
    exp_c <- Cov
  }

  expect_equal(result_obj$M, exp_m)
  expect_equal(result_obj$Cov, exp_c)

}
test_that("matrix content", {
  lapply(results[1:2], test_content, FALSE)
  lapply(results[3:4], test_content, TRUE)
})

exp_dims <- c(7, length(biscuit_bedfiles))
test_dims <- function(result_obj, exp_dims) {
  expect_equal(dim(result_obj$M), exp_dims)
  expect_equal(dim(result_obj$Cov), exp_dims)
}

test_rownames <- function(result_obj) {
  expect_equal(rownames(result_obj$M), row_names)
  expect_equal(rownames(result_obj$Cov), row_names)
  expect_equal(colnames(result_obj$M), col_names)
  expect_equal(colnames(result_obj$Cov), col_names)
  expect_equal(result_obj$sampleNames, col_names)
  expect_equal(result_obj$sampleNames, col_names)
}

test_that("matrix row,colnames,sampleNames", {
  lapply(results, test_dims, exp_dims)
  lapply(results, test_rownames)
})
