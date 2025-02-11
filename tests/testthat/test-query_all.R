library(data.table)
library(Matrix)
extdata <- system.file("extdata", package = "iscream")
biscuit_bedfiles <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = T)
bismark_bedfiles <- list.files(extdata, pattern = "[a|b|c|d].cov.gz$", full.names = T)
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")
regions.dt <- as.data.table(regions)[, tstrsplit(regions, ":|-")]
colnames(regions.dt) <- c("chr", "start", "end")
gr <- GenomicRanges::GRanges(regions)

M <- matrix(
  c(1, 0, 0, 1,
    1, 0, 2, 2,
    0, 2, 0, 0,
    0, 1, 0, 0,
    1, 0, 1, 1,
    2, 0, 0, 0,
    3, 1, 0, 1),
  nrow = 7,
  byrow = TRUE)

Cov <- matrix(
  c(1, 2, 0, 1,
    1, 0, 2, 2,
    2, 2, 0, 0,
    1, 1, 2, 1,
    2, 0, 1, 2,
    2, 2, 0, 0,
    3, 1, 0, 1),
  nrow = 7,
  byrow = TRUE)

beta <- matrix(
  c(1, 0, 0, 1,
    1, 0, 1, 1,
    0, 1, 0, 0,
    0, 1, 0, 0,
    0.5, 0, 1, 0.5,
    1, 0, 0, 0,
    1, 1, 0, 1),
  nrow = 7,
  byrow = TRUE
)

M_sp <- Matrix(M, sparse = T)
Cov_sp <- Matrix(Cov, sparse = T)
beta_sp <- Matrix(beta, sparse = T)

row_names <- c("chr1:2", "chr1:4", "chr1:6", "chr1:8", "chr1:10", "chr1:12", "chr1:14")
col_names <- c("a", "b", "c", "d")
colnames(M) <- col_names
colnames(Cov) <- col_names
colnames(beta) <- col_names
colnames(M_sp) <- col_names
colnames(Cov_sp) <- col_names
colnames(beta_sp) <- col_names

exp_pos <- seq(1, 13, by = 2)
exp_seq <- rep("chr1", 7)

set_threads(1)
biscuit_test <- query_all(biscuit_bedfiles, regions, prealloc = 2)
bismark_test <- query_all(bismark_bedfiles, regions, aligner = "bismark")

biscuit_sparse_test <- query_all(biscuit_bedfiles, regions, sparse = T)
bismark_sparse_test <- query_all(bismark_bedfiles, regions, sparse = T, aligner = "bismark")

biscuit_granges_test <- query_all(biscuit_bedfiles, gr, prealloc = 2)
biscuit_df_test <- query_all(biscuit_bedfiles, regions.dt, prealloc = 2)

biscuit_beta_test <- query_all(biscuit_bedfiles, regions, sparse = F, mval = F)
bismark_beta_test <- query_all(bismark_bedfiles, regions, sparse = F, aligner = "bismark", mval = F)

biscuit_sparse_beta_test <- query_all(biscuit_bedfiles, regions, sparse = T, mval = F)
bismark_sparse_beta_test <- query_all(bismark_bedfiles, regions, sparse = T, aligner = "bismark", mval = F)

results <- list(
  biscuit_test,
  bismark_test,
  biscuit_sparse_test,
  bismark_sparse_test,
  biscuit_granges_test,
  biscuit_df_test,
  biscuit_sparse_beta_test,
  bismark_sparse_beta_test
)

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
  test_types(biscuit_granges_test, sparse = F)
  test_types(biscuit_df_test, sparse = F)
})

test_content <- function(result_obj, sparse, mval) {

  if (sparse) {
    exp_m <- if (mval) {
      M_sp
    } else {
      beta_sp
    }
    exp_c <- Cov_sp
  } else {
    exp_m <- if (mval) {
      M
    } else {
      beta
    }
    exp_c <- Cov
  }


  m_beta <- ifelse(mval, result_obj$M, result_obj$beta)
  if (mval) {
    expect_equal(result_obj$M, exp_m)
  } else {
    expect_equal(result_obj$beta, exp_m)
  }
  expect_equal(result_obj$Cov, exp_c)

}

test_that("matrix content", {
  lapply(results[1:2], test_content, FALSE, TRUE)
  lapply(results[5], test_content, FALSE, TRUE)
  lapply(results[3:4], test_content, TRUE, TRUE)
  lapply(results[7:8], test_content, TRUE, FALSE)
})

exp_dims <- c(7, length(biscuit_bedfiles))
test_dims <- function(result_obj, exp_dims, mval) {
  if ("beta" %in% names(result_obj)) {
    mbeta_name <- "beta"
  } else {
    mbeta_name <- "M"
  }
  expect_equal(dim(result_obj[[mbeta_name]]), exp_dims)
  expect_equal(dim(result_obj$Cov), exp_dims)
}

test_that("matrix row,colnames,sampleNames", {
  lapply(results[1], test_dims, exp_dims)
})

