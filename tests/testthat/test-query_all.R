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
  c(1, 0, 0, 1, 1, 0, 2, 2, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 2, 0, 0, 0, 3, 1, 0, 1),
  nrow = 7,
  byrow = TRUE
)

Cov <- matrix(
  c(1, 2, 0, 1, 1, 0, 2, 2, 2, 2, 0, 0, 1, 1, 2, 1, 2, 0, 1, 2, 2, 2, 0, 0, 3, 1, 0, 1),
  nrow = 7,
  byrow = TRUE
)

beta <- matrix(
  c(1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0.5, 0, 1, 0.5, 1, 0, 0, 0, 1, 1, 0, 1),
  nrow = 7,
  byrow = TRUE
)

M_sp <- Matrix(M, sparse = T)
Cov_sp <- Matrix(Cov, sparse = T)
beta_sp <- Matrix(beta, sparse = T)
M_perc <- beta * 100
M_perc_sp <- Matrix(M_perc, sparse = T)

row_names <- c("chr1:2", "chr1:4", "chr1:6", "chr1:8", "chr1:10", "chr1:12", "chr1:14")
col_names <- c("a", "b", "c", "d")
colnames(M) <- col_names
colnames(Cov) <- col_names
colnames(beta) <- col_names
colnames(M_sp) <- col_names
colnames(Cov_sp) <- col_names
colnames(beta_sp) <- col_names
colnames(M_perc) <- col_names
colnames(M_perc_sp) <- col_names

exp_pos <- seq(1, 13, by = 2)
exp_seq <- rep("chr1", 7)

set_threads(1)

biscuit_bs_MC <- make_bsseq_mat(biscuit_bedfiles, regions, prealloc = 2)
bismark_bs_MC <- make_bsseq_mat(bismark_bedfiles, regions, aligner = "bismark")

biscuit_bs_MC_sp <- make_bsseq_mat(biscuit_bedfiles, regions, sparse = T)
bismark_bs_MC_sp <- make_bsseq_mat(bismark_bedfiles, regions, sparse = T, aligner = "bismark")

biscuit_bs_MC_GR <- make_bsseq_mat(biscuit_bedfiles, gr, prealloc = 2)
biscuit_bs_MC_df <- make_bsseq_mat(biscuit_bedfiles, regions.dt, prealloc = 2)

biscuit_bs_bC <- make_bsseq_mat(biscuit_bedfiles, regions, sparse = F, mval = F)
bismark_bs_bC <- make_bsseq_mat(bismark_bedfiles, regions, sparse = F, aligner = "bismark", mval = F)

biscuit_bs_bC_sp <- make_bsseq_mat(biscuit_bedfiles, regions, sparse = T, mval = F)
bismark_bs_bC_sp <- make_bsseq_mat(bismark_bedfiles, regions, sparse = T, aligner = "bismark", mval = F)

biscuit_b <- make_mat(biscuit_bedfiles, regions, column = 4, mat_name = "beta")
bismark_M <- make_mat(bismark_bedfiles, regions, column = 4, mat_name = "M")

biscuit_b_sp <- make_mat(biscuit_bedfiles, regions, column = 4, sparse = T, mat_name = "beta")
bismark_M_sp <- make_mat(bismark_bedfiles, regions, column = 4, sparse = T, mat_name = "M")

biscuit_C_gr <- make_mat(biscuit_bedfiles, gr, column = 5, prealloc = 2)
biscuit_C_df <- make_mat(biscuit_bedfiles, regions.dt, column = 5, prealloc = 2)

results <- list(
  biscuit_bs_MC,
  bismark_bs_MC,

  biscuit_bs_MC_sp,
  bismark_bs_MC_sp,

  biscuit_bs_MC_GR,
  biscuit_bs_MC_df,

  biscuit_bs_bC,
  biscuit_bs_bC,

  biscuit_bs_bC_sp,
  bismark_bs_bC_sp,

  biscuit_b,
  bismark_M,

  biscuit_b_sp,
  bismark_M_sp,

  biscuit_C_gr
)

MC_results <- list(
  biscuit_bs_MC,
  bismark_bs_MC,
  biscuit_bs_MC_sp,
  bismark_bs_MC_sp,
  biscuit_bs_MC_GR,
  biscuit_bs_MC_df
)

bC_results <- list(
  biscuit_bs_bC,
  biscuit_bs_bC,
  biscuit_bs_bC_sp,
  bismark_bs_bC_sp
)

bs_results <- c(MC_results, bC_results)

mat_results <- list(
  biscuit_b,
  bismark_M,

  biscuit_b_sp,
  bismark_M_sp,

  biscuit_C_gr
)

test_names <- function(result_obj, sparse, mat_names) {
  expected_names <- c(mat_names, c("pos", "chr", "sampleNames"))
  expect_equal(names(result_obj), expected_names)
}


bsseq_names_M <- c("M", "Cov")
bsseq_names_b <- c("beta", "Cov")

bs_MC_results <- list(
  biscuit_bs_MC,
  bismark_bs_MC,
  biscuit_bs_MC_GR,
  biscuit_bs_MC_df
)

bs_MC_sp_results <- list(
  biscuit_bs_MC_sp,
  bismark_bs_MC_sp
)

bs_bC_results <- list(
  biscuit_bs_bC,
  bismark_bs_bC
)

bs_bC_sp_results <- list(
  biscuit_bs_bC_sp,
  bismark_bs_bC_sp
)

test_that("matrix names", {
  lapply(bs_MC_results, test_names, F, bsseq_names_M)
  lapply(bs_MC_sp_results, test_names, T, bsseq_names_M)
  lapply(bs_bC_results, test_names, F, bsseq_names_b)
  lapply(bs_bC_sp_results, test_names, T, bsseq_names_b)

  test_names(biscuit_b, F, "beta")
  test_names(biscuit_b_sp, T, "beta")

  test_names(bismark_M, F, "M")
  test_names(bismark_M_sp, T, "M")
})

test_types <- function(result_obj, sparse, mat_names) {
  dense_class <- "matrix"
  sparse_class <- "dgCMatrix"

  if (sparse) {
    cl <- sparse_class
  } else {
    cl <- dense_class
  }

  for (mat in mat_names) {
    expect_equal(class(result_obj), "list")
    expect_equal(class(result_obj[[mat]])[1], cl)
  }
}

test_that("matrix metadata", {
  lapply(bs_MC_results, test_types, F, bsseq_names_M)
  lapply(bs_MC_sp_results, test_types, T, bsseq_names_M)
  lapply(bs_bC_results, test_types, F, bsseq_names_b)
  lapply(bs_bC_sp_results, test_types, T, bsseq_names_b)

  test_types(biscuit_b, F, "beta")
  test_types(biscuit_b_sp, T, "beta")

  test_types(bismark_M, F, "M")
  test_types(bismark_M_sp, T, "M")
})

test_content <- function(result_obj, sparse, mval, test_mat = NULL, mat_name = NULL) {
  if (!is.null(test_mat) & !is.null(mat_name)) {
    return(expect_equal(result_obj[[mat_name]], test_mat))
  }

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
  lapply(bs_MC_results, test_content, FALSE, TRUE)
  lapply(bs_MC_results, test_content, FALSE, TRUE)
  lapply(bs_MC_sp_results, test_content, TRUE, TRUE)
  lapply(bs_bC_results, test_content, FALSE, FALSE)
  lapply(bs_bC_sp_results, test_content, TRUE, FALSE)

  test_content(biscuit_b, FALSE, test_mat = beta, mat_name = "beta")
  test_content(biscuit_b_sp, T, test_mat = beta_sp, mat_name = "beta")

  test_content(bismark_M, F, test_mat = M_perc, mat_name = "M")
  test_content(bismark_M_sp, T, test_mat = M_perc_sp, mat_name = "M")
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

test_that("SE, GR", {
  grmat <- make_mat_gr(biscuit_bedfiles, regions, column = 4, mat_name = "beta")
  expect_equal(class(grmat)[1], "GRanges")
  semat <- make_mat_se(biscuit_bedfiles, regions, column = 4, mat_name = "beta")
  expect_equal(class(semat)[1], "RangedSummarizedExperiment")
})

test_that("invalid columns", {
  lapply(1:3, function(i) {
    expect_error(
      make_mat(biscuit_bedfiles, regions, column = i, mat_name = "beta"),
      "`col` < 3 - must be a the index of a numeric data column not any of chr, start or end "
    )
  })
})
