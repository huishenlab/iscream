library("stringfish")

# inputs
extdata <- system.file("extdata", package = "iscream")
bedfiles <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")
supported_funcs <- c("sum", "mean", "median", "stddev", "variance", "range")

# utils
get_colnames <- function(mval, funcs) {
  base_colnames <- c("Feature", "Sample")
  values <- c("coverage", ifelse(mval, "M", "beta"))
    function_used <- funcs
  if ("all" %in% funcs) {
    function_used <- supported_funcs
  }
  c(base_colnames, as.vector(outer(values, function_used, paste, sep = ".")))
}

run_test <- function(bedfiles, regions, funcs, mval, nthreads) {
  colnames <- get_colnames(mval, funcs)
  list(
    df = summarize_regions(bedfiles, regions, fun = funcs, mval = mval, nthreads = nthreads),
    colnames = colnames,
    dims = c(length(bedfiles) * length(regions), length(colnames))
  )
}

one_file <- lapply(bedfiles, function(i) {
  run_test(
    i,
    regions,
    funcs = "all",
    mval = FALSE,
    nthreads = 1
  )
})

one_fun <- lapply(supported_funcs, function(i) {
  run_test(
    bedfiles,
    regions,
    funcs = i,
    mval = FALSE,
    nthreads = 1
  )
})

two_fun <- lapply(asplit(combn(supported_funcs, 2), 2), function(i) {
  run_test(
    bedfiles,
    regions,
    funcs = i,
    mval = FALSE,
    nthreads = 1
  )
})

one_reg <- lapply(regions, function(i) {
  run_test(
    bedfiles,
    i,
    funcs = "all",
    mval = FALSE,
    nthreads = 1
  )
})

mbeta <- lapply(c(T, F), function(i) {
  run_test(
    bedfiles,
    regions,
    funcs = "all",
    mval = i,
    nthreads = 1
  )
})

threads <- lapply(c(1, 1), function(i) {
  run_test(
    bedfiles,
    regions,
    funcs = "all",
    mval = F,
    nthreads = i
  )
})

full_set <- list(
  one_file,
  one_fun,
  two_fun,
  one_reg,
  mbeta,
  threads
)

test_dims <- function(set) {
  lapply(set, function(i) {
    expect_equal(
      i$dims,
      dim(i$df)
    )
  })
}

test_colnames <- function(set) {
  lapply(set, function(i) {
    expect_equal(
      i$colnames,
      colnames(i$df)
    )
  })
}

test_that("dims", {
  lapply(full_set, function(i) {
    test_dims(i)
  })
})

test_that("colnames", {
  lapply(full_set, function(i) {
    test_colnames(i)
  })
})

# error on bad 'fun' arguments
test_that("bad fun", {
  expect_error(
    summarize_regions(bedfiles, regions, fun = "none")
  )
  expect_error(
    summarize_regions(bedfiles, regions, fun = c("all", "mean"))
  )
  expect_error(
    summarize_regions(bedfiles, regions, fun = c("none", "mean"))
  )
})

m_sum <- read.csv(file.path(extdata, "summarize_regions_m_sum.test"))
beta_sum <- read.csv(file.path(extdata, "summarize_regions_beta_sum.test"))

m_mean <- read.csv(file.path(extdata, "summarize_regions_m_mean.test"))
beta_mean <- read.csv(file.path(extdata, "summarize_regions_beta_mean.test"))

beta_all <- read.csv(file.path(extdata, "summarize_regions_beta_all.test"))
m_all <- read.csv(file.path(extdata, "summarize_regions_m_all.test"))

test_that("summarize_regions 1 thread sum", {
  expect_equal(
    m_sum,
    summarize_regions(bedfiles, regions, fun = "sum", mval = TRUE, nthreads = 1)
  )
  expect_equal(
    beta_sum,
    summarize_regions(bedfiles, regions, fun = "sum", mval = FALSE, nthreads = 1)
  )
})

test_that("summarize_regions 2 thread sum", {
  expect_equal(
    m_sum,
    summarize_regions(bedfiles, regions, fun = "sum", mval = TRUE, nthreads = 2)
  )
  expect_equal(
    beta_sum,
    summarize_regions(bedfiles, regions, fun = "sum", mval = FALSE, nthreads = 2)
  )
})

test_that("summarize_regions 1 thread sum", {
  expect_equal(
    m_sum,
    summarize_regions(bedfiles, regions, fun = "sum", mval = TRUE, nthreads = 1)
  )
  expect_equal(
    beta_sum,
    summarize_regions(bedfiles, regions, fun = "sum", mval = FALSE, nthreads = 1)
  )
})

test_that("summarize_regions 2 thread sum", {
  expect_equal(
    m_sum,
    summarize_regions(bedfiles, regions, fun = "sum", mval = TRUE, nthreads = 2)
  )
  expect_equal(
    beta_sum,
    summarize_regions(bedfiles, regions, fun = "sum", mval = FALSE, nthreads = 2)
  )
})


test_that("summarize_regions 1 thread all", {
  expect_equal(
    m_all,
    summarize_regions(bedfiles, regions, fun = "all", mval = TRUE, nthreads = 1)
  )
  expect_equal(
    beta_all,
    summarize_regions(bedfiles, regions, fun = "all", mval = FALSE, nthreads = 1)
  )
})

test_that("summarize_regions 2 thread all", {
  expect_equal(
    m_all,
    summarize_regions(bedfiles, regions, fun = "all", mval = TRUE, nthreads = 2)
  )
  expect_equal(
    beta_all,
    summarize_regions(bedfiles, regions, fun = "all", mval = FALSE, nthreads = 2)
  )
})
