library("stringfish")

# inputs
extdata <- system.file("extdata", package = "iscream")
bedfiles <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
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
