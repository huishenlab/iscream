library("stringfish")
library("data.table")

# inputs
extdata <- system.file("extdata", package = "iscream")
biscuit_bedfiles <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
bismark_bedfiles <- list.files(extdata, pattern = "[a|b|c|d].cov.gz$", full.names = TRUE)
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")
regions.dt <- as.data.table(regions)[, tstrsplit(regions, ":|-")][, names := names(regions)]
colnames(regions.dt) <- c("chr", "start", "end", "names")
gr <- GenomicRanges::GRanges(regions)
supported_funcs <- c("sum", "mean", "median", "stddev", "variance", "min", "max", "range")

get_colnames <- function(funcs, col_names) {
  base_colnames <- c("feature", "file")
  values <- c(col_names)
  if ("all" %in% funcs) {
    return(c(base_colnames, as.vector(outer(values, supported_funcs, paste, sep = ".")), "count"))
  }

  if ("count" %in% funcs) {
    index <- which(funcs == "count")
    funcs <- funcs[-index]
  }
  function_used <- funcs
  col_names <- c(
    base_colnames,
    as.vector(outer(values, function_used, paste, sep = "."))
  )

  if ("count" %in% funcs) {
  col_names <- c(col_names, "count")
  }

  return(col_names)
}

run_test <- function(bedfiles, regions, funcs, columns, col_names, nthreads) {
  test_colnames <- get_colnames(funcs, col_names)
  reg_length <- if ("data.frame" %in% class(regions)) {
    nrow(regions)
  } else {
      length(regions)
  }

  list(
    df = summarize_regions(bedfiles, regions, fun = funcs, columns = columns, col_names = col_names, nthreads = nthreads),
    colnames = test_colnames,
    dims = c(length(bedfiles) * reg_length, length(test_colnames))
  )
}

one_file <- lapply(biscuit_bedfiles, function(i) {
  run_test(
    i,
    regions,
    funcs = "all",
    c(4, 5),
    c("beta", "coverage"),
    nthreads = 1
  )
})

one_fun <- lapply(supported_funcs, function(i) {
  run_test(
    biscuit_bedfiles,
    regions,
    funcs = i,
    c(4, 5),
    c("beta", "coverage"),
    nthreads = 1
  )
})

two_fun <- lapply(asplit(combn(supported_funcs, 2), 2), function(i) {
  run_test(
    biscuit_bedfiles,
    regions,
    funcs = i,
    c(4, 5),
    c("cov", "beta"),
    nthreads = 1
  )
})

one_reg <- lapply(regions, function(i) {
  run_test(
    biscuit_bedfiles,
    i,
    funcs = "all",
    c(4, 5),
    c("beta", "coverage"),
    nthreads = 1
  )
})

threads <- lapply(c(1, 1), function(i) {
  run_test(
    biscuit_bedfiles,
    regions,
    funcs = "all",
    c(4, 5),
    c("beta", "coverage"),
    nthreads = i
  )
})

regions_df <- list(run_test(
  biscuit_bedfiles,
  regions.dt,
  funcs = "all",
  c(4, 5),
  c("beta", "coverage"),
  nthreads = 1
))

regions_gr <- list(run_test(
  biscuit_bedfiles,
  gr,
  funcs = "all",
  c(4, 5),
  c("beta", "coverage"),
  nthreads = 1
))


full_set <- list(
  one_file,
  one_fun,
  two_fun,
  one_reg,
  threads,
  regions_df,
  regions_gr
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
  lapply(full_set[6], function(i) {
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
    summarize_regions(biscuit_bedfiles, regions, fun = "none")
  )
  expect_error(
    summarize_regions(biscuit_bedfiles, regions, fun = c("all", "mean"))
  )
  expect_error(
    summarize_regions(biscuit_bedfiles, regions, fun = c("none", "mean"))
  )
})

