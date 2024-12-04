options("iscream.threads" = 1)
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")
gr <- GenomicRanges::GRanges(regions)

with_mocked_bindings(
    {
    expect_error(
      get_granges_string(gr),
      "The 'GenomicRanges' package must be installed for this functionality"
    )
  },
  requireNamespace = function(package, ..., quietly=FALSE) FALSE,
  .package="base"
)

test_that("Test GRanges to string", {
  expect_equal(get_granges_string(gr), regions)
  expect_equal(names(get_granges_string(gr)), names(regions))
})

library(data.table)

regions.dt <- as.data.table(regions, keep.rownames = TRUE)[, tstrsplit(regions, ":|-")]

test_that("Test GRanges to string", {
  expect_error(get_df_string(regions.dt))
})

colnames(regions.dt) <- c("chr", "start", "end")
test_that("Test GRanges to string no names", {
  expect_equal(get_df_string(regions.dt), unname(regions))
})

regions.dt[, names := names(regions)]
colnames(regions.dt) <- c("chr", "start", "end", "names")
test_that("Test GRanges to string with names", {
  expect_equal(get_df_string(regions.dt), regions)
  expect_equal(names(get_df_string(regions.dt)), names(regions))
})

