options("iscream.threads" = 1)
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
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
})

library(data.table)

regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
regions.dt <- as.data.table(regions)[, tstrsplit(regions, ":|-")]

test_that("Test GRanges to string", {
  expect_error(get_df_string(regions.dt))
})

colnames(regions.dt) <- c("chr", "start", "end")
test_that("Test GRanges to string", {
  expect_equal(get_df_string(regions.dt), regions)
})
