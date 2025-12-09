options("iscream.threads" = 1)
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")
gr <- GenomicRanges::GRanges(regions)
feature_col <- c("AA", "BB", "CC")
GenomicRanges::mcols(gr) <- data.frame(feature = feature_col)

regions_single <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14", D = "chr2:5", E = "chrX:3")
regions_single_res <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14", D = "chr2:5-5", E = "chrX:3-3")
gr_single <- GenomicRanges::GRanges(regions_single)

with_mocked_bindings(
  {
    expect_error(
      get_granges_string(gr),
      "The 'GenomicRanges' package must be installed for this functionality"
    )
  },
  requireNamespace = function(package, ..., quietly = FALSE) FALSE,
  .package = "base"
)

test_that("Test GRanges to string", {
  str <- get_granges_string(gr)
  expect_equal(str, regions)
  expect_equal(names(str), names(regions))

  str <- get_granges_string(gr_single)
  expect_equal(str, regions_single_res)
  expect_equal(names(str), names(regions_single))

  str <- get_granges_string(gr_single, feature_col = "feature")
  expect_equal(str, regions_single_res)
  expect_equal(names(str), names(regions_single))

  str <- get_granges_string(gr)
  expect_equal(str, regions)
  expect_equal(names(str), names(regions))

  str <- get_granges_string(gr, feature_col = "feature")
  expect_equal(unname(str), unname(regions))
  expect_equal(names(str), feature_col)
})

library(data.table)

regions.dt <- as.data.table(regions, keep.rownames = TRUE)[, tstrsplit(regions, ":|-")]

test_that("Test GRanges to string", {
  expect_error(get_df_string(regions.dt))
})

colnames(regions.dt) <- c("c", "start", "end")
test_that("Test GRanges to string", {
  expect_error(get_df_string(regions.dt), "colnames must be 'chr', 'start' and 'end'")
})

colnames(regions.dt) <- c("chr", "start", "end")
regions.with_seqnames <- regions.dt
setnames(regions.with_seqnames, "chr", "seqnames")

colnames(regions.dt) <- c("chr", "start", "end")
test_that("Test GRanges to string no names", {
  expect_equal(get_df_string(regions.dt), unname(regions))
  expect_equal(get_df_string(regions.with_seqnames), unname(regions))
})

regions.dt[, names := names(regions)]
colnames(regions.dt) <- c("chr", "start", "end", "names")
test_that("Test GRanges to string with names", {
  expect_equal(get_df_string(regions.dt, feature_col = "names"), regions)
  expect_equal(names(get_df_string(regions.dt, feature_col = "names")), names(regions))
})
