options("iscream.threads" = 1)

regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")
gr <- GenomicRanges::GRanges(regions)

test_that("Test GRanges to string", {
  expect_equal(get_granges_string(gr), regions)
})


