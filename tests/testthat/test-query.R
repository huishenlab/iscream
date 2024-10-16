extdata <- system.file("extdata/test_chroms", package = "iscream")
bedfiles <- list.files(extdata, pattern = "*.bed.gz$", full.names = T)
chrs <- sort(paste0("chr", c(seq(1:22), "M", "X", "Y")))

test_that("query_chroms", {
  expect_equal(
    chrs,
    query_chroms(bedfiles)
  )
})

