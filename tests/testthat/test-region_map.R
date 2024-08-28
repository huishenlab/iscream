extdata <- system.file("extdata", package = "iscream")
bedfiles <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = TRUE)
regions <- c("chr1:1-6", "chr1:7-10", "chr1:11-14")

m_agg <- read.csv(file.path(extdata, "region_map_m_agg.test"))
m_ave <- read.csv(file.path(extdata, "region_map_m_ave.test"))
beta_ave <- read.csv(file.path(extdata, "region_map_beta_ave.test"))
beta_agg <- read.csv(file.path(extdata, "region_map_beta_agg.test"))

test_that("region_map 1 thread", {
  expect_equal(
    beta_agg,
    region_map(bedfiles, regions, fun = "aggregate", mval = FALSE, nthreads = 1)
  )
  expect_equal(
    m_agg,
    region_map(bedfiles, regions, fun = "aggregate", mval = TRUE, nthreads = 1)
  )
  expect_equal(
    beta_ave,
    region_map(bedfiles, regions, fun = "average", mval = FALSE, nthreads = 1)
  )
  expect_equal(
    m_ave,
    region_map(bedfiles, regions, fun = "average", mval = TRUE, nthreads = 1)
  )
})

test_that("region_map 2 thread", {
  expect_equal(
    beta_agg,
    region_map(bedfiles, regions, fun = "aggregate", mval = FALSE, nthreads = 2)
  )
  expect_equal(
    m_agg,
    region_map(bedfiles, regions, fun = "aggregate", mval = TRUE, nthreads = 2)
  )
  expect_equal(
    beta_ave,
    region_map(bedfiles, regions, fun = "average", mval = FALSE, nthreads = 2)
  )
  expect_equal(
    m_ave,
    region_map(bedfiles, regions, fun = "average", mval = TRUE, nthreads = 2)
  )
})
