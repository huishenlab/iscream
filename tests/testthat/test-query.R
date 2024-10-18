library(data.table)
extdata <- system.file("extdata", package = "iscream")
chrom_beds <- list.files(paste0(extdata, "/test_chroms"), pattern = "*.bed.gz$", full.names = T)
tabix_beds <- list.files(extdata, pattern = "*.bed.gz$", full.names = T)
tabix_df_result <- list.files(extdata, pattern = "tabix_dataframe.test", full.names = T)
chrs <- sort(paste0("chr", c(seq(1:22), "M", "X", "Y")))
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")

tabix_raw_res <- list(
    `chr1:1-6` =  c("chr1\t3\t4\t1.000\t2"),
    `chr1:7-10` = c("chr1\t7\t8\t0.000\t2", "chr1\t9\t10\t1.000\t1"),
    `chr1:11-14` = character(0)
)

test_that("query_chroms", {
  expect_equal(
    chrs,
    query_chroms(chrom_beds)
  )
})

test_that("tabix dataframe", {
  expect_equal(
    tabix(tabix_beds[1], regions),
    fread(tabix_df_result, colClasses = c("character", "integer", "integer", "numeric", "integer"))
  )
})

test_that("tabix raw list", {
  expect_equal(
    tabix(tabix_beds[3], regions, raw = T),
    tabix_raw_res
  )
})

test_that("tabix dataframe", {
  expect_warning(
    tabix(chrom_beds[2], regions),
    "No records found"
  )
})
