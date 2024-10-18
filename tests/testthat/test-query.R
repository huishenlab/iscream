library(data.table)
extdata <- system.file("extdata", package = "iscream")
chrom_beds <- list.files(paste0(extdata, "/test_chroms"), pattern = "*.bed.gz$", full.names = T)
biscuit_tabix_beds <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = T)
bismark_tabix_beds <- list.files(extdata, pattern = "[a|b|c|d].cov.gz$", full.names = T)
tabix_df_result <- list.files(extdata, pattern = "tabix_dataframe.test", full.names = T)
tabix_df_result_bismark <- list.files(extdata, pattern = "tabix_dataframe_bismark.test", full.names = T)
mergecg_bed <- list.files(extdata, pattern = "*_mergecg.bed.gz$", full.names = T)
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
    tabix(biscuit_tabix_beds[1], regions),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
  )
  expect_equal(
    tabix(bismark_tabix_beds[1], regions, aligner = 'bismark'),
    fread(tabix_df_result_bismark, colClasses = c("character", "numeric", "numeric", "numeric", "numeric", "numeric"))
  )
})

test_that("tabix raw list", {
  expect_equal(
    tabix(biscuit_tabix_beds[3], regions, raw = T),
    tabix_raw_res
  )
})

test_that("tabix dataframe", {
  expect_warning(
    tabix(chrom_beds[2], regions),
    "No records found"
  )
})

test_that("tabix > 1 file", {
  expect_error(
    tabix(chrom_beds, regions),
    "Cannot tabix multiple files - only single-file queries are currently supported"
  )
})

custom_cols <- paste0("c", 1:5)
test_that("tabix custom colnames", {
  expect_equal(
    colnames(tabix(chrom_beds[1], regions, colnames = custom_cols)),
    custom_cols
  )
})

custom_cols <- paste0("c", 1:6)
test_that("tabix custom colnames too many", {
  expect_warning(
    tabix(chrom_beds[1], regions, colnames = custom_cols),
    "Fewer columns in data than provided colnames"
  )
})

custom_cols <- paste0("c", 1:4)
test_that("tabix custom colnames too few", {
  expect_warning(
    tabix(chrom_beds[1], regions, colnames = custom_cols),
    paste("Did not use input 'colnames' - only", length(custom_cols), "names provided for 5 column data.table")
  )
})


base_colnames <- c("chr", "start", "end")
biscuit_colnames <- c("beta", "coverage")
bismark_colnames <- c("methylation.percentage", "count.methylated", "count.unmethylated")
test_that("tabix mergecg colnames", {
  expect_equal(
    colnames(tabix(mergecg_bed, regions, aligner = "biscuit")),
    c(base_colnames, biscuit_colnames, "mergecg")
  )
})


