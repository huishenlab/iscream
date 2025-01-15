library(data.table)
extdata <- system.file("extdata", package = "iscream")
chrom_beds <- list.files(paste0(extdata, "/test_chroms"), pattern = "*.bed.gz$", full.names = T)
biscuit_tabix_beds <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = T)
bismark_tabix_beds <- list.files(extdata, pattern = "[a|b|c|d].cov.gz$", full.names = T)
tabix_df_result <- list.files(extdata, pattern = "tabix_dataframe.test", full.names = T)
tabix_df_result_bismark <- list.files(extdata, pattern = "tabix_dataframe_bismark.test", full.names = T)
mergecg_bed <- list.files(extdata, pattern = "*_mergecg.bed.gz$", full.names = T)
regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")
regions.dt <- as.data.table(regions)[, tstrsplit(regions, ":|-")]
colnames(regions.dt) <- c("chr", "start", "end")
gr <- GenomicRanges::GRanges(regions)
gr.meta <- GenomicRanges::GRanges(regions)
GenomicRanges::values(gr.meta) <- data.frame(meta = c("s1", "s2", "s3"))

tabix_raw_res <- list(
    `chr1:1-6` =  c("chr1\t3\t4\t1.000\t2"),
    `chr1:7-10` = c("chr1\t7\t8\t0.000\t2", "chr1\t9\t10\t1.000\t1"),
    `chr1:11-14` = character(0)
)

test_that("tabix dataframe", {
  expect_equal(
    tabix(biscuit_tabix_beds[1], regions),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
  )
  expect_equal(
    tabix(biscuit_tabix_beds[1], gr),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric")) |> GenomicRanges::GRanges()
  )
  expect_equal(
    tabix(biscuit_tabix_beds[1], gr.meta),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))[,
      meta := c(rep("s1", 3), rep("s2", 2), rep("s3", 2))
    ] |> GenomicRanges::GRanges()
  )
  expect_equal(
    tabix(biscuit_tabix_beds[1], regions.dt),
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


multi_tabix <- lapply(biscuit_tabix_beds, function(i) {
  tabix(i, regions)[,
  sample := tools::file_path_sans_ext(basename(i), compression = TRUE)
  ]
}) |> rbindlist()

test_that("tabix multi query", {
  expect_equal(
    tabix(biscuit_tabix_beds, regions),
    multi_tabix
  )
})

multi_tabix_raw <- lapply(biscuit_tabix_beds, function(i) {
  tabix(i, regions, raw = T)
})
names(multi_tabix_raw) = tools::file_path_sans_ext(basename(biscuit_tabix_beds), compression = TRUE)

test_that("tabix multi query raw", {
  expect_equal(
    tabix(biscuit_tabix_beds, regions, raw = T),
    multi_tabix_raw
  )
})

