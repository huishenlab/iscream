library(data.table)
tabix_installed <- Sys.which("tabix") != ""
options("tabix.method" = "htslib")

# input data
extdata <- system.file("extdata", package = "iscream")
chrom_beds <- list.files(paste0(extdata, "/test_chroms"), pattern = "*.bed.gz$", full.names = T)
biscuit_tabix_beds <- list.files(extdata, pattern = "[a|b|c|d].bed.gz$", full.names = T)
bismark_tabix_beds <- list.files(extdata, pattern = "[a|b|c|d].cov.gz$", full.names = T)
mergecg_bed <- list.files(extdata, pattern = "*_mergecg.bed.gz$", full.names = T)

regions <- c(A = "chr1:1-6", B = "chr1:7-10", C = "chr1:11-14")
regions.missing_in_3 <- c("chr1:13-14")
regions.dt <- as.data.table(regions)[, tstrsplit(regions, ":|-")]
colnames(regions.dt) <- c("chr", "start", "end")
chrs <- sort(paste0("chr", c(seq(1:22), "M", "X", "Y")))
gr <- GenomicRanges::GRanges(regions)

test_that("query_chroms", {
  expect_equal(
    chrs,
    query_chroms(chrom_beds)
  )
})

# testing data
tabix_df_result <- list.files(extdata, pattern = "tabix_dataframe.test", full.names = T)
tabix_df_result_bismark <- list.files(extdata, pattern = "tabix_dataframe_bismark.test", full.names = T)

# test dataframe output
test_tabix_dataframe <- function(shell = FALSE) {
  if (shell) {
    options("tabix.method" = "shell")
  }
  expect_equal(
    tabix(biscuit_tabix_beds[1], regions, aligner = "biscuit"),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
  )
  expect_equal(
    tabix(biscuit_tabix_beds[1], regions.dt, aligner = "biscuit"),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
  )
  expect_equal(
    tabix(biscuit_tabix_beds[1], gr, aligner = "biscuit"),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
  )
  expect_equal(
    tabix(bismark_tabix_beds[1], regions, aligner = 'bismark'),
    fread(tabix_df_result_bismark, colClasses = c("character", "numeric", "numeric", "numeric", "numeric", "numeric"))
  )
  # empty
  expect_warning(
    tabix(chrom_beds[2], regions),
    "No records found"
  )
}

test_that("tabix dataframe with shell", {
  skip_if_not(tabix_installed)
  test_tabix_dataframe(shell = TRUE)
})

test_that("tabix dataframe with htslib", {
  test_tabix_dataframe(shell = FALSE)
})

test_multi_tabix_dataframe <- function(shell = FALSE) {
  if (shell) {
    options("tabix.method" = "shell")
  }
  multi_tabix <- lapply(biscuit_tabix_beds, function(i) {
    tabix(i, regions)[,
      file := tools::file_path_sans_ext(basename(i), compression = TRUE)
    ]
  }) |>
    rbindlist()
  multi_tabix_biscuit <- lapply(biscuit_tabix_beds, function(i) {
    tabix(i, regions, aligner = "biscuit")[,
      file := tools::file_path_sans_ext(basename(i), compression = TRUE)
    ]
  }) |>
    rbindlist()
  multi_tabix_bismark <- lapply(bismark_tabix_beds, function(i) {
    tabix(i, regions, aligner = "bismark")[,
      file := tools::file_path_sans_ext(basename(i), compression = TRUE)
    ]
  }) |>
    rbindlist()
  multi_tabix_custom <- lapply(biscuit_tabix_beds, function(i) {
    tabix(i, regions, col.names = c("1", "2"))[,
      file := tools::file_path_sans_ext(basename(i), compression = TRUE)
    ]
  }) |>
    rbindlist()

  expect_equal(
    tabix(biscuit_tabix_beds, regions),
    multi_tabix
  )
  expect_equal(
    tabix(biscuit_tabix_beds, regions, aligner = "biscuit"),
    multi_tabix_biscuit
  )
  expect_equal(
    tabix(bismark_tabix_beds, regions, aligner = "bismark"),
    multi_tabix_bismark
  )
  expect_equal(
    tabix(biscuit_tabix_beds, regions, col.names = c("1", "2")),
    multi_tabix_custom
  )
  expect_equal(
    suppressWarnings(tabix(biscuit_tabix_beds, regions.missing_in_3)),
    tabix(biscuit_tabix_beds[-3], regions.missing_in_3)
  )
}

test_that("tabix multi query with shell", {
  skip_if_not(tabix_installed)
  test_multi_tabix_dataframe(shell = TRUE)
})

test_that("tabix multi query with htslib", {
  test_multi_tabix_dataframe(shell = FALSE)
})

gr.meta <- GenomicRanges::GRanges(regions)
GenomicRanges::values(gr.meta) <- data.frame(meta = c("s1", "s2", "s3"))

test_that("tabix_gr", {
  # gr input
  expect_equal(
    tabix_gr(biscuit_tabix_beds[1], gr, aligner = "biscuit"),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric")) |>
      GenomicRanges::makeGRangesFromDataFrame(starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  )
  # string input
  expect_equal(
    tabix_gr(biscuit_tabix_beds[1], regions, aligner = "biscuit"),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric")) |>
      GenomicRanges::makeGRangesFromDataFrame(starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  )
  # dt input
  expect_equal(
    tabix_gr(biscuit_tabix_beds[1], regions.dt, aligner = "biscuit"),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric")) |>
      GenomicRanges::makeGRangesFromDataFrame(starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  )
  # gr input with metadata
  expect_equal(
    tabix_gr(biscuit_tabix_beds[1], gr.meta, aligner = "biscuit"),
    fread(tabix_df_result, colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))[,
      meta := c(rep("s1", 3), rep("s2", 2), rep("s3", 2))
    ] |>
      GenomicRanges::makeGRangesFromDataFrame(starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
  )
  # gr input with metadata, bismark
  expect_equal(
    tabix_gr(bismark_tabix_beds[1], gr.meta, aligner = 'bismark', zero_based = FALSE),
    fread(tabix_df_result_bismark, colClasses = c("character", "numeric", "numeric", "numeric", "numeric", "numeric"))[,
      meta := c(rep("s1", 3), rep("s2", 2), rep("s3", 2))
    ] |>
      GenomicRanges::makeGRangesFromDataFrame(starts.in.df.are.0based = FALSE, keep.extra.columns = TRUE)
  )
  expect_error(suppressWarnings(tabix_gr(biscuit_tabix_beds[1], "chrM:1-10")))
})

test_that("tabix multi GR(List)", {
  grl <- GenomicRanges::GRangesList(
    lapply(seq_len(length(biscuit_tabix_beds)), function(bed) {
      res <- tabix_gr(biscuit_tabix_beds[bed], gr)
      GenomicRanges::mcols(res)$file <- letters[bed]
      res
    })
  )
  names(grl) <- c(letters[1:4])

  res <- tabix_gr(biscuit_tabix_beds, gr)
  expect_s4_class(res, "CompressedGRangesList")
  expect_equal(names(res), letters[1:4])
  expect_equal(res, grl)
})

# test colnames
base_colnames <- c("chr", "start", "end")
custom_cols <- paste0("c", 1:2)
test_that("tabix custom colnames", {
  expect_equal(
    colnames(tabix(chrom_beds[1], regions, col.names = custom_cols)),
    c(base_colnames, custom_cols)
  )
})

custom_cols <- paste0("c", 1:6)
test_that("tabix custom colnames too many", {
  expect_warning(
    tabix(chrom_beds[1], regions, col.names = custom_cols),
    "Fewer columns in data than provided colnames"
  )
})

custom_cols <- paste0("c1")
test_that("tabix custom colnames too few", {
  expect_warning(
    tabix(chrom_beds[1], regions, col.names = custom_cols),
    paste("Did not use input data 'colnames' - only", length(custom_cols), "names provided for 2 data columns")
  )
})


biscuit_colnames <- c("beta", "coverage")
bismark_colnames <- c("methylation.percentage", "count.methylated", "count.unmethylated")
test_that("tabix mergecg colnames", {
  expect_equal(
    colnames(tabix(mergecg_bed, regions, aligner = "biscuit")),
    c(base_colnames, biscuit_colnames, "mergecg")
  )
})

# test raw output
tabix_raw_res <- list(
  `chr1:1-6` = c("chr1\t2\t4\t1.000\t2"),
  `chr1:7-10` = c("chr1\t6\t8\t0.000\t2", "chr1\t8\t10\t1.000\t1"),
  `chr1:11-14` = character(0)
)

test_that("tabix raw list", {
  expect_equal(
    tabix_raw(biscuit_tabix_beds[3], regions),
    tabix_raw_res
  )
})


multi_tabix_raw <- lapply(biscuit_tabix_beds, function(i) {
  tabix_raw(i, regions)
})
names(multi_tabix_raw) = tools::file_path_sans_ext(basename(biscuit_tabix_beds), compression = TRUE)

test_that("tabix multi query raw", {
  expect_equal(
    tabix_raw(biscuit_tabix_beds, regions),
    multi_tabix_raw
  )
})


test_that("check tabix.method", {
  with_mocked_bindings(
    {
      tabix_msg <- "tabix' executable not found in $PATH. tabix() will use htslib to make queries instead which can be slower. See ?tabix for details."
      expect_true(grepl(tabix_msg, iscream:::package_loader(), fixed = TRUE))
      expect_equal(getOption("tabix.method"), "htslib")
    },
    Sys.which = function(input) "",
    .package = "base"
  )
  with_mocked_bindings(
    {
      iscream:::package_loader()
      expect_equal(getOption("tabix.method"), "shell")
    },
    Sys.which = function(input) "tabix",
    .package = "base"
  )
})
