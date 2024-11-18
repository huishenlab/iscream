data_loc <- system.file("extdata", package = "iscream")
bedfiles <- list.files(data_loc, pattern = "*.bed.gz$", full.names = TRUE)
bedfiles_no_path <- list.files(data_loc, pattern = "[a|b].bed.gz$", full.names = FALSE)

test_that("check_files_exist", {
  expect_error(
    iscream:::check_files_exist(bedfiles_no_path[1]),
    paste0("Bedfile: ", bedfiles_no_path[1], " could not be found")
  )
  expect_error(iscream:::verify_files_or_stop(bedfiles_no_path))
  non_tabixed_file <- list.files(data_loc, pattern = "no_tabix.bed.gz", full.names = TRUE)
  expect_error(
    iscream:::verify_files_or_stop(non_tabixed_file),
    paste0("Tabix file: ", non_tabixed_file, ".tbi could not be found")
  )
})

valid_regions <- c("1:1-5", "X:1-5", "chr10:6654-76548", "chr1", "22", "X", "Y", "M")
test_that("verify_regions_or_stop valid", {
  expect_no_error(iscream:::verify_regions_or_stop(valid_regions, nthreads = 1))
  expect_no_error(iscream:::verify_regions_or_stop(valid_regions, nthreads = 2))
})

invalid_regions <- c(
  "ch:1-5",     # malformed chr
  "chri:5-65",  # malformed chr
  "chrU:65-67", # malformed chr
  "chr1-5",     # malformed chr:start
  "chr8:i-5",   # malformed start
  "chr8:5-i",   # malformed end
  "chr:15",     # no end
  "chr1:6-2"    # end > start
)

test_that("verify_regions_or_stop invalid", {
  expect_error(iscream:::verify_regions_or_stop(invalid_regions, 1))
})
