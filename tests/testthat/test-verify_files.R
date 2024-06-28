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

invalid_regions <- c("ch:1-5", "chr1-5", "1:1-5")
test_that("verify_regions_or_stop", {
  expect_error(iscream:::verify_regions_or_stop())
})
