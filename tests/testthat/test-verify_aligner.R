bismark1 <- paste0(letters[1:10], ".cov.gz")
bismark2 <- paste0(letters[1:10], ".cov")
biscuit1 <- paste0(letters[1:10], ".bed.gz")
biscuit2 <- paste0(letters[1:10], ".bed")
biscuit3 <- paste0(letters[1:10], ".mergecg.bed.gz")
bismark <- c(bismark1, bismark2)
biscuit <- c(biscuit1, biscuit2, biscuit3)

test_that("verify_filetype valid bismark", {
  expect_no_error(iscream:::verify_filetype(bismark, "bismark"))
})

test_that("verify_filetype valid biscuit", {
  expect_no_error(iscream:::verify_filetype(bismark, "bsbolt"))
})

test_that("verify_filetype valid", {
  expect_no_error(iscream:::verify_filetype(biscuit, "biscuit"))
})

test_that("verify_filetype warning", {
  expect_warning(
    iscream:::verify_filetype(biscuit, "bismark"),
    "'aligner' set to bismark but no files found with '.cov', extension - verify aligner and data frame colnames"
  )
  expect_warning(
    iscream:::verify_filetype(biscuit, "bsbolt"),
    "'aligner' set to bsbolt but no files found with '.cov', extension - verify aligner and data frame colnames"
  )
  expect_warning(
    iscream:::verify_filetype(bismark, "biscuit"),
    "'aligner' set to 'biscuit' but files found with '.cov', extension - verify aligner and data frame colnames"
  )
})

test_that("verify_filetype error", {
  expect_error(
    iscream:::verify_filetype(biscuit, "bismark", stop_on_error = TRUE),
    "'aligner' set to bismark but no files found with '.cov', extension - verify aligner"
  )
  expect_error(
    iscream:::verify_filetype(biscuit, "bsbolt", stop_on_error = TRUE),
    "'aligner' set to bsbolt but no files found with '.cov', extension - verify aligner"
  )
  expect_error(
    iscream:::verify_filetype(bismark, "biscuit", stop_on_error = TRUE),
    "'aligner' set to 'biscuit' but files found with '.cov', extension - verify aligner"
  )
})

test_that("verify_aligner valid", {
  expect_no_error(iscream:::verify_aligner_or_stop("biscuit"))
  expect_no_error(iscream:::verify_aligner_or_stop("bismark"))
  expect_no_error(iscream:::verify_aligner_or_stop("bsbolt"))
})

test_that("verify_aligner_or_stop invalid", {
  expect_error(iscream:::verify_aligner_or_stop("unsupported"))
})
