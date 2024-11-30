options("iscream.threads" = 2)
library("parallelly")
max_threads <- 2

test_that("set_threads within limit", {
  set_threads(1)
  expect_equal(getOption("iscream.threads"), 1)
  set_threads(max_threads)
  expect_equal(getOption("iscream.threads"), max_threads)
})

test_that("set_threads over limit", {
  expect_error(set_threads(999))
  expect_error(set_threads(availableCores() + 1))
})

test_that("set threads option over limit", {
  expect_error(check_thread_count(999, opt_set = TRUE))
})

options("iscream.threads" = availableCores())
test_that("test thread setting onAttach", {
  expect_no_error(iscream:::package_loader())
})

options("iscream.threads" = NULL)
expected_threads <- c(
  'use_threads' = 1,
  opt_set = FALSE,
  avail_threads = unname(availableCores())
)

test_that("get threads with NULL option", {
  expect_equal(get_threads(), expected_threads)
})
