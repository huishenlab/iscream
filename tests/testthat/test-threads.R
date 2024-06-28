options("iscream.threads" = 2)
library("parallelly")
max_threads <- availableCores()

set_get_threads <- function(nthreads) {
  set_threads(nthreads)
}

test_that("set_threads within limit", {
  set_threads(1)
  expect_equal(getOption("iscream.threads"), 1)
  set_threads(max_threads)
  expect_equal(getOption("iscream.threads"), max_threads)
})

test_that("set_threads over limit", {
  expect_error(set_threads(999))
  expect_error(set_threads(max_threads + 1))
})
