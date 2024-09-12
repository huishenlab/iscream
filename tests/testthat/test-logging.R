options("iscream.threads" = 1)
test_that("set_log_level to debug on single thread", {
  expect_no_error(set_log_level("debug"))
  expect_equal(get_log_level(), "debug")
})

test_that("set_log_level to info on single thread", {
  expect_no_error(set_log_level("info"))
  expect_equal(get_log_level(), "info")
})

options("iscream.threads" = 2)

test_that("set_log_level to info on 2 threads", {
  expect_no_error(set_log_level("info"))
  expect_equal(get_log_level(), "info")
})

test_that("set_log_level to debug on 2 threads", {
  expect_error(set_log_level("debug"))
})
