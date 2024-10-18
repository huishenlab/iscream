options("iscream.threads" = 1)
test_that("set_log_level to trace on single thread", {
  expect_no_error(set_log_level("trace"))
  expect_equal(get_log_level(), "trace")
})

test_that("set_log_level to debug on single thread", {
  expect_no_error(set_log_level("debug"))
  expect_equal(get_log_level(), "debug")
})

test_that("set_log_level to info on single thread", {
  expect_no_error(set_log_level("info"))
  expect_equal(get_log_level(), "info")
})

test_that("set_log_level to debug on single thread", {
  expect_no_error(set_log_level("warn"))
  expect_equal(get_log_level(), "warn")
})

test_that("set_log_level to debug on single thread", {
  expect_no_error(set_log_level("error"))
  expect_equal(get_log_level(), "error")
})

options("iscream.threads" = 2)

test_that("set_log_level to error on 2 threads", {
  expect_no_error(set_log_level("error"))
  expect_equal(get_log_level(), "error")
})

test_that("set_log_level to info on 2 threads", {
  expect_no_error(set_log_level("info"))
  expect_equal(get_log_level(), "info")
})

test_that("set_log_level to info on 2 threads", {
  expect_no_error(set_log_level("off"))
  expect_equal(get_log_level(), "off")
})

test_that("set_log_level to warn on 2 threads", {
  expect_error(set_log_level("warn"))
})

test_that("set_log_level to debug on 2 threads", {
  expect_error(set_log_level("debug"))
})

test_that("set_log_level to warn on 2 threads", {
  expect_error(set_log_level("warn"))
})

test_that("unsupported log level", {
  expect_error(set_log_level("warn"))
})

test_that("invalid log level", {
  expect_error(set_log_level("invalid"))
})

test_that("unsupported log level", {
  expect_error(set_log_level("critical"))
})
