test_that("is_package_loaded", {
  expect_true(iscream:::is_package_loaded("base", "NULL", fail = FALSE))
  expect_false(iscream:::is_package_loaded("does_not_exist", "NULL", fail = FALSE))
  expect_error(
    iscream:::is_package_loaded("does_not_exist", "tester", fail = TRUE),
    "Please load the 'does_not_exist' package to use 'tester' output"
  )
})
