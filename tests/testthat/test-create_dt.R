test_that("create_dt works", {

  dt <- create_dt(dat = mtcars)
  testthat::expect_true(methods::is(dt,"datatables"))
})
