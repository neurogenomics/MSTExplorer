test_that("load_example_results works", {
  res1 <- load_example_results("phenomix_results.tsv.gz")
  testthat::expect_equal(nrow(res1), 2206994)
})
