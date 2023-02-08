test_that("load_example_results works", {

  res1 <- load_example_results("Descartes_All_Results_extras.rds")
  testthat::expect_equal(nrow(res1), 475321)

  res2 <- load_example_results("tabulamuris_merged.rds")
  testthat::expect_equal(nrow(res2), 213028)
})
