test_that("add_symptom_results works", {

  results <- load_example_results(multi_dataset = TRUE)
  results <- map_celltype(results)
  testthat::expect_false("celltype_symptom" %in% names(results))

  res1 <- add_symptom_results(results = results)
  testthat::expect_true("celltype_symptom" %in% names(res1))

  res2 <- add_symptom_results(results = results,
                              celltype_col = "cl_name")
  testthat::expect_true("celltype_symptom" %in% names(res2))
})
