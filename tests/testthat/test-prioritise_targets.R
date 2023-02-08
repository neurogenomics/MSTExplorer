test_that("prioritise_targets works", {

  results <- load_example_results()
  ctd <- load_example_ctd()

  #### Top only ####
  res1 <- prioritise_targets(results = results,
                             ctd = ctd)
  testthat::expect_gte(nrow(res1$top_targets), 300)

  #### All results ####
  res2 <- prioritise_targets(results = results,
                             ctd = ctd,
                             top_n = 2)
  testthat::expect_gte(nrow(res2$top_targets), 200)

})
