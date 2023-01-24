test_that("prioritise_targets works", {

  results <- load_example_results()
  ctd <- load_example_ctd()

  #### Top only ####
  top_targets <- prioritise_targets(results = results,
                                    ctd = ctd)
  testthat::expect_gte(nrow(top_targets), 80)

  #### All results ####
  df <- prioritise_targets(results = results,
                           ctd = ctd,
                           top_n = NULL)
  testthat::expect_gte(nrow(df), 335)

})
