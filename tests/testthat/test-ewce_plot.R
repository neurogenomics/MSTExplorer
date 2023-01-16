test_that("ewce_plot works", {

  full_results <- EWCE::example_bootstrap_results()
  plt <- MultiEWCE::ewce_plot(total_res = full_results$results)
  testthat::expect_true(methods::is(plt$plain,"gg"))
})
