test_that("plot_report works", {

  results <- load_example_results()
  rep_dt <- example_targets$report
  gp <- plot_report(rep_dt=rep_dt, results=results)
  testthat::expect_true(methods::is(gp,"ggplot"))
})
