test_that("report_plot works", {

  results <- load_example_results()
  rep_dt <- example_targets$report
  gp <- report_plot(rep_dt=rep_dt, results=results)
  testthat::expect_true(methods::is(gp,"ggplot"))
})
