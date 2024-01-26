test_that("plot_bar_dendro works", {

  results <- load_example_results(multi_dataset=TRUE)
  out <- plot_bar_dendro(results = results)
  testthat::expect_true(methods::is(out$plot,"gg"))
  testthat::expect_true(methods::is(out$data,"data.table"))
})
