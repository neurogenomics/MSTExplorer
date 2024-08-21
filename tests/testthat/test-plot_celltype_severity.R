test_that("plot_celltype_severity works", {

  results <- load_example_results()
  out <- plot_celltype_severity(results)

  testthat::expect_true(methods::is(out$bar$plot,"gg"))
  testthat::expect_true(methods::is(out$bar$data,"data.table"))
  testthat::expect_true(methods::is(out$dot$plot,"gg"))
  testthat::expect_true(methods::is(out$dot$data,"data.table"))
})
