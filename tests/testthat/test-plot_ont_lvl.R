test_that("plot_ont_lvl works", {

  plts <- plot_ont_lvl()
  testthat::expect_gte(nrow(plts$data),8000)
  testthat::expect_true(methods::is(plts$plot, "ggplot"))
})
