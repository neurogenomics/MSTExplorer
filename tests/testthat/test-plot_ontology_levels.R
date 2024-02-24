test_that("plot_ontology_levels works", {

  plts <- plot_ontology_levels()
  testthat::expect_gte(nrow(plts$data),8000)
  testthat::expect_true(methods::is(plts$plot, "ggplot"))
})
