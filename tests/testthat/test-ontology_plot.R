test_that("ontology_plot works", {

  results = load_example_results()[seq(100000)]
  plt <- ontology_plot(cell_type="ENS_glia", results=results)
  testthat::expect_true(methods::is(plt,"ontology_plot"))
})
