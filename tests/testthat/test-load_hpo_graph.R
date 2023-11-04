test_that("load_hpo_graph works", {

  g <- load_hpo_graph()
  testthat::expect_true(methods::is(g,"igraph"))
})
