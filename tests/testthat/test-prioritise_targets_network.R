test_that("prioritise_targets_network works", {

  res <- prioritise_targets()
  vn <- prioritise_targets_network(top_targets = res$top_targets)
  testthat::expect_true(methods::is(vn$plot,"visNetwork"))
  testthat::expect_true(methods::is(vn$graph,"igraph"))


  # all_targets <- prioritise_targets(keep_tiers = NULL)
  # vn2 <- prioritise_targets_network(top_targets = all_targets,
  #                                   show_plot = FALSE)
  # testthat::expect_true(methods::is(vn2$plot,"visNetwork"))
  # testthat::expect_true(methods::is(vn2$graph,"igraph"))
})
