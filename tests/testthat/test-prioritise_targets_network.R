test_that("prioritise_targets_network works", {

  top_targets <- example_targets$top_targets
  top_targets[,estimate:=fold_change]
  top_targets <- map_celltype(top_targets)
  vn <- prioritise_targets_network(top_targets)
  testthat::expect_true(methods::is(vn$plot,"visNetwork"))
  testthat::expect_true(methods::is(vn$data,"tbl_graph"))


  # all_targets <- prioritise_targets(keep_tiers = NULL)
  # vn2 <- prioritise_targets_network(top_targets = all_targets,
  #                                   show_plot = FALSE)
  # testthat::expect_true(methods::is(vn2$plot,"visNetwork"))
  # testthat::expect_true(methods::is(vn2$graph,"igraph"))
})
