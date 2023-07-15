test_that("correlation_heatmap works", {

  top_targets <- example_targets$top_targets[1:100]
  hm <- correlation_heatmap(top_targets = top_targets)
  testthat::expect_true(methods::is(hm,"Heatmap"))
})
