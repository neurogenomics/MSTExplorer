test_that("ggnetwork_plot_full works", {

  res_list <- ggnetwork_plot_full(cell_type = "Amacrine cells")
  testthat::expect_true(
    all(c("plot","phenos","phenoNet","adjacency") %in% names(res_list))
  )
})
