test_that("ggnetwork_plot_full works", {

  res_list <- ggnetwork_plot_full(filters = list(cell_type = "Cardiomyocytes"),
                                  method = "ggnetwork",
                                  interactive = TRUE)
  testthat::expect_true(methods::is(res_list$plot, "plotly"))
  testthat::expect_true(methods::is(res_list$data, "data.frame"))
})
