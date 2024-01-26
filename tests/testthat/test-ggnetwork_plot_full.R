test_that("ggnetwork_plot_full works", {

  res_list <- ggnetwork_plot_full(cell_type = "Cardiomyocytes")
  testthat::expect_true(methods::is(res_list$plot, "gg"))
  testthat::expect_true(methods::is(res_list$data, "data.frame"))
})
