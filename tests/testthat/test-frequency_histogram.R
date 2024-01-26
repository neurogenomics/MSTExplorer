test_that("frequency_histogram works", {

  results <- load_example_results()[seq(500),]
  fp_res <- frequency_histogram(results=results)
  testthat::expect_true(methods::is(fp_res$data$pheno_plot_df,"data.table"))
  testthat::expect_true(methods::is(fp_res$data$gene_plot,"data.table"))
  testthat::expect_true(methods::is(fp_res$plot,"patchwork"))
})
