test_that("summary_plot works", {

  keep_descendants <- "Neurodevelopmental delay"
  plt_pheno_count <- summary_plot(count_var = "hpo_name",
                                  group_var = "CellType",
                                  keep_descendants = keep_descendants)
  testthat::expect_true(
    methods::is(plt_pheno_count,"plotly")
  )
  plt_cell_count <- summary_plot(count_var = "CellType",
                                  group_var = "hpo_name",
                                 keep_descendants = keep_descendants)
  testthat::expect_true(
    methods::is(plt_pheno_count,"plotly")
  )
})
