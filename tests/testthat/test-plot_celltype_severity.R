test_that("plot_celltype_severity works", {

  results <- load_example_results()
  ct_select <- unique(results$CellType)[seq(10)]
  out <- plot_celltype_severity(results[CellType %in% ct_select],
                                types=c("dot","bar"))

  testthat::expect_true(methods::is(out$bar$plot,"gg"))
  testthat::expect_true(methods::is(out$bar$data,"data.table"))

  testthat::expect_true(methods::is(out$dot$plot,"gg"))
  testthat::expect_true(methods::is(out$dot$data,"data.table"))
})
